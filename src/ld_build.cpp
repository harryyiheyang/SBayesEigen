// Eigen LD builder: per block, correlations straight from PLINK BED or PGEN genotypes
// (bit-plane popcount kernel from CppMatrix, geno_kernel.h), then a partial symmetric
// eigendecomposition that only forms the k leading eigenvectors:
//   dsytrd (tridiagonalise once) -> dsterf (all eigenvalues, O(m^2)) -> k from the
//   variance cut -> dstemr (k eigenvectors of the tridiagonal) -> dormtr (back-transform).
// Blocks run in parallel (OpenMP), each block single-threaded, largest blocks first.
// Output: SBayesRC blockN.eigen.bin (int32 m, int32 k, float sumLambda, float thresh,
// float lambda[k], float U[m*k]) plus per-SNP frequency, N and unbiased LD scores.
#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include "geno_kernel.h"
#include "pgenlib/pgenlib_read.h"
#ifndef FCONE
#define FCONE
#endif

extern "C" void F77_NAME(dstemr)(const char* jobz, const char* range, const int* n, double* d, double* e,
                                 const double* vl, const double* vu, const int* il, const int* iu, int* m,
                                 double* w, double* z, const int* ldz, const int* nzc, int* isuppz,
                                 int* tryrac, double* work, const int* lwork, int* iwork, const int* liwork,
                                 int* info FCLEN FCLEN);

namespace {

// ---- genotype sources; each thread owns one reader ----
struct Source {
  virtual ~Source() {}
  // packed 2-bit genotypes of variant idx (4 samples per byte, first sample in the low bits)
  virtual bool read(std::size_t idx, unsigned char* dst) = 0;
};

struct BedSource : Source {
  std::FILE* fp = nullptr;
  std::size_t stride;
  BedSource(const std::string& path, std::size_t n) : stride((n + 3) / 4) { fp = std::fopen(path.c_str(), "rb"); }
  ~BedSource() { if (fp) std::fclose(fp); }
  bool read(std::size_t idx, unsigned char* dst) override {
    if (!fp) return false;
#ifdef _WIN32
    if (_fseeki64(fp, static_cast<long long>(3 + idx * stride), SEEK_SET) != 0) return false;
#else
    if (fseeko(fp, static_cast<off_t>(3 + idx * stride), SEEK_SET) != 0) return false;
#endif
    return std::fread(dst, 1, stride, fp) == stride;
  }
};

struct PgenShared {
  plink2::PgenFileInfo info;
  unsigned char* alloc = nullptr;
  std::vector<uintptr_t> offsets, nonref;
  uint32_t max_vrec_width = 0;
  uintptr_t pgr_cachelines = 0;
  std::string error;
  PgenShared(const std::string& path, std::size_t n, const std::vector<int>& allele_ct) {
    plink2::PreinitPgfi(&info);
    char err[plink2::kPglErrstrBufBlen];
    err[0] = '\0';
    plink2::PgenHeaderCtrl ctrl;
    uintptr_t cachelines;
    if (plink2::PgfiInitPhase1(path.c_str(), nullptr, UINT32_MAX, UINT32_MAX, &ctrl, &info, &cachelines, err) !=
        plink2::kPglRetSuccess) { error = err; return; }
    if (info.raw_sample_ct != n) { error = "PGEN sample count does not match its .psam"; return; }
    const std::size_t m = info.raw_variant_ct;
    if (allele_ct.size() != m) { error = "PGEN variant count does not match its .pvar"; return; }
    uint32_t mx = 2;
    for (int a : allele_ct) mx = std::max<uint32_t>(mx, a);
    info.max_allele_ct = 2;
    if (mx > 2 || (ctrl & 0x30)) {
      offsets.resize(m + 1);
      offsets[0] = 0;
      for (std::size_t j = 0; j < m; ++j) offsets[j + 1] = offsets[j] + allele_ct[j];
      info.allele_idx_offsets = offsets.data();
      info.max_allele_ct = mx;
    }
    if ((ctrl & 0xc0) == 0xc0) {
      nonref.resize(plink2::DivUp(m, plink2::kBitsPerWord) + 1);
      info.nonref_flags = nonref.data();
    }
    if (plink2::cachealigned_malloc(std::max<uintptr_t>(cachelines, 1) * plink2::kCacheline, &alloc)) {
      error = "out of memory"; return;
    }
    if (plink2::PgfiInitPhase2(ctrl, 1, 0, 0, 0, info.raw_variant_ct, &max_vrec_width, &info, alloc,
                               &pgr_cachelines, err) != plink2::kPglRetSuccess) error = err;
  }
  ~PgenShared() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(&info, &e);
    if (alloc) plink2::aligned_free(alloc);
  }
};

struct PgenSource : Source {
  plink2::PgenReader pgr;
  plink2::PgrSampleSubsetIndex pssi;
  unsigned char* alloc = nullptr;
  unsigned char* genovec = nullptr;
  std::size_t n, stride;
  bool ok = false;
  PgenSource(const std::string& path, PgenShared& sh, std::size_t n_) : n(n_), stride((n_ + 3) / 4) {
    plink2::PreinitPgr(&pgr);
    plink2::PgrSetFreadBuf(nullptr, &pgr);
    if (plink2::cachealigned_malloc(std::max<uintptr_t>(sh.pgr_cachelines, 1) * plink2::kCacheline, &alloc) ||
        plink2::cachealigned_malloc(plink2::NypCtToVecCt(n) * plink2::kBytesPerVec, &genovec)) return;
    ok = plink2::PgrInit(path.c_str(), sh.max_vrec_width, &sh.info, &pgr, alloc) == plink2::kPglRetSuccess;
    plink2::PgrClearSampleSubsetIndex(&pgr, &pssi);
  }
  ~PgenSource() {
    plink2::PglErr e = plink2::kPglRetSuccess;
    plink2::CleanupPgr(&pgr, &e);
    if (alloc) plink2::aligned_free(alloc);
    if (genovec) plink2::aligned_free(genovec);
  }
  bool read(std::size_t idx, unsigned char* dst) override {
    if (!ok) return false;
    if (plink2::PgrGet(nullptr, pssi, static_cast<uint32_t>(n), static_cast<uint32_t>(idx), &pgr,
                       reinterpret_cast<uintptr_t*>(genovec)) != plink2::kPglRetSuccess) return false;
    std::copy(genovec, genovec + stride, dst);
    return true;
  }
};

// dosage of the counted allele per 2-bit code, -1 = missing
// BED: 00 hom A1 (.bim col 5), 01 missing, 10 het, 11 hom A2 -> counts A1
// PGEN genovec: 0 hom REF, 1 het, 2 hom ALT, 3 missing -> counts REF
const std::array<int, 4> kBedCode = {2, -1, 1, 0};
const std::array<int, 4> kPgenCode = {2, 1, 0, -1};

struct BlockOut {
  std::vector<int> keep;             // 1 = polymorphic, used
  std::vector<double> freq, nobs, ld;   // ld: unbiased LD score sum_j [r2 - (1 - r2) / (n - 2)]
  int m = 0, k = 0;
  double sumLambda = 0, trR = 0, relerr = 0;
  std::string error;
};

void eigen_block(Source& src, const std::array<int, 4>& code, std::size_t n, const std::vector<int>& idx,
                 const std::string& out_file, double cut, BlockOut& res) {
  const std::size_t m0 = idx.size(), stride = (n + 3) / 4;
  static const NibbleTable bed_tab(kBedCode), pgen_tab(kPgenCode);
  const NibbleTable& tab = (code == kBedCode) ? bed_tab : pgen_tab;
  // per-byte observed dosage sum and missing count, last byte masked to the n % 4 real samples
  int dsum[256], mct[256], dsum_last[256], mct_last[256];
  const unsigned int tail = n % 4 ? n % 4 : 4;
  for (int b = 0; b < 256; ++b) {
    dsum[b] = mct[b] = dsum_last[b] = mct_last[b] = 0;
    for (unsigned int s = 0; s < 4; ++s) {
      const int x = code[(b >> (2 * s)) & 3];
      const int dd = x < 0 ? 0 : x, mm = x < 0 ? 1 : 0;
      dsum[b] += dd; mct[b] += mm;
      if (s < tail) { dsum_last[b] += dd; mct_last[b] += mm; }
    }
  }
  Planes all;
  reset_planes(all, n, m0);
  res.freq.assign(m0, 0); res.nobs.assign(m0, 0); res.keep.assign(m0, 0);
  std::vector<unsigned char> col(stride);
  for (std::size_t j = 0; j < m0; ++j) {
    if (!src.read(static_cast<std::size_t>(idx[j]), col.data())) { res.error = "genotype read failed"; return; }
    long ds = 0, mc = 0;
    for (std::size_t b = 0; b + 1 < stride; ++b) { ds += dsum[col[b]]; mc += mct[col[b]]; }
    ds += dsum_last[col[stride - 1]]; mc += mct_last[col[stride - 1]];
    const double obs = static_cast<double>(n) - mc;
    res.nobs[j] = obs;
    res.freq[j] = obs > 0 ? ds / (2.0 * obs) : NA_REAL;
    encode_packed(tab, col.data(), n, stride, all.words, all.bits.data() + 2 * j * all.words,
                  all.bits.data() + (2 * j + 1) * all.words, all.s[j], all.v[j]);
    res.keep[j] = all.v[j] > 0;
  }
  // compact to polymorphic variants
  Planes P;
  std::vector<std::size_t> kept;
  for (std::size_t j = 0; j < m0; ++j) if (res.keep[j]) kept.push_back(j);
  const int m = static_cast<int>(kept.size());
  res.m = m;
  if (m == 0) return;
  reset_planes(P, n, m);
  for (int j = 0; j < m; ++j) {
    std::copy(all.lo(kept[j]), all.lo(kept[j]) + 2 * all.words, P.bits.data() + 2 * j * P.words);
    P.s[j] = all.s[kept[j]]; P.v[j] = all.v[kept[j]];
  }
  Planes().bits.swap(all.bits);
  // R (m x m, column-major, both triangles)
  std::vector<double> A(static_cast<std::size_t>(m) * m);
  static const TileKernel kernel = select_kernel();
  cor_block(P, 0, m, P, 0, m, true, static_cast<double>(n), A.data(), m, 1, kernel);
  Planes().bits.swap(P.bits);
  // LD scores, unbiased for the population r^2: sum_j [r2 - (1 - r2) / (n - 2)]
  res.ld.assign(m, 0);
  double sumR2 = 0;
  for (int j = 0; j < m; ++j) {
    double s2 = 0, sa = 0;
    const double* c = A.data() + static_cast<std::size_t>(j) * m;
    for (int i = 0; i < m; ++i) {
      const double r2 = c[i] * c[i];
      s2 += r2;
      sa += r2 - (1 - r2) / (static_cast<double>(n) - 2);
    }
    res.ld[j] = sa; sumR2 += s2;
    res.trR += c[j];
  }
  // tridiagonalise (lower), eigenvalues, cut
  int info = 0, lwork = -1;
  std::vector<double> d(m), e(m), tau(std::max(m - 1, 1));
  double wq;
  F77_CALL(dsytrd)("L", &m, A.data(), &m, d.data(), e.data(), tau.data(), &wq, &lwork, &info FCONE);
  lwork = static_cast<int>(wq);
  std::vector<double> work(std::max(lwork, 1));
  F77_CALL(dsytrd)("L", &m, A.data(), &m, d.data(), e.data(), tau.data(), work.data(), &lwork, &info FCONE);
  if (info) { res.error = "dsytrd failed"; return; }
  std::vector<double> lam(d), e2(e);
  F77_CALL(dsterf)(&m, lam.data(), e2.data(), &info);
  if (info) { res.error = "dsterf failed"; return; }
  std::reverse(lam.begin(), lam.end());   // descending
  int np = 0;
  double sl = 0;
  while (np < m && lam[np] > 0) sl += lam[np++];
  double cs = 0;
  int k = np;
  for (int j = 0; j < np; ++j) { cs += lam[j]; if (cs >= cut * sl) { k = j + 1; break; } }
  if (k < 1) { res.error = "no positive eigenvalues"; return; }
  res.k = k; res.sumLambda = sl;
  double sumD2 = 0;
  for (int j = 0; j < k; ++j) sumD2 += lam[j] * lam[j];
  res.relerr = std::sqrt(std::max(sumR2 - sumD2, 0.0) / sumR2);
  // k leading eigenvectors of the tridiagonal (dstemr returns ascending order)
  const int il = m - k + 1, iu = m, nzc = k;
  int mfound = 0, tryrac = 1, liwork = -1;
  lwork = -1;
  std::vector<double> w(m), Z(static_cast<std::size_t>(m) * k), dd(d), ee(e);
  std::vector<int> isuppz(2 * static_cast<std::size_t>(k));
  int iwq;
  const double vl = 0, vu = 0;
  F77_CALL(dstemr)("V", "I", &m, dd.data(), ee.data(), &vl, &vu, &il, &iu, &mfound, w.data(), Z.data(), &m,
                   &nzc, isuppz.data(), &tryrac, &wq, &lwork, &iwq, &liwork, &info FCONE FCONE);
  lwork = static_cast<int>(wq); liwork = iwq;
  work.assign(std::max(lwork, 1), 0);
  std::vector<int> iwork(std::max(liwork, 1));
  F77_CALL(dstemr)("V", "I", &m, dd.data(), ee.data(), &vl, &vu, &il, &iu, &mfound, w.data(), Z.data(), &m,
                   &nzc, isuppz.data(), &tryrac, work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE);
  if (info || mfound != k) { res.error = "dstemr failed"; return; }
  // back-transform: U = Q Z
  lwork = -1;
  F77_CALL(dormtr)("L", "L", "N", &m, &k, A.data(), &m, tau.data(), Z.data(), &m, &wq, &lwork, &info
                   FCONE FCONE FCONE);
  lwork = static_cast<int>(wq);
  work.assign(std::max(lwork, 1), 0);
  F77_CALL(dormtr)("L", "L", "N", &m, &k, A.data(), &m, tau.data(), Z.data(), &m, work.data(), &lwork, &info
                   FCONE FCONE FCONE);
  if (info) { res.error = "dormtr failed"; return; }
  std::vector<double>().swap(A);
  // write descending: column c of U is column k - 1 - c of Z
  std::vector<float> buf(static_cast<std::size_t>(m) * k);
  for (int c = 0; c < k; ++c) {
    const double* z = Z.data() + static_cast<std::size_t>(k - 1 - c) * m;
    float* u = buf.data() + static_cast<std::size_t>(c) * m;
    for (int i = 0; i < m; ++i) u[i] = static_cast<float>(z[i]);
  }
  std::FILE* fp = std::fopen(out_file.c_str(), "wb");
  if (!fp) { res.error = "cannot write " + out_file; return; }
  const int32_t hdr[2] = {m, k};
  const float fh[2] = {static_cast<float>(sl), static_cast<float>(cut)};
  std::vector<float> lk(lam.begin(), lam.begin() + k);
  bool ok = std::fwrite(hdr, 4, 2, fp) == 2 && std::fwrite(fh, 4, 2, fp) == 2 &&
            std::fwrite(lk.data(), 4, k, fp) == static_cast<std::size_t>(k) &&
            std::fwrite(buf.data(), 4, buf.size(), fp) == buf.size();
  ok = (std::fclose(fp) == 0) && ok;
  if (!ok) res.error = "write failed: " + out_file;
}

} // namespace

// geno: path to .bed or .pgen; blocks: list of 0-based variant indices (file order);
// out_files: one eigen.bin per block; allele_ct: per-variant allele counts (PGEN only).
// [[Rcpp::export]]
Rcpp::List ld_build_cpp(std::string geno, bool pgen, int n_samples, Rcpp::IntegerVector allele_ct,
                        Rcpp::List blocks, Rcpp::CharacterVector out_files, double cut, int threads) {
  const int nb = blocks.size();
  std::vector<std::vector<int>> idx(nb);
  std::vector<std::string> outs(nb);
  for (int b = 0; b < nb; ++b) {
    idx[b] = Rcpp::as<std::vector<int>>(blocks[b]);
    outs[b] = Rcpp::as<std::string>(out_files[b]);
  }
  const std::size_t n = static_cast<std::size_t>(n_samples);
#ifdef _OPENMP
  threads = std::max(1, std::min(threads, nb));
#else
  threads = 1;
#endif
  std::unique_ptr<PgenShared> sh;
  if (pgen) {
    sh.reset(new PgenShared(geno, n, Rcpp::as<std::vector<int>>(allele_ct)));
    if (!sh->error.empty()) Rcpp::stop("Cannot open PGEN " + geno + ": " + sh->error);
  }
  std::vector<std::unique_ptr<Source>> src(threads);
  for (int t = 0; t < threads; ++t) {
    if (pgen) {
      PgenSource* p = new PgenSource(geno, *sh, n);
      src[t].reset(p);
      if (!p->ok) Rcpp::stop("Cannot open PGEN reader for " + geno);
    } else {
      BedSource* p = new BedSource(geno, n);
      src[t].reset(p);
      if (!p->fp) Rcpp::stop("Cannot open " + geno);
    }
  }
  std::vector<int> order(nb);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return idx[a].size() > idx[b].size(); });
  std::vector<BlockOut> res(nb);
  const std::array<int, 4>& code = pgen ? kPgenCode : kBedCode;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int o = 0; o < nb; ++o) {
    const int b = order[o];
#ifdef _OPENMP
    Source& s = *src[omp_get_thread_num()];
#else
    Source& s = *src[0];
#endif
    try {
      eigen_block(s, code, n, idx[b], outs[b], cut, res[b]);
    } catch (const std::exception& ex) {
      res[b].error = ex.what();
    } catch (...) {
      res[b].error = "unknown C++ error";
    }
  }
  Rcpp::List out(nb);
  for (int b = 0; b < nb; ++b) {
    if (!res[b].error.empty()) Rcpp::stop("block file " + outs[b] + ": " + res[b].error);
    const BlockOut& r = res[b];
    out[b] = Rcpp::List::create(
      Rcpp::_["keep"] = Rcpp::LogicalVector(r.keep.begin(), r.keep.end()), Rcpp::_["freq"] = r.freq,
      Rcpp::_["nobs"] = r.nobs, Rcpp::_["ldscore"] = r.ld, Rcpp::_["m"] = r.m, Rcpp::_["k"] = r.k,
      Rcpp::_["sumLambda"] = r.sumLambda, Rcpp::_["trR"] = r.trR, Rcpp::_["relerr_F"] = r.relerr);
  }
  return out;
}
