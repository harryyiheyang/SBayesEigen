// Eigen LD builder: per block, correlations straight from PLINK BED or PGEN genotypes
// (bit-plane popcount kernel from CppMatrix, geno_kernel.h), then a partial symmetric
// eigendecomposition that only forms the k leading eigenvectors:
//   dsytrd (tridiagonalise once) -> dstedc (eigenpairs of the tridiagonal) -> k from the
//   variance cut -> dormtr (back-transform of the k leading vectors).
// Scheduling: the largest blocks, whose single-threaded cost would set the wall time, run one at a time with
// all threads (correlation kernel and BLAS); the rest run in parallel, one thread each, under an optional
// memory cap. Largest blocks first.
// Output: SBayesRC blockN.eigen.bin (int32 m, int32 k, float sumLambda, float thresh,
// float lambda[k], float U[m*k]) plus per-SNP frequency, N and unbiased LD scores.
// With a fixed SNP set A (AB mode) each block is written as blockN.ab.bin instead: R_AA kept, B eigen-
// decomposed after an exact Schur complement with the generalized inverse of R_AA, and R_BA kept so R can be
// rebuilt (layout at write_ab below).
#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#if !defined(_WIN32)
#include <dlfcn.h>
#endif
#include "geno_kernel.h"
#include "pgenlib/pgenlib_read.h"
#ifndef FCONE
#define FCONE
#endif

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
  int m = 0, k = 0, mA = 0, rankA = 0;
  double sumLambda = 0, trR = 0, relerr = 0;
  std::string error;
};

// k leading eigenpairs of the symmetric m x m matrix A (lower triangle used, overwritten): lam = all
// eigenvalues descending, sl = sum of the positive ones, k = first count reaching cut * sl, Z = m x k
// eigenvectors in ascending order (column k - 1 is the largest). The tridiagonal is solved by divide and
// conquer (dstedc, all m vectors, BLAS-3 and so threaded with the BLAS); with cut 0.995 k is close to m, where
// this is far faster than MRRR for k vectors. Peak memory 24 m^2 bytes (A, Z and the dstedc workspace).
bool top_eigen(std::vector<double>& A, int m, double cut, std::vector<double>& lam, std::vector<double>& Z, int& k,
               double& sl, std::string& error) {
  int info = 0, lwork = -1;
  std::vector<double> d(m), e(m), tau(std::max(m - 1, 1));
  double wq;
  F77_CALL(dsytrd)("L", &m, A.data(), &m, d.data(), e.data(), tau.data(), &wq, &lwork, &info FCONE);
  lwork = static_cast<int>(wq);
  std::vector<double> work(std::max(lwork, 1));
  F77_CALL(dsytrd)("L", &m, A.data(), &m, d.data(), e.data(), tau.data(), work.data(), &lwork, &info FCONE);
  if (info) { error = "dsytrd failed"; return false; }
  std::vector<double>().swap(work);
  // all eigenpairs of the tridiagonal, ascending
  Z.assign(static_cast<std::size_t>(m) * m, 0);
  {
    int liwork = -1, iwq;
    lwork = -1;
    F77_CALL(dstedc)("I", &m, d.data(), e.data(), Z.data(), &m, &wq, &lwork, &iwq, &liwork, &info FCONE);
    lwork = static_cast<int>(wq); liwork = iwq;
    std::vector<double> w(std::max(lwork, 1));
    std::vector<int> iw(std::max(liwork, 1));
    F77_CALL(dstedc)("I", &m, d.data(), e.data(), Z.data(), &m, w.data(), &lwork, iw.data(), &liwork, &info FCONE);
    if (info) { error = "dstedc failed"; return false; }
  }
  lam.assign(d.rbegin(), d.rend());   // descending
  int np = 0;
  sl = 0;
  while (np < m && lam[np] > 0) sl += lam[np++];
  double cs = 0;
  k = np;
  for (int j = 0; j < np; ++j) { cs += lam[j]; if (cs >= cut * sl) { k = j + 1; break; } }
  if (k < 1) { error = "no positive eigenvalues"; return false; }
  // keep the k leading vectors (the last k columns), still ascending
  const std::size_t off = static_cast<std::size_t>(m - k) * m;
  if (off) std::copy(Z.begin() + off, Z.end(), Z.begin());
  Z.resize(static_cast<std::size_t>(m) * k);
  // back-transform: Z = Q Z
  lwork = -1;
  F77_CALL(dormtr)("L", "L", "N", &m, &k, A.data(), &m, tau.data(), Z.data(), &m, &wq, &lwork, &info
                   FCONE FCONE FCONE);
  lwork = static_cast<int>(wq);
  work.assign(std::max(lwork, 1), 0);
  F77_CALL(dormtr)("L", "L", "N", &m, &k, A.data(), &m, tau.data(), Z.data(), &m, work.data(), &lwork, &info
                   FCONE FCONE FCONE);
  if (info) { error = "dormtr failed"; return false; }
  return true;
}

// AB mode. R (m x m, both triangles) with a_kept[j] = 1 for SNPs in A. R_AA = Q_A L_A Q_A'; the generalized
// inverse keeps eigenvalues > tolA * max, R_AA^+ = H H' with H = Q_A L_A^-1/2 (rankA columns). G = R_BA H,
// S = R_BB - G G' = Q2 L2 Q2' cut at `cut`; UlB has B rows Q2 and A rows -H G' Q2.
// blockN.ab.bin: int32 m, mA, kB, rankA; float tolA, cutB, sumLambdaB; int32 idxA[mA] (0-based in block);
// float R_AA upper triangle packed column-major (mA(mA+1)/2); float lambdaB[kB]; float UlB[m*kB] column-major;
// float R_BA[(m-mA)*mA] column-major (B rows in block order).
void ab_block(std::vector<double>& R, int m, const std::vector<int>& a_kept, double sumR2, const std::string& out_file,
              double cut, double tolA, BlockOut& res) {
  std::vector<int> iA, iB;
  for (int j = 0; j < m; ++j) (a_kept[j] ? iA : iB).push_back(j);
  const int mA = static_cast<int>(iA.size()), mB = static_cast<int>(iB.size());
  res.mA = mA;
  auto at = [&](int i, int j) { return R[static_cast<std::size_t>(j) * m + i]; };
  std::vector<double> RAA(static_cast<std::size_t>(mA) * mA), RBA(static_cast<std::size_t>(mB) * mA);
  for (int c = 0; c < mA; ++c) {
    for (int r = 0; r < mA; ++r) RAA[static_cast<std::size_t>(c) * mA + r] = at(iA[r], iA[c]);
    for (int r = 0; r < mB; ++r) RBA[static_cast<std::size_t>(c) * mB + r] = at(iB[r], iA[c]);
  }
  // R_BB compacted in place to the leading mB x mB (targets never pass their sources)
  for (int c = 0; c < mB; ++c)
    for (int r = 0; r < mB; ++r) R[static_cast<std::size_t>(c) * mB + r] = at(iB[r], iB[c]);
  R.resize(static_cast<std::size_t>(mB) * mB);
  int info = 0, rankA = 0;
  std::vector<double> H;   // mA x rankA
  if (mA > 0) {
    std::vector<double> Q(RAA), la(mA);
    int lwork = -1;
    double wq;
    F77_CALL(dsyev)("V", "L", &mA, Q.data(), &mA, la.data(), &wq, &lwork, &info FCONE FCONE);
    lwork = static_cast<int>(wq);
    std::vector<double> work(std::max(lwork, 1));
    F77_CALL(dsyev)("V", "L", &mA, Q.data(), &mA, la.data(), work.data(), &lwork, &info FCONE FCONE);
    if (info) { res.error = "dsyev of R_AA failed"; return; }
    const double mx = *std::max_element(la.begin(), la.end());
    for (int j = 0; j < mA; ++j) {
      if (!(la[j] > tolA * mx)) continue;
      const double s = 1 / std::sqrt(la[j]);
      for (int r = 0; r < mA; ++r) H.push_back(Q[static_cast<std::size_t>(j) * mA + r] * s);
      ++rankA;
    }
  }
  res.rankA = rankA;
  std::vector<double> G(static_cast<std::size_t>(mB) * std::max(rankA, 1));
  const double one = 1, zero = 0, mone = -1;
  if (rankA > 0 && mB > 0) {
    F77_CALL(dgemm)("N", "N", &mB, &rankA, &mA, &one, RBA.data(), &mB, H.data(), &mA, &zero, G.data(), &mB
                    FCONE FCONE);
    F77_CALL(dsyrk)("L", "N", &mB, &rankA, &mone, G.data(), &mB, &one, R.data(), &mB FCONE FCONE);
    for (int c = 0; c < mB; ++c)   // symmetric for the S2 sum below
      for (int r = 0; r < c; ++r) R[static_cast<std::size_t>(c) * mB + r] = R[static_cast<std::size_t>(r) * mB + c];
  }
  std::vector<double> lam, Z;
  int k = 0;
  double sl = 0;
  if (mB > 0) {
    double sumS2 = 0;
    for (double x : R) sumS2 += x * x;
    if (!top_eigen(R, mB, cut, lam, Z, k, sl, res.error)) return;
    double sumD2 = 0;
    for (int j = 0; j < k; ++j) sumD2 += lam[j] * lam[j];
    res.relerr = sumS2 > 0 ? std::sqrt(std::max(sumS2 - sumD2, 0.0) / sumS2) : 0;
  }
  std::vector<double>().swap(R);
  res.k = k; res.sumLambda = sl;
  // A rows of UlB: -H G' Q2, computed on the ascending columns of Z (column c of UlB is column k - 1 - c)
  std::vector<double> UA;
  if (rankA > 0 && k > 0) {
    std::vector<double> T(static_cast<std::size_t>(rankA) * k);
    UA.resize(static_cast<std::size_t>(mA) * k);
    F77_CALL(dgemm)("T", "N", &rankA, &k, &mB, &one, G.data(), &mB, Z.data(), &mB, &zero, T.data(), &rankA
                    FCONE FCONE);
    F77_CALL(dgemm)("N", "N", &mA, &k, &rankA, &mone, H.data(), &mA, T.data(), &rankA, &zero, UA.data(), &mA
                    FCONE FCONE);
  }
  std::FILE* fp = std::fopen(out_file.c_str(), "wb");
  if (!fp) { res.error = "cannot write " + out_file; return; }
  const int32_t hdr[4] = {m, mA, k, rankA};
  const float fh[3] = {static_cast<float>(tolA), static_cast<float>(cut), static_cast<float>(sl)};
  std::vector<float> up;
  up.reserve(static_cast<std::size_t>(mA) * (mA + 1) / 2);
  for (int c = 0; c < mA; ++c) for (int r = 0; r <= c; ++r) up.push_back(static_cast<float>(RAA[static_cast<std::size_t>(c) * mA + r]));
  std::vector<float> lk(lam.begin(), lam.begin() + k), rba(RBA.begin(), RBA.end());
  bool ok = std::fwrite(hdr, 4, 4, fp) == 4 && std::fwrite(fh, 4, 3, fp) == 3 &&
            std::fwrite(iA.data(), 4, mA, fp) == static_cast<std::size_t>(mA) &&
            std::fwrite(up.data(), 4, up.size(), fp) == up.size() &&
            std::fwrite(lk.data(), 4, k, fp) == static_cast<std::size_t>(k);
  // UlB, descending columns, one column at a time
  std::vector<float> col(m, 0.0f);
  for (int c = 0; ok && c < k; ++c) {
    const std::size_t zc = static_cast<std::size_t>(k - 1 - c);
    const double* z = Z.data() + zc * mB;
    for (int r = 0; r < mB; ++r) col[iB[r]] = static_cast<float>(z[r]);
    if (!UA.empty()) for (int r = 0; r < mA; ++r) col[iA[r]] = static_cast<float>(UA[zc * mA + r]);
    ok = std::fwrite(col.data(), 4, m, fp) == static_cast<std::size_t>(m);
  }
  ok = ok && std::fwrite(rba.data(), 4, rba.size(), fp) == rba.size();
  ok = (std::fclose(fp) == 0) && ok;
  if (!ok) res.error = "write failed: " + out_file;
}

void eigen_block(Source& src, const std::array<int, 4>& code, std::size_t n, const std::vector<int>& idx,
                 const std::vector<int>& isA, const std::string& out_file, double cut, double tolA, int threads,
                 BlockOut& res) {
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
  cor_block(P, 0, m, P, 0, m, true, static_cast<double>(n), A.data(), m, threads, kernel);
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
  if (!isA.empty()) {
    std::vector<int> a_kept(m);
    for (int j = 0; j < m; ++j) a_kept[j] = isA[kept[j]];
    ab_block(A, m, a_kept, sumR2, out_file, cut, tolA, res);
    return;
  }
  std::vector<double> lam, Z;
  int k = 0;
  double sl = 0;
  if (!top_eigen(A, m, cut, lam, Z, k, sl, res.error)) return;
  res.k = k; res.sumLambda = sl;
  double sumD2 = 0;
  for (int j = 0; j < k; ++j) sumD2 += lam[j] * lam[j];
  res.relerr = std::sqrt(std::max(sumR2 - sumD2, 0.0) / sumR2);
  std::vector<double>().swap(A);
  std::FILE* fp = std::fopen(out_file.c_str(), "wb");
  if (!fp) { res.error = "cannot write " + out_file; return; }
  const int32_t hdr[2] = {m, k};
  const float fh[2] = {static_cast<float>(sl), static_cast<float>(cut)};
  std::vector<float> lk(lam.begin(), lam.begin() + k);
  bool ok = std::fwrite(hdr, 4, 2, fp) == 2 && std::fwrite(fh, 4, 2, fp) == 2 &&
            std::fwrite(lk.data(), 4, k, fp) == static_cast<std::size_t>(k);
  // U descending, one column at a time: column c of U is column k - 1 - c of Z
  std::vector<float> ucol(m);
  for (int c = 0; ok && c < k; ++c) {
    const double* z = Z.data() + static_cast<std::size_t>(k - 1 - c) * m;
    for (int i = 0; i < m; ++i) ucol[i] = static_cast<float>(z[i]);
    ok = std::fwrite(ucol.data(), 4, m, fp) == static_cast<std::size_t>(m);
  }
  ok = (std::fclose(fp) == 0) && ok;
  if (!ok) res.error = "write failed: " + out_file;
}

// BLAS thread count, set at run time through whichever control the loaded BLAS exports (OpenBLAS, MKL,
// FlexiBLAS); without one (reference BLAS) the BLAS stays single-threaded. Restored on exit.
struct BlasThreads {
  void (*set)(int) = nullptr;
  int orig = 0;
  BlasThreads() {
#if !defined(_WIN32)
    const char* sets[] = {"openblas_set_num_threads", "MKL_Set_Num_Threads", "flexiblas_set_num_threads"};
    const char* gets[] = {"openblas_get_num_threads", "MKL_Get_Max_Threads", "flexiblas_get_num_threads"};
    for (int i = 0; i < 3 && !set; ++i) {
      set = reinterpret_cast<void (*)(int)>(dlsym(RTLD_DEFAULT, sets[i]));
      if (set) {
        int (*get)() = reinterpret_cast<int (*)()>(dlsym(RTLD_DEFAULT, gets[i]));
        orig = get ? get() : 0;
      }
    }
#endif
  }
  void operator()(int t) const { if (set) set(t); }
  ~BlasThreads() { if (set && orig > 0) set(orig); }
};

// bytes a block of m variants holds at its peak (R, the tridiagonal's vectors and the dstedc workspace)
double block_bytes(double m) { return 24 * m * m; }

// blocks running in parallel wait here until their peak memory fits under the cap (one always runs)
struct MemGate {
  double cap, used = 0;
  std::mutex mu;
  std::condition_variable cv;
  explicit MemGate(double c) : cap(c) {}
  void acquire(double b) {
    if (cap <= 0) return;
    std::unique_lock<std::mutex> lk(mu);
    cv.wait(lk, [&] { return used == 0 || used + b <= cap; });
    used += b;
  }
  void release(double b) {
    if (cap <= 0) return;
    { std::lock_guard<std::mutex> lk(mu); used -= b; }
    cv.notify_all();
  }
};

} // namespace

// geno: path to .bed or .pgen; blocks: list of 0-based variant indices (file order);
// out_files: one eigen.bin per block; allele_ct: per-variant allele counts (PGEN only).
// a_mask: empty list (eigen.bin) or, per block, 0/1 per variant of `blocks` marking set A (ab.bin).
// mem: memory cap in bytes for the blocks in flight (0 = none).
// [[Rcpp::export]]
Rcpp::List ld_build_cpp(std::string geno, bool pgen, int n_samples, Rcpp::IntegerVector allele_ct,
                        Rcpp::List blocks, Rcpp::CharacterVector out_files, double cut, int threads,
                        Rcpp::List a_mask, double tolA, double mem) {
  const int nb = blocks.size();
  std::vector<std::vector<int>> idx(nb), isA(nb);
  std::vector<std::string> outs(nb);
  const bool ab = a_mask.size() > 0;
  if (ab && a_mask.size() != nb) Rcpp::stop("a_mask must have one entry per block");
  for (int b = 0; b < nb; ++b) {
    idx[b] = Rcpp::as<std::vector<int>>(blocks[b]);
    outs[b] = Rcpp::as<std::string>(out_files[b]);
    if (ab) {
      isA[b] = Rcpp::as<std::vector<int>>(a_mask[b]);
      if (isA[b].size() != idx[b].size()) Rcpp::stop("a_mask and blocks differ in length");
    }
  }
  const std::size_t n = static_cast<std::size_t>(n_samples);
#ifdef _OPENMP
  threads = std::max(1, threads);
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
  // serial part: a block runs alone with all threads while its O(m^3) cost exceeds an equal share of what is
  // left (only when the BLAS can be threaded), or while `threads` copies of it would not fit under the memory cap
  const BlasThreads blas;
  double left = 0;
  for (int b = 0; b < nb; ++b) left += std::pow(static_cast<double>(idx[b].size()), 3);
  int nbig = 0;
  if (threads > 1)
    for (; nbig < nb; ++nbig) {
      const double m = static_cast<double>(idx[order[nbig]].size()), c = m * m * m;
      if (!((blas.set && c > left / threads) || (mem > 0 && threads * block_bytes(m) > mem))) break;
      left -= c;
    }
  std::vector<BlockOut> res(nb);
  const std::array<int, 4>& code = pgen ? kPgenCode : kBedCode;
  auto run = [&](int b, Source& s, int nt) {
    try {
      eigen_block(s, code, n, idx[b], isA[b], outs[b], cut, tolA, nt, res[b]);
    } catch (const std::exception& ex) {
      res[b].error = ex.what();
    } catch (...) {
      res[b].error = "unknown C++ error";
    }
  };
  blas(threads);
  for (int o = 0; o < nbig; ++o) run(order[o], *src[0], threads);
  blas(1);
  MemGate gate(mem);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(std::max(1, std::min(threads, nb - nbig))) schedule(dynamic, 1)
#endif
  for (int o = nbig; o < nb; ++o) {
    const int b = order[o];
#ifdef _OPENMP
    Source& s = *src[omp_get_thread_num()];
#else
    Source& s = *src[0];
#endif
    const double need = block_bytes(static_cast<double>(idx[b].size()));
    gate.acquire(need);
    run(b, s, 1);
    gate.release(need);
  }
  Rcpp::List out(nb);
  for (int b = 0; b < nb; ++b) {
    if (!res[b].error.empty()) Rcpp::stop("block file " + outs[b] + ": " + res[b].error);
    const BlockOut& r = res[b];
    out[b] = Rcpp::List::create(
      Rcpp::_["keep"] = Rcpp::LogicalVector(r.keep.begin(), r.keep.end()), Rcpp::_["freq"] = r.freq,
      Rcpp::_["nobs"] = r.nobs, Rcpp::_["ldscore"] = r.ld, Rcpp::_["m"] = r.m, Rcpp::_["k"] = r.k,
      Rcpp::_["sumLambda"] = r.sumLambda, Rcpp::_["trR"] = r.trR, Rcpp::_["relerr_F"] = r.relerr,
      Rcpp::_["mA"] = r.mA, Rcpp::_["rankA"] = r.rankA);
  }
  return out;
}
