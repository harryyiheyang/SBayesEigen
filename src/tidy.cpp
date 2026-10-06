// Fast tidy of COJO summary data against snp.info, same QC as SBayesRC::tidy:
// whole files read into memory, lines split in parallel, SNPs matched with a hash table on
// string views (no R strings until the result), output lines copied from the input text.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <cstdint>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

std::string slurp(const std::string& path) {
  std::FILE* fp = std::fopen(path.c_str(), "rb");
  if (!fp) Rcpp::stop("cannot open " + path);
  std::string s;
  if (std::fseek(fp, 0, SEEK_END) == 0) {
    const long sz = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);
    if (sz > 0) {
      s.resize(static_cast<size_t>(sz) + 1);
      const size_t r = std::fread(&s[0], 1, static_cast<size_t>(sz), fp);
      s.resize(r);
    }
  }
  char buf[1 << 16];
  size_t r;
  while ((r = std::fread(buf, 1, sizeof buf, fp)) > 0) s.append(buf, r);
  std::fclose(fp);
  if (!s.empty() && s.back() != '\n') s.push_back('\n');
  return s;
}

constexpr int kMaxF = 256;
inline bool is_sep(char c) { return c == ' ' || c == '\t' || c == '\r'; }

inline std::string_view unquote(const char* p, const char* q) {
  while (q > p && (q[-1] == '\r' || q[-1] == ' ')) --q;
  while (p < q && *p == ' ') ++p;
  if (q - p >= 2 && *p == '"' && q[-1] == '"') { ++p; --q; }
  return std::string_view(p, q - p);
}

// split one line: on tabs when tab_only (header has a tab), else on runs of whitespace;
// surrounding double quotes are stripped
inline int split(const char* p, const char* end, std::string_view* f, int maxf, bool tab_only) {
  int n = 0;
  if (tab_only) {
    while (n < maxf) {
      const char* q = p;
      while (q < end && *q != '\t') ++q;
      f[n++] = unquote(p, q);
      if (q >= end) break;
      p = q + 1;
    }
    if (n == 1 && f[0].empty()) n = 0;
    return n;
  }
  while (p < end && n < maxf) {
    while (p < end && is_sep(*p)) ++p;
    if (p >= end) break;
    const char* q = p;
    while (q < end && !is_sep(*q)) ++q;
    f[n++] = unquote(p, q);
    p = q;
  }
  return n;
}

std::vector<size_t> line_starts(const std::string& s, size_t from) {
  std::vector<size_t> st;
  st.reserve(s.size() / 64);
  size_t p = from;
  while (p < s.size()) {
    const void* nl = std::memchr(s.data() + p, '\n', s.size() - p);
    const size_t e = nl ? static_cast<const char*>(nl) - s.data() : s.size();
    size_t q = p;
    while (q < e && is_sep(s[q])) ++q;
    if (q < e) st.push_back(p);
    p = e + 1;
  }
  return st;
}

inline const char* line_end(const std::string& s, size_t start) {
  const void* nl = std::memchr(s.data() + start, '\n', s.size() - start);
  return nl ? static_cast<const char*>(nl) : s.data() + s.size();
}

inline double num(std::string_view v) {
  char buf[64];
  if (v.empty() || v.size() >= sizeof buf) return NAN;
  std::memcpy(buf, v.data(), v.size());
  buf[v.size()] = '\0';
  char* e;
  const double x = std::strtod(buf, &e);
  return (e == buf + v.size()) ? x : NAN;
}

int find_col(const std::vector<std::string_view>& h, std::initializer_list<const char*> names) {
  for (const char* nm : names)
    for (size_t i = 0; i < h.size(); ++i) if (h[i] == nm) return static_cast<int>(i);
  return -1;
}

std::vector<std::string_view> header(const std::string& s, bool& tab_only) {
  const char* e = line_end(s, 0);
  tab_only = std::memchr(s.data(), '\t', e - s.data()) != nullptr;
  std::string_view f[kMaxF];
  const int n = split(s.data(), e, f, kMaxF, tab_only);
  return std::vector<std::string_view>(f, f + n);
}

double median(std::vector<double> x) {
  const size_t n = x.size(), h = n / 2;
  std::nth_element(x.begin(), x.begin() + h, x.end());
  const double hi = x[h];
  if (n % 2) return hi;
  return 0.5 * (hi + *std::max_element(x.begin(), x.begin() + h));
}

inline uint64_t fnv(std::string_view v) {
  uint64_t h = 1469598103934665603ULL;
  for (unsigned char ch : v) { h ^= ch; h *= 1099511628211ULL; }
  return h ^ (h >> 29);
}

// open addressing, linear probing; slots hold row index + 1 (0 = empty); first copy of an ID wins
struct FlatMap {
  const std::vector<std::string_view>& key;
  std::vector<int> slot;
  uint64_t mask;
  FlatMap(const std::vector<std::string_view>& k, int threads) : key(k) {
    const size_t n = k.size();
    size_t cap = 1;
    while (cap < 2 * n + 2) cap <<= 1;
    slot.assign(cap, 0);
    mask = cap - 1;
    std::vector<uint64_t> h(n);
#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (long i = 0; i < static_cast<long>(n); ++i) h[i] = fnv(k[i]);
    for (size_t i = 0; i < n; ++i) {
      uint64_t j = h[i] & mask;
      while (slot[j] && key[slot[j] - 1] != k[i]) j = (j + 1) & mask;
      if (!slot[j]) slot[j] = static_cast<int>(i) + 1;
    }
  }
  int find(std::string_view v) const {
    uint64_t j = fnv(v) & mask;
    while (slot[j]) {
      if (key[slot[j] - 1] == v) return slot[j] - 1;
      j = (j + 1) & mask;
    }
    return -1;
  }
};

struct Row { int si; size_t line; double freq, b, se, p, N; bool flip, allele_ok, freq_ok; };

} // namespace

// [[Rcpp::export]]
Rcpp::List tidy_cpp(std::string mafile, std::string snpinfo, std::string output, double freq_thresh,
                    double N_sd_range, double rate2pq, bool want_strings,
                    int threads) {
#ifdef _OPENMP
  threads = std::max(1, threads);
#else
  threads = 1;
#endif
  // ---- snp.info: ID, A1, A2, A1Freq ----
  const std::string si = slurp(snpinfo);
  bool st;
  const auto sh = header(si, st);
  int cID = find_col(sh, {"ID", "SNP"}), cA1 = find_col(sh, {"A1"}), cA2 = find_col(sh, {"A2"}),
      cF = find_col(sh, {"A1Freq", "freq", "A1Frq"});
  if (cID < 0 || cA1 < 0 || cA2 < 0 || cF < 0) Rcpp::stop("snp.info needs columns ID, A1, A2, A1Freq");
  const std::vector<size_t> sl = line_starts(si, line_end(si, 0) - si.data() + 1);
  const int M = static_cast<int>(sl.size());
  const int needs = std::max(std::max(cID, cA1), std::max(cA2, cF)) + 1;
  if (needs > kMaxF) Rcpp::stop("snp.info: needed columns beyond column 256");
  std::vector<std::string_view> sid(M), sa1(M), sa2(M);
  std::vector<double> sfreq(M);
  int bad = 0;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:bad)
#endif
  for (int i = 0; i < M; ++i) {
    std::string_view f[kMaxF];
    const int n = split(si.data() + sl[i], line_end(si, sl[i]), f, needs, st);
    if (n < needs) { ++bad; continue; }
    sid[i] = f[cID]; sa1[i] = f[cA1]; sa2[i] = f[cA2]; sfreq[i] = num(f[cF]);
  }
  if (bad) Rcpp::stop("snp.info has short lines");
  FlatMap map(sid, threads);

  // ---- summary data ----
  const std::string ma = slurp(mafile);
  bool mt;
  const auto mh = header(ma, mt);
  const int c[8] = {find_col(mh, {"SNP"}), find_col(mh, {"A1"}), find_col(mh, {"A2"}), find_col(mh, {"freq"}),
                    find_col(mh, {"b"}), find_col(mh, {"se"}), find_col(mh, {"p"}), find_col(mh, {"N"})};
  for (int k = 0; k < 8; ++k) if (c[k] < 0) Rcpp::stop("The summary data is not a valid COJO format (SNP A1 A2 freq b se p N)");
  const int mneeds = *std::max_element(c, c + 8) + 1;
  if (mneeds > kMaxF) Rcpp::stop("summary data: needed columns beyond column 256");
  const std::vector<size_t> ml = line_starts(ma, line_end(ma, 0) - ma.data() + 1);
  const long nma = static_cast<long>(ml.size());
  // per line: -3 short/invalid, -2 valid but not in LD, else the matched snp.info row
  std::vector<int> hit(nma);
  std::vector<Row> rows(nma);
  long n_valid = 0;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:n_valid)
#endif
  for (long i = 0; i < nma; ++i) {
    std::string_view f[kMaxF];
    const int n = split(ma.data() + ml[i], line_end(ma, ml[i]), f, mneeds, mt);
    hit[i] = -3;
    if (n < mneeds) continue;
    Row r;
    r.freq = num(f[c[3]]); r.b = num(f[c[4]]); r.se = num(f[c[5]]); r.p = num(f[c[6]]); r.N = num(f[c[7]]);
    if (!(std::isfinite(r.N) && std::isfinite(r.b) && std::isfinite(r.se) && std::isfinite(r.freq) &&
          r.p >= 0 && r.p <= 1)) continue;
    ++n_valid;
    const int s = map.find(f[c[0]]);
    if (s < 0) { hit[i] = -2; continue; }
    const bool same = f[c[1]] == sa1[s] && f[c[2]] == sa2[s];
    const bool flip = f[c[1]] == sa2[s] && f[c[2]] == sa1[s];
    r.si = s; r.line = ml[i];
    r.allele_ok = same || flip;
    r.flip = flip;   // as SBayesRC: flipped when A1/A2 match the reference A2/A1
    const double fref = r.flip ? 1 - sfreq[s] : sfreq[s];
    r.freq_ok = r.allele_ok && std::fabs(fref - r.freq) <= freq_thresh;
    rows[i] = r;
    hit[i] = s;
  }
  // as SBayesRC (intersect + match): the first valid row of each SNP ID is kept, then checked
  std::vector<long> owner(M, -1);
  for (long i = 0; i < nma; ++i) if (hit[i] >= 0 && owner[hit[i]] < 0) owner[hit[i]] = i;
  long n_common = 0, n_allele = 0, n_dup = 0;
  for (long i = 0; i < nma; ++i) if (hit[i] >= 0 && owner[hit[i]] != i) ++n_dup;
  std::vector<long> ord;   // snp.info order
  ord.reserve(nma);
  for (int s = 0; s < M; ++s) {
    if (owner[s] < 0) continue;
    ++n_common;
    const Row& r = rows[owner[s]];
    if (!r.allele_ok) continue;
    ++n_allele;
    if (r.freq_ok) ord.push_back(owner[s]);
  }
  const long n_freq = static_cast<long>(ord.size());

  auto write_rows = [&](const std::string& path, const std::vector<long>& which) {
    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) Rcpp::stop("cannot write " + path);
    const char* hd = "SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n";
    std::fwrite(hd, 1, std::strlen(hd), fp);
    const size_t chunk = 1 << 18;
    for (size_t c0 = 0; c0 < which.size(); c0 += chunk * threads) {
      std::vector<std::string> buf(threads);
#ifdef _OPENMP
      #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
      for (int t = 0; t < threads; ++t) {
        const size_t a = std::min(which.size(), c0 + t * chunk), z = std::min(which.size(), a + chunk);
        std::string& out = buf[t];
        out.reserve((z - a) * 96);
        for (size_t q = a; q < z; ++q) {
          std::string_view f[kMaxF];
          const size_t ln = rows[which[q]].line;
          split(ma.data() + ln, line_end(ma, ln), f, mneeds, mt);
          for (int k = 0; k < 8; ++k) { out.append(f[c[k]].data(), f[c[k]].size()); out.push_back(k < 7 ? '\t' : '\n'); }
        }
      }
      for (int t = 0; t < threads; ++t) std::fwrite(buf[t].data(), 1, buf[t].size(), fp);
    }
    if (std::fclose(fp) != 0) Rcpp::stop("write failed: " + path);
  };

  // N within mean +- N_sd_range SD
  long double sumN = 0, ssN = 0;
  for (long i : ord) sumN += rows[i].N;
  const double mN = static_cast<double>(sumN / n_freq);
  for (long i : ord) ssN += (static_cast<long double>(rows[i].N) - mN) * (static_cast<long double>(rows[i].N) - mN);
  const double sN = std::sqrt(static_cast<double>(ssN / (n_freq - 1)));
  std::vector<long> ord2;
  ord2.reserve(n_freq);
  for (long i : ord) if (rows[i].N >= mN - N_sd_range * sN && rows[i].N <= mN + N_sd_range * sN) ord2.push_back(i);
  std::vector<double> vp(ord2.size());
  std::vector<double> Ns(ord2.size());
  for (size_t j = 0; j < ord2.size(); ++j) {
    const Row& r = rows[ord2[j]];
    vp[j] = 2 * r.freq * (1 - r.freq) * (r.N * r.se * r.se + r.b * r.b);
    Ns[j] = r.N;
  }
  const double medN = ord2.empty() ? NAN : median(Ns);
  const double vary = vp.empty() ? NAN : median(vp);
  std::vector<long> ord3;
  ord3.reserve(ord2.size());
  for (size_t j = 0; j < ord2.size(); ++j) {
    const double ind = std::sqrt(vp[j] / vary);
    if (ind > 1 - rate2pq && ind < 1 + rate2pq) ord3.push_back(ord2[j]);
  }
  if (!output.empty()) write_rows(output, ord3);

  Rcpp::NumericVector counts = Rcpp::NumericVector::create(
    Rcpp::_["n_ld"] = M, Rcpp::_["n_ma"] = static_cast<double>(nma), Rcpp::_["n_valid"] = n_valid,
    Rcpp::_["n_common"] = n_common, Rcpp::_["n_dup"] = n_dup, Rcpp::_["n_allele"] = n_allele,
    Rcpp::_["n_freq"] = n_freq, Rcpp::_["mean_N"] = mN, Rcpp::_["sd_N"] = sN,
    Rcpp::_["n_N"] = static_cast<double>(ord2.size()), Rcpp::_["median_N"] = medN, Rcpp::_["vary"] = vary,
    Rcpp::_["n_out"] = static_cast<double>(ord3.size()));
  Rcpp::List res = Rcpp::List::create(Rcpp::_["counts"] = counts);
  {
    const R_xlen_t n = ord3.size();
    Rcpp::IntegerVector idx(n);
    Rcpp::LogicalVector flp(n);
    Rcpp::NumericVector fq(n), b(n), se(n), p(n), N(n);
    for (R_xlen_t j = 0; j < n; ++j) {
      const Row& r = rows[ord3[j]];
      idx[j] = r.si + 1; flp[j] = r.flip;
      fq[j] = r.freq; b[j] = r.b; se[j] = r.se; p[j] = r.p; N[j] = r.N;
    }
    Rcpp::List tab = Rcpp::List::create(Rcpp::_["idx"] = idx, Rcpp::_["flip"] = flp, Rcpp::_["freq"] = fq,
                                        Rcpp::_["b"] = b, Rcpp::_["se"] = se, Rcpp::_["p"] = p, Rcpp::_["N"] = N);
    if (want_strings) {
      Rcpp::CharacterVector snp(n), a1(n), a2(n);
      for (R_xlen_t j = 0; j < n; ++j) {
        const Row& r = rows[ord3[j]];
        const std::string_view id = sid[r.si];
        snp[j] = Rf_mkCharLenCE(id.data(), id.size(), CE_NATIVE);
        const std::string_view x1 = r.flip ? sa2[r.si] : sa1[r.si], x2 = r.flip ? sa1[r.si] : sa2[r.si];
        a1[j] = Rf_mkCharLenCE(x1.data(), x1.size(), CE_NATIVE);
        a2[j] = Rf_mkCharLenCE(x2.data(), x2.size(), CE_NATIVE);
      }
      tab["SNP"] = snp; tab["A1"] = a1; tab["A2"] = a2;
    }
    res["table"] = tab;
  }
  return res;
}
