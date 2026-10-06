// Genotypes are held as two bit planes per variant, lo = [x == 1] and
// hi = [x == 2], so that
//   sum_k x_k y_k = pc(lo & lo') + 2 pc((lo & hi') | (hi & lo')) + 4 pc(hi & hi').
// All sums are exact integers; r = (n G - s_a s_b) / sqrt(v_a v_b), with
// v = n ss - s^2. Missing calls are filled with the variant's lower median.
#ifndef CPPMATRIX_GENO_KERNEL_H
#define CPPMATRIX_GENO_KERNEL_H

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define GENO_KERNEL_X86 1
#endif

namespace {

using u64 = std::uint64_t;
constexpr std::size_t kTile = 64;          // variants per tile side
constexpr std::size_t kChunk = 512;        // 64-bit words per cache chunk

struct Planes {
  std::size_t m = 0, words = 0;       // words is a multiple of 4
  std::vector<u64> bits;              // variant j: lo at 2jW, hi at (2j+1)W
  std::vector<std::int64_t> s, v;     // sum x and n*sum x^2 - (sum x)^2
  const u64* lo(std::size_t j) const { return bits.data() + 2 * j * words; }
  const u64* hi(std::size_t j) const { return bits.data() + (2 * j + 1) * words; }
};

inline std::size_t plane_words(std::size_t n) { return ((n + 63) / 64 + 3) / 4 * 4; }

inline int popcount64(u64 x) { return __builtin_popcountll(x); }

inline void fill_median(std::size_t n, std::size_t words, u64* lo, u64* hi,
                        const u64* miss, std::int64_t& s, std::int64_t& v) {
  std::int64_t c1 = 0, c2 = 0, cm = 0;
  for (std::size_t w = 0; w < words; ++w) {
    c1 += popcount64(lo[w]);
    c2 += popcount64(hi[w]);
    cm += popcount64(miss[w]);
  }
  const std::int64_t observed = static_cast<std::int64_t>(n) - cm;
  const std::int64_t c0 = observed - c1 - c2;
  if (cm && observed) {
    u64* fill = 2 * c0 >= observed ? nullptr : (2 * (c0 + c1) >= observed ? lo : hi);
    if (fill) {
      for (std::size_t w = 0; w < words; ++w) fill[w] |= miss[w];
      (fill == lo ? c1 : c2) += cm;
    }
  }
  s = c1 + 2 * c2;
  v = static_cast<std::int64_t>(n) * (c1 + 4 * c2) - s * s;
}

// Packed 2-bit genotypes, 4 samples per byte with the first sample in the low
// bits (PLINK 1 BED columns and pgenlib genovecs). dosage[code] is the
// genotype x in {0, 1, 2} for that code, or -1 for missing.
struct NibbleTable {
  std::array<unsigned char, 256> lo{}, hi{}, miss{};
  explicit NibbleTable(const std::array<int, 4>& dosage) {
    for (unsigned int b = 0; b < 256; ++b)
      for (unsigned int k = 0; k < 4; ++k) {
        const int x = dosage[(b >> (2 * k)) & 3u];
        if (x == 2) hi[b] |= 1u << k;
        else if (x == 1) lo[b] |= 1u << k;
        else if (x < 0) miss[b] |= 1u << k;
      }
  }
};

inline void encode_packed(const NibbleTable& table, const unsigned char* column,
                          std::size_t n, std::size_t stride, std::size_t words,
                          u64* lo, u64* hi, std::int64_t& s, std::int64_t& v) {
  std::vector<u64> miss(words, 0);
  for (std::size_t w = 0; w < words; ++w) {
    u64 l = 0, h = 0, m = 0;
    const std::size_t b0 = 16 * w, b1 = std::min(stride, b0 + 16);
    for (std::size_t b = b0; b < b1; ++b) {
      const unsigned char byte = column[b];
      const unsigned int shift = 4 * (b - b0);
      l |= static_cast<u64>(table.lo[byte]) << shift;
      h |= static_cast<u64>(table.hi[byte]) << shift;
      m |= static_cast<u64>(table.miss[byte]) << shift;
    }
    lo[w] = l;
    hi[w] = h;
    miss[w] = m;
  }
  const std::size_t full = n / 64;
  if (n % 64) {
    const u64 mask = (u64(1) << (n % 64)) - 1;
    lo[full] &= mask;
    hi[full] &= mask;
    miss[full] &= mask;
  }
  fill_median(n, words, lo, hi, miss.data(), s, v);
}

inline void reset_planes(Planes& out, std::size_t n, std::size_t m) {
  out.m = m;
  out.words = plane_words(n);
  out.bits.assign(2 * m * out.words, 0);
  out.s.assign(m, 0);
  out.v.assign(m, 0);
}

// Cross products of one tile: acc[a * kTile + b] = sum x_{i0+a} y_{j0+b}.
using TileKernel = void (*)(const Planes&, std::size_t, std::size_t,
                            const Planes&, std::size_t, std::size_t,
                            std::int64_t*);

void tile_generic(const Planes& A, std::size_t i0, std::size_t ni,
                  const Planes& B, std::size_t j0, std::size_t nj,
                  std::int64_t* acc) {
  const std::size_t W = A.words;
  for (std::size_t w0 = 0; w0 < W; w0 += kChunk) {
    const std::size_t w1 = std::min(W, w0 + kChunk);
    for (std::size_t a0 = 0; a0 < ni; a0 += 4)
      for (std::size_t b0 = 0; b0 < nj; b0 += 4) {
        const u64 *al[4], *ah[4], *bl[4], *bh[4];
        for (std::size_t t = 0; t < 4; ++t) {
          const std::size_t i = i0 + std::min(a0 + t, ni - 1);
          const std::size_t j = j0 + std::min(b0 + t, nj - 1);
          al[t] = A.lo(i); ah[t] = A.hi(i); bl[t] = B.lo(j); bh[t] = B.hi(j);
        }
        std::int64_t c[4][4] = {};
        for (std::size_t w = w0; w < w1; ++w) {
          u64 xl[4], xh[4], yl[4], yh[4];
          for (int t = 0; t < 4; ++t) {
            xl[t] = al[t][w]; xh[t] = ah[t][w]; yl[t] = bl[t][w]; yh[t] = bh[t][w];
          }
          for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
              c[a][b] += popcount64(xl[a] & yl[b]) +
                2 * popcount64((xl[a] & yh[b]) | (xh[a] & yl[b])) +
                4 * popcount64(xh[a] & yh[b]);
        }
        for (std::size_t a = 0; a < 4 && a0 + a < ni; ++a)
          for (std::size_t b = 0; b < 4 && b0 + b < nj; ++b)
            acc[(a0 + a) * kTile + b0 + b] += c[a][b];
      }
  }
}

#ifdef GENO_KERNEL_X86
// Nibble-lookup popcount on 256-bit vectors. Per byte, one vector adds at
// most 8 + 2*8 + 4*8 = 56, so byte counters are flushed every 4 vectors.
// GCC lambdas do not inherit target("avx2"), hence always_inline helpers.
#define GENO_KERNEL_AVX2 __attribute__((target("avx2"), always_inline)) inline

GENO_KERNEL_AVX2 __m256i count_avx2(__m256i x) {
  const __m256i lut = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                       0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
  const __m256i low = _mm256_set1_epi8(0x0f);
  return _mm256_add_epi8(
    _mm256_shuffle_epi8(lut, _mm256_and_si256(x, low)),
    _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(x, 4), low)));
}

GENO_KERNEL_AVX2 __m256i dot_avx2(__m256i xl, __m256i xh, __m256i yl, __m256i yh) {
  __m256i c = count_avx2(_mm256_and_si256(xl, yl));
  c = _mm256_add_epi8(c, _mm256_slli_epi16(count_avx2(_mm256_or_si256(
    _mm256_and_si256(xl, yh), _mm256_and_si256(xh, yl))), 1));
  return _mm256_add_epi8(c, _mm256_slli_epi16(count_avx2(_mm256_and_si256(xh, yh)), 2));
}

__attribute__((target("avx2")))
void tile_avx2(const Planes& A, std::size_t i0, std::size_t ni,
               const Planes& B, std::size_t j0, std::size_t nj,
               std::int64_t* acc) {
  const std::size_t V = A.words / 4;
  const __m256i zero = _mm256_setzero_si256();
  const std::size_t chunk = kChunk / 4;
  for (std::size_t v0 = 0; v0 < V; v0 += chunk) {
    const std::size_t v1 = std::min(V, v0 + chunk);
    for (std::size_t a0 = 0; a0 < ni; a0 += 2)
      for (std::size_t b0 = 0; b0 < nj; b0 += 4) {
        const __m256i *al[2], *ah[2], *bl[4], *bh[4];
        #pragma GCC unroll 4
        for (std::size_t t = 0; t < 2; ++t) {
          const std::size_t i = i0 + std::min(a0 + t, ni - 1);
          al[t] = reinterpret_cast<const __m256i*>(A.lo(i));
          ah[t] = reinterpret_cast<const __m256i*>(A.hi(i));
        }
        #pragma GCC unroll 4
        for (std::size_t t = 0; t < 4; ++t) {
          const std::size_t j = j0 + std::min(b0 + t, nj - 1);
          bl[t] = reinterpret_cast<const __m256i*>(B.lo(j));
          bh[t] = reinterpret_cast<const __m256i*>(B.hi(j));
        }
        __m256i total[2][4], bytes[2][4];
        #pragma GCC unroll 4
        for (int a = 0; a < 2; ++a)
          #pragma GCC unroll 4
          for (int b = 0; b < 4; ++b) total[a][b] = zero;
        for (std::size_t v = v0; v < v1;) {
          const std::size_t end = std::min(v1, v + 4);
          #pragma GCC unroll 4
          for (int a = 0; a < 2; ++a)
            #pragma GCC unroll 4
            for (int b = 0; b < 4; ++b) bytes[a][b] = zero;
          for (; v < end; ++v) {
            const __m256i xl0 = _mm256_loadu_si256(al[0] + v);
            const __m256i xh0 = _mm256_loadu_si256(ah[0] + v);
            const __m256i xl1 = _mm256_loadu_si256(al[1] + v);
            const __m256i xh1 = _mm256_loadu_si256(ah[1] + v);
            #pragma GCC unroll 4
            for (int b = 0; b < 4; ++b) {
              const __m256i yl = _mm256_loadu_si256(bl[b] + v);
              const __m256i yh = _mm256_loadu_si256(bh[b] + v);
              bytes[0][b] = _mm256_add_epi8(bytes[0][b], dot_avx2(xl0, xh0, yl, yh));
              bytes[1][b] = _mm256_add_epi8(bytes[1][b], dot_avx2(xl1, xh1, yl, yh));
            }
          }
          #pragma GCC unroll 4
          for (int a = 0; a < 2; ++a)
            #pragma GCC unroll 4
            for (int b = 0; b < 4; ++b)
              total[a][b] = _mm256_add_epi64(total[a][b], _mm256_sad_epu8(bytes[a][b], zero));
        }
        #pragma GCC unroll 4
        for (std::size_t a = 0; a < 2 && a0 + a < ni; ++a)
          #pragma GCC unroll 4
          for (std::size_t b = 0; b < 4 && b0 + b < nj; ++b) {
            alignas(32) std::int64_t lanes[4];
            _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), total[a][b]);
            acc[(a0 + a) * kTile + b0 + b] += lanes[0] + lanes[1] + lanes[2] + lanes[3];
          }
      }
  }
}
#endif

TileKernel select_kernel() {
#ifdef GENO_KERNEL_X86
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx2")) return tile_avx2;
#endif
  return tile_generic;
}

void cor_block(const Planes& L, std::size_t ib, std::size_t pb,
               const Planes& R, std::size_t jb, std::size_t qb, bool self,
               double n, double* out, std::size_t ld, int threads,
               TileKernel kernel) {
  const std::size_t ti = (pb + kTile - 1) / kTile, tj = (qb + kTile - 1) / kTile;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<std::int64_t> acc(kTile * kTile);
    #pragma omp for schedule(dynamic, 1)
    for (std::ptrdiff_t t = 0; t < static_cast<std::ptrdiff_t>(ti * tj); ++t) {
      const std::size_t i0 = (t / tj) * kTile, j0 = (t % tj) * kTile;
      if (self && ib + i0 > jb + j0 + kTile - 1) continue;
      const std::size_t ni = std::min(kTile, pb - i0), nj = std::min(kTile, qb - j0);
      std::fill(acc.begin(), acc.end(), 0);
      kernel(L, i0, ni, R, j0, nj, acc.data());
      for (std::size_t x = 0; x < ni; ++x)
        for (std::size_t y = 0; y < nj; ++y) {
          const std::size_t i = i0 + x, j = j0 + y;
          const double va = static_cast<double>(L.v[i]);
          const double vb = static_cast<double>(R.v[j]);
          double r = NA_REAL;
          if (va > 0 && vb > 0) {
            r = (n * static_cast<double>(acc[x * kTile + y]) -
                 static_cast<double>(L.s[i]) * static_cast<double>(R.s[j])) /
              std::sqrt(va * vb);
            r = std::min(1.0, std::max(-1.0, r));
          }
          out[(ib + i) + (jb + j) * ld] = r;
          if (self) out[(jb + j) + (ib + i) * ld] = r;
        }
    }
  }
}

constexpr double kPlaneBudget = 268435456; // bytes of bit planes per side

template <class LoadA, class LoadB>
Rcpp::NumericMatrix cor_driver(std::size_t p, std::size_t q, std::size_t n_samples,
                               bool self, int threads, LoadA load_a, LoadB load_b) {
  const TileKernel kernel = select_kernel();
  const double n = static_cast<double>(n_samples);
  const double plane_bytes = 2.0 * 8.0 * plane_words(n_samples);
  const std::size_t block = std::max<std::size_t>(
    kTile, static_cast<std::size_t>(kPlaneBudget / plane_bytes) / kTile * kTile);

  Rcpp::NumericMatrix result(p, q);
  double* out = result.begin();
  Planes PA, PB;
  for (std::size_t jb = 0; jb < q; jb += block) {
    const std::size_t qb = std::min(block, q - jb);
    load_b(jb, qb, PB);
    for (std::size_t ib = 0; ib < (self ? jb + qb : p); ib += block) {
      const std::size_t pb = std::min(block, p - ib);
      const Planes* left = &PA;
      if (self && ib == jb) left = &PB;
      else if (self) load_b(ib, pb, PA);
      else load_a(ib, pb, PA);
      cor_block(*left, ib, pb, PB, jb, qb, self, n, out, p, threads, kernel);
      Rcpp::checkUserInterrupt();
    }
  }
  return result;
}

} // namespace

#endif
