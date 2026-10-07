// blockN.ab.bin from LDbuild(A = ...), mapped read-only (layout at .read_ab in R/ld_build.R):
// int32 m, mA, kB, rankA; float tolA, cutB, sumLambdaB; int32 idxA[mA]; float R_AA upper triangle packed
// column-major; float lambdaB[kB]; float UlB[m*kB] column-major; float R_BA[(m-mA)*mA] column-major.
#ifndef SBAYESEIGEN_AB_FILE_H
#define SBAYESEIGEN_AB_FILE_H
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

struct AbFile {
  int m = 0, mA = 0, kB = 0, rankA = 0;
  float tolA = 0, cutB = 0, sumLambda = 0;
  const int32_t* idxA = nullptr;   // 0-based rows of A in the block
  const float* RAA = nullptr;      // packed upper: (r, c), r <= c, at c (c + 1) / 2 + r
  const float* lam = nullptr;
  const float* UlB = nullptr;      // column j at UlB + j * m
  const float* RBA = nullptr;      // (m - mA) x mA, column c at RBA + c * (m - mA)
  float raa(int r, int c) const { if (r > c) std::swap(r, c); return RAA[static_cast<size_t>(c) * (c + 1) / 2 + r]; }
  bool open(const std::string& f) {
#ifdef _WIN32
    std::FILE* fp = std::fopen(f.c_str(), "rb");
    if (!fp) return false;
    std::fseek(fp, 0, SEEK_END);
    len_ = static_cast<size_t>(std::ftell(fp));
    std::fseek(fp, 0, SEEK_SET);
    buf_.resize(len_);
    const bool ok = std::fread(buf_.data(), 1, len_, fp) == len_;
    std::fclose(fp);
    if (!ok) return false;
    base_ = buf_.data();
#else
    const int fd = ::open(f.c_str(), O_RDONLY);
    if (fd < 0) return false;
    struct stat st;
    if (fstat(fd, &st) != 0) { ::close(fd); return false; }
    len_ = static_cast<size_t>(st.st_size);
    void* p = len_ >= 28 ? mmap(nullptr, len_, PROT_READ, MAP_PRIVATE, fd, 0) : MAP_FAILED;
    ::close(fd);
    if (p == MAP_FAILED) return false;
    base_ = static_cast<const char*>(p);
    mapped_ = true;
    madvise(const_cast<char*>(base_), len_, MADV_SEQUENTIAL);
    madvise(const_cast<char*>(base_), len_, MADV_WILLNEED);
#endif
    if (len_ < 28) return false;
    int32_t h[4];
    float hf[3];
    std::memcpy(h, base_, 16);
    std::memcpy(hf, base_ + 16, 12);
    m = h[0]; mA = h[1]; kB = h[2]; rankA = h[3];
    tolA = hf[0]; cutB = hf[1]; sumLambda = hf[2];
    if (m <= 0 || mA < 0 || mA > m || kB < 0) return false;
    const size_t mB = static_cast<size_t>(m - mA);
    const size_t need = 28 + 4 * (static_cast<size_t>(mA) + static_cast<size_t>(mA) * (mA + 1) / 2 + kB +
                                  static_cast<size_t>(m) * kB + mB * mA);
    if (len_ != need) return false;
    const char* p = base_ + 28;
    idxA = reinterpret_cast<const int32_t*>(p); p += 4 * static_cast<size_t>(mA);
    RAA = reinterpret_cast<const float*>(p); p += 4 * (static_cast<size_t>(mA) * (mA + 1) / 2);
    lam = reinterpret_cast<const float*>(p); p += 4 * static_cast<size_t>(kB);
    UlB = reinterpret_cast<const float*>(p); p += 4 * static_cast<size_t>(m) * kB;
    RBA = reinterpret_cast<const float*>(p);
    return true;
  }
  ~AbFile() {
#ifndef _WIN32
    if (mapped_) munmap(const_cast<char*>(base_), len_);
#endif
  }

 private:
  const char* base_ = nullptr;
  size_t len_ = 0;
  bool mapped_ = false;
  std::vector<char> buf_;
};
#endif
