// SBayesRC blockN.eigen.bin: int32 m, int32 k, float sumLambda, float thresh, float lambda[k], float U[m*k]
// (column-major). Mapped read-only (mmap on POSIX, one read on Windows): pass 2 re-maps the same file and
// the page cache serves it when RAM allows.
#ifndef SBAYESEIGEN_EIG_FILE_H
#define SBAYESEIGEN_EIG_FILE_H
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

// number of leading components kept: the first k with cumsum(lambda) >= thresh * sumLambda (as
// SBayesRC read1LD); all stored components when thresh is <= 0 or not below the file's own threshold
inline int eig_cut(const float* lam, int k, float sumL, float file_thresh, double thresh) {
  if (thresh <= 0 || thresh >= static_cast<double>(file_thresh) - 1e-6) return k;
  double cs = 0;
  for (int j = 0; j < k; ++j) {
    cs += lam[j];
    if (cs >= thresh * sumL) return j + 1;
  }
  return k;
}

struct EigFile {
  int m = 0, k = 0;              // k after the cut
  const float* lam = nullptr;
  const float* U = nullptr;      // column j at U + j * m
  bool open(const std::string& f, double thresh) {
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
    void* p = len_ >= 16 ? mmap(nullptr, len_, PROT_READ, MAP_PRIVATE, fd, 0) : MAP_FAILED;
    ::close(fd);
    if (p == MAP_FAILED) return false;
    base_ = static_cast<const char*>(p);
    mapped_ = true;
#endif
    if (len_ < 16) return false;
    int32_t h[2];
    float hf[2];
    std::memcpy(h, base_, 8);
    std::memcpy(hf, base_ + 8, 8);
    m = h[0];
    const int kf = h[1];
    if (m <= 0 || kf <= 0 || len_ != 16 + 4 * static_cast<size_t>(kf) * (1 + static_cast<size_t>(m))) return false;
    lam = reinterpret_cast<const float*>(base_ + 16);
    U = lam + kf;
    k = eig_cut(lam, kf, hf[0], hf[1], thresh);
#ifndef _WIN32
    // ask the kernel to read the used prefix (header, lambda, first k columns) ahead in large requests;
    // matters on network file systems, where 4 KB page faults leave the CPUs waiting
    const size_t need = 16 + 4 * static_cast<size_t>(kf) + 4 * static_cast<size_t>(m) * k;
    madvise(const_cast<char*>(base_), len_, MADV_SEQUENTIAL);
    madvise(const_cast<char*>(base_), need, MADV_WILLNEED);
#endif
    return true;
  }
  ~EigFile() {
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
