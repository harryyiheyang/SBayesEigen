// Summary-statistics imputation in the eigen space of SBayesRC block*.eigen.bin (Yihe Yang,
// impute.sbrc.fast 0.1.3). For each block the smallest exact solve is used: a missing-SNP system of
// size n_missing, a Woodbury eigen-space system of size k, or a typed-SNP system of size n_typed, with
// R_OO = U_O diag(lambda) U_O' + 0.1 I. Float32, blocks in parallel, Eigen single-threaded per block.
// With return_w the block's pass 1 is done while U is in memory: bhat = z / sqrt(N + z^2) on every SNP
// (typed: own N, imputed: N_median), w = Lambda^{-1/2} U' bhat.
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "eig_file.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using Eigen::LLT;
using Eigen::MatrixXf;
using Eigen::Ref;
using Eigen::VectorXf;

struct EigenBlock {
  int m;
  int k;
  Eigen::Map<const VectorXf> lambda;
  Eigen::Map<const MatrixXf> U;
  explicit EigenBlock(const EigFile& e) : m(e.m), k(e.k), lambda(e.lam, e.k), U(e.U, e.m, e.k) {}
};

static void validate_inputs(const EigenBlock& ld, const std::vector<int>& typed, const std::vector<float>& z) {
  if (typed.size() != z.size()) throw std::runtime_error("typed_index and z have different lengths");
  int previous = -1;
  for (std::size_t j = 0; j < typed.size(); ++j) {
    if (typed[j] < 0 || typed[j] >= ld.m) throw std::runtime_error("typed_index is outside the block");
    if (typed[j] <= previous) throw std::runtime_error("typed_index must be strictly increasing within each block");
    if (!std::isfinite(z[j])) throw std::runtime_error("z contains a non-finite value");
    previous = typed[j];
  }
}

static std::vector<int> missing_indices(int m, const std::vector<int>& typed) {
  std::vector<int> missing;
  missing.reserve(static_cast<std::size_t>(m) - typed.size());
  std::size_t j = 0;
  for (int i = 0; i < m; ++i) {
    if (j < typed.size() && typed[j] == i) ++j;
    else missing.push_back(i);
  }
  return missing;
}

static VectorXf solve_typed(const MatrixXf& UO, const VectorXf& lambda, const VectorXf& z, float diag_mod) {
  MatrixXf ROO = UO * lambda.asDiagonal() * UO.transpose();
  ROO.diagonal().array() += diag_mod;
  LLT<Ref<MatrixXf>, Eigen::Lower> llt(ROO);
  if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky factorization of R_OO failed");
  VectorXf solved = llt.solve(z);
  if (llt.info() != Eigen::Success || !solved.allFinite()) throw std::runtime_error("Solve of R_OO failed");
  return solved;
}

static VectorXf direct_eigen_coefficients(MatrixXf& UO, const VectorXf& lambda, const VectorXf& z, float diag_mod) {
  const int n_typed = UO.rows();
  const int k = UO.cols();
  if (k < n_typed) {
    const VectorXf sqrt_lambda = lambda.array().sqrt();
    for (int j = 0; j < k; ++j) UO.col(j) *= sqrt_lambda[j];
    MatrixXf system = (UO.transpose() * UO) / diag_mod;
    system.diagonal().array() += 1.0f;
    LLT<Ref<MatrixXf>, Eigen::Lower> llt(system);
    if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky factorization of the eigen-space system failed");
    const VectorXf alpha = llt.solve((UO.transpose() * z) / diag_mod);
    if (llt.info() != Eigen::Success || !alpha.allFinite()) throw std::runtime_error("Solve of the eigen-space system failed");
    return sqrt_lambda.array() * alpha.array();
  }
  const VectorXf solved = solve_typed(UO, lambda, z, diag_mod);
  return lambda.array() * (UO.transpose() * solved).array();
}

static std::vector<float> impute_missing_space(const EigenBlock& ld, const std::vector<int>& typed, const std::vector<float>& z, const std::vector<int>& missing, float diag_mod) {
  VectorXf z_full = VectorXf::Zero(ld.m);
  for (std::size_t i = 0; i < typed.size(); ++i) z_full[typed[i]] = z[i];
  const VectorXf dO = ld.U.transpose() * z_full;
  const int n_missing = static_cast<int>(missing.size());
  MatrixXf UM(n_missing, ld.k);
  for (int i = 0; i < n_missing; ++i) UM.row(i) = ld.U.row(missing[static_cast<std::size_t>(i)]);
  const VectorXf weight = ld.lambda.array() / (ld.lambda.array() + diag_mod);
  const VectorXf weighted_dO = weight.array() * dO.array();
  const VectorXf rhs = UM * weighted_dO;
  const VectorXf sqrt_weight = weight.array().sqrt();
  for (int j = 0; j < ld.k; ++j) UM.col(j) *= sqrt_weight[j];
  MatrixXf system = MatrixXf::Identity(n_missing, n_missing);
  system.noalias() -= UM * UM.transpose();
  UM.resize(0, 0);
  LLT<Ref<MatrixXf>, Eigen::Lower> llt(system);
  if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky factorization of the missing-SNP system failed");
  const VectorXf solution = llt.solve(rhs);
  if (llt.info() != Eigen::Success || !solution.allFinite()) throw std::runtime_error("Solve of the missing-SNP system failed");
  return std::vector<float>(solution.data(), solution.data() + solution.size());
}

static std::string solver_method(int n_typed, int n_missing, int k) {
  if (n_missing <= n_typed && n_missing <= k) return "missing";
  if (k < n_typed) return "eigen";
  return "typed";
}

static std::vector<float> impute_eigenspace(const EigenBlock& ld, const std::vector<int>& typed, const std::vector<float>& z, float diag_mod) {
  validate_inputs(ld, typed, z);
  const std::vector<int> missing = missing_indices(ld.m, typed);
  if (missing.empty()) return std::vector<float>();
  if (typed.empty()) throw std::runtime_error("Cannot impute a block with no typed SNPs");
  const int n_typed = static_cast<int>(typed.size());
  if (solver_method(n_typed, static_cast<int>(missing.size()), ld.k) == "missing") return impute_missing_space(ld, typed, z, missing, diag_mod);
  MatrixXf UO(n_typed, ld.k);
  VectorXf zf(n_typed);
  for (int i = 0; i < n_typed; ++i) {
    UO.row(i) = ld.U.row(typed[static_cast<std::size_t>(i)]);
    zf[i] = z[static_cast<std::size_t>(i)];
  }
  const VectorXf eig_coef = direct_eigen_coefficients(UO, ld.lambda, zf, diag_mod);
  std::vector<float> answer(missing.size());
  for (std::size_t i = 0; i < missing.size(); ++i) answer[i] = ld.U.row(missing[i]).dot(eig_coef);
  return answer;
}


struct BlockResult {
  std::vector<std::vector<float>> z_missing;   // per trait
  std::vector<std::vector<double>> w;          // per trait; empty when the trait has no typed SNP here
  std::vector<double> lam, ld;
  std::vector<int> method;                     // per trait: 0 missing, 1 eigen, 2 typed, -1 nothing to impute
  std::string error;
};

// K traits, each block mapped once. typed_index[[t]][[b]]: zero-based, strictly increasing typed SNPs of
// trait t in block b (empty: trait t skips the block); z, n_typed: per typed SNP; n_missing[t]: N used
// for imputed SNPs. return_w: pass 1 per trait while U is in memory, bhat = z / sqrt(N + z^2) on every
// SNP, w = Lambda^{-1/2} U' bhat. want_ld: eigen LD score sum_j U_ij^2 lambda_j^2 per SNP.
// [[Rcpp::export]]
Rcpp::List impute_blocks_eigen_cpp(Rcpp::CharacterVector files, Rcpp::List typed_index, Rcpp::List z,
                                   Rcpp::List n_typed, Rcpp::NumericVector n_missing, double thresh, int threads,
                                   bool return_w, bool want_ld) {
  const float diag_mod = 0.1f;
  const int nb = files.size(), K = typed_index.size();
  std::vector<std::string> fs(nb);
  for (int b = 0; b < nb; ++b) fs[b] = Rcpp::as<std::string>(files[b]);
  std::vector<std::vector<std::vector<int>>> typed(K, std::vector<std::vector<int>>(nb));
  std::vector<std::vector<std::vector<float>>> zt(K, std::vector<std::vector<float>>(nb));
  std::vector<std::vector<std::vector<double>>> nt(K, std::vector<std::vector<double>>(nb));
  for (int t = 0; t < K; ++t) {
    Rcpp::List ti = typed_index[t], zi = z[t];
    for (int b = 0; b < nb; ++b) {
      typed[t][b] = Rcpp::as<std::vector<int>>(ti[b]);
      zt[t][b] = Rcpp::as<std::vector<float>>(zi[b]);
    }
    if (return_w) {
      Rcpp::List ni = n_typed[t];
      for (int b = 0; b < nb; ++b) nt[t][b] = Rcpp::as<std::vector<double>>(ni[b]);
    }
  }
  std::vector<double> nmiss(n_missing.begin(), n_missing.end());
  std::vector<BlockResult> res(nb);
  Eigen::setNbThreads(1);
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
#endif
  for (int b = 0; b < nb; ++b) {
    BlockResult& r = res[b];
    r.z_missing.resize(K); r.w.resize(K); r.method.assign(K, -1);
    try {
      EigFile f;
      if (!f.open(fs[b], thresh)) throw std::runtime_error("cannot read " + fs[b]);
      const EigenBlock ld(f);
      std::vector<std::vector<double>> bh(K);
      for (int t = 0; t < K; ++t) {
        if (typed[t][b].empty()) continue;
        const std::vector<int> missing = missing_indices(ld.m, typed[t][b]);
        if (!missing.empty()) {
          const std::string meth = solver_method(static_cast<int>(typed[t][b].size()), static_cast<int>(missing.size()), ld.k);
          r.method[t] = meth == "missing" ? 0 : (meth == "eigen" ? 1 : 2);
        }
        r.z_missing[t] = impute_eigenspace(ld, typed[t][b], zt[t][b], diag_mod);
        if (!return_w) continue;
        bh[t].resize(ld.m);
        for (size_t i = 0; i < typed[t][b].size(); ++i) {
          const double zz = zt[t][b][i];
          bh[t][typed[t][b][i]] = zz / std::sqrt(nt[t][b][i] + zz * zz);
        }
        for (size_t i = 0; i < missing.size(); ++i) {
          const double zz = r.z_missing[t][i];
          bh[t][missing[i]] = zz / std::sqrt(nmiss[t] + zz * zz);
        }
        r.w[t].resize(ld.k);
      }
      if (return_w || want_ld) {
        r.lam.assign(f.lam, f.lam + ld.k);
        if (want_ld) r.ld.assign(ld.m, 0.0);
        for (int j = 0; j < ld.k; ++j) {   // one sweep over U for all traits
          const float* u = f.U + static_cast<size_t>(j) * ld.m;
          const double isl = 1.0 / std::sqrt(static_cast<double>(f.lam[j])), l2 = static_cast<double>(f.lam[j]) * f.lam[j];
          for (int t = 0; t < K; ++t) {
            if (bh[t].empty()) continue;
            double s = 0;
            for (int i = 0; i < ld.m; ++i) s += u[i] * bh[t][i];
            r.w[t][j] = s * isl;
          }
          if (want_ld) for (int i = 0; i < ld.m; ++i) r.ld[i] += static_cast<double>(u[i]) * u[i] * l2;
        }
      }
    } catch (const std::exception& e) {
      r.error = e.what();
    } catch (...) {
      r.error = "unknown C++ error";
    }
  }
  Rcpp::List out(nb);
  Rcpp::IntegerMatrix method(nb, K);
  for (int b = 0; b < nb; ++b) {
    if (!res[b].error.empty()) Rcpp::stop(fs[b] + ": " + res[b].error);
    Rcpp::List zl(K), wl(K);
    for (int t = 0; t < K; ++t) {
      method(b, t) = res[b].method[t];
      zl[t] = res[b].z_missing[t];
      wl[t] = res[b].w[t];
    }
    out[b] = Rcpp::List::create(Rcpp::_["z"] = zl, Rcpp::_["w"] = wl, Rcpp::_["lam"] = res[b].lam,
                                Rcpp::_["ld"] = res[b].ld);
  }
  out.attr("method") = method;
  return out;
}
