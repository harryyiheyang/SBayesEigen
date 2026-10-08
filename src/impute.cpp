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
#include "ab_file.h"
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
  std::vector<double> urow;                    // rows[[b]] of U (|rows| x k, column-major), ABC candidates
  std::vector<int> method;                     // per trait: 0 missing, 1 eigen, 2 typed, -1 nothing to impute
  std::string error;
};

// K traits, each block mapped once. typed_index[[t]][[b]]: zero-based, strictly increasing typed SNPs of
// trait t in block b (empty: trait t skips the block); z, n_typed: per typed SNP; n_missing[t]: N used
// for imputed SNPs. return_w: pass 1 per trait while U is in memory, bhat = z / sqrt(N + z^2) on every
// SNP, w = Lambda^{-1/2} U' bhat. want_ld: eigen LD score sum_j U_ij^2 lambda_j^2 per SNP. rows (optional, per
// block, zero-based): SNPs whose U rows on the kept components are returned as urow (|rows| x k), the ABC candidates.
// [[Rcpp::export]]
Rcpp::List impute_blocks_eigen_cpp(Rcpp::CharacterVector files, Rcpp::List typed_index, Rcpp::List z,
                                   Rcpp::List n_typed, Rcpp::NumericVector n_missing, double thresh, int threads,
                                   bool return_w, bool want_ld, Rcpp::List rows) {
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
  std::vector<std::vector<int>> rw(nb);
  if (rows.size() == nb) for (int b = 0; b < nb; ++b) rw[b] = Rcpp::as<std::vector<int>>(rows[b]);
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
      if (!rw[b].empty()) {
        const size_t s = rw[b].size();
        for (size_t i = 0; i < s; ++i) if (rw[b][i] < 0 || rw[b][i] >= ld.m) throw std::runtime_error("rows outside the block");
        r.urow.resize(s * ld.k);
        for (int j = 0; j < ld.k; ++j) {
          const float* u = f.U + static_cast<size_t>(j) * ld.m;
          for (size_t i = 0; i < s; ++i) r.urow[j * s + i] = u[rw[b][i]];
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
    const int s = static_cast<int>(rw[b].size());
    Rcpp::NumericMatrix ur(s, s ? static_cast<int>(res[b].urow.size() / s) : 0);
    std::copy(res[b].urow.begin(), res[b].urow.end(), ur.begin());
    out[b] = Rcpp::List::create(Rcpp::_["z"] = zl, Rcpp::_["w"] = wl, Rcpp::_["lam"] = res[b].lam,
                                Rcpp::_["ld"] = res[b].ld, Rcpp::_["urow"] = ur);
  }
  out.attr("method") = method;
  return out;
}

// ---- ab.bin (LDbuild with A): pass 1 for the joint ABC fit ----
// R ~ F F', F = [R_.A H, (0; Q2 Lambda2^{1/2})], H = V_A Lambda_A^{-1/2} (R_AA = V_A Lambda_A V_A', kept rankA),
// Q2 = B rows of UlB. Imputation uses F as a general low-rank factor (lambda = 1; Woodbury when rank < typed,
// else the typed system). Per trait: wA = H' bhat_A, wB = Lambda2^{-1/2} UlB' bhat. Per block (trait-free):
// lamA, lamB, XA = H' R_AA (rankA x mA), Cm = H' R_AB Q2 (rankA x kB), and for the candidate rows in B
// (rows, 0-based in the block; A rows dropped) XCa = H' R_A,i, XCb = Lambda2^{1/2} Q2[i, ] and their rows.
struct AbResult {
  std::vector<std::vector<float>> z_missing;
  std::vector<std::vector<double>> wA, wB;
  std::vector<double> lamA, lamB, XA, Cm, XCa, XCb;
  std::vector<int> idxA, crow;
  int rankA = 0, kB = 0, mA = 0;
  std::string error;
};

// [[Rcpp::export]]
Rcpp::List ab_pass1_cpp(Rcpp::CharacterVector files, Rcpp::List typed_index, Rcpp::List z, Rcpp::List n_typed,
                        Rcpp::NumericVector n_missing, Rcpp::List rows, int threads) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  const float diag_mod = 0.1f;
  const int nb = files.size(), K = typed_index.size();
  std::vector<std::string> fs(nb);
  for (int b = 0; b < nb; ++b) fs[b] = Rcpp::as<std::string>(files[b]);
  std::vector<std::vector<std::vector<int>>> typed(K, std::vector<std::vector<int>>(nb));
  std::vector<std::vector<std::vector<float>>> zt(K, std::vector<std::vector<float>>(nb));
  std::vector<std::vector<std::vector<double>>> nt(K, std::vector<std::vector<double>>(nb));
  for (int t = 0; t < K; ++t) {
    Rcpp::List ti = typed_index[t], zi = z[t], ni = n_typed[t];
    for (int b = 0; b < nb; ++b) {
      typed[t][b] = Rcpp::as<std::vector<int>>(ti[b]);
      zt[t][b] = Rcpp::as<std::vector<float>>(zi[b]);
      nt[t][b] = Rcpp::as<std::vector<double>>(ni[b]);
    }
  }
  std::vector<std::vector<int>> rw(nb);
  for (int b = 0; b < nb; ++b) rw[b] = Rcpp::as<std::vector<int>>(rows[b]);
  std::vector<double> nmiss(n_missing.begin(), n_missing.end());
  std::vector<AbResult> res(nb);
  Eigen::setNbThreads(1);
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
#endif
  for (int b = 0; b < nb; ++b) {
    AbResult& r = res[b];
    r.z_missing.resize(K); r.wA.resize(K); r.wB.resize(K);
    try {
      AbFile f;
      if (!f.open(fs[b])) throw std::runtime_error("cannot read " + fs[b]);
      const int m = f.m, mA = f.mA, mB = m - mA, kB = f.kB;
      r.mA = mA; r.kB = kB;
      std::vector<int> isA(m, -1), posB(m, -1);
      for (int a = 0; a < mA; ++a) { isA[f.idxA[a]] = a; r.idxA.push_back(f.idxA[a]); }
      for (int i = 0, q = 0; i < m; ++i) if (isA[i] < 0) posB[i] = q++;
      // A: eigen of R_AA, H, R_AA H = V Lambda^{1/2}, G = R_BA H
      int rk = 0;
      MatrixXd H, RAH, G(mB, 0);
      if (mA > 0) {
        MatrixXd RAA(mA, mA);
        for (int c = 0; c < mA; ++c) for (int q = 0; q < mA; ++q) RAA(q, c) = f.raa(q, c);
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(RAA);
        const VectorXd ev = es.eigenvalues();   // ascending
        const double mx = ev[mA - 1];
        // at most the top rankA components (LDbuild's cut), and only eigenvalues > 1e-6 max: R_AA is stored as float,
        // smaller ones are rounding noise whose Lambda^{-1/2} would blow up w_A (Yihe 2026-10-08)
        std::vector<int> keep;
        for (int j = mA - 1; j >= 0 && static_cast<int>(keep.size()) < f.rankA; --j) if (ev[j] > 1e-6 * mx) keep.push_back(j);
        rk = static_cast<int>(keep.size());
        H.resize(mA, rk); RAH.resize(mA, rk);
        for (int j = 0; j < rk; ++j) {
          const double l = ev[keep[j]];
          r.lamA.push_back(l);
          H.col(j) = es.eigenvectors().col(keep[j]) / std::sqrt(l);
          RAH.col(j) = es.eigenvectors().col(keep[j]) * std::sqrt(l);
        }
        Eigen::Map<const Eigen::MatrixXf> RBA(f.RBA, mB, mA);
        G = RBA.cast<double>() * H;
      }
      r.rankA = rk;
      // Q2 (B rows of UlB), lamB
      MatrixXd Q2(mB, kB);
      for (int j = 0; j < kB; ++j) {
        const float* u = f.UlB + static_cast<size_t>(j) * m;
        for (int i = 0; i < m; ++i) if (posB[i] >= 0) Q2(posB[i], j) = u[i];
        r.lamB.push_back(f.lam[j]);
      }
      // per trait: impute, bhat, wA, wB
      const int nr = rk + kB;
      for (int t = 0; t < K; ++t) {
        const std::vector<int>& ty = typed[t][b];
        if (ty.empty()) continue;
        std::vector<float> zz(zt[t][b]);
        for (size_t i = 0; i < ty.size(); ++i) if ((i > 0 && ty[i] <= ty[i - 1]) || ty[i] < 0 || ty[i] >= m) throw std::runtime_error("typed_index must be increasing within the block");
        const std::vector<int> missing = missing_indices(m, ty);
        std::vector<double> bh(m, 0.0);
        auto frow = [&](int i, float* out) {   // row i of F
          if (isA[i] >= 0) { for (int j = 0; j < rk; ++j) out[j] = static_cast<float>(RAH(isA[i], j)); for (int j = 0; j < kB; ++j) out[rk + j] = 0; }
          else { const int q = posB[i]; for (int j = 0; j < rk; ++j) out[j] = static_cast<float>(G(q, j));
                 for (int j = 0; j < kB; ++j) out[rk + j] = static_cast<float>(Q2(q, j) * std::sqrt(static_cast<double>(f.lam[j]))); }
        };
        if (!missing.empty()) {
          const int no = static_cast<int>(ty.size());
          Eigen::MatrixXf FO(no, nr); Eigen::VectorXf zf(no);
          std::vector<float> row(nr);
          for (int i = 0; i < no; ++i) { frow(ty[i], row.data()); FO.row(i) = Eigen::Map<Eigen::RowVectorXf>(row.data(), nr); zf[i] = zz[i]; }
          const Eigen::VectorXf coef = direct_eigen_coefficients(FO, Eigen::VectorXf::Ones(nr), zf, diag_mod);
          r.z_missing[t].resize(missing.size());
          for (size_t i = 0; i < missing.size(); ++i) {
            frow(missing[i], row.data());
            r.z_missing[t][i] = Eigen::Map<Eigen::RowVectorXf>(row.data(), nr).dot(coef);
          }
          for (size_t i = 0; i < missing.size(); ++i) { const double v = r.z_missing[t][i]; bh[missing[i]] = v / std::sqrt(nmiss[t] + v * v); }
        }
        for (size_t i = 0; i < ty.size(); ++i) { const double v = zz[i]; bh[ty[i]] = v / std::sqrt(nt[t][b][i] + v * v); }
        r.wA[t].assign(rk, 0.0);
        for (int j = 0; j < rk; ++j) { double s = 0; for (int a = 0; a < mA; ++a) s += H(a, j) * bh[r.idxA[a]]; r.wA[t][j] = s; }
        r.wB[t].assign(kB, 0.0);
        for (int j = 0; j < kB; ++j) {
          const float* u = f.UlB + static_cast<size_t>(j) * m; double s = 0;
          for (int i = 0; i < m; ++i) s += u[i] * bh[i];
          r.wB[t][j] = s / std::sqrt(static_cast<double>(f.lam[j]));
        }
      }
      // design pieces
      r.XA.assign(static_cast<size_t>(rk) * mA, 0.0);
      for (int a = 0; a < mA; ++a) for (int j = 0; j < rk; ++j) r.XA[static_cast<size_t>(a) * rk + j] = RAH(a, j);   // H'R_AA = (R_AA H)'
      r.Cm.assign(static_cast<size_t>(rk) * kB, 0.0);
      if (rk > 0) { const MatrixXd Cm = G.transpose() * Q2; std::copy(Cm.data(), Cm.data() + Cm.size(), r.Cm.begin()); }
      for (int i : rw[b]) if (i >= 0 && i < m && isA[i] < 0) r.crow.push_back(i);
      const size_t s = r.crow.size();
      r.XCa.assign(static_cast<size_t>(rk) * s, 0.0); r.XCb.assign(static_cast<size_t>(kB) * s, 0.0);
      for (size_t c = 0; c < s; ++c) {
        const int q = posB[r.crow[c]];
        for (int j = 0; j < rk; ++j) r.XCa[c * rk + j] = G(q, j);
        for (int j = 0; j < kB; ++j) r.XCb[c * kB + j] = Q2(q, j) * std::sqrt(static_cast<double>(f.lam[j]));
      }
    } catch (const std::exception& e) {
      r.error = e.what();
    } catch (...) {
      r.error = "unknown C++ error";
    }
  }
  Rcpp::List out(nb);
  for (int b = 0; b < nb; ++b) {
    const AbResult& r = res[b];
    if (!r.error.empty()) Rcpp::stop(fs[b] + ": " + r.error);
    Rcpp::List zl(K), wl(K);
    for (int t = 0; t < K; ++t) {
      zl[t] = r.z_missing[t];
      std::vector<double> w(r.wA[t]); w.insert(w.end(), r.wB[t].begin(), r.wB[t].end());
      wl[t] = r.wA[t].empty() && r.wB[t].empty() ? std::vector<double>() : w;
    }
    const int rk = r.rankA, kB = r.kB, s = static_cast<int>(r.crow.size());
    Rcpp::NumericMatrix XA(rk, r.mA), Cm(rk, kB), XCa(rk, s), XCb(kB, s);
    std::copy(r.XA.begin(), r.XA.end(), XA.begin()); std::copy(r.Cm.begin(), r.Cm.end(), Cm.begin());
    std::copy(r.XCa.begin(), r.XCa.end(), XCa.begin()); std::copy(r.XCb.begin(), r.XCb.end(), XCb.begin());
    std::vector<double> lam(r.lamA); lam.insert(lam.end(), r.lamB.begin(), r.lamB.end());
    std::vector<int> ia(r.idxA), cr(r.crow);
    for (int& v : ia) ++v;
    for (int& v : cr) ++v;
    out[b] = Rcpp::List::create(Rcpp::_["z"] = zl, Rcpp::_["w"] = wl, Rcpp::_["lam"] = lam, Rcpp::_["rankA"] = rk,
                                Rcpp::_["iA"] = ia, Rcpp::_["XA"] = XA, Rcpp::_["Cm"] = Cm,
                                Rcpp::_["crow"] = cr, Rcpp::_["XCa"] = XCa, Rcpp::_["XCb"] = XCb);
  }
  return out;
}
