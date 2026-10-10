// ABC model (Yihe 2026-10-07): beta = beta_A + U alpha + gamma, joint on one eigen.bin.
// In component space (blocks pooled): w = sum_{i in A} x_i beta_i + c alpha + sum_{i in C} x_i gamma_i + e,
// x_i = sqrt(lambda) * U[i, ] (k_b, zero outside its block), e_j ~ N(0, 1 / p_j), p_j = n_j / (ve0 + kappa / lambda_j).
// One call does one coordinate sweep over the A candidates (mr.ash: exact zero + gammaA * s2A) and/or the
// C candidates (MCP on the standardised scale h = gamma sqrt(D), D = sum_j p_j x_j^2, so the threshold is in
// units of the residual z noise sd: tau sqrt(ve0) under constant noise) on the residual r, updated in place.
#include <Rcpp.h>
#include <string>
#include "ab_file.h"
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static inline double mcp_prox(double r, double tau, double a) {
  const double ar = std::fabs(r);
  if (ar <= tau) return 0;
  if (ar <= a * tau) return (r > 0 ? 1 : -1) * (ar - tau) / (1 - 1 / a);
  return r;
}

// X[[b]]: k_b x s_b candidate columns of block b; off[b]: first component of block b in r and p;
// role[[b]]: 1 = A, 2 = C, per column; coef[[b]]: beta_i (A) or gamma_i (C) per column (updated in place).
// [[Rcpp::export]]
List abc_sweep_cpp(List X, IntegerVector off, NumericVector p, NumericVector r, List role, List coef,
                   NumericVector gammaA, NumericVector piA, double s2A, double tau, double a,
                   bool doA, bool doC, int threads) {
  const int nb = X.size(), K = gammaA.size();
  std::vector<NumericMatrix> Xb(nb); std::vector<IntegerVector> rl(nb); std::vector<NumericVector> cf(nb);
  for (int b = 0; b < nb; b++) { Xb[b] = as<NumericMatrix>(X[b]); rl[b] = as<IntegerVector>(role[b]); cf[b] = as<NumericVector>(coef[b]); }
  double* pr = r.begin(); const double* pp = p.begin();
  std::vector<double> lpi(K);
  for (int k = 0; k < K; k++) lpi[k] = std::log(piA[k]);
  std::vector<double> sphi(K, 0.0); double sa = 0, snz = 0, mx = 0;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> lw(K), mu(K), v(K), sphil(K, 0.0); double sal = 0, snzl = 0, mxl = 0;
    #pragma omp for schedule(dynamic)
    for (int b = 0; b < nb; b++) {
      const NumericMatrix& x = Xb[b]; const int kb = x.nrow(), s = x.ncol();
      double* rb = pr + off[b]; const double* pb = pp + off[b];
      for (int i = 0; i < s; i++) {
        const int ro = rl[b][i];
        if (!((ro == 1 && doA) || (ro == 2 && doC))) continue;
        const double* xi = &x(0, i);
        double D = 0, rho = 0;
        for (int j = 0; j < kb; j++) { const double px = pb[j] * xi[j]; D += px * xi[j]; rho += px * rb[j]; }
        if (D <= 0) continue;
        const double old = cf[b][i];
        rho += D * old;                                  // x_i' P (r + x_i coef_i)
        double nw;
        if (ro == 1) {
          double m = -INFINITY;
          for (int k = 0; k < K; k++) {
            if (gammaA[k] == 0) { lw[k] = lpi[k]; mu[k] = 0; v[k] = 0; }
            else {
              const double g = gammaA[k] * s2A;
              v[k] = 1 / (D + 1 / g); mu[k] = v[k] * rho;
              lw[k] = lpi[k] + 0.5 * std::log(v[k] / g) + 0.5 * mu[k] * mu[k] / v[k];
            }
            if (lw[k] > m) m = lw[k];
          }
          double sw = 0; for (int k = 0; k < K; k++) { lw[k] = std::exp(lw[k] - m); sw += lw[k]; }
          nw = 0;
          for (int k = 0; k < K; k++) {
            const double phi = lw[k] / sw; sphil[k] += phi; nw += phi * mu[k];
            if (gammaA[k] > 0) { sal += phi * (mu[k] * mu[k] + v[k]) / gammaA[k]; snzl += phi; }
          }
        } else {
          const double sd = std::sqrt(D);
          nw = mcp_prox(rho / sd, tau, a) / sd;
          mxl = std::max(mxl, std::fabs(nw - old) * sd);
        }
        const double dm = nw - old;
        if (dm != 0) { for (int j = 0; j < kb; j++) rb[j] -= xi[j] * dm; cf[b][i] = nw; }
      }
    }
    #pragma omp critical
    {
      for (int k = 0; k < K; k++) sphi[k] += sphil[k];
      sa += sal; snz += snzl; mx = std::max(mx, mxl);
    }
  }
  return List::create(_["sphi"] = NumericVector(sphi.begin(), sphi.end()), _["sa"] = sa, _["snz"] = snz, _["dz"] = mx);
}

// ---- ABC on ab.bin (LD-build thread AB_JOINT_FIT.md): per block, with rk = rankA,
//   w_A = XA beta_A + Cm alpha + XCa gamma + e_A,   w_B = sqrt(lamB) alpha + XCb gamma + e_B,
// e ~ N(0, 1/p) per component (p = n / (ve0 + kappa / lambda)). One call: one coordinate sweep of A (mr.ash),
// B (mr.ash per component; design column (Cm[, j]; sqrt(lamB_j) e_j)) and C (MCP on the standardised scale)
// on the residuals rA = r[off + 0..rk), rB = r[off + rk..rk+kB), updated in place. Blocks in parallel.
// blk[[b]]: list(XA, Cm, sl, XCa, XCb); mA, al, gC: per block coefficients (updated in place).
// [[Rcpp::export]]
List abj_sweep_cpp(List blk, IntegerVector off, NumericVector p, NumericVector r, List mA, List al, List gC,
                   NumericVector gammaA, NumericVector piA, double s2A, NumericVector gammaB, NumericVector piB,
                   double s2B, double tau, double a, int threads) {
  const int nb = blk.size(), KA = gammaA.size(), KB = gammaB.size();
  std::vector<NumericMatrix> XA(nb), Cm(nb), XCa(nb), XCb(nb);
  std::vector<NumericVector> sl(nb), ma(nb), aa(nb), gc(nb);
  for (int b = 0; b < nb; b++) {
    List bk = blk[b];
    XA[b] = as<NumericMatrix>(bk["XA"]); Cm[b] = as<NumericMatrix>(bk["Cm"]); XCa[b] = as<NumericMatrix>(bk["XCa"]);
    XCb[b] = as<NumericMatrix>(bk["XCb"]); sl[b] = as<NumericVector>(bk["sl"]);
    ma[b] = as<NumericVector>(mA[b]); aa[b] = as<NumericVector>(al[b]); gc[b] = as<NumericVector>(gC[b]);
  }
  double* pr = r.begin(); const double* pp = p.begin();
  std::vector<double> lpA(KA), lpB(KB);
  for (int k = 0; k < KA; k++) lpA[k] = std::log(piA[k]);
  for (int k = 0; k < KB; k++) lpB[k] = std::log(piB[k]);
  std::vector<double> sphiA(KA, 0.0), sphiB(KB, 0.0); double saA = 0, snzA = 0, saB = 0, snzB = 0, mx = 0;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> lw(std::max(KA, KB)), mu(lw.size()), v(lw.size()), shA(KA, 0.0), shB(KB, 0.0);
    double sAl = 0, nAl = 0, sBl = 0, nBl = 0, mxl = 0;
    // mr.ash posterior mean for precision D and score rho (own term included); accumulates EM sums
    auto mix = [&](double D, double rho, const NumericVector& g, const std::vector<double>& lp, double s2,
                   std::vector<double>& sh, double& sa, double& snz) {
      const int K = g.size(); double m = -INFINITY;
      for (int k = 0; k < K; k++) {
        if (g[k] == 0) { lw[k] = lp[k]; mu[k] = 0; v[k] = 0; }
        else { const double gg = g[k] * s2; v[k] = 1 / (D + 1 / gg); mu[k] = v[k] * rho;
               lw[k] = lp[k] + 0.5 * std::log(v[k] / gg) + 0.5 * mu[k] * mu[k] / v[k]; }
        if (lw[k] > m) m = lw[k];
      }
      double sw = 0; for (int k = 0; k < K; k++) { lw[k] = std::exp(lw[k] - m); sw += lw[k]; }
      double nw = 0;
      for (int k = 0; k < K; k++) { const double phi = lw[k] / sw; sh[k] += phi; nw += phi * mu[k];
        if (g[k] > 0) { sa += phi * (mu[k] * mu[k] + v[k]) / g[k]; snz += phi; } }
      return nw;
    };
    #pragma omp for schedule(dynamic)
    for (int b = 0; b < nb; b++) {
      const int rk = XA[b].nrow(), kB = sl[b].size();
      double* ra = pr + off[b]; double* rb = ra + rk; const double* pa = pp + off[b]; const double* pb = pa + rk;
      for (int i = 0; i < XA[b].ncol(); i++) {             // A
        const double* x = &XA[b](0, i); double D = 0, rho = 0;
        for (int j = 0; j < rk; j++) { const double px = pa[j] * x[j]; D += px * x[j]; rho += px * ra[j]; }
        if (D <= 0) continue;
        const double nw = mix(D, rho + D * ma[b][i], gammaA, lpA, s2A, shA, sAl, nAl), dm = nw - ma[b][i];
        if (dm != 0) { for (int j = 0; j < rk; j++) ra[j] -= x[j] * dm; ma[b][i] = nw; }
      }
      for (int j = 0; j < kB; j++) {                        // B components
        const double* c = rk ? &Cm[b](0, j) : nullptr; const double s = sl[b][j];
        double D = pb[j] * s * s, rho = pb[j] * s * rb[j];
        for (int q = 0; q < rk; q++) { const double pc = pa[q] * c[q]; D += pc * c[q]; rho += pc * ra[q]; }
        const double nw = mix(D, rho + D * aa[b][j], gammaB, lpB, s2B, shB, sBl, nBl), dm = nw - aa[b][j];
        if (dm != 0) { rb[j] -= s * dm; for (int q = 0; q < rk; q++) ra[q] -= c[q] * dm; aa[b][j] = nw; }
      }
      for (int i = 0; i < XCb[b].ncol(); i++) {            // C (gamma in beta units)
        const double* xa = rk ? &XCa[b](0, i) : nullptr; const double* xb = &XCb[b](0, i); double D = 0, rho = 0;
        for (int q = 0; q < rk; q++) { const double px = pa[q] * xa[q]; D += px * xa[q]; rho += px * ra[q]; }
        for (int j = 0; j < kB; j++) { const double px = pb[j] * xb[j]; D += px * xb[j]; rho += px * rb[j]; }
        if (D <= 0) continue;
        const double sd = std::sqrt(D), old = gc[b][i];
        const double nw = mcp_prox((rho + D * old) / sd, tau, a) / sd, dm = nw - old;
        if (dm != 0) { mxl = std::max(mxl, std::fabs(dm) * sd);
          for (int q = 0; q < rk; q++) ra[q] -= xa[q] * dm; for (int j = 0; j < kB; j++) rb[j] -= xb[j] * dm; gc[b][i] = nw; }
      }
    }
    #pragma omp critical
    {
      for (int k = 0; k < KA; k++) sphiA[k] += shA[k];
      for (int k = 0; k < KB; k++) sphiB[k] += shB[k];
      saA += sAl; snzA += nAl; saB += sBl; snzB += nBl; mx = std::max(mx, mxl);
    }
  }
  return List::create(_["sphiA"] = NumericVector(sphiA.begin(), sphiA.end()), _["saA"] = saA, _["snzA"] = snzA,
                      _["sphiB"] = NumericVector(sphiB.begin(), sphiB.end()), _["saB"] = saB, _["snzB"] = snzB, _["dz"] = mx);
}

// pass 2 on ab.bin: beta_B = Q2 alpha on the B rows (A rows 0; the caller adds beta_A and gamma), K traits per sweep
// [[Rcpp::export]]
List ab_beta_cpp(CharacterVector files, List alpha, double thresh, int threads) {
  const int nb = files.size(), K = alpha.size();
  std::vector<std::string> fs(nb);
  for (int b = 0; b < nb; b++) fs[b] = as<std::string>(files[b]);
  std::vector<std::vector<std::vector<double>>> al(K, std::vector<std::vector<double>>(nb)), beta = al;
  for (int t = 0; t < K; t++) { List x = alpha[t]; for (int b = 0; b < nb; b++) al[t][b] = as<std::vector<double>>(x[b]); }
  std::vector<int> err(nb, 0);
  #pragma omp parallel for num_threads(threads) schedule(dynamic)
  for (int b = 0; b < nb; b++) {
    AbFile f;
    if (!f.open(fs[b])) { err[b] = 1; continue; }
    f.cut(thresh);
    std::vector<char> isA(f.m, 0);
    for (int q = 0; q < f.mA; q++) isA[f.idxA[q]] = 1;
    for (int t = 0; t < K; t++) {
      if (al[t][b].empty()) continue;
      if ((int)al[t][b].size() != f.kB) { err[b] = 1; break; }
      beta[t][b].assign(f.m, 0.0);
    }
    if (err[b]) continue;
    for (int j = 0; j < f.kB; j++) {
      const float* u = f.UlB + (size_t)j * f.m;
      for (int t = 0; t < K; t++) {
        if (al[t][b].empty() || al[t][b][j] == 0) continue;
        const double x = al[t][b][j]; double* bt = beta[t][b].data();
        for (int i = 0; i < f.m; i++) if (!isA[i]) bt[i] += u[i] * x;
      }
    }
  }
  for (int b = 0; b < nb; b++) if (err[b]) stop("cannot read " + fs[b] + " (or its components differ from pass 1)");
  List out(K);
  for (int t = 0; t < K; t++) { List o(nb); for (int b = 0; b < nb; b++) o[b] = beta[t][b]; out[t] = o; }
  return out;
}
