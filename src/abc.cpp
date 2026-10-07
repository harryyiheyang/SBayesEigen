// ABC model (Yihe 2026-10-07): beta = beta_A + U alpha + gamma, joint on one eigen.bin.
// In component space (blocks pooled): w = sum_{i in A} x_i beta_i + c alpha + sum_{i in C} x_i gamma_i + e,
// x_i = sqrt(lambda) * U[i, ] (k_b, zero outside its block), e_j ~ N(0, 1 / p_j), p_j = n_j / (ve0 + kappa / lambda_j).
// One call does one coordinate sweep over the A candidates (mr.ash: exact zero + gammaA * s2A) and/or the
// C candidates (MCP on the standardised scale h = gamma sqrt(D), D = sum_j p_j x_j^2, so the threshold is in
// units of the residual z noise sd: tau sqrt(ve0) under constant noise) on the residual r, updated in place.
#include <Rcpp.h>
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
