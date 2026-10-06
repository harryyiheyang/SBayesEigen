// Eigen-space VI for SBayesEigen.
// Model per component j (blocks pooled): w_j = c_j alpha_j + e_j, e_j ~ N(0, ve / n_j),
// w = D^{-1/2} U' bhat, c = sqrt(lambda), beta = U alpha, Vg = sum c_j^2 alpha_j^2;
// alpha_j ~ pi_0 delta_0 + sum_k pi_k N(0, gamma_k sigma2). The likelihood is diagonal, so the
// mean-field posterior is exact and variational EM is EM on the marginal likelihood.
#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "eig_file.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// pass 2: beta = U alpha per block for K traits in one sweep over U; alpha[[t]][[b]] empty: trait t
// skips block b (beta empty)
// [[Rcpp::export]]
List eig_beta_cpp(CharacterVector files, List alpha, double thresh, int threads) {
  const int nb = files.size(), K = alpha.size();
  std::vector<std::string> fs(nb);
  for (int b = 0; b < nb; b++) fs[b] = as<std::string>(files[b]);
  std::vector<std::vector<std::vector<double>>> al(K, std::vector<std::vector<double>>(nb)), beta = al;
  for (int t = 0; t < K; t++) {
    List a = alpha[t];
    for (int b = 0; b < nb; b++) al[t][b] = as<std::vector<double>>(a[b]);
  }
  std::vector<int> err(nb, 0);
  #pragma omp parallel for num_threads(threads) schedule(dynamic)
  for (int b = 0; b < nb; b++) {
    EigFile e;
    if (!e.open(fs[b], thresh)) { err[b] = 1; continue; }
    for (int t = 0; t < K; t++) {
      if (al[t][b].empty()) continue;
      if ((int)al[t][b].size() != e.k) { err[b] = 1; break; }
      beta[t][b].assign(e.m, 0.0);
    }
    if (err[b]) continue;
    for (int j = 0; j < e.k; j++) {
      const float* u = e.U + (size_t)j * e.m;
      for (int t = 0; t < K; t++) {
        if (al[t][b].empty()) continue;
        const double a = al[t][b][j];
        if (a == 0) continue;
        double* bt = beta[t][b].data();
        for (int i = 0; i < e.m; i++) bt[i] += u[i] * a;
      }
    }
  }
  for (int b = 0; b < nb; b++) if (err[b]) stop("cannot read " + fs[b] + " (or its components differ from pass 1)");
  List out(K);
  for (int t = 0; t < K; t++) {
    List o(nb);
    for (int b = 0; b < nb; b++) o[b] = beta[t][b];
    out[t] = o;
  }
  return out;
}

// one EM map at theta = (pi, sigma2), ve fixed: single pass over components, no J x K storage;
// returns ll, obj, vg at theta and the updated (pi, sigma2)
// [[Rcpp::export]]
List em_step_cpp(NumericVector w, NumericVector c, NumericVector n, NumericVector gamma, NumericVector pi,
                 double sigma2, double ve, double s2p, double nu, double A0, int threads) {
  int J = w.size(), K = gamma.size();
  const double *pw = w.begin(), *pc = c.begin(), *pn = n.begin();
  std::vector<double> g(gamma.begin(), gamma.end()), lpi(K), gs(K);
  for (int k = 0; k < K; k++) { lpi[k] = std::log(pi[k]); gs[k] = g[k] * sigma2; }
  double ll = 0, vg = 0, sa = 0, snz = 0;
  std::vector<double> sr(K, 0.0);
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> lp(K), M(K), S(K), srl(K, 0.0);
    double lll = 0, vgl = 0, sal = 0, snzl = 0;
    #pragma omp for schedule(static)
    for (int j = 0; j < J; j++) {
      double c2 = pc[j] * pc[j], tau = ve / pn[j], w2 = pw[j] * pw[j], prec = pn[j] * c2 / ve, mx = -1e300;
      for (int k = 0; k < K; k++) {
        double V = c2 * gs[k] + tau;
        lp[k] = lpi[k] - 0.5 * std::log(2 * M_PI * V) - 0.5 * w2 / V;
        if (lp[k] > mx) mx = lp[k];
        if (g[k] > 0) { S[k] = 1.0 / (prec + 1.0 / gs[k]); M[k] = S[k] * pn[j] * pc[j] * pw[j] / ve; }
        else { S[k] = 0; M[k] = 0; }
      }
      double se = 0; for (int k = 0; k < K; k++) { lp[k] = std::exp(lp[k] - mx); se += lp[k]; }
      lll += mx + std::log(se);
      double ea2 = 0;
      for (int k = 0; k < K; k++) {
        double r = lp[k] / se, e2 = M[k] * M[k] + S[k];
        srl[k] += r; ea2 += r * e2;
        if (g[k] > 0) { sal += r * e2 / g[k]; snzl += r; }
      }
      vgl += c2 * ea2;
    }
    #pragma omp critical
    {
      ll += lll; vg += vgl; sa += sal; snz += snzl;
      for (int k = 0; k < K; k++) sr[k] += srl[k];
    }
  }
  double tau0 = sigma2 * A0;
  double obj = ll - (nu / 2 + 1) * std::log(tau0) - s2p / (2 * tau0);
  NumericVector pin(K);
  for (int k = 0; k < K; k++) pin[k] = std::max(sr[k] / J, 1e-300);
  double s2n = (sa + s2p / A0) / (snz + nu + 2);
  return List::create(_["ll"] = ll, _["obj"] = obj, _["vg"] = vg, _["pi"] = pin, _["sigma2"] = s2n);
}

// posterior summaries at theta: E[alpha], Vg, Vg_sd
// [[Rcpp::export]]
List post_cpp(NumericVector w, NumericVector c, NumericVector n, NumericVector gamma, NumericVector pi,
              double sigma2, double ve, int threads) {
  int J = w.size(), K = gamma.size();
  std::vector<double> g(gamma.begin(), gamma.end()), lpi(K);
  for (int k = 0; k < K; k++) lpi[k] = std::log(pi[k]);
  NumericVector a(J); double *pa = a.begin();
  const double *pw = w.begin(), *pc = c.begin(), *pn = n.begin();
  double vg = 0, vv = 0;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> lp(K), M(K), S(K); double vgl = 0, vvl = 0;
    #pragma omp for schedule(static)
    for (int j = 0; j < J; j++) {
      double c2 = pc[j] * pc[j], tau = ve / pn[j], w2 = pw[j] * pw[j], prec = pn[j] * c2 / ve, mx = -1e300;
      for (int k = 0; k < K; k++) {
        double gs = g[k] * sigma2, V = c2 * gs + tau;
        lp[k] = lpi[k] - 0.5 * std::log(V) - 0.5 * w2 / V;
        if (lp[k] > mx) mx = lp[k];
        if (g[k] > 0) { S[k] = 1.0 / (prec + 1.0 / gs); M[k] = S[k] * pn[j] * pc[j] * pw[j] / ve; } else { S[k] = 0; M[k] = 0; }
      }
      double se = 0; for (int k = 0; k < K; k++) { lp[k] = std::exp(lp[k] - mx); se += lp[k]; }
      double e1 = 0, e2 = 0, e4 = 0;
      for (int k = 0; k < K; k++) {
        double r = lp[k] / se, m2 = M[k] * M[k];
        e1 += r * M[k]; e2 += r * (m2 + S[k]); e4 += r * (m2 * m2 + 6 * m2 * S[k] + 3 * S[k] * S[k]);
      }
      pa[j] = e1; vgl += c2 * e2; vvl += c2 * c2 * (e4 - e2 * e2);
    }
    #pragma omp critical
    { vg += vgl; vv += vvl; }
  }
  return List::create(_["alpha"] = a, _["Vg"] = vg, _["Vg_sd"] = std::sqrt(vv));
}
