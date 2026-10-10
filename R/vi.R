# Model per eigen component j (blocks pooled):
#   w_j = c_j alpha_j + e_j,  e_j ~ N(0, ve / n_j),  w = Lambda^{-1/2} U' bhat,  c = sqrt(lambda)
#   alpha_j ~ pi_0 delta_0 + sum_k pi_k N(0, gamma_k sigma2);  beta = U alpha,  Vg = sum c_j^2 alpha_j^2
# The likelihood is diagonal, so the mean-field posterior q(alpha_j, z_j) is exact and variational
# EM is EM on the marginal likelihood. ve is fixed. Prior on tau = sigma2 A0 (A0 = sum(pi_start gamma)
# sum(c^2), fixed): scaled-inv-chi2(nu, s2) with E[tau] = h2p, i.e. s2p = nu s2 = (nu - 2) h2p.

# w, c, n: per component. h2p: prior centre (LDSC h2, floored by the caller). Returns pi, sigma2,
# Vg, Vg_sd, alpha (posterior means), traces, iter, converged.
vi_eigen <- function(w, c, n, ve, h2p, threads = 4, gamma = c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1),
                     tol = 1e-4, maxit = 2000) {
  K <- length(gamma); nu <- 4
  c2 <- c^2; sc2 <- sum(c2)
  pi <- rep(1 / K, K)
  # moment start E[n w^2] = n c^2 s2a + ve, floored
  s2a <- max(sum((n * w^2 - ve) * n * c2) / sum((n * c2)^2), 1e-8 / sc2)
  A0 <- sum(pi * gamma) * sc2
  s2p <- (nu - 2) * h2p
  em <- function(th) {
    p <- exp(th[1:K] - max(th[1:K])); p <- p / sum(p)
    e <- em_step_cpp(w, c, n, gamma, p, exp(th[K + 1]), ve, s2p, nu, A0, threads)
    e$th <- c(log(e$pi), log(e$sigma2)); e$th0 <- th
    e
  }
  th <- c(log(pi), log(s2a / sum(pi * gamma)))
  e <- em(th)
  ll <- e$ll; vg <- e$vg; obj <- e$obj
  for (it in 1:maxit) {
    # SQUAREM (Varadhan & Roland 2008), falling back to two plain EM steps when it does not ascend
    e1 <- em(e$th)
    rr <- e$th - th; v <- e1$th - e$th - rr
    live <- c(e$th[1:K] > log(1e-8), TRUE)   # step length from sigma2 and live classes only
    a <- min(-sqrt(sum(rr[live]^2) / max(sum(v[live]^2), 1e-300)), -1)
    tp <- th - 2 * a * rr + a^2 * v
    tp[1:K] <- pmax(tp[1:K], -700); tp[K + 1] <- min(max(tp[K + 1], -50), 50)
    en <- em(em(tp)$th)
    e2 <- em(e1$th)
    if (!is.finite(en$obj) || en$obj < e2$obj) en <- e2
    th <- en$th0; e <- en
    ll[it + 1] <- e$ll; vg[it + 1] <- e$vg; obj[it + 1] <- e$obj
    # stop when posterior Vg changes < tol (relative) twice in a row: pi/sigma2 crawl along a flat
    # ridge for hundreds of iterations while Vg and alpha barely move
    if (it > 4 && abs(vg[it + 1] - vg[it]) < tol * vg[it + 1] && abs(vg[it] - vg[it - 1]) < tol * vg[it]) break
  }
  p <- exp(th[1:K] - max(th[1:K])); p <- p / sum(p)
  ps <- post_cpp(w, c, n, gamma, p, exp(th[K + 1]), ve, threads)
  list(pi = p, sigma2 = exp(th[K + 1]), gamma = gamma, Vg = ps$Vg, Vg_sd = ps$Vg_sd, alpha = ps$alpha,
       loglik = ll, objective = obj, vg_trace = vg, iter = it, converged = it < maxit)
}

# LD score regression E[n bhat^2] = a + n h2 l / M, two reweighting steps (LDSC weights)
ldsc_eigen <- function(bhat, l, n, M = length(bhat)) {
  if (length(n) == 1) n <- rep(n, length(bhat))
  chi2 <- n * bhat^2
  X <- cbind(1, n * l / M)
  wt <- 1 / pmax(l, 1)
  for (it in 1:2) {
    f <- stats::lm.wfit(X, chi2, wt)
    cf <- pmax(f$coefficients, c(1e-3, 0))
    wt <- 1 / (pmax(l, 1) * (cf[1] + n * cf[2] * l / M)^2)
  }
  f <- stats::lm.wfit(X, chi2, wt)
  list(intercept = unname(f$coefficients[1]), h2 = unname(f$coefficients[2]))
}

# LD-mismatch noise per eigen component (Yihe 2026-10-07): var(e_j) = (ve0 + kappa / lambda_j) / n_j.
# MoM: E[n w^2] = ve0 + kappa / lambda + n lambda tau, tau = h2 / sum(lambda) from LDSC. The signal term is
# subtracted, components with lambda < lam0 are put in nbin equal-count bins by lambda, and the bin means
# are regressed on (1, 1/lambda) by weighted least squares (weights 1/fit^2, as LDSC), with
# ve0 in ve_range and kappa >= 0. Returns ve0, kappa and the number of components used.
noise_mom <- function(w, lam, n, h2, lam0 = 1, nbin = 100, ve_range = c(0.5, 3)) {
  tau <- h2 / sum(lam)
  s <- lam < lam0
  if (sum(s) < 2 * nbin) return(list(ve0 = 1, kappa = 0, n_comp = sum(s)))
  y <- (n * w^2 - n * lam * tau)[s]; x <- 1 / lam[s]
  bin <- cut(rank(-x, ties.method = "first"), nbin, labels = FALSE)
  my <- vapply(split(y, bin), mean, 0); mx <- vapply(split(x, bin), mean, 0)
  # box-constrained weighted LS for (ve0, kappa): best of the unconstrained solution and the four edges
  fit2 <- function(wt) {
    sse <- function(a, k) sum(wt * (my - a - k * mx)^2)
    kfix <- function(a) max(sum(wt * mx * (my - a)) / sum(wt * mx^2), 0)
    afix <- function(k) min(max(sum(wt * (my - k * mx)) / sum(wt), ve_range[1]), ve_range[2])
    f <- stats::lm.wfit(cbind(1, mx), my, wt)$coefficients
    cand <- list(c(ve_range[1], kfix(ve_range[1])), c(ve_range[2], kfix(ve_range[2])), c(afix(0), 0))
    if (f[1] >= ve_range[1] && f[1] <= ve_range[2] && f[2] >= 0) cand <- c(cand, list(unname(f)))
    cand[[which.min(vapply(cand, function(p) sse(p[1], p[2]), 0))]]
  }
  p <- fit2(rep(1, nbin))
  for (it in 1:2) p <- fit2(1 / pmax(p[1] + p[2] * mx, 1e-3)^2)
  list(ve0 = p[1], kappa = p[2], n_comp = sum(s))
}
