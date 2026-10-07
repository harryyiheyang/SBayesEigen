# ABC (Yihe 2026-10-07): beta = beta_A + U alpha + gamma on one eigen.bin, fitted jointly.
#   A: annotation SNPs with |z| > zA, thinned to greedy r2 < r2A leads (largest |z| first); mr.ash in beta space,
#      grid exact 0 / 1 / 10.
#   B: eigen-space mixture alpha (exact 0 / 1 / 100 / 500), one EM map per iteration, no prior on its variance.
#   C: SNPs with |z| > zC not in A, MCP(tau, a) on the standardised scale (threshold tau in residual z noise sd,
#      i.e. tau sqrt(ve0) under constant noise; kappa enters through the per-component precision).
# One iteration: A sweep, B EM map, C sweep, all on the shared component-space residual (no outer loop).
# Stops when the relative change of beta'R beta < tol and no C effect moved by more than stopz (z units).
.abc_const <- list(zA = 4.5, zC = 4, r2A = 0.9, gA = c(0, 1, 10), gB = c(0, 1, 100, 500))

# candidate rows (zero-based, per block) of one trait: typed (not imputed) SNPs with |z| above the lower of zA and zC
.abc_rows <- function(x) Map(function(ti, z, ty) ti[ty & abs(z) > min(.abc_const$zA, .abc_const$zC)], x$ti, x$z, x$typ)

# roles of one trait's candidates in one block. U: candidate rows x k, lam: k, z: candidate z, ann: in annotation
.abc_roles <- function(U, lam, z, ann) {
  role <- integer(length(z))
  a <- which(ann & abs(z) > .abc_const$zA)
  if (length(a)) {
    Ua <- U[a, , drop = FALSE]
    R <- Ua %*% (lam * t(Ua)); d <- sqrt(pmax(diag(R), 1e-12)); R2 <- (R / outer(d, d))^2
    keep <- logical(length(a))
    for (i in order(-abs(z[a]))) if (!any(keep & R2[i, ] > .abc_const$r2A)) keep[i] <- TRUE
    role[a[keep]] <- 1L
  }
  role[role == 0L & abs(z) > .abc_const$zC] <- 2L
  role
}

# w, c (= sqrt(lambda)), p (= n / residual variance): per component, blocks pooled; kb: components per block;
# X: per block k_b x s_b candidate columns sqrt(lambda) * U[i, ]; role: per block, 1 = A, 2 = C, 0 = unused
abc_vi <- function(w, c, p, kb, X, role, tau = 5, a = 2.5, gA = .abc_const$gA, gB = .abc_const$gB,
                   tol = 5e-4, stopz = 0.05, maxit = 1000, threads = 4) {
  off <- as.integer(c(0, cumsum(kb))[seq_along(kb)])
  coef <- lapply(role, function(r) numeric(length(r)))
  nA <- sum(unlist(role) == 1L); nC <- sum(unlist(role) == 2L)
  piA <- rep(1 / length(gA), length(gA)); piB <- rep(1 / length(gB), length(gB))
  s2a <- max(sum((p * w^2 - 1) * p * c^2) / sum((p * c^2)^2), 1e-8 / sum(c^2))
  s2B <- s2a / sum(piB * gB); s2A <- 1e-3 / max(nA, 1) / sum(piA * gA) * 10
  Ea <- numeric(length(w)); r <- w + 0; vg_old <- Inf; dz <- 0
  for (it in 1:maxit) {
    if (nA > 0) {
      sw <- abc_sweep_cpp(X, off, p, r, role, coef, gA, piA, s2A, tau, a, TRUE, FALSE, threads)
      piA <- pmax(sw$sphi / nA, 1e-300); s2A <- max(sw$sa / max(sw$snz, 1e-12), 1e-12)
    }
    rB <- r + c * Ea
    e <- em_step_cpp(rB, c, p, gB, piB, s2B, 1, 0, -2, 1, threads)
    piB <- e$pi; s2B <- e$sigma2
    Ea <- post_cpp(rB, c, p, gB, piB, s2B, 1, threads)$alpha
    r <- rB - c * Ea
    if (nC > 0) dz <- abc_sweep_cpp(X, off, p, r, role, coef, gA, piA, s2A, tau, a, FALSE, TRUE, threads)$dz
    vg <- sum((w - r)^2)
    if (it > 3 && abs(vg - vg_old) < tol * vg && dz < stopz) break
    vg_old <- vg
  }
  list(alpha = Ea, coef = coef, Vg = vg, Vg_B = sum((c * Ea)^2), piA = piA, s2A = s2A, piB = piB, s2B = s2B,
       nA = nA, nC_cand = nC, nC = sum(unlist(coef)[unlist(role) == 2L] != 0), iter = it, converged = it < maxit)
}
