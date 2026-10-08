# ABC (Yihe 2026-10-07): beta = beta_A + U alpha + gamma on one eigen.bin, fitted jointly.
#   A: annotation SNPs with |z| > zA, thinned to greedy r2 < r2A leads (largest |z| first); mr.ash in beta space,
#      grid exact 0 / 1 / 10.
#   B: eigen-space mixture alpha (exact 0 / 1 / 100 / 500), one EM map per iteration, no prior on its variance.
#   C: SNPs with |z| > zC not in A, MCP(tau, a) on the standardised scale (threshold tau in residual z noise sd,
#      i.e. tau sqrt(ve0) under constant noise; kappa enters through the per-component precision).
# One iteration: A sweep, B EM map, C sweep, all on the shared component-space residual (no outer loop).
# Stops when the relative change of beta'R beta < tol and no C effect moved by more than stopz (z units).
.abc_const <- list(zA = 4.5, zC = 4, r2A = 0.9, gA = c(0, 1, 10), gB = c(0, 1, 100, 500))
# Diagnostic switch (not exported; real-data chr1 ablation 2026-10-07): options(SBayesEigen.bprior = "eigen") gives B
# the eigen-VI prior (grid 0/1e-4/.../1, scaled-inv-chi2 on tau centred at the LDSC h2) and
# options(SBayesEigen.nemB = k) runs k EM maps of B's hyperparameters per iteration (whole-U ABC only).
.abc_bprior <- function(h2p, sc2) {
  if (!identical(getOption("SBayesEigen.bprior", "flat"), "eigen") || is.null(h2p))
    return(list(gB = .abc_const$gB, s2p = 0, nu = -2, A0 = 1))
  gB <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1)
  list(gB = gB, s2p = 2 * h2p, nu = 4, A0 = mean(gB) * sc2)
}

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
                   tol = 5e-4, stopz = 0.05, maxit = 1000, threads = 4, h2p = NULL) {
  bp <- .abc_bprior(h2p, sum(c^2)); if (!is.null(h2p)) gB <- bp$gB; nem <- getOption("SBayesEigen.nemB", 1L)
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
    for (q in seq_len(nem)) { e <- em_step_cpp(rB, c, p, gB, piB, s2B, 1, bp$s2p, bp$nu, bp$A0, threads)
      piB <- e$pi; s2B <- e$sigma2 }
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

# ABC on ab.bin (LDbuild with A; LD-build thread AB_JOINT_FIT.md): A = the stored annotation set, all in beta
# space; B = the Schur-complement eigen components; C = typed B SNPs with |z| > zC. w = (w_A, w_B) per block,
# p = per-component precision; blk[[b]] = list(XA, Cm, sl = sqrt(lamB), XCa, XCb). One CAVI sweep of A, B and C
# per iteration (abj_sweep_cpp); stops as abc_vi.
abj_vi <- function(w, p, blk, tau = 5, a = 2.5, gA = .abc_const$gA, gB = .abc_const$gB,
                   tol = 5e-4, stopz = 0.05, maxit = 1000, threads = 4, h2p = NULL) {
  if (!all(is.finite(w)) || !all(is.finite(p) & p > 0)) stop("abj_vi: non-finite w or non-positive precision p")
  rk <- vapply(blk, function(x) nrow(x$XA), 0L); kB <- lengths(lapply(blk, `[[`, "sl"))
  off <- as.integer(c(0, cumsum(rk + kB))[seq_along(blk)])
  iB <- unlist(lapply(seq_along(blk), function(b) off[b] + rk[b] + seq_len(kB[b])))
  sl <- unlist(lapply(blk, `[[`, "sl")); nA <- sum(vapply(blk, function(x) ncol(x$XA), 0L)); nCc <- sum(vapply(blk, function(x) ncol(x$XCb), 0L))
  mA <- lapply(blk, function(x) numeric(ncol(x$XA))); al <- lapply(kB, numeric); gC <- lapply(blk, function(x) numeric(ncol(x$XCb)))
  piA <- rep(1 / length(gA), length(gA)); piB <- rep(1 / length(gB), length(gB))
  pB <- p[iB]; wB <- w[iB]
  bp <- .abc_bprior(h2p, sum(sl^2)); if (!is.null(h2p)) gB <- bp$gB; piB <- rep(1 / length(gB), length(gB))
  s2a <- max(sum((pB * wB^2 - 1) * pB * sl^2) / sum((pB * sl^2)^2), 1e-8 / sum(sl^2))
  s2B <- s2a / sum(piB * gB); s2A <- 1e-3 / max(nA, 1) / sum(piA * gA) * 10
  r <- w + 0; vg_old <- Inf
  for (it in 1:maxit) {
    s <- abj_sweep_cpp(blk, off, p, r, mA, al, gC, gA, piA, s2A, gB, piB, s2B, tau, a, threads)
    if (nA > 0) { piA <- pmax(s$sphiA / nA, 1e-300); s2A <- max(s$saA / max(s$snzA, 1e-12), 1e-12) }
    piB <- pmax(s$sphiB / length(sl), 1e-300); s2B <- max((s$saB + bp$s2p / bp$A0) / max(s$snzB + bp$nu + 2, 1e-12), 1e-12)
    vg <- sum((w - r)^2)   # whitened fit: (w - r)'(w - r) = beta' R beta
    if (it > 3 && abs(vg - vg_old) < tol * vg && s$dz < stopz) break
    vg_old <- vg
  }
  list(alpha = al, mA = mA, gC = gC, Vg = vg, Vg_B = sum((sl * unlist(al))^2), piA = piA, s2A = s2A, piB = piB, s2B = s2B,
       nA = nA, nC_cand = nCc, nC = sum(unlist(gC) != 0), iter = it, converged = it < maxit)
}
