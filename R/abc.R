# ABC (Yihe 2026-10-07): beta = beta_A + U alpha + gamma on one eigen.bin, fitted jointly.
#   A: annotation SNPs with |z| > zA, thinned to greedy r2 < r2A leads (largest |z| first); mr.ash in beta space,
#      grid exact 0 / 1 / 10.
#   B: eigen-space mixture alpha (exact 0 / 1 / 100 / 500), SQUAREM-EM of its hyperparameters (full in iterations 1-3, then <= 3 steps), no prior on its variance.
#   C: SNPs with |z| > zC not in A, MCP(tau, a) on the standardised scale (threshold tau in residual z noise sd,
#      i.e. tau sqrt(ve0) under constant noise; kappa enters through the per-component precision).
# One iteration: A sweep, B EM map, C sweep, all on the shared component-space residual (no outer loop).
# Stops when the relative change of beta'R beta < tol and no C effect moved by more than stopz (z units).
.abc_const <- list(zA = 4.5, zC = 4, r2A = 0.9, gA = c(0, 1, 10), gB = c(0, 1, 100, 500), tolB = 1e-3)

# start of E[alpha^2] for B: the moment estimate, but at least h2p / sum(lambda), i.e. Vg_B = LDSC h2 at the start.
# A start at ~0 is an absorbing state of the EM: with s2B -> 0 every class has the same likelihood, pi_B stays at
# its start (uniform) and s2B at its value (real 7M fits 2026-10-09, where overstated noise made the moment negative).
# B hyperparameters (pi_B, s2B) by SQUAREM-EM to convergence on the independent-normal-means problem
# y_j = c_j alpha_j + e_j, e_j ~ N(0, 1 / p_j), with the current residual (other parts fixed), as vi_eigen does.
# One plain EM map per outer iteration crawls (2026-10-09 real 7M fits); this is the B update of each iteration.
# Run to convergence only in the first iterations; after that, warm-started, at most 3 SQUAREM steps with tolerance
# 0.1 x the last outer change (0ea0e9a ran to 1e-6 every iteration: 18-trait HPC fit > 4.5x slower, 2026-10-10).
.abc_emB_ctl <- function(it, dh) if (it <= 3) list(maxit = 200, tol = 1e-6) else list(maxit = 3, tol = max(1e-6, 0.1 * dh))
.abc_emB <- function(y, c, p, gB, piB, s2B, bp, threads, maxit = 200, tol = 1e-6) {
  K <- length(gB)
  em <- function(th) {
    q <- exp(th[1:K] - max(th[1:K])); q <- q / sum(q)
    e <- em_step_cpp(y, c, p, gB, q, exp(th[K + 1]), 1, bp$s2p, bp$nu, bp$A0, threads)
    e$th <- c(log(pmax(e$pi, 1e-300)), log(max(e$sigma2, 1e-300))); e
  }
  th <- c(log(pmax(piB, 1e-300)), log(s2B)); e <- em(th)
  for (it in seq_len(maxit)) {
    e1 <- em(e$th); rr <- e$th - th; v <- e1$th - e$th - rr
    live <- c(e$th[1:K] > log(1e-8), TRUE)
    a <- min(-sqrt(sum(rr[live]^2) / max(sum(v[live]^2), 1e-300)), -1)
    tp <- th - 2 * a * rr + a^2 * v; tp[1:K] <- pmax(tp[1:K], -700); tp[K + 1] <- min(max(tp[K + 1], -50), 50)
    en <- em(em(tp)$th)
    if (!is.finite(en$obj) || en$obj < e1$obj) en <- em(e1$th)
    dth <- max(abs(en$th - th)); th <- en$th; e <- en
    if (dth < tol) break
  }
  q <- exp(th[1:K] - max(th[1:K])); list(pi = q / sum(q), s2 = exp(th[K + 1]))
}

.abc_s2B0 <- function(p, w, c, h2p) {
  mom <- sum((p * w^2 - 1) * p * c^2) / sum((p * c^2)^2)
  max(mom, (if (is.null(h2p)) 1e-8 else h2p) / sum(c^2))
}
# B prior (default since 2026-10-10, HPC: the eigen-VI prior recovered Vg on every trait and was insensitive to its centre
# x0.25 / x4, while 0/1/100/500 left pi_B uniform): grid 0/1e-4/.../1, scaled-inv-chi2 (nu 4) on tau centred at the
# LDSC h2. options(SBayesEigen.bprior = "flat") restores the old 0/1/100/500 grid without a prior;
# options(SBayesEigen.nemB = k) runs k plain EM maps of B's hyperparameters per iteration instead of SQUAREM-EM to convergence (whole-U ABC only; diagnostic).
.abc_bprior <- function(h2p, sc2) {
  if (identical(getOption("SBayesEigen.bprior", "eigen"), "flat") || is.null(h2p))
    return(list(gB = .abc_const$gB, s2p = 0, nu = -2, A0 = 1))
  gB <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1); nu <- .abc_nu("B")
  list(gB = gB, s2p = (nu - 2) * h2p, nu = nu, A0 = mean(gB) * sc2)
}
# Priors on the A and B variance scales (Yihe 2026-10-10 06:20: both centred at h2, A concentrated, B wide).
# For a part with grid g and scale s2, tau = s2 * mean(g) * S (S = sum(lambda_B) for B, the number of A SNPs for A,
# i.e. tau ~ Vg of that part with equal class weights; A in beta space with standardised SNPs and thinned leads,
# so beta_A' R_AA beta_A ~ sum beta_A^2) has a scaled-inv-chi2(nu, s^2) prior with mean nu s^2 / (nu - 2) = LDSC h2.
# M-step: s2 = (sum_k phi (mu^2 + v) / g_k + (nu - 2) h2 / A0) / (sum phi_nonzero + nu + 2), A0 = mean(g) * S.
# Both parts are centred at h2, so the prior expectation of the total Vg is about 2 h2 (the data decide the split).
# nu: options(SBayesEigen.nuA = 50, SBayesEigen.nuB = 4) (> 2); nuA = 0 turns the A prior off (old plain M-step).
.abc_nu <- function(part) {
  nu <- getOption(paste0("SBayesEigen.nu", part), if (part == "A") 50 else 4)
  if (!(part == "A" && nu == 0) && !(is.numeric(nu) && nu > 2)) stop("SBayesEigen.nu", part, " must be > 2", if (part == "A") " (or 0)")
  nu
}
.abc_aprior <- function(h2p, nA, gA) {
  nu <- .abc_nu("A")
  if (is.null(h2p) || nu == 0 || nA == 0) return(list(s2p = 0, nu = -2, A0 = 1, s2_0 = NULL))
  A0 <- mean(gA) * nA; list(s2p = (nu - 2) * h2p, nu = nu, A0 = A0, s2_0 = h2p / A0)
}
.abc_mA <- function(sa, snz, ap) max((sa + ap$s2p / ap$A0) / max(snz + ap$nu + 2, 1e-12), 1e-12)

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
  bp <- .abc_bprior(h2p, sum(c^2)); if (!is.null(h2p)) gB <- bp$gB; nem <- getOption("SBayesEigen.nemB", 0L)
  off <- as.integer(c(0, cumsum(kb))[seq_along(kb)])
  coef <- lapply(role, function(r) numeric(length(r)))
  nA <- sum(unlist(role) == 1L); nC <- sum(unlist(role) == 2L)
  piA <- rep(1 / length(gA), length(gA)); piB <- rep(1 / length(gB), length(gB))
  s2B <- .abc_s2B0(p, w, c, h2p) / sum(piB * gB); s2A <- 1e-3 / max(nA, 1) / sum(piA * gA) * 10
  ap <- .abc_aprior(h2p, nA, gA); if (!is.null(ap$s2_0)) s2A <- ap$s2_0
  Ea <- numeric(length(w)); r <- w + 0; vg_old <- Inf; dz <- 0; dh <- Inf; hy_old <- c(piB, log(s2B))
  for (it in 1:maxit) {
    if (nA > 0) {
      sw <- abc_sweep_cpp(X, off, p, r, role, coef, gA, piA, s2A, tau, a, TRUE, FALSE, threads)
      piA <- pmax(sw$sphi / nA, 1e-300); s2A <- .abc_mA(sw$sa, sw$snz, ap)
    }
    rB <- r + c * Ea
    if (nem > 0) for (q in seq_len(nem)) { e <- em_step_cpp(rB, c, p, gB, piB, s2B, 1, bp$s2p, bp$nu, bp$A0, threads)
      piB <- e$pi; s2B <- e$sigma2 }
    else { ct <- .abc_emB_ctl(it, dh); h <- .abc_emB(rB, c, p, gB, piB, s2B, bp, threads, ct$maxit, ct$tol); piB <- h$pi; s2B <- h$s2 }
    Ea <- post_cpp(rB, c, p, gB, piB, s2B, 1, threads)$alpha
    r <- rB - c * Ea
    if (nC > 0) dz <- abc_sweep_cpp(X, off, p, r, role, coef, gA, piA, s2A, tau, a, FALSE, TRUE, threads)$dz
    vg <- sum((w - r)^2)
    hy <- c(piB, log(s2B)); dh <- max(abs(hy - hy_old)); hy_old <- hy
    if (it > 3 && abs(vg - vg_old) < tol * vg && dz < stopz && dh < .abc_const$tolB) break
    vg_old <- vg
  }
  list(alpha = Ea, coef = coef, Vg = vg, Vg_B = sum((c * Ea)^2), piA = piA, s2A = s2A, piB = piB, s2B = s2B,
       nA = nA, nC_cand = nC, nC = sum(unlist(coef)[unlist(role) == 2L] != 0), iter = it, converged = it < maxit)
}

# ABC on ab.bin (LDbuild with A; LD-build thread AB_JOINT_FIT.md): A = the stored annotation set, all in beta
# space; B = the Schur-complement eigen components; C = typed B SNPs with |z| > zC. w = (w_A, w_B) per block,
# p = per-component precision; blk[[b]] = list(XA, Cm, sl = sqrt(lamB), XCa, XCb). One CAVI sweep of A, B and C
# per iteration (abj_sweep_cpp); stops as abc_vi.
# conditional B scores after a sweep: y_j = (rho_j + D_j alpha_j) / D_j with
# rho_j = p_Bj s_j r_Bj + sum_q p_Aq Cm_qj r_Aq (the score the sweep's mixture update used)
.abj_yD <- function(blk, off, p, r) {   # B residual scores rho; conditional score y = rho / D + alpha
  unlist(lapply(seq_along(blk), function(b) { x <- blk[[b]]; rk <- nrow(x$XA); kB <- length(x$sl)
    ia <- off[b] + seq_len(rk); ib <- off[b] + rk + seq_len(kB)
    rho <- p[ib] * x$sl * r[ib] + if (rk) drop(crossprod(x$Cm, p[ia] * r[ia])) else 0
    rho }))
}

# whitened fit of each part per block (A rows | B rows): A = (XA m_A | 0), B = (Cm alpha | sl alpha),
# C = (XCa gamma | XCb gamma); Vg = |fA + fB + fC|^2 = Vg_A + Vg_B + Vg_C + 2 (fA'fB + fA'fC + fB'fC) (crosses reported x2)
.abj_parts <- function(blk, mA, al, gC, w = NULL) {
  rk0 <- vapply(blk, function(x) nrow(x$XA), 0L); off <- c(0, cumsum(rk0 + lengths(lapply(blk, `[[`, "sl"))))
  o <- vapply(seq_along(blk), function(b) { x <- blk[[b]]; rk <- nrow(x$XA)
    fA <- c(if (rk) drop(x$XA %*% mA[[b]]) else numeric(0), numeric(length(x$sl)))
    fB <- c(if (rk) drop(x$Cm %*% al[[b]]) else numeric(0), x$sl * al[[b]])
    fC <- if (ncol(x$XCb)) c(if (rk) drop(x$XCa %*% gC[[b]]) else numeric(0), drop(x$XCb %*% gC[[b]])) else 0 * fB
    wb <- if (is.null(w)) 0 * fA else w[off[b] + seq_along(fA)]
    c(A = sum(fA^2), B = sum(fB^2), C = sum(fC^2), AB = 2 * sum(fA * fB), AC = 2 * sum(fA * fC), BC = 2 * sum(fB * fC),
      wA = sum(fA * wb), wB = sum(fB * wb), wC = sum(fC * wb)) },
    numeric(9))
  rowSums(o)
}

abj_vi <- function(w, p, blk, tau = 5, a = 2.5, gA = .abc_const$gA, gB = .abc_const$gB,
                   tol = 5e-4, stopz = 0.05, maxit = 1000, threads = 4, h2p = NULL) {
  if (!all(is.finite(w)) || !all(is.finite(p) & p > 0)) stop("abj_vi: non-finite w or non-positive precision p")
  rk <- vapply(blk, function(x) nrow(x$XA), 0L); kB <- lengths(lapply(blk, `[[`, "sl"))
  off <- as.integer(c(0, cumsum(rk + kB))[seq_along(blk)])
  iB <- unlist(lapply(seq_along(blk), function(b) off[b] + rk[b] + seq_len(kB[b])))
  sl <- unlist(lapply(blk, `[[`, "sl")); nA <- sum(vapply(blk, function(x) ncol(x$XA), 0L)); nCc <- sum(vapply(blk, function(x) ncol(x$XCb), 0L))
  mA <- lapply(blk, function(x) numeric(ncol(x$XA))); al <- lapply(kB, numeric); gC <- lapply(blk, function(x) numeric(ncol(x$XCb)))
  piA <- rep(1 / length(gA), length(gA)); piB <- rep(1 / length(gB), length(gB))
  # fitted B components: the first nB of each block (threshB); the rest keep alpha = 0 but stay in the data (A, C see them)
  fB <- unlist(lapply(blk, function(x) seq_along(x$sl) <= (if (is.null(x$nB)) length(x$sl) else x$nB)))
  pB <- p[iB][fB]; wB <- w[iB][fB]; slf <- sl[fB]
  bp <- .abc_bprior(h2p, sum(slf^2)); if (!is.null(h2p)) gB <- bp$gB; piB <- rep(1 / length(gB), length(gB))
  s2B <- .abc_s2B0(pB, wB, slf, h2p) / sum(piB * gB); s2A <- 1e-3 / max(nA, 1) / sum(piA * gA) * 10
  ap <- .abc_aprior(h2p, nA, gA); if (!is.null(ap$s2_0)) s2A <- ap$s2_0
  DB <- unlist(lapply(seq_along(blk), function(b) { x <- blk[[b]]; pa <- p[off[b] + seq_len(rk[b])]
    p[off[b] + rk[b] + seq_len(kB[b])] * x$sl^2 + if (rk[b]) colSums(x$Cm^2 * pa) else 0 }))[fB]
  # s2B is per unit of alpha; with y = score / D the design is 1, so the prior centre A0 uses sum(sl^2) as before
  r <- w + 0; vg_old <- Inf; dh <- Inf; hy_old <- c(piB, log(s2B))
  # options(SBayesEigen.blockve = "rho") (benchmark thread 2026-10-10, after SBayesRC's allMixVe): in a block with
  # rho_b = beta'beta / beta'R beta > 1.1 (here (|alpha|^2 + |beta_A|^2 + |gamma|^2) / |fit_b|^2, Q2 orthonormal, cross
  # terms ignored) the noise multiplier is s_b = (sum p0 r^2 + nu) / (k_b + nu), nu = 4, snapped to 1/1.1/1.25/1.5/2/3
  # (residual only: no posterior-variance term). Off by default; for ablation.
  brho <- identical(getOption("SBayesEigen.blockve", FALSE), "rho"); p0 <- p; sb <- rep(1, length(blk))
  bid <- rep(seq_along(blk), rk + kB); dsb <- 0
  for (it in 1:maxit) {
    s <- abj_sweep_cpp(blk, off, p, r, mA, al, gC, gA, piA, s2A, gB, piB, s2B, tau, a, threads)
    if (brho && it > 1) {
      fw2 <- vapply(split((w - r)^2, bid), sum, 0)
      bb2 <- vapply(seq_along(blk), function(b) sum(al[[b]]^2) + sum(mA[[b]]^2) + sum(gC[[b]]^2), 0)
      sv <- (vapply(split(p0 * r^2, bid), sum, 0) + 4) / (rk + kB + 4)
      gr <- c(1, 1.1, 1.25, 1.5, 2, 3); sn <- ifelse(bb2 > 1.1 * fw2 & fw2 > 0, gr[findInterval(sv, (gr[-1] + gr[-6]) / 2) + 1], sb)
      dsb <- max(abs(log(sn / sb))); sb <- sn; p <- p0 / sb[bid]
      DB <- unlist(lapply(seq_along(blk), function(b) { x <- blk[[b]]; pa <- p[off[b] + seq_len(rk[b])]
        p[off[b] + rk[b] + seq_len(kB[b])] * x$sl^2 + if (rk[b]) colSums(x$Cm^2 * pa) else 0 }))[fB]
    }
    if (nA > 0) { piA <- pmax(s$sphiA / nA, 1e-300); s2A <- .abc_mA(s$saA, s$snzA, ap) }
    # B: conditional score of each component given A, C and the other components, (y, D); SQUAREM-EM to convergence
    yd <- .abj_yD(blk, off, p, r)[fB] / DB + unlist(al)[fB]
    ct <- .abc_emB_ctl(it, dh); h <- .abc_emB(yd, rep(1, length(DB)), DB, gB, piB, s2B, bp, threads, ct$maxit, ct$tol); piB <- h$pi; s2B <- max(h$s2, 1e-12)
    vg <- sum((w - r)^2)   # whitened fit: (w - r)'(w - r) = beta' R beta
    hy <- c(piB, log(s2B)); dh <- max(abs(hy - hy_old), dsb); hy_old <- hy
    if (it > 3 && abs(vg - vg_old) < tol * vg && s$dz < stopz && dh < .abc_const$tolB) break
    vg_old <- vg
  }
  # Pratt split (Yihe 2026-10-10 06:33): moment version with the GWAS as y, V_k = beta_k' bhat = f_k' w (no LD needed),
  # share = V_k / beta' bhat. Vg_parts (debug): own squares, crosses (x2), and the LD version beta_k' R beta (sums to Vg).
  pt <- .abj_parts(blk, mA, al, gC, w)
  list(alpha = al, mA = mA, gC = gC, sb = if (brho) sb, Vg = vg, Vg_A = pt[["wA"]], Vg_B = pt[["wB"]], Vg_C = pt[["wC"]],
       Vg_parts = c(pt[1:6], RA = pt[["A"]] + (pt[["AB"]] + pt[["AC"]]) / 2, RB = pt[["B"]] + (pt[["AB"]] + pt[["BC"]]) / 2,
                    RC = pt[["C"]] + (pt[["AC"]] + pt[["BC"]]) / 2), piA = piA, s2A = s2A, piB = piB, s2B = s2B,
       fitw = w - r, nB_fit = sum(fB), nA = nA, nC_cand = nCc, nC = sum(unlist(gC) != 0), iter = it, converged = it < maxit)
}
