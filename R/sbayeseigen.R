#' Genome-wide PRS from eigen LD and GWAS summary statistics
#'
#' All-in-one: \code{\link{tidy}} and \code{\link{impute}} in memory, LD score regression for
#' the heritability prior, variational EM in the eigen space and per-SNP effects. Each
#' \code{block<b>.eigen.bin} is read once for imputation and the rotation
#' \eqn{w = \Lambda^{-1/2} U' \hat b} (\eqn{\hat b = z / \sqrt{N + z^2}}, as SBayesRC), and once
#' more for \eqn{\beta = U E[\alpha]}. Summary data that are already imputed (an \code{r2}
#' column and one row per \code{snp.info} SNP) skip tidy and impute.
#'
#' Several traits run together read each eigen file the same two times as one trait: imputation,
#' the rotation and \eqn{\beta = U E[\alpha]} are done for all traits while a block is in memory.
#' LDSC and VI are fitted per trait, so every trait gets the same result as on its own.
#'
#' @param ma GWAS summary statistics in COJO format (SNP A1 A2 freq b se p N) or the simplified format
#'   SNP A1 A2 A1freq Z N (see \code{\link{tidy}}; output beta is then the per-dosage effect in phenotype
#'   SD with the LD reference freq), raw or already
#'   imputed: a path or a data.frame, or for several traits a character vector of paths or a list
#'   of data.frames. Trait names are \code{names(ma)}, else the file names without extension.
#' @param ld LD folder with \code{snp.info} and \code{block<b>.eigen.bin} (from
#'   \code{\link{LDbuild}} or SBayesRC). Without \code{ldscore.txt} the LD scores are computed
#'   from the eigen files during the first pass.
#' @param out One trait: output prefix, writes \code{<out>.snpRes}, \code{<out>.par.rds} and
#'   \code{<out>.log}. Several traits: output directory, writes
#'   \code{<out>/<trait>/<trait>_sbeigen.snpRes}, \code{.par.rds} and \code{.log} (trait folders are
#'   created when missing; next to SBayesRC's \code{<trait>_sbrc} files when \code{out} is the same
#'   folder). The log holds all messages and timings; with several traits the reading and pass
#'   timings are for the whole batch. \code{NULL} writes nothing.
#' @param threads Number of OpenMP threads.
#' @param ve Residual variance: a number (default 1) or \code{"ldsc"} for the LDSC intercept
#'   clamped to [0.9, 2]. Not used when \code{kappa = "mom"}, which estimates it.
#' @param kappa LD-mismatch noise: the residual variance of eigen component j is
#'   \eqn{ve_0 + \kappa / \lambda_j}. \code{"mom"} (default) estimates \eqn{ve_0} (within [0.5, 3])
#'   and \eqn{\kappa \ge 0} by moments on the components with \eqn{\lambda < 1} (100 bins, LDSC
#'   signal subtracted); a number fixes \eqn{\kappa} with \eqn{ve_0} = \code{ve}; 0 gives the
#'   constant residual variance \code{ve}.
#' @param thresh Proportion of eigenvalue mass kept per block. On ab.bin LD it truncates B at read time to the leading
#'   components reaching \code{thresh} of the Schur complement's positive eigenvalue mass (A is not affected);
#'   values at or above the build threshold keep every stored component.
#' @param threshB ab.bin only: B (eigen-space) effects are fitted only on the leading components reaching
#'   \code{threshB} of the Schur complement's positive eigenvalue mass; the remaining stored components keep
#'   \eqn{\alpha = 0} but stay in the data, so A and C are fitted in the whole \code{thresh} space (Yihe
#'   2026-10-10: e.g. \code{threshB = 0.95}, the rest left to C). \code{NULL} fits B on every component;
#'   \code{"auto"} picks it from 0.995/0.99/0.95/0.9 by pseudo-validation (90/10 split of the summary data in
#'   component space; a lower value only when its score is > 25\% higher), at about 5 times the fitting time.
#' @param tol VI stops when the posterior genetic variance changes by less than \code{tol}
#'   (relative) in two consecutive iterations (\code{method = "eigen"}; ABC uses 5e-4 together with
#'   a 0.05 z-unit change of the C effects).
#' @param method \code{"abc"} (default): \eqn{\beta = \beta_A + U\alpha + \gamma}, fitted jointly on the
#'   same eigen files. A: annotation SNPs with |z| > 4.5, thinned to r2 < 0.9 leads (largest |z| first),
#'   mr.ash in \eqn{\beta} space (exact 0 / 1 / 10); B: eigen-space mixture (exact 0 / 1 / 100 / 500);
#'   C: other SNPs with |z| > 4, MCP on the residual z scale. One iteration sweeps A, updates B and
#'   sweeps C on a shared residual. The U rows of the candidates are read in pass 1, so no extra pass.
#'   \code{"eigen"}: eigen-space VI only (six-class grid, LDSC prior on Vg).
#' @param annot Annotation SNPs for A (e.g. coding and xQTL): a character vector of SNP IDs or a file
#'   with a \code{SNP} column (else its first column). \code{NULL}: A is empty.
#' @param mcp MCP threshold \code{tau} and concavity \code{a} for C (Yihe 2026-10-07: 5 and 2.5). The
#'   threshold is in units of the residual z noise sd, \eqn{\tau\sqrt{ve_0}} under constant noise; with
#'   \code{kappa} the per-component noise \eqn{ve_0 + \kappa/\lambda} enters through the precision.
#' @param beta_ref \code{TRUE}: output beta = beta_std \eqn{\sqrt{V_y / 2p(1-p)}} with p the LD reference
#'   A1 frequency for COJO input as well (always so for the simplified format).
#' @return Invisibly. One trait: a list with \code{snpRes} (SNP, A1, A2, Block, beta, beta_std,
#'   every \code{snp.info} SNP; 0 in blocks without typed SNPs) and \code{par} (Vg, Vg_sd, ve, pi,
#'   sigma2, LDSC fit, timings, and \code{comp}: per eigen component its Block, lambda, n, w and
#'   posterior mean alpha, so \eqn{\beta' R \beta = \sum \lambda \alpha^2} for the B part). ABC adds
#'   \code{abc}: the A SNPs and C candidates (SNP, Block, set, z, beta_std; beta_std of a C SNP is 0 unless
#'   selected) and \code{abc_fit} (Vg_B, A mixture, counts); Vg is then that of the whole fit. On ab.bin LD the
#'   \code{comp} rows of the A components (eigenvalues of \eqn{R_{AA}}) have \code{alpha = NA}. Several traits: a list with \code{beta} (a SNP x trait matrix in
#'   \code{snp.info} order) and \code{par} (one entry per trait).
#' @examples
#' \dontrun{
#' fit <- sbayeseigen("trait.ma", "ukbEUR_LD", out = "trait", threads = 8)
#' fit$par$Vg
#' sbayeseigen(c(LDL = "ldl.ma", HDL = "hdl.ma"), "ukbEUR_LD", out = "prs", threads = 8)
#' }
#' @export
sbayeseigen <- function(ma, ld, out = NULL, threads = 4, ve = 1, kappa = "mom", thresh = 0.995, tol = 1e-4,
                        method = c("abc", "eigen"), annot = NULL, mcp = c(tau = 5, a = 2.5), beta_ref = FALSE, threshB = NULL) {
  method <- match.arg(method)
  if (!is.null(threshB) && !identical(threshB, "auto") && !(is.numeric(threshB) && length(threshB) == 1 && threshB > 0 && threshB <= 1))
    stop("threshB must be NULL, \"auto\" or in (0, 1]")
  if (!is.null(names(mcp)) && all(c("tau", "a") %in% names(mcp))) mcp <- mcp[c("tau", "a")]
  if (length(mcp) != 2 || mcp[[1]] <= 0 || mcp[[2]] <= 1) stop("mcp must be c(tau = > 0, a = > 1)")
  # messages still go to the console; with out they are also written to <prefix>.log per trait
  msg <- character(0)
  r <- withCallingHandlers(.sbayeseigen(ma, ld, out, threads, ve, kappa, thresh, tol, method, annot, mcp, beta_ref, threshB),
                           message = function(m) msg <<- c(msg, sub("\n$", "", conditionMessage(m))))
  for (p in attr(r, "prefix")) writeLines(c(format(Sys.time()), msg), paste0(p, ".log"))
  attr(r, "prefix") <- NULL
  invisible(r)
}

.sbayeseigen <- function(ma, ld, out, threads, ve, kappa, thresh, tol, method = "eigen", annot = NULL, mcp = c(5, 2.5),
                         beta_ref = FALSE, threshB = NULL) {
  if (!identical(kappa, "mom") && !(is.numeric(kappa) && length(kappa) == 1 && kappa >= 0))
    stop("kappa must be \"mom\" or a number >= 0")
  tm <- c(start = proc.time()[[3]])
  if (is.data.frame(ma)) ma <- list(ma)
  K <- length(ma)
  tn <- .trait_names(ma)
  si <- .read_snpinfo(ld)
  rows <- .block_rows(si)
  abc <- method == "abc"
  abld <- length(list.files(ld, "^block.*\\.ab\\.bin$")) > 0   # LDbuild with A: joint ABC on ab.bin
  if (abld && !abc) stop(ld, " holds ab.bin LD (LDbuild with A); use method = \"abc\"")
  if (abld && !is.null(annot)) message("annot ignored: A is the annotation set stored in the ab.bin files")
  isann <- if (abc && !abld) .read_annot(annot, si) else NULL

  # ---- summary data per trait (one full table in memory at a time) ----
  inp <- lapply(seq_len(K), function(t) {
    if (K > 1) message("---- ", tn[t])
    .trait_input(ma[[t]], ld, si, rows, beta_ref)
  })
  tm["read"] <- proc.time()[[3]]

  # ---- pass 1, all traits per block: impute, w = Lambda^{-1/2} U' bhat (and LD scores if needed) ----
  nty <- matrix(vapply(inp, `[[`, numeric(length(rows)), "nty"), ncol = K)
  run <- which(rowSums(nty > 0) > 0)
  files <- if (abld) file.path(ld, paste0("block", names(rows)[run], ".ab.bin")) else .eig_files(ld, names(rows)[run])
  # ABC: U rows of every trait's candidates (union over traits), read in pass 1
  if (abc) for (t in seq_len(K)) inp[[t]]$cand <- .abc_rows(inp[[t]])
  crow <- if (abc) lapply(run, function(j) sort(unique(unlist(lapply(inp, function(x) x$cand[[j]])))))
  lf <- file.path(ld, "ldscore.txt")
  if (abld && !file.exists(lf)) stop("ab.bin LD needs ldscore.txt (written by LDbuild)")
  if (!file.exists(lf)) message("ldscore.txt not found; LD scores computed from the eigen files in pass 1")
  p1 <- if (abld) ab_pass1_cpp(files, lapply(inp, function(x) x$ti[run]), lapply(inp, function(x) x$z[run]),
                               lapply(inp, function(x) x$n[run]), vapply(inp, `[[`, 0, "nmiss"), lapply(crow, as.integer), thresh, threads) else
        impute_blocks_eigen_cpp(files, lapply(inp, function(x) x$ti[run]), lapply(inp, function(x) x$z[run]),
                                lapply(inp, function(x) x$n[run]), vapply(inp, `[[`, 0, "nmiss"), thresh, threads,
                                TRUE, !file.exists(lf), if (abc) lapply(crow, as.integer) else list())
  lds <- if (file.exists(lf)) .ldscore_file(lf, si) else .ldscore_eig(p1, rows[run], si)
  tm["pass1"] <- proc.time()[[3]]

  # ---- per trait: LDSC on typed SNPs, VI ----
  alpha <- vector("list", K); par <- vector("list", K); add <- vector("list", K)
  for (t in seq_len(K)) {
    t0 <- proc.time()[[3]]
    x <- inp[[t]]
    use <- nty[run, t] > 0
    w <- unlist(lapply(p1[use], function(b) b$w[[t]]), use.names = FALSE)
    lam <- unlist(lapply(p1[use], `[[`, "lam"), use.names = FALSE)
    kb <- lengths(lapply(p1[use], `[[`, "lam"))
    ldsc <- ldsc_eigen(x$bh_ty, lds[x$ty], x$n_ty, M = nrow(si))
    h2p <- max(ldsc$h2, 0.01)
    nc <- rep(x$nbar[run][use], kb)
    # ab.bin: the first rankA components of a block are R_AA's (A); LD-mismatch noise is modelled on B only
    isB <- if (abld) unlist(lapply(p1[use], function(b) rep(c(FALSE, TRUE), c(b$rankA, length(b$lam) - b$rankA)))) else rep(TRUE, length(w))
    bad <- !is.finite(w) | !is.finite(nc) | !(lam > 0)
    if (any(bad)) stop(sprintf("%s%d of %d components have a non-finite w or N, or lambda <= 0 (blocks %s); check b, se, N (and r2) in the summary data",
                               if (K > 1) paste0(tn[t], ": ") else "", sum(bad), length(bad),
                               paste(utils::head(unique(rep(names(rows)[run][use], kb)[bad]), 10), collapse = ", ")))
    if (identical(kappa, "mom")) {
      nm <- noise_mom(w[isB], lam[isB], nc[isB], h2p); vt <- nm$ve0; kp <- nm$kappa
    } else {
      vt <- if (identical(ve, "ldsc")) min(max(ldsc$intercept, 0.9), 2) else as.numeric(ve); kp <- kappa
    }
    vj <- vt + kp / lam * isB   # residual variance per component; VI sees n / vj with ve = 1 (A components: ve0)
    # experimental, off by default (Yihe 2026-10-10, after SBayesRC's per-block ve and outlier removal). Under the model a
    # B component has E[n w^2] = vj + n lam tau0, tau0 = h2p / sum(lam_B) (Vg_B = LDSC h2), so signal shrinks with lam.
    # options(SBayesEigen.blockve = TRUE): fixed per-block noise multiplier s_b = 1 + max(0, mean_b(n w^2 / E) - 1)
    #   k_b / (k_b + 50) from the block's B components (only inflation), applied to all its components.
    # options(SBayesEigen.dropT = 30): B components with lam < 1 and n w^2 / E > dropT get precision ~0 (out of the fit;
    #   also the components beyond threshB, so C does not chase them).
    bb <- rep(seq_along(kb), kb); tau0 <- h2p / sum(lam[isB]); ndrop <- 0L; sbk <- NULL
    if (isTRUE(getOption("SBayesEigen.blockve", FALSE))) {
      e <- nc * w^2 / (vj + nc * lam * tau0); fb <- factor(bb[isB], levels = seq_along(kb))
      mb <- tapply(e[isB], fb, mean); kbB <- tabulate(bb[isB], length(kb)); mb[is.na(mb)] <- 1
      sbk <- 1 + pmax(0, mb - 1) * kbB / (kbB + 50); vj <- vj * sbk[bb]
      message(sprintf("%sper-block noise: %d of %d blocks inflated, median %.3f, max %.3f", if (K > 1) paste0(tn[t], ": ") else "",
                      sum(sbk > 1.001), length(sbk), stats::median(sbk), max(sbk)))
    }
    dT <- getOption("SBayesEigen.dropT", NULL)
    if (!is.null(dT)) {
      # options(SBayesEigen.dropLam = "median"): only components below their block's median B lambda (default lambda < 1)
      dl <- getOption("SBayesEigen.dropLam", 1)
      small <- if (identical(dl, "median")) lam < ave(ifelse(isB, lam, NA), bb, FUN = function(v) stats::median(v, na.rm = TRUE)) else lam < dl
      drop <- isB & small & nc * w^2 / (vj + nc * lam * tau0) > dT; ndrop <- sum(drop); vj[drop] <- vj[drop] * 1e8
      message(sprintf("%sdropped %d B components (lambda < %s, n w^2 / E > %g) in %d blocks", if (K > 1) paste0(tn[t], ": ") else "",
                      ndrop, if (identical(dl, "median")) "block median" else format(dl), dT, length(unique(bb[drop]))))
    }
    if (abld) {
      ab <- .abj_setup(x, p1[use], rows[run][use], t)
      pv <- NULL; thB <- threshB
      if (identical(threshB, "auto")) {
        pv <- .abj_pseudo_threshB(w, nc / vj, ab$blk, mcp, threads, h2p); thB <- pv$threshB
        message(sprintf("%spseudo-validation of threshB: %s; chosen %s", if (K > 1) paste0(tn[t], ": ") else "",
                        paste(sprintf("%g: %.4f", pv$grid, pv$score), collapse = ", "), if (is.null(thB)) "all" else format(thB)))
      }
      ab$blk <- .abj_nB(ab$blk, thB)
      fit <- abj_vi(w, nc / vj, ab$blk, tau = mcp[[1]], a = mcp[[2]], threads = threads, h2p = h2p)
      fit$threshB <- if (is.null(thB)) NA_real_ else thB; fit$pv <- pv
      if (!is.null(fit$sb)) message(sprintf("%sper-block noise (rho > 1.1): %d of %d blocks inflated, max %.2f", if (K > 1) paste0(tn[t], ": ") else "",
                                            sum(fit$sb > 1), length(fit$sb), max(fit$sb)))
      if (!fit$converged) warning("ABC did not converge in ", fit$iter, " iterations")
      fit$Vg_sd <- NA_real_; fit$gamma <- .abc_bprior(h2p, 1)$gB; fit$pi <- fit$piB; fit$sigma2 <- fit$s2B
      abc_tab <- data.table(SNP = si$SNP[c(unlist(ab$snpA), unlist(ab$snpC))],
                            Block = c(rep(names(rows)[run][use], lengths(ab$snpA)), rep(names(rows)[run][use], lengths(ab$snpC))),
                            set = rep(c("A", "C"), c(length(unlist(ab$snpA)), length(unlist(ab$snpC)))),
                            z = c(unlist(ab$zA), unlist(ab$zC)), beta_std = c(unlist(fit$mA), unlist(fit$gC)),
                            frac_trunc = c(rep(NA_real_, length(unlist(ab$snpA))), .abj_ctrunc(ab$blk, nc / vj)))
      fit$alpha <- unlist(lapply(seq_along(kb), function(j) c(rep(NA_real_, p1[use][[j]]$rankA), fit$alpha[[j]])))
      alB <- fit$alpha
    } else if (abc) {
      ab <- .abc_setup(x, p1[use], crow[use], rows[run][use], isann)
      fit <- abc_vi(w, sqrt(lam), nc / vj, kb, ab$X, ab$role, tau = mcp[[1]], a = mcp[[2]], threads = threads, h2p = h2p)
      if (!fit$converged) warning("ABC did not converge in ", fit$iter, " iterations")
      fit$Vg_sd <- NA_real_; fit$gamma <- .abc_bprior(h2p, 1)$gB; fit$pi <- fit$piB; fit$sigma2 <- fit$s2B
      cset <- unlist(ab$role) > 0
      abc_tab <- data.table(SNP = si$SNP[unlist(ab$snp)[cset]], Block = rep(names(rows)[run][use], lengths(ab$role))[cset],
                            set = c("A", "C")[unlist(ab$role)[cset]], z = unlist(ab$z)[cset], beta_std = unlist(fit$coef)[cset])
    } else {
      fit <- vi_eigen(w, sqrt(lam), nc / vj, ve = 1, h2p = h2p, threads = threads, tol = tol)
    }
    message(sprintf("%sLDSC on %d typed SNPs: h2 = %.4f, intercept = %.3f; ve0 = %.3f, kappa = %.4f%s; VI: %d iterations, Vg = %.4f (sd %.4f)",
                    if (K > 1) paste0(tn[t], ": ") else "", length(x$ty), ldsc$h2, ldsc$intercept, vt, kp,
                    if (identical(kappa, "mom")) " (MoM)" else "", fit$iter, fit$Vg, fit$Vg_sd))
    tl <- if (K > 1) paste0(tn[t], ": ") else ""
    if (identical(kappa, "mom") && vt < 1) message(sprintf("%snote: MoM ve0 = %.3f < 1 (N overstated or sample overlap?)", tl, vt))
    if (abld) message(sprintf("%sABC on ab.bin: %d A SNPs (stored annotation set), %d of %d C candidates selected (MCP tau = %g, a = %g); B fitted on %d of %d components; Pratt beta_k' bhat: A %.4f, B %.4f, C %.4f (shares %.3f, %.3f, %.3f)",
                              if (K > 1) paste0(tn[t], ": ") else "", fit$nA, fit$nC, fit$nC_cand, mcp[[1]], mcp[[2]],
                              fit$nB_fit, sum(isB), fit$Vg_A, fit$Vg_B, fit$Vg_C, fit$Vg_A / (fit$Vg_A + fit$Vg_B + fit$Vg_C),
                              fit$Vg_B / (fit$Vg_A + fit$Vg_B + fit$Vg_C), fit$Vg_C / (fit$Vg_A + fit$Vg_B + fit$Vg_C)))
    else if (abc) message(sprintf("%sABC: %d A SNPs (annotation, |z| > %g, r2 < %g leads), %d of %d C candidates selected (MCP tau = %g, a = %g); Vg_B = %.4f",
                             if (K > 1) paste0(tn[t], ": ") else "", fit$nA, .abc_const$zA, .abc_const$r2A, fit$nC, fit$nC_cand,
                             mcp[[1]], mcp[[2]], fit$Vg_B))
    alpha[[t]] <- rep(list(numeric(0)), length(run))
    alpha[[t]][use] <- if (abld) lapply(split(alB, rep(seq_along(kb), kb)), function(v) v[!is.na(v)]) else split(fit$alpha, rep(seq_along(kb), kb))
    par[[t]] <- list(Vg = fit$Vg, Vg_sd = fit$Vg_sd, ve = vt, kappa = kp, pi = fit$pi, sigma2 = fit$sigma2, gamma = fit$gamma,
                     iter = fit$iter, converged = fit$converged, ldsc = ldsc, h2_prior = h2p, n_drop = ndrop,
                     s_block = if (!is.null(fit$sb)) fit$sb else sbk,
                     n_snp = nrow(si), n_typed = length(x$ty), n_comp = length(w), n_block = sum(use),
                     time_fit = proc.time()[[3]] - t0,
                     # eigen-space fit: refit VI or get beta' R beta = sum(lam * alpha^2) without reading LD (eigen.bin; on ab.bin the A rows have alpha = NA)
                     comp = data.table(Block = rep(names(rows)[run][use], kb), lam = lam, n = nc,
                                       ve = vj, w = w, alpha = fit$alpha))
    if (abc) {
      par[[t]]$abc <- abc_tab
      par[[t]]$abc_fit <- fit[intersect(c("Vg_A", "Vg_B", "Vg_C", "Vg_parts", "piA", "s2A", "nA", "nC_cand", "nC", "nB_fit", "threshB", "pv"), names(fit))]
      add[[t]] <- abc_tab[beta_std != 0]
    }
  }
  tm["fit"] <- proc.time()[[3]]

  # ---- pass 2, all traits per block: beta = U E[alpha] ----
  p2 <- if (abld) ab_beta_cpp(files, alpha, thresh, threads) else eig_beta_cpp(files, alpha, thresh, threads)
  tm["pass2"] <- proc.time()[[3]]
  message("Time (s", if (K > 1) paste0(", whole batch of ", K, " traits"), "): ",
          paste(names(diff(tm)), sprintf("%.1f", diff(tm)), sep = " ", collapse = ", "),
          sprintf("; total %.1f", tm[["pass2"]] - tm[["start"]]))
  rr <- rows[run]
  bstd <- function(t) {
    b <- numeric(nrow(si))
    use <- nty[run, t] > 0
    pb <- p2[[t]][use]
    e <- lengths(pb) == 0   # ab.bin block with no B components (all SNPs in A): no U alpha part
    pb[e] <- lapply(rr[use][e], function(r) numeric(length(r)))
    stopifnot(lengths(pb) == lengths(rr[use]))
    b[unlist(rr[use], use.names = FALSE)] <- unlist(pb, use.names = FALSE)
    if (!is.null(add[[t]])) { i <- match(add[[t]]$SNP, si$SNP); b[i] <- b[i] + add[[t]]$beta_std }   # beta_A and gamma
    b
  }
  if (K == 1) {
    beta_std <- bstd(1)
    snpRes <- data.table(SNP = si$SNP, A1 = si$A1, A2 = si$A2, Block = si$Block,
                         beta = beta_std * inp[[1]]$scale, beta_std = beta_std)
    par[[1]]$time <- diff(tm)
    if (!is.null(out)) .write_fit(snpRes, par[[1]], out)
    return(structure(list(snpRes = snpRes, par = par[[1]]), prefix = out))
  }
  beta <- matrix(0, nrow(si), K, dimnames = list(NULL, tn))
  if (!is.null(out)) dir.create(out, showWarnings = FALSE, recursive = TRUE)
  for (t in seq_len(K)) {
    beta_std <- bstd(t)
    beta[, t] <- beta_std * inp[[t]]$scale
    par[[t]]$time <- diff(tm)
    if (!is.null(out)) {
      d <- file.path(out, tn[t])
      if (!dir.exists(d)) dir.create(d)
      .write_fit(data.table(SNP = si$SNP, A1 = si$A1, A2 = si$A2, Block = si$Block, beta = beta[, t],
                            beta_std = beta_std), par[[t]], file.path(d, paste0(tn[t], "_sbeigen")))
    }
  }
  names(par) <- tn
  structure(list(beta = beta, par = par), prefix = if (!is.null(out)) file.path(out, tn, paste0(tn, "_sbeigen")))
}

.write_fit <- function(snpRes, par, out) {
  fwrite(snpRes, paste0(out, ".snpRes"), sep = "\t")
  saveRDS(par, paste0(out, ".par.rds"))
  message("wrote ", out, ".snpRes, .par.rds and .log")
}

.trait_names <- function(ma) {
  tn <- names(ma)
  if (is.null(tn)) tn <- character(length(ma))
  if (is.character(ma)) {
    fn <- sub("\\.(gz|bz2|zip)$", "", basename(ma))
    tn[tn == ""] <- sub("\\.[^.]*$", "", fn)[tn == ""]
  }
  tn[tn == ""] <- paste0("trait", seq_along(ma))[tn == ""]
  if (anyDuplicated(tn)) stop("trait names must be unique: ", paste(unique(tn[duplicated(tn)]), collapse = ", "))
  tn
}

# one trait, slim: per block typed index (0-based), z, N; scale = sqrt(N se^2 + b^2) on every SNP (imputed
# SNPs: sqrt(var_y / 2pq), as impute() writes them); typed SNPs for LDSC; mean typed N per block
.trait_input <- function(x, ld, si, rows, beta_ref = FALSE) {
  if (is.character(x) && "r2" %in% names(fread(x, nrows = 0, showProgress = FALSE))) x <- fread(x, showProgress = FALSE)
  if (is.data.frame(x) && .is_imputed(x, si)) {
    message("Summary data already imputed")
    if (!identical(as.character(x$SNP), si$SNP) || !identical(x$A1, si$A1))
      stop("imputed summary data are not in snp.info order and allele coding; run impute() on the tidied data")
    # no second imputation: every row is used as given; rows with a non-finite b/se/N (or se <= 0) get z = 0
    b <- x$b; se <- x$se; N <- x$N; f <- x$freq
    ok <- is.finite(b) & is.finite(se) & se > 0 & is.finite(N) & N > 0
    if (!all(ok)) message(sum(!ok), " rows with a non-finite b, se or N (or se <= 0) set to z = 0")
    nmiss <- stats::median(N[ok])
    vp <- stats::median((2 * f * (1 - f) * (N * se^2 + b^2))[ok])
    obs <- rep(TRUE, length(b)); typed <- ok & !is.na(x$r2) & x$r2 == 1
    scale <- ifelse(ok, sqrt(N * se^2 + b^2), sqrt(vp / (2 * f * (1 - f))))
    b[!ok] <- 0; se[!ok] <- 1; N[!ok] <- nmiss
    if (beta_ref) scale <- sqrt(vp / (2 * si$freq * (1 - si$freq)))
  } else {
    zin <- .is_zinput(x)
    a <- .align(.tidy(x, ld, si = si, idx = TRUE), si)
    b <- a$res$b; se <- a$res$se; N <- a$res$N; f <- a$res$freq
    obs <- typed <- is.finite(b) & is.finite(se) & se > 0 & is.finite(N) & N > 0
    scale <- ifelse(obs, sqrt(N * se^2 + b^2), sqrt(a$vp / (2 * f * (1 - f))))
    # Z input (Var_y = 1) or beta_ref: per-dosage effect with the LD reference freq
    if (zin || beta_ref) scale <- sqrt(a$vp / (2 * si$freq * (1 - si$freq)))
    nmiss <- a$Nmed
  }
  z <- b / se
  ty <- which(typed)
  list(ti = lapply(rows, function(r) which(obs[r]) - 1L), z = lapply(rows, function(r) z[r][obs[r]]),
       n = lapply(rows, function(r) N[r][obs[r]]), typ = lapply(rows, function(r) typed[r][obs[r]]), nmiss = nmiss, nty = vapply(rows, function(r) sum(obs[r]), 0),
       nbar = vapply(rows, function(r) if (any(typed[r])) mean(N[r][typed[r]]) else if (any(obs[r])) mean(N[r][obs[r]]) else nmiss, 0),
       scale = scale,
       ty = ty, bh_ty = b[ty] / scale[ty], n_ty = N[ty])
}

.eig_files <- function(ld, blocks) file.path(ld, paste0("block", blocks, ".eigen.bin"))

.ldscore_file <- function(f, si) {
  l <- fread(f, showProgress = FALSE)
  if (!"ldscore" %in% names(l)) stop(f, " has no ldscore column")
  v <- l$ldscore[match(si$SNP, l$SNP)]
  if (anyNA(v)) stop(f, " does not cover every snp.info SNP")
  v
}

# unbiased LD score from the eigen LD score of pass 1 (blocks in rr only; NA elsewhere):
# sum_j r2_adj = (l_eig (n - 1) - m_block) / (n - 2), n the LD reference sample size
.ldscore_eig <- function(p1, rr, si) {
  v <- rep(NA_real_, nrow(si))
  i <- unlist(rr, use.names = FALSE)
  n <- si$N[i]
  v[i] <- (unlist(lapply(p1, `[[`, "ld"), use.names = FALSE) * (n - 1) - rep(lengths(rr), lengths(rr))) / (n - 2)
  v
}

# annotation (fixed A set): SNP IDs, or a file whose SNP column (else first column) lists them
.read_annot <- function(annot, si) {
  if (is.null(annot)) return(logical(nrow(si)))
  if (length(annot) == 1 && file.exists(annot)) {
    a <- fread(annot, header = FALSE, showProgress = FALSE)
    annot <- if (any(a[1] == "SNP")) a[[which(unlist(a[1]) == "SNP")[1]]][-1] else a[[1]]
  }
  isa <- si$SNP %in% annot
  message(sprintf("Annotation: %d of %d snp.info SNPs", sum(isa), nrow(si)))
  isa
}

# one trait's ABC design in the blocks it uses: p1b pass-1 results, cr the union candidate rows (0-based), rr the
# blocks' snp.info rows. X: k x s candidate columns sqrt(lambda) U[i, ], role, snp (snp.info index), z per column
.abc_setup <- function(x, p1b, cr, rr, isann) {
  bi <- match(names(rr), names(x$ti))
  out <- list(X = list(), role = list(), snp = list(), z = list())
  for (j in seq_along(p1b)) {
    tr <- x$cand[[bi[j]]]; q <- match(tr, cr[[j]])
    zz <- x$z[[bi[j]]][match(tr, x$ti[[bi[j]]])]
    lam <- p1b[[j]]$lam
    U <- if (length(q)) p1b[[j]]$urow[q, , drop = FALSE] else matrix(0, 0, length(lam))
    role <- .abc_roles(U, lam, zz, isann[rr[[j]][tr + 1L]])
    out$X[[j]] <- t(U) * sqrt(lam); out$role[[j]] <- as.integer(role)
    out$snp[[j]] <- rr[[j]][tr + 1L]; out$z[[j]] <- zz
  }
  out
}

# one trait's ab.bin ABC design in the blocks it uses: the trait's C candidates among the union rows of pass 1
# C diagnostic: share of each C candidate's precision-weighted signal (its D) on the B components beyond threshB
# (alpha fixed at 0 there), i.e. how much of what C can fit lies in the truncated directions
.abj_ctrunc <- function(blk, p) {
  rk <- vapply(blk, function(x) nrow(x$XA), 0L); kB <- lengths(lapply(blk, `[[`, "sl")); off <- c(0, cumsum(rk + kB))
  unlist(lapply(seq_along(blk), function(b) { x <- blk[[b]]; if (!ncol(x$XCb)) return(numeric(0))
    pa <- p[off[b] + seq_len(rk[b])]; pb <- p[off[b] + rk[b] + seq_len(kB[b])]; nB <- if (is.null(x$nB)) kB[b] else x$nB
    e <- x$XCb^2 * pb; D <- colSums(e) + (if (rk[b]) colSums(x$XCa^2 * pa) else 0)
    (if (nB < kB[b]) colSums(e[-seq_len(nB), , drop = FALSE]) else 0) / D }))
}

# B fitted on the leading components reaching th of the block's Schur-complement mass (NULL: all)
.abj_nB <- function(blk, th) lapply(blk, function(x) { lb <- x$sl^2
  x$nB <- if (is.null(th)) NULL else min(length(lb), which(cumsum(lb) >= th * x$sumL)[1], na.rm = TRUE); x })

# threshB = "auto" (benchmark thread 2026-10-10, after SBayesRC's pseudo-validation): split w into a 90% training and
# a 10% validation sample in component space (w_t = w + sqrt((1/0.9 - 1) / p) e, w_v = (w - 0.9 w_t) / 0.1), fit on
# w_t for each threshB, score r = fit' w_v / ||fit|| (= beta' bhat_v / sqrt(beta' R beta)); leave 0.995 (all) only
# when the best is more than 25% higher. Fixed seed, so a rerun gives the same choice.
.abj_pseudo_threshB <- function(w, p, blk, mcp, threads, h2p, grid = c(0.995, 0.99, 0.95, 0.9), gain = 1.25) {
  if (exists(".Random.seed", globalenv())) { rs <- get(".Random.seed", globalenv()); on.exit(assign(".Random.seed", rs, globalenv())) }
  set.seed(20261010); wt <- w + sqrt((1 / 0.9 - 1) / p) * stats::rnorm(length(w)); wv <- (w - 0.9 * wt) / 0.1
  sc <- vapply(grid, function(th) {
    f <- abj_vi(wt, 0.9 * p, .abj_nB(blk, th), tau = mcp[[1]], a = mcp[[2]], threads = threads, h2p = h2p)
    sum(f$fitw * wv) / sqrt(sum(f$fitw^2)) }, 0)
  b <- which.max(sc); better <- if (sc[1] > 0) sc[b] > gain * sc[1] else b != 1
  list(grid = grid, score = sc, threshB = if (better) grid[b] else NULL)
}

.abj_setup <- function(x, p1b, rr, t) {
  bi <- match(names(rr), names(x$ti))
  out <- list(blk = list(), snpA = list(), snpC = list(), zA = list(), zC = list())
  for (j in seq_along(p1b)) {
    pb <- p1b[[j]]; ti <- x$ti[[bi[j]]] + 1L; zz <- x$z[[bi[j]]]
    q <- which(pb$crow %in% (x$cand[[bi[j]]] + 1L))
    lb <- pb$lam[pb$rankA + seq_len(ncol(pb$Cm))]
    # Yihe 2026-10-10: B fitted only up to threshB of the Schur complement's mass (sumL; .abj_nB); A and C use the whole stored space
    out$blk[[j]] <- list(XA = pb$XA, Cm = pb$Cm, sl = sqrt(lb), XCa = pb$XCa[, q, drop = FALSE], XCb = pb$XCb[, q, drop = FALSE],
                         sumL = pb$sumLambdaB)
    out$snpA[[j]] <- rr[[j]][pb$iA]; out$zA[[j]] <- zz[match(pb$iA, ti)]
    out$snpC[[j]] <- rr[[j]][pb$crow[q]]; out$zC[[j]] <- zz[match(pb$crow[q], ti)]
  }
  out
}
