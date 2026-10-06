#' Genome-wide PRS from eigen LD and GWAS summary statistics
#'
#' All-in-one: \code{\link{tidy}} and \code{\link{impute}} in memory, LD score regression for
#' the heritability prior, variational EM in the eigen space and per-SNP effects. Each
#' \code{block<b>.eigen.bin} is read once for imputation and the rotation
#' \eqn{w = \Lambda^{-1/2} U' \hat b} (\eqn{\hat b = z / \sqrt{N + z^2}}, as SBayesRC), and once
#' more for \eqn{\beta = U E[\alpha]}. Summary data that are already imputed (an \code{r2}
#' column and one row per \code{snp.info} SNP) skip tidy and impute.
#'
#' @param ma GWAS summary statistics in COJO format (SNP A1 A2 freq b se p N): a path or a
#'   data.frame, raw or already imputed.
#' @param ld LD folder with \code{snp.info} and \code{block<b>.eigen.bin} (from
#'   \code{\link{LDbuild}} or SBayesRC). Without \code{ldscore.txt} the LD scores are computed
#'   from the eigen files (and saved there when the folder is writable).
#' @param out Output prefix; writes \code{<out>.snpRes} and \code{<out>.par.rds}.
#'   \code{NULL} writes nothing.
#' @param threads Number of OpenMP threads.
#' @param ve Residual variance: a number (default 1) or \code{"ldsc"} for the LDSC intercept
#'   clamped to [0.9, 2].
#' @param thresh Proportion of eigenvalue mass kept per block.
#' @return Invisibly, a list with \code{snpRes} (SNP, A1, A2, Block, beta, beta_std, every
#'   \code{snp.info} SNP; 0 in blocks without typed SNPs) and \code{par} (Vg, Vg_sd, ve, pi,
#'   sigma2, LDSC fit, timings).
#' @examples
#' \dontrun{
#' fit <- sbayeseigen("trait.ma", "ukbEUR_LD", out = "trait", threads = 8)
#' fit$par$Vg
#' }
#' @export
sbayeseigen <- function(ma, ld, out = NULL, threads = 4, ve = 1, thresh = 0.995) {
  tm <- c(start = proc.time()[[3]])
  si <- .read_snpinfo(ld)

  # ---- summary data -> bhat (all SNPs) and pass 1 w = Lambda^{-1/2} U' bhat ----
  if (is.character(ma) && "r2" %in% names(fread(ma, nrows = 0, showProgress = FALSE))) ma <- fread(ma, showProgress = FALSE)
  if (is.data.frame(ma) && "r2" %in% names(ma) && nrow(ma) == nrow(si)) {
    message("Summary data already imputed")
    if (!identical(as.character(ma$SNP), si$SNP) || !identical(ma$A1, si$A1))
      stop("imputed summary data are not in snp.info order and allele coding; run impute() on the tidied data")
    gw <- ma
    bh <- gw$b / sqrt(gw$N * gw$se^2 + gw$b^2)
    rows <- split(seq_len(nrow(si)), factor(si$Block, levels = unique(si$Block)))
    run <- names(rows)[vapply(rows, function(r) any(gw$r2[r] >= 0), TRUE)]
    tm["impute"] <- proc.time()[[3]]
    blk <- eig_w_cpp(.eig_files(ld, run), lapply(rows[run], function(r) bh[r]), thresh, threads)
  } else {
    gw <- .tidy(ma, ld, si = si)
    tm["tidy"] <- proc.time()[[3]]
    r <- .impute(gw, ld, si, thresh, threads, return_w = TRUE)
    gw <- r$ma; run <- r$run; blk <- r$blocks
    rows <- split(seq_len(nrow(si)), factor(si$Block, levels = unique(si$Block)))
    tm["impute"] <- proc.time()[[3]]
  }
  scale <- sqrt(gw$N * gw$se^2 + gw$b^2)
  kb <- lengths(lapply(blk, `[[`, "w"))
  w <- unlist(lapply(blk, `[[`, "w"), use.names = FALSE)
  lam <- unlist(lapply(blk, `[[`, "lam"), use.names = FALSE)
  nn <- rep(vapply(rows[run], function(r) mean(gw$N[r][gw$r2[r] == 1]), 0), kb)   # typed SNPs only
  tm["pass1"] <- proc.time()[[3]]

  # ---- LDSC on typed SNPs: slope -> h2 prior centre, intercept -> optional ve ----
  lds <- .ldscore(ld, si)
  ty <- which(gw$r2 == 1)
  ldsc <- ldsc_eigen(gw$b[ty] / scale[ty], lds[ty], gw$N[ty], M = nrow(si))
  ve <- if (identical(ve, "ldsc")) min(max(ldsc$intercept, 0.9), 2) else as.numeric(ve)
  h2p <- max(ldsc$h2, 0.01)
  message(sprintf("LDSC on %d typed SNPs: h2 = %.4f, intercept = %.3f; ve = %.3f", length(ty), ldsc$h2, ldsc$intercept, ve))
  tm["ldsc"] <- proc.time()[[3]]

  # ---- VI (hyperparameters shared genome-wide) ----
  fit <- vi_eigen(w, sqrt(lam), nn, ve = ve, h2p = h2p, threads = threads)
  message(sprintf("VI: %d iterations, Vg = %.4f (sd %.4f)", fit$iter, fit$Vg, fit$Vg_sd))
  tm["vi"] <- proc.time()[[3]]

  # ---- pass 2: beta = U E[alpha] in the same basis ----
  beta_std <- numeric(nrow(si))
  b2 <- eig_beta_cpp(.eig_files(ld, run), split(fit$alpha, rep(seq_along(kb), kb)), thresh, threads)
  beta_std[unlist(rows[run], use.names = FALSE)] <- unlist(b2, use.names = FALSE)
  tm["pass2"] <- proc.time()[[3]]

  snpRes <- data.table(SNP = si$SNP, A1 = si$A1, A2 = si$A2, Block = si$Block,
                       beta = beta_std * scale, beta_std = beta_std)
  snpRes[!is.finite(beta), beta := 0]
  par <- list(Vg = fit$Vg, Vg_sd = fit$Vg_sd, ve = ve, pi = fit$pi, sigma2 = fit$sigma2, gamma = fit$gamma,
              iter = fit$iter, converged = fit$converged, ldsc = ldsc, h2_prior = h2p,
              n_snp = nrow(si), n_typed = length(ty), n_comp = length(w), n_block = length(run),
              time = diff(tm))
  if (!is.null(out)) {
    fwrite(snpRes, paste0(out, ".snpRes"), sep = "\t")
    saveRDS(par, paste0(out, ".par.rds"))
    message("wrote ", out, ".snpRes and ", out, ".par.rds")
  }
  invisible(list(snpRes = snpRes, par = par))
}

.eig_files <- function(ld, blocks) file.path(ld, paste0("block", blocks, ".eigen.bin"))

# unbiased LD score per snp.info SNP: ldscore.txt (LDbuild) or, if absent, from the eigen files,
# sum_j r2_adj = (l_eig (n - 1) - m_block) / (n - 2) with n the LD reference sample size
.ldscore <- function(ld, si) {
  f <- file.path(ld, "ldscore.txt")
  if (file.exists(f)) {
    l <- fread(f, showProgress = FALSE)
    if (!"ldscore" %in% names(l)) stop(f, " has no ldscore column")
    v <- l$ldscore[match(si$SNP, l$SNP)]
    if (anyNA(v)) stop(f, " does not cover every snp.info SNP")
    return(v)
  }
  message("ldscore.txt not found; computing LD scores from the eigen files")
  blocks <- unique(si$Block)
  le <- eig_ldscore_cpp(.eig_files(ld, blocks))
  mb <- rep(lengths(le), lengths(le))
  n <- si$N[order(match(si$Block, blocks))]
  v <- (unlist(le, use.names = FALSE) * (n - 1) - mb) / (n - 2)
  v <- v[order(order(match(si$Block, blocks)))]
  if (file.access(ld, 2) == 0)
    fwrite(data.table(SNP = si$SNP, Block = si$Block, ldscore = v), f, sep = "\t")
  v
}
