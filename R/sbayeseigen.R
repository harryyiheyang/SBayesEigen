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
#' @param ma GWAS summary statistics in COJO format (SNP A1 A2 freq b se p N), raw or already
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
#'   clamped to [0.9, 2].
#' @param thresh Proportion of eigenvalue mass kept per block.
#' @param tol VI stops when the posterior genetic variance changes by less than \code{tol}
#'   (relative) in two consecutive iterations.
#' @return Invisibly. One trait: a list with \code{snpRes} (SNP, A1, A2, Block, beta, beta_std,
#'   every \code{snp.info} SNP; 0 in blocks without typed SNPs) and \code{par} (Vg, Vg_sd, ve, pi,
#'   sigma2, LDSC fit, timings, and \code{comp}: per eigen component its Block, lambda, n, w and
#'   posterior mean alpha, so \eqn{\beta' R \beta = \sum \lambda \alpha^2}). Several traits: a list with \code{beta} (a SNP x trait matrix in
#'   \code{snp.info} order) and \code{par} (one entry per trait).
#' @examples
#' \dontrun{
#' fit <- sbayeseigen("trait.ma", "ukbEUR_LD", out = "trait", threads = 8)
#' fit$par$Vg
#' sbayeseigen(c(LDL = "ldl.ma", HDL = "hdl.ma"), "ukbEUR_LD", out = "prs", threads = 8)
#' }
#' @export
sbayeseigen <- function(ma, ld, out = NULL, threads = 4, ve = 1, thresh = 0.995, tol = 1e-4) {
  # messages still go to the console; with out they are also written to <prefix>.log per trait
  msg <- character(0)
  r <- withCallingHandlers(.sbayeseigen(ma, ld, out, threads, ve, thresh, tol),
                           message = function(m) msg <<- c(msg, sub("\n$", "", conditionMessage(m))))
  for (p in attr(r, "prefix")) writeLines(c(format(Sys.time()), msg), paste0(p, ".log"))
  attr(r, "prefix") <- NULL
  invisible(r)
}

.sbayeseigen <- function(ma, ld, out, threads, ve, thresh, tol) {
  tm <- c(start = proc.time()[[3]])
  if (is.data.frame(ma)) ma <- list(ma)
  K <- length(ma)
  tn <- .trait_names(ma)
  si <- .read_snpinfo(ld)
  rows <- .block_rows(si)

  # ---- summary data per trait (one full table in memory at a time) ----
  inp <- lapply(seq_len(K), function(t) {
    if (K > 1) message("---- ", tn[t])
    .trait_input(ma[[t]], ld, si, rows)
  })
  tm["read"] <- proc.time()[[3]]

  # ---- pass 1, all traits per block: impute, w = Lambda^{-1/2} U' bhat (and LD scores if needed) ----
  nty <- matrix(vapply(inp, `[[`, numeric(length(rows)), "nty"), ncol = K)
  run <- which(rowSums(nty > 0) > 0)
  files <- .eig_files(ld, names(rows)[run])
  lf <- file.path(ld, "ldscore.txt")
  if (!file.exists(lf)) message("ldscore.txt not found; LD scores computed from the eigen files in pass 1")
  p1 <- impute_blocks_eigen_cpp(files, lapply(inp, function(x) x$ti[run]), lapply(inp, function(x) x$z[run]),
                                lapply(inp, function(x) x$n[run]), vapply(inp, `[[`, 0, "nmiss"), thresh, threads,
                                TRUE, !file.exists(lf))
  lds <- if (file.exists(lf)) .ldscore_file(lf, si) else .ldscore_eig(p1, rows[run], si)
  tm["pass1"] <- proc.time()[[3]]

  # ---- per trait: LDSC on typed SNPs, VI ----
  alpha <- vector("list", K); par <- vector("list", K)
  for (t in seq_len(K)) {
    t0 <- proc.time()[[3]]
    x <- inp[[t]]
    use <- nty[run, t] > 0
    w <- unlist(lapply(p1[use], function(b) b$w[[t]]), use.names = FALSE)
    lam <- unlist(lapply(p1[use], `[[`, "lam"), use.names = FALSE)
    kb <- lengths(lapply(p1[use], `[[`, "lam"))
    ldsc <- ldsc_eigen(x$bh_ty, lds[x$ty], x$n_ty, M = nrow(si))
    vt <- if (identical(ve, "ldsc")) min(max(ldsc$intercept, 0.9), 2) else as.numeric(ve)
    h2p <- max(ldsc$h2, 0.01)
    fit <- vi_eigen(w, sqrt(lam), rep(x$nbar[run][use], kb), ve = vt, h2p = h2p, threads = threads, tol = tol)
    message(sprintf("%sLDSC on %d typed SNPs: h2 = %.4f, intercept = %.3f; ve = %.3f; VI: %d iterations, Vg = %.4f (sd %.4f)",
                    if (K > 1) paste0(tn[t], ": ") else "", length(x$ty), ldsc$h2, ldsc$intercept, vt, fit$iter,
                    fit$Vg, fit$Vg_sd))
    alpha[[t]] <- rep(list(numeric(0)), length(run))
    alpha[[t]][use] <- split(fit$alpha, rep(seq_along(kb), kb))
    par[[t]] <- list(Vg = fit$Vg, Vg_sd = fit$Vg_sd, ve = vt, pi = fit$pi, sigma2 = fit$sigma2, gamma = fit$gamma,
                     iter = fit$iter, converged = fit$converged, ldsc = ldsc, h2_prior = h2p,
                     n_snp = nrow(si), n_typed = length(x$ty), n_comp = length(w), n_block = sum(use),
                     time_fit = proc.time()[[3]] - t0,
                     # eigen-space fit: refit VI or get beta' R beta = sum(lam * alpha^2) without reading LD
                     comp = data.table(Block = rep(names(rows)[run][use], kb), lam = lam, n = rep(x$nbar[run][use], kb),
                                       w = w, alpha = fit$alpha))
  }
  tm["fit"] <- proc.time()[[3]]

  # ---- pass 2, all traits per block: beta = U E[alpha] ----
  p2 <- eig_beta_cpp(files, alpha, thresh, threads)
  tm["pass2"] <- proc.time()[[3]]
  message("Time (s", if (K > 1) paste0(", whole batch of ", K, " traits"), "): ",
          paste(names(diff(tm)), sprintf("%.1f", diff(tm)), sep = " ", collapse = ", "),
          sprintf("; total %.1f", tm[["pass2"]] - tm[["start"]]))
  rr <- rows[run]
  bstd <- function(t) {
    b <- numeric(nrow(si))
    use <- nty[run, t] > 0
    b[unlist(rr[use], use.names = FALSE)] <- unlist(p2[[t]][use], use.names = FALSE)
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
.trait_input <- function(x, ld, si, rows) {
  if (is.character(x) && "r2" %in% names(fread(x, nrows = 0, showProgress = FALSE))) x <- fread(x, showProgress = FALSE)
  if (is.data.frame(x) && .is_imputed(x, si)) {
    message("Summary data already imputed")
    if (!identical(as.character(x$SNP), si$SNP) || !identical(x$A1, si$A1))
      stop("imputed summary data are not in snp.info order and allele coding; run impute() on the tidied data")
    obs <- x$r2 >= 0; typed <- x$r2 == 1
    b <- x$b; se <- x$se; N <- x$N
    scale <- sqrt(N * se^2 + b^2)
    nmiss <- stats::median(N)
  } else {
    a <- .align(.tidy(x, ld, si = si), si)
    b <- a$res$b; se <- a$res$se; N <- a$res$N; f <- a$res$freq
    obs <- typed <- is.finite(b)
    scale <- ifelse(obs, sqrt(N * se^2 + b^2), sqrt(a$vp / (2 * f * (1 - f))))
    nmiss <- a$Nmed
  }
  z <- b / se
  ty <- which(typed)
  list(ti = lapply(rows, function(r) which(obs[r]) - 1L), z = lapply(rows, function(r) z[r][obs[r]]),
       n = lapply(rows, function(r) N[r][obs[r]]), nmiss = nmiss, nty = vapply(rows, function(r) sum(obs[r]), 0),
       nbar = vapply(rows, function(r) mean(N[r][typed[r]]), 0), scale = scale,
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
