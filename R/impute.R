#' Impute missing summary statistics from eigen LD
#'
#' Yihe Yang's block-parallel imputation working directly in the eigen space of
#' \code{block*.eigen.bin}: per block the smallest exact solve among the missing-SNP,
#' eigen-space (Woodbury) and typed-SNP systems, in float32, with
#' \eqn{R_{OO} = U_O \Lambda U_O' + 0.1 I}. Typed SNPs keep their values (\code{r2 = 1}),
#' imputed SNPs get \code{r2 = 0}, SNPs in blocks without typed SNPs get \code{b = 0},
#' \code{se = 1}, \code{r2 = -1}. Input that is already imputed (an \code{r2} column and one
#' row per \code{snp.info} SNP) is returned unchanged.
#'
#' @param ma Summary data (SNP A1 A2 freq b se p N): a path or a data.frame, e.g. the output
#'   of \code{\link{tidy}}.
#' @param ld LD folder with \code{snp.info} and \code{block*.eigen.bin}; the leading eigen
#'   components holding 99.5\% of the eigenvalue mass are used, as in \code{\link{sbayeseigen}}.
#' @param out Output path (may equal \code{ma}; replaced atomically); \code{NULL} writes nothing.
#' @param threads Blocks imputed in parallel.
#' @return Invisibly, the summary data for every SNP in \code{snp.info} order, with an
#'   \code{r2} column.
#' @export
impute <- function(ma, ld, out = NULL, threads = 4) {
  r <- .impute(ma, ld, .read_snpinfo(ld), thresh = 0.995, threads = threads)
  if (!is.null(out) && !r$done) {
    tmp <- tempfile(pattern = paste0(".", basename(out), ".tmp-"), tmpdir = dirname(out))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    fwrite(r$ma, tmp, sep = "\t", quote = FALSE, na = "NA")
    if (!file.rename(tmp, out)) stop("imputation finished, but replacing ", out, " failed")
  }
  invisible(r$ma)
}

# si: .std_snpinfo table. thresh: eigen cut (0 = all stored components). return_w: also return
# pass 1 (w, lambda) for every block with a typed SNP. Returns list(ma, done, blocks, run).
.impute <- function(ma, ld, si, thresh, threads, return_w = FALSE) {
  t0 <- proc.time()[[3]]
  ma <- if (is.data.frame(ma)) as.data.table(ma) else fread(ma, showProgress = FALSE)
  cols <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  if (!all(cols %in% names(ma))) stop("missing columns in the summary data: ", paste(setdiff(cols, names(ma)), collapse = ", "))
  if ("r2" %in% names(ma) && nrow(ma) == nrow(si)) {
    message("Already imputed: r2 exists and the summary data has one row per snp.info SNP")
    return(list(ma = ma, done = TRUE))
  }
  vp <- stats::median(2 * ma$freq * (1 - ma$freq) * (ma$N * ma$se^2 + ma$b^2), na.rm = TRUE)
  Nmed <- stats::median(ma$N, na.rm = TRUE)
  if (!is.finite(vp) || !is.finite(Nmed)) stop("cannot compute a finite median var(y) or N from the summary data")

  # align to snp.info
  m <- match(ma$SNP, si$SNP)
  i <- which(!is.na(m)); j <- m[i]
  same <- ma$A1[i] == si$A1[j] & ma$A2[i] == si$A2[j]
  flip <- ma$A1[i] == si$A2[j] & ma$A2[i] == si$A1[j]
  message(length(i), " SNPs in common between the summary data and LD; ", sum(same), " as is, ",
          sum(flip), " with flipped alleles")
  res <- data.table(SNP = si$SNP, A1 = si$A1, A2 = si$A2, freq = si$freq, b = NA_real_, se = NA_real_,
                    p = NA_real_, N = NA_real_)
  k <- i[same | flip]; s <- ifelse(flip[same | flip], -1, 1)
  set(res, j[same | flip], c("freq", "b", "se", "p", "N"),
      list(ifelse(s < 0, 1 - ma$freq[k], ma$freq[k]), s * ma$b[k], ma$se[k], ma$p[k], ma$N[k]))
  res[is.na(N), N := Nmed]

  # blocks: typed index (0-based) and z per block
  typed <- is.finite(res$b)
  rows <- split(seq_len(nrow(si)), factor(si$Block, levels = unique(si$Block)))
  nty <- vapply(rows, function(r) sum(typed[r]), 0)
  run <- if (return_w) nty > 0 else nty > 0 & nty < lengths(rows)
  tr <- lapply(rows[run], function(r) r[typed[r]])
  out <- impute_blocks_eigen_cpp(file.path(ld, paste0("block", names(rows)[run], ".eigen.bin")),
                                 Map(function(r, t) match(t, r) - 1L, rows[run], tr),
                                 lapply(tr, function(t) res$b[t] / res$se[t]),
                                 if (return_w) lapply(tr, function(t) res$N[t]) else list(),
                                 Nmed, thresh, threads, return_w)

  # fill imputed SNPs (same scale as SBayesRC: b = z sqrt(var_y) / sqrt(2pq (N + z^2)))
  res[, r2 := 1]
  # set(i = NULL) would touch every row, hence as.integer()
  miss <- as.integer(unlist(lapply(rows[run], function(r) r[!typed[r]]), use.names = FALSE))
  z <- unlist(lapply(out, `[[`, "z"), use.names = FALSE)
  base <- sqrt(2 * res$freq[miss] * (1 - res$freq[miss]) * (res$N[miss] + z^2))
  set(res, miss, c("b", "se", "r2", "p"), list(z * sqrt(vp) / base, sqrt(vp) / base, 0, stats::pchisq(z^2, 1, lower.tail = FALSE)))
  none <- as.integer(unlist(rows[nty == 0], use.names = FALSE))
  set(res, none, c("b", "se", "r2", "p"), list(0, 1, -1, 1))
  res[!is.finite(b), b := 0]
  res[!is.finite(se), se := 1]
  res[!is.finite(p), p := 1]
  meth <- attr(out, "method")
  message(sprintf("Imputed %d SNPs in %d blocks (solver missing/eigen/typed: %d/%d/%d); %d blocks without typed SNPs; %.1f s",
                  length(miss), sum(meth >= 0), sum(meth == 0), sum(meth == 1), sum(meth == 2), sum(nty == 0),
                  proc.time()[[3]] - t0))
  list(ma = res, done = FALSE, blocks = if (return_w) out, run = names(rows)[run])
}
