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
#' @param ma Summary data (SNP A1 A2 freq b se p N, or SNP A1 A2 A1freq Z N as in \code{\link{tidy}}):
#'   a path or a data.frame, e.g. the output of \code{\link{tidy}}.
#' @param ld LD folder with \code{snp.info} and \code{block*.eigen.bin}; the leading eigen
#'   components holding 99.5\% of the eigenvalue mass are used, as in \code{\link{sbayeseigen}}.
#' @param out Output path (may equal \code{ma}; replaced atomically); \code{NULL} writes nothing.
#' @param threads Blocks imputed in parallel.
#' @return Invisibly, the summary data for every SNP in \code{snp.info} order, with an
#'   \code{r2} column.
#' @export
impute <- function(ma, ld, out = NULL, threads = 4) {
  t0 <- proc.time()[[3]]
  si <- .read_snpinfo(ld)
  ma <- if (.is_zinput(ma)) .z_to_cojo(ma, si) else if (is.data.frame(ma)) as.data.table(ma) else fread(ma, showProgress = FALSE)
  if (.is_imputed(ma, si)) {
    message("Already imputed: r2 exists and the summary data has one row per snp.info SNP")
    return(invisible(ma))
  }
  a <- .align(ma, si); res <- a$res
  rows <- .block_rows(si)
  obs <- is.finite(res$b)
  nty <- vapply(rows, function(r) sum(obs[r]), 0)
  run <- nty > 0 & nty < lengths(rows)
  tr <- lapply(rows[run], function(r) r[obs[r]])
  imp <- impute_blocks_eigen_cpp(.eig_files(ld, names(rows)[run]), list(Map(function(r, t) match(t, r) - 1L, rows[run], tr)),
                                 list(lapply(tr, function(t) res$b[t] / res$se[t])), list(), a$Nmed, 0.995, threads,
                                 FALSE, FALSE, list())
  # fill imputed SNPs on the SBayesRC scale: b = z sqrt(var_y) / sqrt(2pq (N + z^2)); as.integer() because
  # set(i = NULL) would touch every row
  res[, r2 := 1]
  miss <- as.integer(unlist(lapply(rows[run], function(r) r[!obs[r]]), use.names = FALSE))
  z <- unlist(lapply(imp, function(x) x$z[[1]]), use.names = FALSE)
  base <- sqrt(2 * res$freq[miss] * (1 - res$freq[miss]) * (res$N[miss] + z^2))
  set(res, miss, c("b", "se", "r2", "p"), list(z * sqrt(a$vp) / base, sqrt(a$vp) / base, 0, stats::pchisq(z^2, 1, lower.tail = FALSE)))
  none <- as.integer(unlist(rows[nty == 0], use.names = FALSE))
  set(res, none, c("b", "se", "r2", "p"), list(0, 1, -1, 1))
  meth <- attr(imp, "method")
  message(sprintf("Imputed %d SNPs in %d blocks (solver missing/eigen/typed: %d/%d/%d); %d blocks without typed SNPs; %.1f s",
                  length(miss), sum(meth >= 0), sum(meth == 0), sum(meth == 1), sum(meth == 2), sum(nty == 0),
                  proc.time()[[3]] - t0))
  if (!is.null(out)) {
    tmp <- tempfile(pattern = paste0(".", basename(out), ".tmp-"), tmpdir = dirname(out))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    fwrite(res, tmp, sep = "\t", quote = FALSE, na = "NA")
    if (!file.rename(tmp, out)) stop("imputation finished, but replacing ", out, " failed")
  }
  invisible(res)
}

.is_imputed <- function(ma, si) "r2" %in% names(ma) && nrow(ma) == nrow(si)

# snp.info row indices per block, blocks in snp.info order
.block_rows <- function(si) {
  r <- rle(si$Block)
  if (anyDuplicated(r$values)) return(split(seq_len(nrow(si)), factor(si$Block, levels = unique(si$Block))))
  e <- cumsum(r$lengths)
  stats::setNames(Map(seq.int, e - r$lengths + 1L, e), as.character(r$values))
}

# summary data aligned to snp.info (alleles flipped to snp.info A1): one row per snp.info SNP, b/se/p NA
# where missing, N filled with the median; vp = median var(y) estimate, Nmed = median N
.align <- function(ma, si) {
  cols <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  if (!is.null(ma$idx) && !is.null(ma$flip)) return(.align_idx(ma, si))
  if (!all(cols %in% names(ma))) stop("missing columns in the summary data: ", paste(setdiff(cols, names(ma)), collapse = ", "))
  m <- chmatch(as.character(ma$SNP), si$SNP)
  i <- which(!is.na(m)); j <- m[i]
  same <- ma$A1[i] == si$A1[j] & ma$A2[i] == si$A2[j]
  flip <- ma$A1[i] == si$A2[j] & ma$A2[i] == si$A1[j]
  message(length(i), " SNPs in common between the summary data and LD; ", sum(same), " as is, ",
          sum(flip), " with flipped alleles")
  k <- same | flip
  .align_rows(si, i[k], j[k], flip[k], ma$freq, ma$b, ma$se, ma$p, ma$N)
}

# from .tidy(idx = TRUE): snp.info rows and flips are known, no ID matching
.align_idx <- function(tb, si) {
  message(nrow(tb), " SNPs in common between the summary data and LD; ", sum(!tb$flip), " as is, ",
          sum(tb$flip), " with flipped alleles")
  .align_rows(si, seq_len(nrow(tb)), tb$idx, tb$flip, tb$freq, tb$b, tb$se, tb$p, tb$N)
}

# summary rows k -> snp.info rows j (flip: alleles swapped); vp and Nmed from all input rows
.align_rows <- function(si, k, j, flip, freq, b, se, p, N) {
  vp <- stats::median(2 * freq * (1 - freq) * (N * se^2 + b^2), na.rm = TRUE)
  Nmed <- stats::median(N, na.rm = TRUE)
  if (!is.finite(vp) || !is.finite(Nmed)) stop("cannot compute a finite median var(y) or N from the summary data")
  n <- nrow(si); s <- 1 - 2 * flip
  fb <- si$freq; bb <- rep(NA_real_, n); sb <- bb; pb <- bb; Nb <- rep(Nmed, n)
  fb[j] <- flip + s * freq[k]; bb[j] <- s * b[k]; sb[j] <- se[k]; pb[j] <- p[k]; Nb[j] <- N[k]
  bad <- !is.finite(bb) | !is.finite(sb); bb[bad] <- NA_real_; sb[bad] <- NA_real_
  Nb[is.na(Nb)] <- Nmed
  res <- setDT(list(SNP = si$SNP, A1 = si$A1, A2 = si$A2, freq = fb, b = bb, se = sb, p = pb, N = Nb))
  list(res = res, vp = vp, Nmed = Nmed)
}
