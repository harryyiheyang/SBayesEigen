#' Tidy GWAS summary data (fast)
#'
#' Same quality control as \code{SBayesRC::tidy}, in the same order and with the same
#' thresholds: finite values and \eqn{0 \le p \le 1}; the first valid row of each SNP ID;
#' SNPs in the LD reference; alleles consistent in either orientation; allele frequency within
#' \code{freq_thresh} of the reference; per-SNP N within mean \eqn{\pm} \code{N_sd_range} SD;
#' and the rate2pq check \eqn{|\sqrt{2pq(N se^2 + b^2)/V_y} - 1| <} \code{rate2pq}, with
#' \eqn{V_y} the median of \eqn{2pq(N se^2 + b^2)}. A plain-text file is processed in C++
#' (files read whole, lines split in parallel, SNPs matched by hashing, output lines copied
#' from the input); a data.frame or a compressed file goes through R.
#'
#' @param ma Summary data, a path or a data.frame: COJO format (SNP A1 A2 freq b se p N), or the
#'   simplified format SNP A1 A2 A1freq Z N (a Z / Zscore column and no b, se; names case-insensitive;
#'   freq may be missing, then the LD reference's is used; rows with a non-finite Z are dropped).
#'   Z input is converted to COJO on the phenotype-SD scale (\eqn{V_y = 1}).
#' @param ld LD folder with \code{snp.info}.
#' @param out Output path; \code{NULL} writes nothing.
#' @param freq_thresh Maximum allele frequency difference to the LD reference.
#' @param N_sd_range Keep SNPs with N within mean +- \code{N_sd_range} SD.
#' @param rate2pq Tolerance of the rate2pq check.
#' @return Invisibly, the tidied summary data in \code{snp.info} order.
#' @export
tidy <- function(ma, ld, out = NULL, freq_thresh = 0.2, N_sd_range = 3, rate2pq = 0.5) {
  .tidy(ma, ld, out, freq_thresh, N_sd_range, rate2pq)
}

# si: standardised snp.info (.std_snpinfo) already in memory; its strings are reused. idx = TRUE (C++ path
# only) also returns the snp.info row (idx) and allele flip of each SNP, so .align can skip matching IDs.
.tidy <- function(ma, ld, out = NULL, freq_thresh = 0.2, N_sd_range = 3, rate2pq = 0.5, si = NULL, idx = FALSE) {
  cols <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  if (.is_zinput(ma)) ma <- .z_to_cojo(ma, if (is.null(si)) .read_snpinfo(ld) else si)
  if (is.character(ma) && !grepl("\\.(gz|bz2|zip)$", ma)) {
    r <- tidy_cpp(path.expand(ma), path.expand(file.path(ld, "snp.info")), if (is.null(out)) "" else path.expand(out),
                  freq_thresh, N_sd_range, rate2pq, is.null(si), getDTthreads())
    k <- r$counts
    message(k[["n_ma"]], " SNPs in summary data")
    message(k[["n_valid"]], " valid SNPs in summary data")
    message(k[["n_ld"]], " SNPs in LD information")
    message(k[["n_common"]], " SNPs in common with LD information")
    message(k[["n_allele"]], " SNPs have consistent alleles (A1, A2) between the summary data and LD")
    if (k[["n_dup"]] > 0) message(k[["n_dup"]], " later rows of duplicated SNP IDs ignored (first valid row kept)")
    message(k[["n_freq"]], " SNPs passed the allele frequency checking with threshold ", freq_thresh)
    message("Mean N = ", k[["mean_N"]], ", SD = ", k[["sd_N"]])
    message(k[["n_N"]], " SNPs have sample size within mean +- ", N_sd_range, "SD")
    message("After QC, Median sample size: ", k[["median_N"]])
    message("Var_y: ", k[["vary"]])
    message(k[["n_out"]], " SNPs remained after QC by rate2pq ", rate2pq)
    if (k[["n_out"]] / k[["n_ld"]] < 0.7)
      warning("Too many SNPs (>30%) were missing in the summary data after QC. The results may be unreliable.")
    tb <- r$table
    if (idx && !is.null(si)) return(invisible(setDT(tb)))
    if (is.null(tb$SNP)) {
      if (nrow(si) != k[["n_ld"]]) stop("snp.info in memory does not match ", file.path(ld, "snp.info"))
      i <- tb$idx
      tb$SNP <- si$SNP[i]
      tb$A1 <- ifelse(tb$flip, si$A2[i], si$A1[i])
      tb$A2 <- ifelse(tb$flip, si$A1[i], si$A2[i])
    }
    return(invisible(setDT(tb[cols])))
  }
  ma <- if (is.data.frame(ma)) as.data.table(ma) else fread(ma, showProgress = FALSE)
  if (!all(cols %in% names(ma))) stop("The summary data is not a valid COJO format (SNP A1 A2 freq b se p N)")
  ma <- ma[, ..cols]
  for (x in cols[4:8]) if (!is.double(ma[[x]])) set(ma, j = x, value = suppressWarnings(as.numeric(ma[[x]])))
  message(nrow(ma), " SNPs in summary data")
  ma <- ma[is.finite(N) & is.finite(b) & is.finite(se) & is.finite(freq) & p >= 0 & p <= 1]
  message(nrow(ma), " valid SNPs in summary data")
  ma <- ma[!duplicated(SNP)]   # as SBayesRC: first valid row of each SNP ID
  if (is.null(si)) si <- .read_snpinfo(ld)
  message(nrow(si), " SNPs in LD information")
  # match() is much faster than a keyed join on 7M strings; result in snp.info order
  m <- match(ma$SNP, si$SNP)
  d <- ma[!is.na(m)]
  m <- m[!is.na(m)]
  d[, `:=`(rA1 = si$A1[m], rA2 = si$A2[m], rfreq = si$freq[m], ord = m)]
  message(nrow(d), " SNPs in common with LD information")
  same <- d$A1 == d$rA1 & d$A2 == d$rA2
  flip <- d$A1 == d$rA2 & d$A2 == d$rA1
  d <- d[same | flip]
  message(nrow(d), " SNPs have consistent alleles (A1, A2) between the summary data and LD")
  d[, frq_ref := ifelse(A1 == rA1, rfreq, 1 - rfreq)]
  d <- d[abs(frq_ref - freq) <= freq_thresh]
  message(nrow(d), " SNPs passed the allele frequency checking with threshold ", freq_thresh)
  setorder(d, ord)
  mN <- mean(d$N); sN <- stats::sd(d$N)
  message("Mean N = ", mN, ", SD = ", sN)
  d <- d[N >= mN - N_sd_range * sN & N <= mN + N_sd_range * sN]
  message(nrow(d), " SNPs have sample size within mean +- ", N_sd_range, "SD")
  message("After QC, Median sample size: ", stats::median(d$N))
  vp <- 2 * d$freq * (1 - d$freq) * (d$N * d$se^2 + d$b^2)
  vary <- stats::median(vp)
  indic <- sqrt(vp / vary)
  message("Var_y: ", vary)
  d <- d[indic > 1 - rate2pq & indic < 1 + rate2pq, ..cols]
  message(nrow(d), " SNPs remained after QC by rate2pq ", rate2pq)
  if (nrow(d) / nrow(si) < 0.7)
    warning("Too many SNPs (>30%) were missing in the summary data after QC. The results may be unreliable.")
  if (!is.null(out)) fwrite(d, out, sep = "\t", quote = FALSE, na = "NA")
  invisible(d)
}

# snp.info with standard names SNP, A1, A2, freq, N, Block (SBayesRC: Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block)
.read_snpinfo <- function(ld) {
  f <- file.path(ld, "snp.info")
  h <- names(fread(f, nrows = 0, showProgress = FALSE))
  need <- c("ID", "A1", "A2", "A1Freq", "N", "Block")
  if (all(need %in% h)) {   # read only the needed columns (2-3x faster on 7M SNPs)
    si <- fread(f, select = need, colClasses = list(character = c("ID", "A1", "A2")), showProgress = FALSE)
    setnames(si, c("ID", "A1Freq"), c("SNP", "freq"))
    setcolorder(si, c("SNP", "A1", "A2", "freq", "N", "Block"))   # in place, no copy of 7M strings
    return(si)
  }
  .std_snpinfo(fread(f, showProgress = FALSE))
}

.std_snpinfo <- function(si) {
  si <- as.data.table(si)
  if ("ID" %in% names(si)) {
    out <- si[, .(SNP = ID, A1, A2, freq = A1Freq, N, Block)]
  } else {
    nm <- switch(as.character(ncol(si)), "8" = c(2, 4:8), "9" = c(2, 5:9), "10" = c(2, 6:10),
                 stop("the LD information looks odd"))
    out <- si[, nm, with = FALSE]
    setnames(out, c("SNP", "A1", "A2", "freq", "N", "Block"))
  }
  out[, SNP := as.character(SNP)]
  out
}

# Z-score input (SNP A1 A2 [A1freq] Z N, column names case-insensitive; Yihe 2026-10-08): the model only needs z
# and N, so a Z column (when there are no b/se columns) is converted to COJO on the phenotype-SD scale, b = Z / sqrt(2pq (N + Z^2)), se = b / Z, Var_y = 1.
# A missing freq column (or value) takes the LD reference's frequency of A1.
.zcols <- c("z", "zscore")
.fcols <- c("freq", "a1freq", "frq", "af", "eaf")
# simplified format: a Z column (case-insensitive) and no b/se columns
.is_zinput <- function(ma) {
  h <- tolower(if (is.character(ma)) names(fread(ma, nrows = 0, showProgress = FALSE)) else names(ma))
  any(.zcols %in% h) && !all(c("b", "se") %in% h)
}
.z_to_cojo <- function(ma, si) {
  ma <- if (is.data.frame(ma)) as.data.table(ma) else fread(ma, showProgress = FALSE)
  lc <- tolower(names(ma)); col <- function(x) names(ma)[match(x, lc)]
  std <- c(snp = "SNP", a1 = "A1", a2 = "A2", n = "N")
  for (k in names(std)) if (k %in% lc && !std[[k]] %in% names(ma)) setnames(ma, col(k), std[[k]])
  lc <- tolower(names(ma))
  miss <- setdiff(c("SNP", "A1", "A2", "N"), names(ma))
  if (length(miss)) stop("Z-score summary data needs SNP A1 A2 Z N; missing: ", paste(miss, collapse = ", "))
  z <- suppressWarnings(as.numeric(ma[[col(intersect(.zcols, lc)[1])]]))
  fc <- intersect(.fcols, lc)
  f <- if (length(fc)) suppressWarnings(as.numeric(ma[[col(fc[1])]])) else rep(NA_real_, nrow(ma))
  if (anyNA(f)) {
    m <- chmatch(as.character(ma$SNP), si$SNP); fr <- ifelse(ma$A1 == si$A1[m], si$freq[m], 1 - si$freq[m])
    message(sum(is.na(f)), " SNPs without freq take the LD reference frequency")
    f[is.na(f)] <- fr[is.na(f)]
  }
  N <- suppressWarnings(as.numeric(ma$N)); s <- 1 / sqrt(2 * f * (1 - f) * (N + z^2))
  message("Z-score input: b and se on the phenotype-SD scale (Var_y = 1)")
  data.table(SNP = as.character(ma$SNP), A1 = as.character(ma$A1), A2 = as.character(ma$A2), freq = f,
             b = z * s, se = s, p = stats::pchisq(z^2, 1, lower.tail = FALSE), N = N)
}
