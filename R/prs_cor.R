# Correlation of two PRS weight sets under the LD of an ab.bin reference (Yihe 2026-10-10 06:21: after a prior change,
# compare new and old PRS directly; traits with r >= 0.99 are not rescored).
#   r = b1' R b2 / sqrt(b1' R b1 * b2' R b2), b = standardised effects (beta_std), R ~ F F' per block (as in pass 1:
#   F = [R_AA H, 0; R_BA H, Q2 Lambda^{1/2}], H = V_A Lambda_A^{-1/2} over R_AA eigenvalues > tolA * max), so
#   b' R c = sum_blocks (F'b)'(F'c) with F'b = (H'(R_AA b_A + R_AB b_B), Lambda^{1/2} Q2' b_B). Streams block by block.
# x1, x2: snpRes tables (or paths) with SNP, A1 and beta_std (or beta, scaled by sqrt(2 p q) of the reference);
# matched to snp.info by SNP, sign flipped where A1 is the reference A2; SNPs missing from one set count as 0.
# Returns r overall and per chromosome, with b'Rb of each set (= Vg of each PRS under the reference LD).
# Internal (SBayesEigen:::prs_cor), not exported.
prs_cor <- function(x1, x2, ld, threads = 1) {
  si <- fread(file.path(ld, "snp.info"))
  b1 <- .prs_std(x1, si); b2 <- .prs_std(x2, si)
  blk <- unique(si$Block)
  one <- function(b) {
    i <- which(si$Block == b); f <- file.path(ld, paste0("block", b, ".ab.bin"))
    if (!file.exists(f)) stop("prs_cor needs ab.bin LD: ", f, " not found")
    fb <- .read_ab(f); u1 <- .ab_Ft(fb, b1[i]); u2 <- .ab_Ft(fb, b2[i])
    c(s12 = sum(u1 * u2), s11 = sum(u1^2), s22 = sum(u2^2))
  }
  s <- if (threads > 1 && .Platform$OS.type == "unix") parallel::mclapply(blk, one, mc.cores = threads) else lapply(blk, one)
  s <- data.table(Block = blk, Chrom = si$Chrom[match(blk, si$Block)], do.call(rbind, s))
  sm <- function(d) d[, .(r = sum(s12) / sqrt(sum(s11) * sum(s22)), bRb_1 = sum(s11), bRb_2 = sum(s22))]
  list(r = sm(s)$r, total = sm(s), by_chr = s[, sm(.SD), by = Chrom], by_block = s)
}

.ab_Ft <- function(f, b) {
  iA <- f$iA; bB <- b[-iA]; bA <- b[iA]
  uB <- sqrt(pmax(f$lambda, 0)) * drop(crossprod(f$UlB[-iA, , drop = FALSE], bB))
  if (!f$mA) return(uB)
  e <- eigen(f$RAA, symmetric = TRUE); k <- e$values > f$tolA * e$values[1]
  H <- e$vectors[, k, drop = FALSE] %*% diag(1 / sqrt(e$values[k]), sum(k))
  c(drop(crossprod(H, f$RAA %*% bA + crossprod(f$RBA, bB))), uB)
}

.prs_std <- function(x, si) {
  if (is.character(x)) x <- fread(x)
  x <- as.data.table(x); j <- match(si$ID, x$SNP)
  b <- if ("beta_std" %in% names(x)) x$beta_std[j] else x$beta[j] * sqrt(2 * si$A1Freq * (1 - si$A1Freq))
  s <- ifelse(x$A1[j] == si$A1, 1, ifelse(x$A1[j] == si$A2, -1, NA))
  b <- b * s; b[is.na(b)] <- 0; b
}
