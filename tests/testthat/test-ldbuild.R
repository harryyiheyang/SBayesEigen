test_that("LDbuild eigen files reproduce the sample correlation", {
  fx <- make_ld()
  si <- data.table::fread(file.path(fx$ld, "snp.info"))
  expect_equal(nrow(si), 300)
  for (b in unique(si$Block)) {
    e <- read_eig(file.path(fx$ld, paste0("block", b, ".eigen.bin")))
    j <- match(si$ID[si$Block == b], paste0("rs", 1:300))
    R <- cor(fx$G[, j])
    ev <- eigen(R, symmetric = TRUE)
    k <- length(e$lam)
    expect_true(sum(e$lam) >= 0.995 * sum(ev$values) - 1e-3)
    expect_equal(crossprod(e$U), diag(k), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(e$lam, ev$values[1:k], tolerance = 1e-5)
    Rk <- ev$vectors[, 1:k] %*% (ev$values[1:k] * t(ev$vectors[, 1:k]))
    expect_lt(max(abs(e$U %*% (e$lam * t(e$U)) - Rk)), 1e-4)
  }
  l <- data.table::fread(file.path(fx$ld, "ldscore.txt"))
  expect_named(l, c("SNP", "Block", "ldscore"))
  R1 <- cor(fx$G[, match(l$SNP[l$Block == 1], paste0("rs", 1:300))]); n <- 600
  expect_equal(l$ldscore[l$Block == 1], colSums(R1^2 - (1 - R1^2) / (n - 2)), tolerance = 1e-6)
  expect_error(suppressMessages(LDbuild(file.path(tempdir(), "geno.bed"), fx$ld)), "already")
})

test_that("LDbuild with A writes ab.bin: exact R_AA, R_BA and a whitening B basis orthogonal to A", {
  G <- sim_geno(600, 300, seed = 3)
  pre <- file.path(tempdir(), "genoab"); write_bed(G, pre)
  br <- tempfile(fileext = ".txt")   # LDetect-style intervals; the 40-SNP middle block is merged
  writeLines(c("chr \t start \t stop", "chr1 \t 1 \t 130500", "chr1 \t 130500 \t 170500", "chr1 \t 170500 \t 1000000"), br)
  A <- paste0("rs", c(seq(5, 300, 9), 50, 51))   # rs50/rs51 strongly correlated
  dir <- tempfile("ldab")
  suppressMessages(LDbuild(paste0(pre, ".bed"), dir, threads = 2, blockRef = br, minsnp = 100, A = A))
  si <- data.table::fread(file.path(dir, "snp.info"))
  expect_equal(as.integer(table(si$Block)), c(170L, 130L))
  e <- data.table::fread(file.path(dir, "eigen.info"))
  expect_true(all(c("mA", "rankA") %in% names(e)))
  for (b in unique(si$Block)) {
    ids <- si$ID[si$Block == b]
    R <- cor(G[, match(ids, paste0("rs", 1:300))])
    f <- SBayesEigen:::.read_ab(file.path(dir, paste0("block", b, ".ab.bin")))
    expect_identical(ids[f$iA], ids[ids %in% A])
    expect_equal(f$RAA, R[f$iA, f$iA], tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(f$RBA, R[-f$iA, f$iA], tolerance = 1e-6, ignore_attr = TRUE)
    P <- sweep(f$UlB, 2, sqrt(f$lambda), "/")
    expect_lt(max(abs(crossprod(P, R %*% P) - diag(f$kB))), 1e-4)
    expect_lt(max(abs(R[f$iA, ] %*% P)), 1e-4)
    pinv <- function(X) { e <- eigen(X, symmetric = TRUE); k <- e$values > 1e-12 * max(e$values)
      e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) / e$values[k]) }
    S <- R[-f$iA, -f$iA] - R[-f$iA, f$iA] %*% pinv(R[f$iA, f$iA]) %*% R[f$iA, -f$iA]
    expect_equal(f$lambda, eigen(S, symmetric = TRUE, only.values = TRUE)$values[1:f$kB], tolerance = 1e-4)
  }
})
