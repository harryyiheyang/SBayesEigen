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
