sim_ma <- function(si, keep = 0.7, seed = 2) {
  set.seed(seed); m <- nrow(si)
  z <- rnorm(m)
  f <- si$A1Freq; N <- 1e5 + round(rnorm(m, 0, 1e3))
  se <- 1 / sqrt(2 * f * (1 - f) * (N + z^2)); b <- z * se
  ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = b, se = se,
                               p = 2 * pnorm(-abs(z)), N = N)
  fl <- sample(m, 20); ma[fl, `:=`(A1 = si$A2[fl], A2 = si$A1[fl], freq = 1 - freq, b = -b)]
  ma[sort(sample(m, round(keep * m)))]
}

test_that("tidy: C++ path equals R path", {
  fx <- make_ld(); si <- data.table::fread(file.path(fx$ld, "snp.info"))
  ma <- sim_ma(si); f <- tempfile(); data.table::fwrite(ma, f, sep = "\t")
  a <- suppressWarnings(suppressMessages(tidy(f, fx$ld)))
  b <- suppressWarnings(suppressMessages(tidy(ma, fx$ld)))
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("impute equals the direct R reconstruction", {
  fx <- make_ld(); si <- data.table::fread(file.path(fx$ld, "snp.info"))
  ma <- sim_ma(si)
  imp <- suppressWarnings(suppressMessages(impute(suppressMessages(tidy(ma, fx$ld)), fx$ld)))
  expect_equal(imp$SNP, si$ID)
  expect_equal(sum(imp$r2 == 1), nrow(ma))
  z <- imp$b / imp$se
  for (b in unique(si$Block)) {
    e <- read_eig(file.path(fx$ld, paste0("block", b, ".eigen.bin")))
    i <- which(si$Block == b); o <- imp$r2[i] == 1
    UO <- e$U[o, ]; UM <- e$U[!o, ]
    zm <- UM %*% (e$lam * t(UO)) %*% solve(UO %*% (e$lam * t(UO)) + 0.1 * diag(sum(o)), z[i][o])
    expect_equal(z[i][!o], drop(zm), tolerance = 1e-3)
  }
})

test_that("sbayeseigen: fused path equals tidy -> impute -> fit (causal fraction 1%)", {
  fx <- make_ld(); si <- data.table::fread(file.path(fx$ld, "snp.info"))
  set.seed(3); m <- nrow(si); n <- 2e4
  X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1)
  beta <- numeric(m); cs <- sample(m, 3); beta[cs] <- rnorm(3, 0, sqrt(0.05 / 3))
  R <- crossprod(X); bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
  f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s,
                               se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n)
  set.seed(4); ma <- ma[sort(sample(m, 0.8 * m))]
  f1 <- tempfile(); data.table::fwrite(ma, f1, sep = "\t")
  a <- suppressMessages(sbayeseigen(f1, fx$ld, threads = 2))
  f2 <- tempfile(); suppressMessages(impute(suppressMessages(tidy(f1, fx$ld)), fx$ld, f2))
  b <- suppressMessages(sbayeseigen(f2, fx$ld, threads = 2))
  expect_equal(a$snpRes$beta_std, b$snpRes$beta_std, tolerance = 1e-6)
  expect_equal(nrow(a$snpRes), m)
  bs <- a$snpRes$beta_std
  expect_gt(sum(bs * (R %*% beta)) / sqrt(sum(bs * (R %*% bs)) * sum(beta * (R %*% beta))), 0.7)
})
