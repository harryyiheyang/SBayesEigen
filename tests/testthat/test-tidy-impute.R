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

test_that("several traits together equal each trait on its own", {
  fx <- make_ld(); si <- data.table::fread(file.path(fx$ld, "snp.info"))
  m <- nrow(si); X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1); R <- crossprod(X)
  f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  mk <- function(seed, n, keep) {
    set.seed(seed); beta <- numeric(m); cs <- sample(m, 3); beta[cs] <- rnorm(3, 0, sqrt(0.05 / 3))
    bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
    ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s,
                                 se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n + round(rnorm(m, 0, 100)))
    ma[sort(sample(m, keep * m))]
  }
  d <- tempfile(); dir.create(d)
  fs <- file.path(d, c("A.ma", "B.ma", "C.ma"))
  data.table::fwrite(mk(11, 2e4, 0.8), fs[1], sep = "\t")
  data.table::fwrite(mk(12, 5e4, 0.6), fs[2], sep = "\t")
  m3 <- mk(13, 1e4, 0.9); m3 <- m3[!si$Block[match(m3$SNP, si$ID)] == 2]   # trait C: no SNP in block 2
  data.table::fwrite(m3, fs[3], sep = "\t")
  fi <- file.path(d, "B.imp")
  suppressWarnings(suppressMessages(impute(suppressWarnings(suppressMessages(tidy(fs[2], fx$ld))), fx$ld, fi)))
  one <- lapply(c(fs[1], fi, fs[3]), function(x) suppressWarnings(suppressMessages(sbayeseigen(x, fx$ld, threads = 2))))
  o <- file.path(d, "out"); dir.create(file.path(o, "A"), recursive = TRUE)   # existing trait folder is reused
  all <- suppressWarnings(suppressMessages(sbayeseigen(c(A = fs[1], B = fi, C = fs[3]), fx$ld, out = o, threads = 2)))
  expect_equal(colnames(all$beta), c("A", "B", "C"))
  for (t in 1:3) {
    expect_equal(all$beta[, t], one[[t]]$snpRes$beta, tolerance = 1e-10)
    expect_equal(all$par[[t]]$Vg, one[[t]]$par$Vg, tolerance = 1e-10)
    r <- data.table::fread(file.path(o, c("A", "B", "C")[t], paste0(c("A", "B", "C")[t], "_sbeigen.snpRes")))
    expect_true(file.exists(file.path(o, c("A", "B", "C")[t], paste0(c("A", "B", "C")[t], "_sbeigen.log"))))
    expect_equal(r$beta, one[[t]]$snpRes$beta, tolerance = 1e-6)
  }
  expect_true(all(all$beta[si$Block == 2, 3] == 0))
  # stored eigen-space fit: beta' R beta from alpha equals the SNP-space value
  cp <- all$par$A$comp; bs <- one[[1]]$snpRes$beta_std
  q <- sum(sapply(unique(si$Block), function(b) {
    e <- read_eig(file.path(fx$ld, paste0("block", b, ".eigen.bin"))); x <- bs[si$Block == b]
    sum(e$lam * crossprod(e$U, x)^2)
  }))
  expect_equal(sum(cp$lam * cp$alpha^2), q, tolerance = 1e-5)
  expect_equal(SBayesEigen:::.trait_names(c("x/APOA1.ma", "y/LDL.txt.gz")), c("APOA1", "LDL"))
  expect_error(SBayesEigen:::.trait_names(c("a/T.ma", "b/T.ma")), "unique")
})
