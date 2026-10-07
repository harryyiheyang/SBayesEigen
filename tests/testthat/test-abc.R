test_that("ABC: A = annotation z-leads, C selected by MCP, beta = U alpha + A/C effects", {
  fx <- make_ld(); si <- data.table::fread(file.path(fx$ld, "snp.info"))
  set.seed(21); m <- nrow(si); n <- 5e4
  X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1); R <- crossprod(X)
  beta <- numeric(m); cs <- c(20, 90, 200); beta[cs] <- c(0.12, -0.1, 0.08)
  bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
  f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s,
                               se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n)
  af <- tempfile(); writeLines(si$ID[c(15:25, 195:205)], af)       # annotation covers causal 20 and 200
  fit <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, annot = af))
  ab <- fit$par$abc
  expect_true(all(ab[set == "A"]$SNP %in% si$ID[c(15:25, 195:205)]))
  expect_true(all(abs(ab$z) > 4))
  expect_true(any(ab$set == "A") && any(ab$set == "C"))
  # A leads are not in r2 > 0.9 with each other (within a block)
  ia <- match(ab[set == "A"]$SNP, si$ID)
  if (length(ia) > 1) { r2 <- (cov2cor(R[ia, ia])^2)[upper.tri(diag(length(ia)))]; same <- outer(si$Block[ia], si$Block[ia], "==")[upper.tri(diag(length(ia)))]
    expect_true(all(r2[same] <= 0.9 + 0.02)) }   # eigen-truncated R vs sample R
  # beta_std = U alpha + A/C effects, alpha from par$comp
  cp <- fit$par$comp; bs <- numeric(m)
  for (b in unique(si$Block)) {
    e <- read_eig(file.path(fx$ld, paste0("block", b, ".eigen.bin"))); a <- cp[Block == b]$alpha
    bs[si$Block == b] <- drop(e$U[, seq_along(a)] %*% a)
  }
  i <- match(ab$SNP, si$ID); bs[i] <- bs[i] + ab$beta_std
  expect_equal(fit$snpRes$beta_std, bs, tolerance = 1e-5)
  acc <- function(x) sum(x * (R %*% beta)) / sqrt(sum(x * (R %*% x)) * sum(beta * (R %*% beta)))
  expect_gt(acc(fit$snpRes$beta_std), 0.9)
  # method = "eigen" ignores annot and has no ABC table
  e <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, method = "eigen", annot = af))
  expect_null(e$par$abc)
})

test_that("annotation file with or without a SNP header", {
  si <- data.frame(SNP = paste0("rs", 1:10))
  f1 <- tempfile(); writeLines(paste0("rs", 3:5), f1)
  f2 <- tempfile(); writeLines(c("CHR\tSNP", paste0("1\trs", 3:5)), f2)
  for (f in c(f1, f2)) expect_equal(which(suppressMessages(SBayesEigen:::.read_annot(f, si))), 3:5)
})
