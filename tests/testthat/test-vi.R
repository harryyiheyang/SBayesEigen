sim_w <- function(J = 3000, seed = 5) {
  set.seed(seed); n <- rep(5e4, J)
  c <- sqrt(rexp(J)); a <- numeric(J); nz <- sample(J, 30); a[nz] <- rnorm(30, 0, sqrt(0.1 / sum(c[nz]^2)))
  list(w = c * a + rnorm(J, 0, sqrt(1 / n)), c = c, n = n, a = a)
}

test_that("em_step_cpp equals the R EM map", {
  s <- sim_w(); g <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1); p <- c(0.9, rep(0.02, 5))
  a <- SBayesEigen:::em_step_cpp(s$w, s$c, s$n, g, p, 0.01, 1, 0.2, 4, 0.5, 2)
  b <- em_step_r(s$w, s$c, s$n, g, p, 0.01, 1, 0.2, 4, 0.5)
  for (x in names(b)) expect_equal(a[[x]], b[[x]], tolerance = 1e-10)
})

test_that("vi_eigen ascends, reaches the plain-EM optimum and recovers Vg", {
  s <- sim_w(); g <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1); K <- 6
  f <- SBayesEigen:::vi_eigen(s$w, s$c, s$n, ve = 1, h2p = 0.1, threads = 2)
  expect_true(f$converged)
  expect_true(all(diff(f$objective) > -1e-6))
  # 5000 plain EM steps from the same start
  sc2 <- sum(s$c^2); p <- rep(1 / K, K); A0 <- sum(p * g) * sc2
  s2 <- max(sum((s$n * s$w^2 - 1) * s$n * s$c^2) / sum((s$n * s$c^2)^2), 1e-8 / sc2) / sum(p * g)
  for (it in 1:5000) { e <- em_step_r(s$w, s$c, s$n, g, p, s2, 1, 0.2, 4, A0); p <- e$pi; s2 <- e$sigma2 }
  # the Vg stop rule leaves pi/sigma2 on a flat ridge: objective within 0.002 per component
  expect_gte(max(f$objective), e$obj - 2e-3 * length(s$w))
  ps <- SBayesEigen:::post_cpp(s$w, s$c, s$n, g, p, s2, 1, 2)
  expect_gt(cor(f$alpha, ps$alpha), 0.999)
  expect_equal(f$Vg, e$vg, tolerance = 0.02)
  expect_equal(f$Vg, sum(s$c^2 * s$a^2), tolerance = 0.2)
})

test_that("ldsc_eigen recovers h2 and intercept", {
  set.seed(6); M <- 2e4; l <- rgamma(M, 2, 0.05); n <- 1e5; h2 <- 0.3
  chi2 <- (1 + n * h2 * l / M) * rchisq(M, 1)
  r <- SBayesEigen:::ldsc_eigen(sqrt(chi2 / n), l, n, M)
  expect_equal(r$h2, h2, tolerance = 0.15)
  expect_equal(r$intercept, 1, tolerance = 0.05)
})

test_that("noise_mom recovers ve0 and kappa and respects the bounds", {
  set.seed(7)
  K <- 2e5; lam <- exp(runif(K, log(0.02), log(50))); n <- rep(3e5, K)
  h2 <- 0.2; tau <- h2 / sum(lam)
  a <- rnorm(K, 0, sqrt(tau))
  sim <- function(ve0, kappa) sqrt(lam) * a + rnorm(K, 0, sqrt((ve0 + kappa / lam) / n))
  m <- SBayesEigen:::noise_mom(sim(1.1, 0.05), lam, n, h2)
  expect_equal(m$ve0, 1.1, tolerance = 0.03)
  expect_equal(m$kappa, 0.05, tolerance = 0.15)
  m0 <- SBayesEigen:::noise_mom(sim(1, 0), lam, n, h2)          # matched LD: kappa ~ 0
  expect_lt(m0$kappa, 0.005); expect_equal(m0$ve0, 1, tolerance = 0.03)
  mb <- SBayesEigen:::noise_mom(sim(2, 0.05), lam, n, h2)        # ve0 2 inside [0.5, 3] (MVP: 1.3-1.9)
  expect_equal(mb$ve0, 2, tolerance = 0.03); expect_gte(mb$kappa, 0)
  mc <- SBayesEigen:::noise_mom(sim(4, 0.05), lam, n, h2)        # clamped to 3, kappa absorbs some
  expect_equal(mc$ve0, 3); expect_gte(mc$kappa, 0)
  ml <- SBayesEigen:::noise_mom(sim(0.55, 0), lam, n, h2)       # overstated N (QuickDraws-like): ve0 < 1 allowed
  expect_equal(ml$ve0, 0.55, tolerance = 0.03)
})
