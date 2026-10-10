ab_fixture <- function() {
  G <- sim_geno(600, 300, seed = 3)
  pre <- file.path(tempdir(), "genoabj"); write_bed(G, pre)
  br <- tempfile(fileext = ".txt")
  writeLines(c("chr \t start \t stop", "chr1 \t 1 \t 150500", "chr1 \t 150500 \t 1000000"), br)
  A <- paste0("rs", c(seq(5, 300, 9), 50, 51))
  dir <- tempfile("ldabj")
  suppressMessages(LDbuild(paste0(pre, ".bed"), dir, threads = 2, blockRef = br, minsnp = 100, A = A))
  list(G = G, ld = dir, A = A, si = data.table::fread(file.path(dir, "snp.info")))
}
# R pieces of one ab block: H, F (R ~ F F'), XA, Cm, Q2
ab_pieces <- function(f) {
  e <- eigen(f$RAA, symmetric = TRUE); k <- e$values > f$tolA * e$values[1]
  H <- e$vectors[, k, drop = FALSE] %*% diag(1 / sqrt(e$values[k]), sum(k))
  Q2 <- f$UlB[-f$iA, , drop = FALSE]; G <- f$RBA %*% H
  Fm <- matrix(0, f$m, ncol(H) + f$kB)
  Fm[f$iA, seq_len(ncol(H))] <- f$RAA %*% H
  Fm[-f$iA, ] <- cbind(G, sweep(Q2, 2, sqrt(f$lambda), "*"))
  list(H = H, Q2 = Q2, G = G, F = Fm, XA = crossprod(H, f$RAA), Cm = crossprod(G, Q2), lamA = e$values[k])
}

test_that("ab.bin pass 1: imputation with R ~ F F', w_A = H' bhat_A, w_B = Lambda^-1/2 UlB' bhat, design pieces", {
  fx <- ab_fixture(); si <- fx$si
  set.seed(5); z <- rnorm(nrow(si)) * 2; N <- 1e4
  for (b in unique(si$Block)) {
    f <- SBayesEigen:::.read_ab(file.path(fx$ld, paste0("block", b, ".ab.bin"))); P <- ab_pieces(f)
    zb <- z[si$Block == b]; ty <- sort(sample(f$m, 0.8 * f$m)); cr <- c(0L, 3L, 7L, 40L)
    o <- SBayesEigen:::ab_pass1_cpp(file.path(fx$ld, paste0("block", b, ".ab.bin")), list(list(ty - 1L)), list(list(zb[ty])),
                                    list(list(rep(N, length(ty)))), N, list(cr), 1, 1)[[1]]
    R <- tcrossprod(P$F)
    zm <- R[-ty, ty] %*% solve(R[ty, ty] + 0.1 * diag(length(ty)), zb[ty])
    expect_equal(o$z[[1]], drop(zm), tolerance = 1e-3)
    zf <- numeric(f$m); zf[ty] <- zb[ty]; zf[-ty] <- o$z[[1]]; bh <- zf / sqrt(N + zf^2)
    # H is unique up to the signs and order of R_AA's eigenvectors: compare basis-free products
    rk <- ncol(P$H); wA <- o$w[[1]][seq_len(rk)]
    expect_equal(o$w[[1]][-seq_len(rk)], drop(crossprod(f$UlB, bh)) / sqrt(f$lambda), tolerance = 1e-4)
    expect_equal(drop(crossprod(o$XA, wA)), drop(crossprod(P$XA, crossprod(P$H, bh[f$iA]))), tolerance = 1e-4)
    expect_equal(crossprod(o$XA), crossprod(P$XA), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(crossprod(o$XA, o$Cm), crossprod(P$XA, P$Cm), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(crossprod(o$Cm), crossprod(P$Cm), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(o$iA, f$iA)
    cb <- setdiff(cr + 1L, f$iA); expect_equal(o$crow, cb)
    qb <- match(cb, seq_len(f$m)[-f$iA])
    expect_equal(crossprod(o$XA, o$XCa), crossprod(P$XA, t(P$G[qb, , drop = FALSE])), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(o$XCb, t(sweep(P$Q2[qb, , drop = FALSE], 2, sqrt(f$lambda), "*")), tolerance = 1e-5, ignore_attr = TRUE)
    # whitening: the stacked map T = [H'; Lambda^-1/2 UlB'] gives T R T' = I, so (w - r)'(w - r) = beta' R beta
    Tm <- rbind(t(P$H) %*% diag(f$m)[f$iA, , drop = FALSE], t(sweep(f$UlB, 2, sqrt(f$lambda), "/")))
    expect_lt(max(abs(Tm %*% R %*% t(Tm) - diag(nrow(Tm)))), 1e-3)
  }
})

test_that("abj_sweep_cpp equals a pure-R CAVI sweep", {
  set.seed(8); rk <- 4; kB <- 6; nA <- 3; nC <- 2
  blk <- list(list(XA = matrix(rnorm(rk * nA), rk), Cm = matrix(rnorm(rk * kB, 0, 0.3), rk), sl = sqrt(runif(kB, 0.2, 3)),
                   XCa = matrix(rnorm(rk * nC, 0, 0.3), rk), XCb = matrix(rnorm(kB * nC), kB)))
  w <- rnorm(rk + kB, 0, 0.05); p <- runif(rk + kB, 800, 1200)
  gA <- c(0, 1, 10); gB <- c(0, 1, 100, 500); piA <- c(0.5, 0.3, 0.2); piB <- c(0.4, 0.3, 0.2, 0.1)
  mA <- list(rnorm(nA, 0, 0.01)); al <- list(rnorm(kB, 0, 0.01)); gC <- list(c(0, 0.02))
  r <- w - c(blk[[1]]$XA %*% mA[[1]] + blk[[1]]$Cm %*% al[[1]] + blk[[1]]$XCa %*% gC[[1]],
             blk[[1]]$sl * al[[1]] + blk[[1]]$XCb %*% gC[[1]])
  ref <- function(r, mA, al, gC) {
    X <- rbind(cbind(blk[[1]]$XA, blk[[1]]$Cm, blk[[1]]$XCa), cbind(matrix(0, kB, nA), diag(blk[[1]]$sl, kB), blk[[1]]$XCb))
    co <- c(mA, al, gC); grid <- list(gA, gB); pis <- list(piA, piB); s2 <- c(1e-5, 2e-5)
    for (i in seq_along(co)) {
      x <- X[, i]; D <- sum(p * x^2); rho <- sum(p * x * r) + D * co[i]
      if (i <= nA + kB) {
        g <- if (i <= nA) 1 else 2; gg <- grid[[g]] * s2[g]; v <- ifelse(gg > 0, 1 / (D + 1 / gg), 0); mu <- v * rho
        lw <- log(pis[[g]]) + ifelse(gg > 0, 0.5 * log(v / gg) + 0.5 * mu^2 / v, 0); ph <- exp(lw - max(lw)); ph <- ph / sum(ph)
        nw <- sum(ph * mu)
      } else { h <- rho / sqrt(D); a <- 2.5; tau <- 5
        nw <- (if (abs(h) <= tau) 0 else if (abs(h) <= a * tau) sign(h) * (abs(h) - tau) / (1 - 1 / a) else h) / sqrt(D) }
      r <- r - x * (nw - co[i]); co[i] <- nw
    }
    list(r = r, co = co)
  }
  e <- ref(r, mA[[1]], al[[1]], gC[[1]])
  r2 <- r + 0
  SBayesEigen:::abj_sweep_cpp(blk, 0L, p, r2, mA, al, gC, gA, piA, 1e-5, gB, piB, 2e-5, 5, 2.5, 1)
  expect_equal(r2, e$r, tolerance = 1e-10)
  expect_equal(c(mA[[1]], al[[1]], gC[[1]]), e$co, tolerance = 1e-10)
})

test_that("sbayeseigen on ab.bin: joint ABC, beta = A effects + Q2 alpha + gamma, good accuracy", {
  fx <- ab_fixture(); si <- fx$si; m <- nrow(si); n <- 5e4
  X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1); R <- crossprod(X)
  set.seed(9); beta <- numeric(m); cs <- c(match(c("rs50", "rs140"), si$ID), 30, 220); beta[cs] <- c(0.1, -0.08, 0.09, 0.06)
  bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
  f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s,
                               se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n)[sort(sample(m, 0.9 * m))]
  fit <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2))
  ab <- fit$par$abc
  expect_true(all(ab[set == "A"]$SNP %in% fx$A)); expect_equal(sum(ab$set == "A"), sum(si$ID %in% fx$A))
  expect_true(all(!ab[set == "C"]$SNP %in% fx$A))
  # beta_std = Q2 alpha on B rows + A and C effects
  cp <- fit$par$comp; bs <- numeric(m)
  for (b in unique(si$Block)) {
    fb <- SBayesEigen:::.read_ab(file.path(fx$ld, paste0("block", b, ".ab.bin"))); a <- cp[Block == b & !is.na(alpha)]$alpha
    x <- drop(fb$UlB %*% a); x[fb$iA] <- 0; bs[si$Block == b] <- x
  }
  i <- match(ab$SNP, si$ID); bs[i] <- bs[i] + ab$beta_std
  expect_equal(fit$snpRes$beta_std, bs, tolerance = 1e-5)
  acc <- function(x) sum(x * (R %*% beta)) / sqrt(sum(x * (R %*% x)) * sum(beta * (R %*% beta)))
  expect_gt(acc(fit$snpRes$beta_std), 0.89)   # 0.9002 at 0ea0e9a, 0.8999 with the capped B EM
  expect_error(suppressMessages(sbayeseigen(ma, fx$ld, method = "eigen")), "ab.bin")
  # thresh truncates B at read time (leading components reaching thresh of the Schur complement's mass)
  f9 <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, thresh = 0.9)); cp9 <- f9$par$comp; bs9 <- numeric(m)
  expect_lt(sum(!is.na(cp9$alpha)), sum(!is.na(cp$alpha)))
  for (b in unique(si$Block)) {
    fb <- SBayesEigen:::.read_ab(file.path(fx$ld, paste0("block", b, ".ab.bin"))); a <- cp9[Block == b & !is.na(alpha)]$alpha
    k <- which(cumsum(fb$lambda) >= 0.9 * fb$sumLambda)[1]; expect_length(a, k)
    x <- drop(fb$UlB[, seq_len(k), drop = FALSE] %*% a); x[fb$iA] <- 0; bs9[si$Block == b] <- x
  }
  i <- match(f9$par$abc$SNP, si$ID); bs9[i] <- bs9[i] + f9$par$abc$beta_std
  expect_equal(f9$snpRes$beta_std, bs9, tolerance = 1e-5)
  expect_equal(suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, thresh = 1))$snpRes$beta_std, fit$snpRes$beta_std, tolerance = 1e-6)
  # threshB: B fitted on the leading components only, the rest alpha = 0 but kept in the data (A and C see them)
  expect_equal(suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, threshB = 1))$snpRes$beta_std, fit$snpRes$beta_std, tolerance = 1e-6)
  fb9 <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, threshB = 0.9)); cb <- fb9$par$comp[!is.na(alpha)]
  expect_equal(nrow(cb), sum(!is.na(cp$alpha)))
  expect_equal(fb9$par$abc_fit$nB_fit, f9$par$abc_fit$nB_fit)
  for (b in unique(si$Block)) {
    fb <- SBayesEigen:::.read_ab(file.path(fx$ld, paste0("block", b, ".ab.bin")))
    k <- which(cumsum(fb$lambda) >= 0.9 * fb$sumLambda)[1]; a <- cb[Block == b]$alpha
    expect_true(all(a[-seq_len(k)] == 0)); expect_true(any(a[seq_len(k)] != 0))
  }
  expect_gt(acc(fb9$snpRes$beta_std), 0.85)
  ct <- fb9$par$abc[set == "C"]$frac_trunc; expect_true(all(ct >= 0 & ct <= 1))
  expect_true(all(fit$par$abc[set == "C"]$frac_trunc == 0)); expect_true(all(is.na(fit$par$abc[set == "A"]$frac_trunc)))
  # threshB = "auto": pseudo-validation over 0.995/0.99/0.95/0.9, reproducible, chosen value used
  fa <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, threshB = "auto")); pv <- fa$par$abc_fit$pv
  expect_length(pv$score, 4); expect_true(all(is.finite(pv$score)))
  fa2 <- suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, threshB = "auto")); expect_equal(fa2$snpRes$beta_std, fa$snpRes$beta_std)
  th <- fa$par$abc_fit$threshB
  expect_equal(fa$snpRes$beta_std, suppressMessages(sbayeseigen(ma, fx$ld, threads = 2, threshB = if (is.na(th)) NULL else th))$snpRes$beta_std)
})

test_that("ab.bin: a block whose SNPs are all in A (no B components) gives zero U alpha there, not recycled values", {
  G <- sim_geno(600, 300, seed = 3)
  pre <- file.path(tempdir(), "genoabj2"); write_bed(G, pre)
  br <- tempfile(fileext = ".txt")
  writeLines(c("chr \t start \t stop", "chr1 \t 1 \t 150500", "chr1 \t 150500 \t 1000000"), br)
  A <- paste0("rs", c(seq(5, 140, 9), 151:300)); dir <- tempfile("ldabj2")
  suppressMessages(LDbuild(paste0(pre, ".bed"), dir, threads = 2, blockRef = br, minsnp = 100, A = A))
  si <- data.table::fread(file.path(dir, "snp.info")); m <- nrow(si); n <- 5e4
  X <- scale(G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(599); R <- crossprod(X)
  set.seed(4); beta <- numeric(m); beta[c(30, 200)] <- c(0.1, 0.08)
  bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(600))) / sqrt(n); f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  ma <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s, se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n)
  fit <- suppressMessages(sbayeseigen(ma, dir, threads = 2))
  b2 <- si$Block == unique(si$Block)[2]
  ab <- fit$par$abc[set == "A"]
  expect_equal(fit$snpRes$beta_std[b2], ab$beta_std[match(si$ID[b2], ab$SNP)], tolerance = 1e-12)
})

test_that("ab.bin: imputed input with se = 0, r2 = NA and a block without typed SNPs gives a finite fit", {
  fx <- ab_fixture(); si <- fx$si; m <- nrow(si); n <- 5e4
  X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1); R <- crossprod(X)
  set.seed(9); beta <- numeric(m); beta[c(30, 220)] <- c(0.09, 0.06)
  bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
  f <- si$A1Freq; s <- sqrt(2 * f * (1 - f))
  # already-imputed input (snp.info order, r2 column; r2 = 1 typed)
  im <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = bh / s,
                               se = sqrt((1 - bh^2) / n) / s, p = 0.5, N = n, r2 = ifelse(runif(m) < 0.9, 1, 0.9))
  b2 <- si$Block == unique(si$Block)[2]
  im$r2[b2 & im$r2 == 1] <- 0.95          # block 2: no typed SNP
  im$se[3] <- 0; im$r2[5] <- NA; im$N[7] <- NA; im$b[9] <- Inf
  expect_message(fit <- sbayeseigen(im, fx$ld, threads = 2), "3 rows with a non-finite")
  expect_true(is.finite(fit$par$Vg)); expect_true(all(is.finite(fit$snpRes$beta_std)))
})

test_that("Z-score input (SNP A1 A2 freq Z N) fits like the equivalent COJO input; beta is per dosage with reference freq", {
  fx <- ab_fixture(); si <- fx$si; m <- nrow(si); n <- 5e4
  X <- scale(fx$G[, match(si$ID, paste0("rs", 1:300))]) / sqrt(nrow(fx$G) - 1); R <- crossprod(X)
  set.seed(9); beta <- numeric(m); beta[c(30, 220)] <- c(0.09, 0.06)
  bh <- drop(R %*% beta) + drop(crossprod(X, rnorm(nrow(X)))) / sqrt(n)
  z <- bh * sqrt(n) / sqrt(1 - bh^2); f <- si$A1Freq; k <- sort(sample(m, 0.9 * m))
  zd <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, A1freq = f, Z = z, N = n)[k]
  s <- 1 / sqrt(2 * f * (1 - f) * (n + z^2))
  cj <- data.table::data.table(SNP = si$ID, A1 = si$A1, A2 = si$A2, freq = f, b = z * s, se = s,
                               p = 0.5, N = n)[k]
  fz <- suppressMessages(sbayeseigen(zd, fx$ld, threads = 2)); fc <- suppressMessages(sbayeseigen(cj, fx$ld, threads = 2))
  expect_equal(fz$snpRes$beta_std, fc$snpRes$beta_std, tolerance = 1e-6)
  expect_equal(fz$snpRes$beta, fz$snpRes$beta_std / sqrt(2 * f * (1 - f)), tolerance = 1e-8)
  # no freq column: the reference's is used; flipped alleles keep the sign right
  zf <- data.table::copy(zd)[, A1freq := NULL]; fl <- 1:20
  zf[fl, `:=`(A1 = A2, A2 = A1, Z = -Z)]
  # names case-insensitive; Z together with b/se: COJO is used
  zl <- data.table::copy(zd); data.table::setnames(zl, c("snp", "a1", "a2", "A1FREQ", "zscore", "n"))
  expect_equal(suppressMessages(sbayeseigen(zl, fx$ld, threads = 2))$snpRes$beta, fz$snpRes$beta, tolerance = 1e-8)
  expect_false(SBayesEigen:::.is_zinput(cbind(cj, Z = 1)))
  # beta_ref: COJO output on the reference-freq dosage scale too (Var_y = 1 here)
  fr <- suppressMessages(sbayeseigen(cj, fx$ld, threads = 2, beta_ref = TRUE))
  expect_equal(fr$snpRes$beta, fz$snpRes$beta, tolerance = 1e-6)
  tz <- suppressMessages(tidy(zf, fx$ld)); td <- suppressMessages(tidy(zd, fx$ld)); i <- match(td$SNP, tz$SNP)
  sg <- ifelse(td$SNP %in% zd$SNP[fl], -1, 1)
  expect_equal(tz$b[i], sg * td$b, tolerance = 1e-4)
})

test_that("abj_vi: B hyperparameters learn when the moment start of s2B is negative (no absorbing s2B ~ 0)", {
  # old start (moment, floored at ~0) gave pi_B exactly uniform, Vg_B ~ 5e-11 and a stop at iteration 4
  set.seed(1); nb <- 50; k <- 200; n <- 1e4
  blk <- lapply(1:nb, function(b) list(XA = matrix(0, 0, 0), Cm = matrix(0, 0, k), sl = sqrt(rexp(k, 1 / 3) + 1e-3),
                                        XCa = matrix(0, 0, 0), XCb = matrix(0, k, 0)))
  sl <- unlist(lapply(blk, `[[`, "sl")); M <- length(sl)
  al <- rnorm(M) * (runif(M) < 0.1); al <- al * sqrt(0.1 / sum((sl * al)^2))
  w <- sl * al + rnorm(M) * sqrt(0.55 / n); p <- rep(n / 0.9, M)   # noise overstated: moment start < 0
  expect_lt(sum((p * w^2 - 1) * p * sl^2), 0)
  f <- SBayesEigen:::abj_vi(w, p, blk, threads = 2, h2p = 0.1)
  expect_gt(max(abs(f$piB - 0.25)), 0.05)
  expect_gt(f$Vg_B, 0.005)
})
