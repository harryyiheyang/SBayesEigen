# synthetic genotypes with LD: thresholded AR(1) latent normals, n samples x m SNPs (dosage of A1)
sim_geno <- function(n, m, rho = 0.9, seed = 1) {
  set.seed(seed)
  hap <- function() {
    Z <- matrix(0, n, m); Z[, 1] <- rnorm(n)
    for (j in 2:m) Z[, j] <- rho * Z[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
    sweep(Z, 2, qnorm(runif(m, 0.1, 0.9)), ">") * 1L
  }
  hap() + hap()
}

# PLINK 1 .bed/.bim/.fam; A1 (bim column 5) is the counted allele
write_bed <- function(G, prefix, chr = 1, pos = seq_len(ncol(G)) * 1000) {
  n <- nrow(G); m <- ncol(G); nb <- ceiling(n / 4)
  code <- c(3L, 2L, 0L)   # dosage 0, 1, 2 -> 11, 10, 00
  raw <- unlist(lapply(seq_len(m), function(j) {
    x <- c(code[G[, j] + 1L], rep(0L, nb * 4 - n))
    x <- matrix(x, 4)
    as.raw(x[1, ] + 4L * x[2, ] + 16L * x[3, ] + 64L * x[4, ])
  }))
  writeBin(c(as.raw(c(0x6c, 0x1b, 0x01)), raw), paste0(prefix, ".bed"))
  data.table::fwrite(data.table::data.table(chr, paste0("rs", seq_len(m)), 0, pos, "A", "G"),
                     paste0(prefix, ".bim"), sep = "\t", col.names = FALSE)
  data.table::fwrite(data.table::data.table(seq_len(n), seq_len(n), 0, 0, 0, -9),
                     paste0(prefix, ".fam"), sep = " ", col.names = FALSE)
}

read_eig <- function(f) {
  con <- file(f, "rb"); on.exit(close(con))
  mk <- readBin(con, "integer", 2, size = 4)
  hd <- readBin(con, "numeric", 2, size = 4)
  lam <- readBin(con, "numeric", mk[2], size = 4)
  U <- matrix(readBin(con, "numeric", mk[1] * mk[2], size = 4), mk[1])
  list(lam = lam, U = U)
}

# shared fixture: 2 blocks of 150 SNPs, 600 samples
make_ld <- function(dir = tempfile("ld")) {
  G <- sim_geno(600, 300)
  pre <- file.path(tempdir(), "geno"); write_bed(G, pre)
  br <- tempfile(fileext = ".pos")
  data.table::fwrite(data.table::data.table(Block = 1:2, Chrom = 1, StartBP = c(1, 150500), EndBP = c(150500, 1e6)), br, sep = " ")
  suppressMessages(LDbuild(paste0(pre, ".bed"), dir, threads = 2, blockRef = br))
  list(G = G, ld = dir)
}

# pure-R EM map for the same model (ve fixed, prior on tau = sigma2 A0 always on): reference for em_step_cpp
em_step_r <- function(w, c, n, gamma, pi, sigma2, ve, s2p, nu, A0) {
  K <- length(gamma); nz <- gamma > 0; c2 <- c^2; J <- length(w)
  V <- outer(c2, gamma * sigma2) + ve / n
  lp <- sweep(-0.5 * log(2 * base::pi * V) - 0.5 * w^2 / V, 2, log(pi), "+")
  mx <- apply(lp, 1, max); lse <- mx + log(rowSums(exp(lp - mx)))
  r <- exp(lp - lse)
  S <- M <- matrix(0, J, K)
  S[, nz] <- 1 / (outer(n * c2 / ve, rep(1, sum(nz))) + outer(rep(1, J), 1 / (gamma[nz] * sigma2)))
  M[, nz] <- S[, nz] * (n * c * w / ve)
  tau <- sigma2 * A0
  list(ll = sum(lse), obj = sum(lse) - (nu / 2 + 1) * log(tau) - s2p / (2 * tau),
       vg = sum(c2 * rowSums(r * (M^2 + S))), pi = pmax(colSums(r) / J, 1e-300),
       sigma2 = (sum(sweep(r[, nz] * (M[, nz]^2 + S[, nz]), 2, gamma[nz], "/")) + s2p / A0) / (sum(r[, nz]) + nu + 2))
}
