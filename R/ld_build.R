#' Build an eigen-decomposed LD reference from PLINK genotypes
#'
#' Replaces the GCTB route of SBayesRC (\code{LDstep1}-\code{LDstep4}). Correlations are
#' computed per block straight from the genotype file with a bit-plane popcount kernel,
#' and only the \eqn{k} leading eigenvectors that explain \code{thresh} of the variance
#' are formed (one tridiagonalisation, all eigenvalues, then \eqn{k} eigenvectors). Blocks
#' run in parallel. Autosomes only; chromosome X and other non-numeric chromosomes are skipped.
#'
#' @param geno Genotypes: a PLINK 1 \code{.bed} or PLINK 2 \code{.pgen} file or prefix, a
#'   prefix containing \code{{CHR}} (expanded to 1-22), or a directory holding one file per
#'   chromosome. A single file named like \code{..._chr1} whose siblings \code{chr2},
#'   \code{chr3}, ... exist is treated as the per-chromosome set. BED is used when both
#'   formats exist. PGEN needs a plain-text \code{.pvar}.
#' @param out Output LD folder (created).
#' @param blockRef Block file: columns Block, Chrom, StartBP, EndBP, or an interval file with a header and
#'   three columns chr, start, stop such as LDetect (\code{chr1}/\code{1} both accepted; blocks numbered
#'   1, 2, ... in genome order). Intervals are [start, stop) in the genome build of \code{geno}. Default is
#'   the SBayesRC 4cM blocks on GRCh37 (\code{ref4cM_v37.pos}).
#' @param snps Optional SNP filter: a character vector of IDs, or a file with an ID per
#'   line or with columns SNP, A1, A2 (e.g. a GWAS \code{.ma} file; alleles must then
#'   match in either orientation, as in \code{SBayesRC::LDstep1}).
#' @param thresh Proportion of the positive eigenvalue mass kept per block (of the Schur complement of A in
#'   AB mode).
#' @param minsnp Blocks with fewer variants are merged into the smaller neighbouring block on the same
#'   chromosome (keeping the first block's number) until every block has at least \code{minsnp}; 0 turns
#'   merging off.
#' @param A Optional fixed SNP set A (character vector of IDs or a file with one ID per line): every block is
#'   then written as \code{block<b>.ab.bin}: \eqn{R_{AA}} as is, the B part eigen-decomposed after an exact
#'   Schur complement \eqn{S = R_{BB} - R_{BA} R_{AA}^+ R_{AB}} (generalized inverse keeping eigenvalues
#'   above \code{tolA} times the largest), and \eqn{R_{BA}}, so that R can be rebuilt; see
#'   \code{read_ab()} in the package source for the layout.
#' @param tolA Relative eigenvalue cut of the generalized inverse of \eqn{R_{AA}} (R_AA is stored as float; sbayeseigen also drops A components below 1e-6 times the largest).
#' @param threads Threads. The largest blocks, whose single-threaded cost would set the wall time, run one at a
#'   time on all threads (correlation kernel and, with OpenBLAS, MKL or FlexiBLAS, the LAPACK steps; the BLAS
#'   thread count is set internally and restored); the other blocks run in parallel, one thread each.
#' @param mem Memory cap in GB for the blocks in flight; a block of \eqn{m} SNPs needs about \eqn{24 m^2}
#'   bytes at its peak. Default: 80\% of the smallest of physical memory, the cgroup limit (e.g. a Slurm
#'   allocation) and \code{SLURM_MEM_PER_NODE}. \code{Inf} turns the cap off.
#' @return Invisibly, the \code{eigen.info} table. Writes \code{block<b>.eigen.bin},
#'   \code{snp.info}, \code{ldm.info}, \code{ldscore.txt} and \code{eigen.info}.
#'   \code{ldscore.txt} holds one column, the LD score unbiased for the population
#'   \eqn{r^2}: \eqn{\sum_j [r^2_{ij} - (1 - r^2_{ij})/(n - 2)]} with \eqn{n} the reference size.
#' @details The correlations are of the counted allele, recorded as A1 in \code{snp.info}:
#'   .bim column 5 for BED, REF for PGEN. Missing genotypes take the variant's median;
#'   monomorphic variants are dropped.
#' @examples
#' \dontrun{
#' LDbuild("ukb_imp_chr{CHR}.bed", "ukbEUR_LD", threads = 16)
#' LDbuild("/data/ukb_pgen/", "ukbEUR_LD", snps = "trait.ma", threads = 16)
#' }
#' @export
LDbuild <- function(geno, out, thresh = 0.995, threads = 4, snps = NULL, blockRef = NULL, minsnp = 100,
                    A = NULL, tolA = 1e-6, mem = NULL) {
  t0 <- proc.time()[[3]]
  if (is.null(blockRef)) blockRef <- system.file("extdata", "ref4cM_v37.pos", package = "SBayesEigen")
  pos <- .read_blocks(blockRef)
  pos <- pos[Chrom %in% 1:22]
  setorder(pos, Chrom, StartBP)
  gf <- .geno_files(geno)
  message(nrow(gf), " genotype file(s), ", toupper(gf$type[1]), ": ", paste(basename(gf$prefix), collapse = ", "))

  # ---- variants: autosomes, biallelic, unique IDs, optional filter, block assignment ----
  vl <- lapply(seq_len(nrow(gf)), function(f) .read_variants(gf$prefix[f], gf$type[f])[, file := f])
  v <- rbindlist(vl)
  message(nrow(v), " variants in genotype files")
  v <- v[Chrom %in% 1:22 & nallele == 2]
  dup <- v$ID[duplicated(v$ID)]
  if (length(dup)) {
    v <- v[!ID %in% dup]
    message(length(unique(dup)), " duplicated variant IDs removed")
  }
  if (!is.null(snps)) v <- .filter_snps(v, snps)
  v[, Block := NA_integer_]
  for (ch in unique(v$Chrom)) {
    p <- pos[Chrom == ch]
    if (!nrow(p)) next
    ii <- v$Chrom == ch
    j <- findInterval(v$PhysPos[ii], p$StartBP)
    ok <- j > 0
    ok[ok] <- v$PhysPos[ii][ok] < p$EndBP[j[ok]]
    v[which(ii), Block := ifelse(ok, p$Block[pmax(j, 1L)], NA_integer_)]
  }
  v <- v[!is.na(Block)]
  if (!nrow(v)) stop("no variants fall inside the blocks of blockRef")
  setorder(v, Block, file, fidx)
  if (minsnp > 0) v[, Block := .merge_blocks(Block, Chrom, minsnp)]
  ab <- !is.null(A)
  if (ab) {
    if (length(A) == 1 && file.exists(A)) A <- fread(A, header = FALSE)[[1]]
    v[, isA := ID %in% A]
    message(sum(v$isA), " of ", length(unique(A)), " A SNPs found; blocks written as block<b>.ab.bin")
  }
  bl <- v[, .(m = .N), by = .(Block, file)]
  if (anyDuplicated(bl$Block)) stop("a block spans several genotype files")
  mx <- max(bl$m)
  if (is.null(mem)) mem <- 0.8 * .mem_limit() / 2^30
  peak <- 24 * mx^2 / 2^30
  message(nrow(v), " variants in ", nrow(bl), " blocks (largest m = ", mx, "), ", threads, " threads, ",
          "memory cap ", if (is.finite(mem)) paste(format(signif(mem, 3)), "GB") else "none",
          " (largest block about ", format(signif(peak, 2)), " GB)")
  if (is.finite(mem) && peak > mem) warning("the largest block needs about ", signif(peak, 2), " GB, above mem")

  # ---- per file: LD, eigen, write blocks ----
  if (length(list.files(out, pattern = "^block\\d+\\.(eigen|ab)\\.bin$")))
    stop(out, " already holds block*.eigen.bin or block*.ab.bin files; use an empty folder")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  res <- vector("list", nrow(bl))
  for (f in unique(bl$file)) {
    ib <- which(bl$file == f)
    fb <- factor(v$Block, levels = bl$Block)
    idx <- split(v$fidx, fb)[ib]
    am <- if (ab) lapply(split(v$isA, fb)[ib], as.integer) else list()
    outs <- file.path(out, paste0("block", bl$Block[ib], if (ab) ".ab.bin" else ".eigen.bin"))
    tf <- proc.time()[[3]]
    res[ib] <- ld_build_cpp(gf$path[f], gf$type[f] == "pgen", gf$n[f], vl[[f]]$nallele, idx, outs,
                            thresh, threads, am, tolA, if (is.finite(mem)) mem * 2^30 else 0)
    message(sprintf("  %s: %d blocks, %.1f s", basename(gf$prefix[f]), length(ib), proc.time()[[3]] - tf))
  }

  # ---- snp.info, ldm.info, ldscore.txt, eigen.info ----
  keep <- unlist(lapply(res, `[[`, "keep"))
  v[, `:=`(A1Freq = unlist(lapply(res, `[[`, "freq")), N = unlist(lapply(res, `[[`, "nobs")))]
  if (any(!keep)) message(sum(!keep), " monomorphic variants dropped")
  v <- v[keep]
  v[, Index := 0:(.N - 1)]
  fwrite(v[, .(Chrom, ID, Index, GenPos, PhysPos, A1, A2, A1Freq, N, Block)], file.path(out, "snp.info"), sep = "\t")
  fwrite(v[, .(Chrom = Chrom[1], StartSnpIdx = Index[1], StartSnpID = ID[1], EndSnpIdx = Index[.N],
               EndSnpID = ID[.N], NumSnps = .N), by = Block], file.path(out, "ldm.info"), sep = "\t")
  fwrite(data.table(SNP = v$ID, Block = v$Block, ldscore = unlist(lapply(res, `[[`, "ldscore"))),
         file.path(out, "ldscore.txt"), sep = "\t")
  einfo <- data.table(Block = bl$Block, m = vapply(res, `[[`, 0L, "m"), k = vapply(res, `[[`, 0L, "k"),
                      trR = vapply(res, `[[`, 0, "trR"), sumLambda = vapply(res, `[[`, 0, "sumLambda"),
                      relerr_F = vapply(res, `[[`, 0, "relerr_F"))
  if (ab) einfo[, `:=`(mA = vapply(res, `[[`, 0L, "mA"), rankA = vapply(res, `[[`, 0L, "rankA"))]
  einfo <- einfo[m > 0]   # blocks left without polymorphic variants have no eigen file
  fwrite(einfo, file.path(out, "eigen.info"), sep = "\t")
  message(sprintf("LD reference written to %s (%d SNPs, %d blocks) in %.1f s", out, nrow(v), nrow(einfo),
                  proc.time()[[3]] - t0))
  invisible(einfo)
}

# memory available to this process in bytes: physical memory, the cgroup (v2 or v1) limit and a Slurm
# allocation, whichever is smallest; Inf when none can be read
.mem_limit <- function() {
  rd <- function(f) {
    x <- suppressWarnings(tryCatch(as.numeric(readLines(f, n = 1, warn = FALSE)), error = function(e) NA))
    if (length(x) == 1 && is.finite(x) && x > 0 && x < 2^60) x else Inf
  }
  lim <- Inf
  if (file.exists("/proc/meminfo")) {
    mi <- readLines("/proc/meminfo", warn = FALSE)
    tot <- suppressWarnings(as.numeric(sub("^MemTotal:\\s+(\\d+) kB.*", "\\1", grep("^MemTotal:", mi, value = TRUE))))
    if (length(tot) == 1 && is.finite(tot)) lim <- tot * 1024
  }
  if (file.exists("/proc/self/cgroup")) {
    for (l in readLines("/proc/self/cgroup", warn = FALSE)) {
      f <- strsplit(l, ":", fixed = TRUE)[[1]]
      if (length(f) < 3) next
      path <- paste(f[-(1:2)], collapse = ":")
      if (f[1] == "0" && f[2] == "") lim <- min(lim, rd(file.path("/sys/fs/cgroup", path, "memory.max")))
      if ("memory" %in% strsplit(f[2], ",", fixed = TRUE)[[1]])
        lim <- min(lim, rd(file.path("/sys/fs/cgroup/memory", path, "memory.limit_in_bytes")))
    }
  }
  sl <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")))
  if (is.finite(sl) && sl > 0) lim <- min(lim, sl * 2^20)
  lim
}

# block definitions: Block, Chrom, StartBP, EndBP (autosomes, sorted); also a chr/start/stop interval file
.read_blocks <- function(f) {
  pos <- fread(f, strip.white = TRUE)
  if (!all(c("Block", "Chrom", "StartBP", "EndBP") %in% names(pos))) {
    if (ncol(pos) < 3) stop("blockRef needs columns Block, Chrom, StartBP, EndBP or chr, start, stop")
    pos <- pos[, 1:3]
    setnames(pos, c("Chrom", "StartBP", "EndBP"))
    pos[, Chrom := sub("^chr", "", trimws(as.character(Chrom)), ignore.case = TRUE)]
    pos <- pos[Chrom %in% as.character(1:22)]
    pos[, `:=`(Chrom = as.integer(Chrom), StartBP = as.numeric(StartBP), EndBP = as.numeric(EndBP))]
    setorder(pos, Chrom, StartBP)
    pos[, Block := seq_len(.N)]
  }
  pos[, Chrom := as.integer(as.character(Chrom))]
  pos <- pos[Chrom %in% 1:22]
  setorder(pos, Chrom, StartBP)
  if (pos[, any(utils::head(EndBP, -1) > utils::tail(StartBP, -1)), by = Chrom][, any(V1)])
    stop("blockRef has overlapping blocks")
  pos
}

# per variant (sorted by block, blocks in genome order) the merged block: a block with fewer than minsnp
# variants joins its smaller neighbour on the same chromosome, until all reach minsnp or one block is left
.merge_blocks <- function(block, chrom, minsnp) {
  r <- rle(block)
  ch <- chrom[cumsum(r$lengths)]
  id <- r$values; cnt <- r$lengths
  grp <- seq_along(id)
  for (c0 in unique(ch)) {
    w <- which(ch == c0)
    g <- seq_along(w); n <- cnt[w]
    while (length(n) > 1 && min(n) < minsnp) {
      s <- which.min(n)
      nb <- if (s == 1) 2 else if (s == length(n)) s - 1 else c(s - 1, s + 1)[which.min(n[c(s - 1, s + 1)])]
      lo <- min(s, nb)
      g[g == lo + 1] <- lo; g[g > lo + 1] <- g[g > lo + 1] - 1L
      n[lo] <- n[lo] + n[lo + 1]; n <- n[-(lo + 1)]
    }
    grp[w] <- id[w][match(g, g)]   # first original block of each group
  }
  nm <- sum(grp != id)
  if (nm) message(nm, " blocks with fewer than ", minsnp, " variants merged into neighbours")
  rep(grp, r$lengths)
}

# genotype files: prefix, type (bed/pgen), path to .bed/.pgen, sample count
.geno_files <- function(geno) {
  ext <- "\\.(bed|bim|fam|pgen|pvar|psam)$"
  has <- function(p, t) all(file.exists(paste0(p, if (t == "bed") c(".bed", ".bim", ".fam") else c(".pgen", ".pvar", ".psam"))))
  type_of <- function(p) if (has(p, "bed")) "bed" else if (has(p, "pgen")) "pgen" else NA_character_
  if (dir.exists(geno)) {
    fs <- sub(ext, "", list.files(geno, pattern = "\\.(bed|pgen)$", full.names = TRUE))
    bn <- basename(fs)
    ch <- ifelse(grepl("chr\\d+", bn, ignore.case = TRUE),
                 sub(".*?chr(\\d+).*$", "\\1", bn, ignore.case = TRUE), sub(".*?(\\d+)\\D*$", "\\1", bn))
    ch <- suppressWarnings(as.integer(ch))
    pre <- unique(fs[ch %in% 1:22])
    if (!length(pre)) stop("no per-chromosome .bed/.pgen files for chromosomes 1-22 in ", geno)
  } else {
    p <- sub(ext, "", geno)
    if (grepl("{CHR}", p, fixed = TRUE)) {
      pre <- vapply(1:22, function(c) gsub("{CHR}", c, p, fixed = TRUE), "")
    } else {
      pre <- p
      b <- basename(p)
      if (grepl("chr\\d+", b, ignore.case = TRUE)) {
        tmpl <- file.path(dirname(p), sub("(chr)\\d+", "\\1{CHR}", b, ignore.case = TRUE))
        cand <- vapply(1:22, function(c) gsub("{CHR}", c, tmpl, fixed = TRUE), "")
        if (sum(!is.na(vapply(cand, type_of, ""))) > 1) pre <- cand
      }
    }
  }
  ty <- vapply(pre, type_of, "", USE.NAMES = FALSE)
  if (all(is.na(ty))) stop("no PLINK .bed/.bim/.fam or .pgen/.pvar/.psam found for ", geno)
  pre <- pre[!is.na(ty)]; ty <- ty[!is.na(ty)]
  if (length(unique(ty)) > 1) ty[vapply(pre, has, TRUE, t = "bed")] <- "bed"
  n <- vapply(seq_along(pre), function(i) {
    if (ty[i] == "bed") return(nrow(fread(paste0(pre[i], ".fam"), header = FALSE, select = 1L)))
    s <- readLines(paste0(pre[i], ".psam"))
    sum(!startsWith(s, "#") & nzchar(s))
  }, 0)
  if (length(unique(n)) > 1) stop("genotype files have different sample counts")
  for (i in which(ty == "bed")) {
    con <- file(paste0(pre[i], ".bed"), "rb")
    magic <- readBin(con, "raw", 3L)
    close(con)
    nv <- length(utils::count.fields(paste0(pre[i], ".bim"), quote = "", comment.char = ""))
    if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01))))
      stop(pre[i], ".bed is not a SNP-major PLINK 1 BED file")
    if (file.size(paste0(pre[i], ".bed")) != 3 + nv * ceiling(n[i] / 4))
      stop(pre[i], ".bed size does not match its .bim and .fam")
  }
  data.table(prefix = pre, type = ty, path = paste0(pre, ifelse(ty == "bed", ".bed", ".pgen")), n = n)
}

# variant table in file order; A1 = counted allele (.bim col 5 for BED, REF for PGEN)
.read_variants <- function(prefix, type) {
  if (type == "bed") {
    b <- fread(paste0(prefix, ".bim"), header = FALSE, colClasses = list(character = c(1, 2, 5, 6)))
    v <- data.table(Chrom = suppressWarnings(as.integer(sub("^chr", "", b$V1, ignore.case = TRUE))), ID = b$V2,
                    GenPos = b$V3, PhysPos = b$V4, A1 = b$V5, A2 = b$V6, nallele = 2L)
  } else {
    pv <- paste0(prefix, ".pvar")
    hdr <- 0L
    con <- file(pv, "r")
    repeat {
      l <- readLines(con, n = 1L)
      if (!length(l) || !startsWith(l, "##")) break
      hdr <- hdr + 1L
    }
    close(con)
    if (!length(l) || !startsWith(l, "#CHROM")) stop(pv, ": expected a #CHROM header line")
    b <- fread(pv, skip = hdr, select = c("#CHROM", "POS", "ID", "REF", "ALT"),
               colClasses = list(character = c("#CHROM", "ID", "REF", "ALT")))
    setnames(b, c("CHR", "POS", "ID", "REF", "ALT"))
    v <- data.table(Chrom = suppressWarnings(as.integer(sub("^chr", "", b$CHR, ignore.case = TRUE))), ID = b$ID,
                    GenPos = 0, PhysPos = b$POS, A1 = b$REF, A2 = b$ALT,
                    nallele = 1L + lengths(regmatches(b$ALT, gregexpr(",", b$ALT, fixed = TRUE))) + 1L)
  }
  v[, fidx := 0:(.N - 1)]
  v
}

.filter_snps <- function(v, snps) {
  if (length(snps) == 1 && file.exists(snps)) {
    first <- strsplit(readLines(snps, n = 1L), "[[:space:],]+")[[1]]
    s <- fread(snps, header = "SNP" %in% first)
    if (all(c("SNP", "A1", "A2") %in% names(s))) {
      s <- s[, .(SNP, sA1 = A1, sA2 = A2)]
      v <- merge(v, s[!duplicated(SNP)], by.x = "ID", by.y = "SNP")
      v <- v[(A1 == sA1 & A2 == sA2) | (A1 == sA2 & A2 == sA1)][, c("sA1", "sA2") := NULL]
      message(nrow(v), " variants left after matching SNP and alleles to ", basename(snps))
      return(v)
    }
    snps <- s[[1]]
  }
  v <- v[ID %in% snps]
  message(nrow(v), " variants left after the SNP filter")
  v
}

# blockN.ab.bin (LDbuild with A): int32 m, mA, kB, rankA; float tolA, cutB, sumLambdaB; int32 idxA[mA]
# (0-based within the block); float R_AA upper triangle, column-major packed; float lambdaB[kB];
# float UlB[m * kB] column-major (B rows: eigenvectors Q2 of S = R_BB - R_BA R_AA^+ R_AB; A rows:
# -R_AA^+ R_AB Q2); float R_BA[(m - mA) * mA] column-major (B rows in block order).
# P = UlB diag(lambdaB)^-1/2 gives P'RP = I and R[A, ] P = 0 in sample; R is rebuilt as
# [R_AA, R_AB; R_BA, Q2 diag(lambdaB) Q2' + R_BA R_AA^+ R_AB].
.read_ab <- function(file) {
  h <- file(file, "rb"); on.exit(close(h))
  mk <- readBin(h, integer(), n = 4, size = 4); m <- mk[1]; mA <- mk[2]; kB <- mk[3]
  ct <- readBin(h, numeric(), n = 3, size = 4)
  iA <- readBin(h, integer(), n = mA, size = 4) + 1L
  up <- readBin(h, numeric(), n = mA * (mA + 1) / 2, size = 4)
  lambda <- readBin(h, numeric(), n = kB, size = 4)
  UlB <- matrix(readBin(h, numeric(), n = m * kB, size = 4), m, kB)
  RBA <- matrix(readBin(h, numeric(), n = (m - mA) * mA, size = 4), m - mA, mA)
  RAA <- matrix(0, mA, mA); RAA[upper.tri(RAA, diag = TRUE)] <- up; RAA[lower.tri(RAA)] <- t(RAA)[lower.tri(RAA)]
  list(m = m, mA = mA, kB = kB, rankA = mk[4], tolA = ct[1], cutB = ct[2], sumLambda = ct[3],
       iA = iA, RAA = RAA, lambda = lambda, UlB = UlB, RBA = RBA)
}
