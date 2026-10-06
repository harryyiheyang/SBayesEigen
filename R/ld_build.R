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
#' @param blockRef Block file with columns Block, Chrom, StartBP, EndBP; default is the
#'   SBayesRC 4cM blocks on GRCh37 (\code{ref4cM_v37.pos}).
#' @param snps Optional SNP filter: a character vector of IDs, or a file with an ID per
#'   line or with columns SNP, A1, A2 (e.g. a GWAS \code{.ma} file; alleles must then
#'   match in either orientation, as in \code{SBayesRC::LDstep1}).
#' @param thresh Proportion of the positive eigenvalue mass kept per block.
#' @param threads Blocks processed in parallel (each block single-threaded; with a threaded
#'   OpenBLAS set \code{OPENBLAS_NUM_THREADS=1}). Memory per thread is about
#'   \eqn{8 m^2 (1 + k/m)} bytes for the largest block of \eqn{m} SNPs.
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
LDbuild <- function(geno, out, thresh = 0.995, threads = 4, snps = NULL, blockRef = NULL) {
  t0 <- proc.time()[[3]]
  if (is.null(blockRef)) blockRef <- system.file("extdata", "ref4cM_v37.pos", package = "SBayesEigen")
  pos <- fread(blockRef)
  if (!all(c("Block", "Chrom", "StartBP", "EndBP") %in% names(pos)))
    stop("blockRef needs columns Block, Chrom, StartBP, EndBP")
  pos[, Chrom := as.integer(as.character(Chrom))]
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
  bl <- v[, .(m = .N), by = .(Block, file)]
  if (anyDuplicated(bl$Block)) stop("a block spans several genotype files")
  mx <- max(bl$m)
  message(nrow(v), " variants in ", nrow(bl), " blocks (largest m = ", mx, "), ", threads, " threads, ",
          "about ", format(signif(threads * 16 * mx^2 / 2^30, 2)), " GB peak memory")

  # ---- per file: LD, eigen, write blocks ----
  if (length(list.files(out, pattern = "^block\\d+\\.eigen\\.bin$")))
    stop(out, " already holds block*.eigen.bin files; use an empty folder")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  res <- vector("list", nrow(bl))
  for (f in unique(bl$file)) {
    ib <- which(bl$file == f)
    idx <- split(v$fidx, factor(v$Block, levels = bl$Block))[ib]
    outs <- file.path(out, paste0("block", bl$Block[ib], ".eigen.bin"))
    tf <- proc.time()[[3]]
    res[ib] <- ld_build_cpp(gf$path[f], gf$type[f] == "pgen", gf$n[f], vl[[f]]$nallele, idx, outs,
                            thresh, threads)
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
  einfo <- einfo[m > 0]   # blocks left without polymorphic variants have no eigen file
  fwrite(einfo, file.path(out, "eigen.info"), sep = "\t")
  message(sprintf("LD reference written to %s (%d SNPs, %d blocks) in %.1f s", out, nrow(v), nrow(einfo),
                  proc.time()[[3]] - t0))
  invisible(einfo)
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
