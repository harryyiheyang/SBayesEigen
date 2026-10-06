# Rscript ld_build.R geno=<.bed/.pgen, prefix with {CHR}, or directory> out=<LD dir> [threads=4] [thresh=0.995]
#   [snps=<file>] [blockRef=<pos file>]
a <- commandArgs(TRUE); arg <- setNames(sub(".*?=", "", a), sub("=.*", "", a))
opt <- function(k, def) if (is.na(arg[k])) def else arg[[k]]
SBayesEigen::LDbuild(geno = arg[["geno"]], out = arg[["out"]], blockRef = opt("blockRef", NULL),
                     snps = opt("snps", NULL), thresh = as.numeric(opt("thresh", 0.995)),
                     threads = as.integer(opt("threads", 4)))
