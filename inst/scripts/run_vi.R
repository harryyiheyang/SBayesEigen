# Command-line wrapper:
# Rscript $(Rscript -e 'cat(system.file("scripts/run_vi.R", package = "SBayesEigen"))') \
#   ma=<COJO .ma, raw or imputed> ld=<LD dir> out=<prefix> [threads=4] [ve=1|ldsc] [thresh=0.995]
a <- commandArgs(TRUE); arg <- setNames(sub(".*?=", "", a), sub("=.*", "", a))
opt <- function(k, def) if (is.na(arg[k])) def else arg[[k]]
ve <- opt("ve", "1")
SBayesEigen::sbayeseigen(ma = arg[["ma"]], ld = arg[["ld"]], out = opt("out", "sbayeseigen"),
                         threads = as.integer(opt("threads", 4)), ve = if (ve == "ldsc") ve else as.numeric(ve),
                         thresh = as.numeric(opt("thresh", 0.995)))
