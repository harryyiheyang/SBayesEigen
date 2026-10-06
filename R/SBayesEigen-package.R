#' SBayesEigen: fast PRS from eigen-decomposed LD by variational EM
#'
#' Pipeline: \code{\link{LDbuild}} (LD reference from PLINK genotypes, once),
#' \code{\link{tidy}} and \code{\link{impute}} (summary data), \code{\link{sbayeseigen}}
#' (fit and per-SNP effects). LD folders use the SBayesRC layout (\code{snp.info},
#' \code{block<b>.eigen.bin}), so official SBayesRC LD references work as well.
#'
#' @keywords internal
#' @useDynLib SBayesEigen, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import data.table
"_PACKAGE"

utils::globalVariables(c(".", "..cols", "A1", "A1Freq", "A2", "Block", "ID", "N", "SNP", "b", "beta", "Chrom", "GenPos", "Index", "PhysPos",
                         "StartBP", "fidx", "m", "nallele", "sA1", "sA2",
                         "freq", "frq_ref", "ord", "p", "r2", "rA1", "rA2", "rfreq", "se"))
