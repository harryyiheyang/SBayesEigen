# SBayesEigen

Fast SBayes-style polygenic risk scores from GWAS summary statistics and eigen-decomposed LD,
fitted by variational EM instead of MCMC. LD folders use the SBayesRC layout
(`snp.info`, `block<b>.eigen.bin`), so official SBayesRC LD references work as well.

## Installation

```r
# install.packages("remotes")
remotes::install_github("harryyiheyang/SBayesEigen")
```

Building from source needs a C++17 compiler and LAPACK/BLAS development files (on
Debian/Ubuntu: `liblapack-dev libblas-dev gfortran`). A multithreaded BLAS (OpenBLAS or MKL)
is not required, because `LDbuild()` parallelises over blocks.

## Pipeline

```r
library(SBayesEigen)

# 1. LD reference from PLINK genotypes (once). geno can be one .bed/.pgen, a prefix with
#    {CHR}, a directory of per-chromosome files, or one ..._chr1 file whose chr2..chr22
#    siblings exist. Autosomes only.
LDbuild("ukb_imp_chr{CHR}.pgen", "ukbEUR_LD", threads = 16)

# 2. tidy -> impute -> LDSC -> VI -> beta, in memory; each eigen file is read for imputation
#    and the rotation together, then once more for beta
fit <- sbayeseigen("trait.ma", "ukbEUR_LD", out = "trait", threads = 8)
fit$par$Vg
head(fit$snpRes)   # SNP A1 A2 Block beta beta_std, every snp.info SNP

# several traits: each eigen file is still read twice in total; results equal single-trait runs.
# out is a directory; writes prs/<trait>/<trait>_sbeigen.snpRes, .par.rds and .log
fit <- sbayeseigen(c(LDL = "ldl.ma", HDL = "hdl.ma", TG = "tg.ma"), "ukbEUR_LD", out = "prs", threads = 8)
```

The steps also run on their own, like SBayesRC, and `sbayeseigen()` accepts an imputed file:

| SBayesRC | SBayesEigen |
|---|---|
| `tidy(mafile, LDdir, output)` | `tidy(ma, ld, out)` |
| `impute(mafile, LDdir, output)` | `impute(ma, ld, out)` |
| `sbrc(mafile, LDdir, outPrefix, annot)` | `sbayeseigen(ma, ld, out)` (no annotations) |

`ldscore.txt` is written by `LDbuild()`. For an SBayesRC LD folder without it, `sbayeseigen()`
computes the LD scores from the eigen files during its first pass, at no extra read.
The output has posterior mean effects only; there are no PIPs.

Command line:

```sh
S=$(Rscript -e 'cat(system.file("scripts", package = "SBayesEigen"))')
Rscript $S/ld_build.R geno=ukb_imp_chr{CHR}.pgen out=ukbEUR_LD threads=16
Rscript $S/run_vi.R ma=trait.ma ld=ukbEUR_LD out=trait threads=8
Rscript $S/run_vi.R ma=ldl.ma,hdl.ma,tg.ma ld=ukbEUR_LD out=prs threads=8
```

## Methods

- **LDbuild:** for each block, reads the variants straight from the BED/PGEN file and
  computes correlations with a bit-plane popcount kernel (from CppMatrix). It then
  tridiagonalises once, computes all eigenvalues, and forms only the k eigenvectors needed
  to keep `thresh` (0.995) of the variance (LAPACK `dsytrd`/`dsterf`/`dstemr`/`dormtr`).
  Blocks run in parallel, largest first. It also writes `snp.info`, `ldm.info`,
  `eigen.info` and `ldscore.txt`, which holds one LD score per SNP,
  sum_j [r2_ij - (1 - r2_ij)/(n - 2)], unbiased for the population LD score. The counted allele is A1: .bim column 5 for BED, REF for
  PGEN.
- **tidy:** applies the same QC as `SBayesRC::tidy`, and its output matches SBayesRC's row
  for row. It reads the files whole, splits lines in parallel and matches SNPs with a hash
  table.
- **impute:** Yihe Yang's block-parallel imputation, which works directly in the eigen space
  (impute.sbrc.fast).
- **Model:** bhat = z / sqrt(N + z^2) as in SBayesRC; w = D^{-1/2} U' bhat = D^{1/2} alpha + e and beta = U alpha. alpha has a
  six-class prior: an exact zero plus variance ratios 1:10:100:1000:10000.
  - LD score regression on the typed (not imputed) SNPs gives the centre of a
    scaled-inverse-chi-squared prior (4 df) on Vg; h2 is floored at 0.01.
  - Variational EM with SQUAREM. The likelihood is diagonal in the eigen basis, so the
    mean-field posterior is exact. It stops when Vg changes by less than `tol` (1e-4) twice.
  - Residual variance per eigen component is ve0 + kappa / lambda. A reference LD panel that does
    not match the GWAS adds noise that grows as lambda shrinks; kappa absorbs it.
  - With `kappa = "mom"` (the default), ve0 (within [0.9, 1.2]) and kappa (>= 0) come from a moment fit.
    The fit uses components with lambda < 1 in 100 equal-count bins, with the LDSC signal subtracted.
    When the LD matches, kappa is about 0.
  - A numeric `kappa` is used as given, with ve0 = `ve`. `kappa = 0` gives a constant `ve`
    (default 1, or `"ldsc"` for the LDSC intercept clamped to [0.9, 2]).

## License

GPL (>= 3). `src/pgenlib/` holds unmodified pgenlib sources from plink-ng
(Christopher Chang, LGPL >= 3). `src/simde/` holds MIT-licensed SIMDe headers.
