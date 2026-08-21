suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(betterChromVAR)
})

# Load shared example data once for all test files. Note this must go into
# this file's own environment (which testthat shares with the test files):
# under R CMD check the parent is the locked package namespace.
data(NR3C1example, package="weightedMotifAccess", envir=environment())
