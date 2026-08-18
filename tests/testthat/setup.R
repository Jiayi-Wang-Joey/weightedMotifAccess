suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(betterChromVAR)
})

# Load shared example data once for all test files
load(system.file("data", "NR3C1example.RData", package="weightedMotifAccess"),
     envir=parent.env(environment()))
