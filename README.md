# weightedMotifAccess

<!-- badges: start -->
![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)
<!-- badges: end -->

`weightedMotifAccess` provides two complementary approaches for
differential motif accessibility analysis from ATAC-seq data:

- **The weighted model** corrects GC-content and fragment-length biases by
  assigning weights to individual fragments, followed by cyclic loess
  normalization of peak weights to correct for enrichment bias. This
  produces an unbiased peak-by-sample accessibility matrix that can be used
  both for differential motif accessibility analysis (via
  [betterChromVAR](https://github.com/ETHZ-INS/betterChromVAR)) and for
  broader downstream analyses at the peak level, such as differential
  accessibility testing or bias-corrected coverage profiles for
  visualization.
- **The insertion model** incorporates Tn5 footprint patterns around motif
  sites, assigning bias-corrected, weighted insertion scores to each motif
  match. This enables discovery and characterization of interactions
  between pairs of motifs — testing whether the accessibility response at
  co-occurring motif pairs differs from that at either motif alone, and how
  this depends on the distance between the motifs.

The package integrates with standard Bioconductor classes (`GRanges`,
`SummarizedExperiment`) throughout.

## Installation

`weightedMotifAccess` depends on `betterChromVAR` and `epiwraps`, which are
not on Bioconductor and must be installed from GitHub first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("remotes")
remotes::install_github("ETHZ-INS/betterChromVAR")
remotes::install_github("ETHZ-INS/epiwraps")

remotes::install_github("Jiayi-Wang-Joey/weightedMotifAccess")
```

## Overview

The package covers two main workflows:

- **Bias-corrected accessibility quantification**: `getWeightedCounts()`
  counts (weighted) ATAC-seq fragments in peak regions, correcting for
  fragment-length and GC-content biases; `weightedInsertions()` computes
  bias-corrected, weighted Tn5-insertion counts around motif matches to
  derive per-motif activity scores.
- **Motif interaction analysis**: `getInteractionsDeviations()` and
  `discoverMotifInteractions()` test whether the accessibility response at
  a "bait" motif differs depending on co-occurrence with other motifs, and
  `exploreMotifInteraction()` characterizes how the interaction between a
  given pair of motifs depends on the distance between them.

Additionally, `genomicRangesMapping()` is a general-purpose helper for
mapping and aggregating other genomic-coordinate data (e.g. motif scores,
fragment counts) onto a set of reference ranges.

## Quick example

```r
library(weightedMotifAccess)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)

# example peaks bundled with the package
data(NR3C1example, package = "weightedMotifAccess")
peaks <- rowRanges(peakSE)

# ATAC-seq fragments per sample, as a data.table with
# seqnames/start/end columns (or paths to fragment/bam files)
atacFrag <- list(sample1 = myFragmentsDataTable)

se <- getWeightedCounts(
    files = NULL,
    atacFrag = atacFrag,
    ranges = peaks,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    fragWeight = TRUE,
    peakWeight = TRUE
)
se
```

## Getting help

Questions and bug reports are welcome via
[GitHub issues](https://github.com/Jiayi-Wang-Joey/weightedMotifAccess/issues).

## Authors

- Jiayi Wang (jiayi.wang2@uzh.ch)
- Emanuel Sonder
- Pierre-Luc Germain

## License

Artistic-2.0
