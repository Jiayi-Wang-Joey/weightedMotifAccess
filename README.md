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
  both for differential motif accessibility analysis and for
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

`weightedMotifAccess` is not on Bioconductor yet and must be installed from GitHub first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("remotes")
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

se <- getWeightedCounts(
    # one entry per sample; see below for the other accepted forms
    atac = c(sample1 = "sample1.fragments.tsv.gz",
             sample2 = "sample2.fragments.tsv.gz"),
    ranges = peaks,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    fragWeight = TRUE,
    peakWeight = TRUE
)
se
```

### Fragment input

The `atac` argument takes one entry per sample, in whichever form is most
convenient:

- a named character vector of **file paths**: `.bam` (paired-end; fragments
  are assembled from the read pairs and ATAC-shifted), `.bed`, `.tsv` or
  `.txt` (optionally `.gz`-compressed; the first three columns are read as
  chromosome, start and end, and any leading `#` comment lines are skipped),
  or `.rds`;
- a named list of **`data.table`s** (or `data.frame`s) with
  `seqnames`/`start`/`end` columns (`chr` is accepted instead of `seqnames`);
- a named list of **`GRanges`**, or a `GRangesList`.

The names are used as sample names; if you don't provide any, they are
derived from the file names. If you only have one sample, you can skip the
list and pass the object directly:

```r
se <- getWeightedCounts(atac = myFragmentsDataTable, ranges = peaks,
                        genome = BSgenome.Hsapiens.UCSC.hg38)
```

Only the coordinate columns are used; any further columns are ignored.

`weightedInsertions()` also takes one entry per sample, through its `lf`
argument, but is stricter: because it streams over the genome chromosome by
chromosome, it only accepts **paths to indexed files** — either indexed
`.bam` files, or tabix-indexed fragment files (see `Rsamtools::bgzip()` and
`Rsamtools::indexTabix()`).

## Getting help

Questions and bug reports are welcome via
[GitHub issues](https://github.com/Jiayi-Wang-Joey/weightedMotifAccess/issues).

## Authors

- Jiayi Wang (jiayi.wang2@uzh.ch)
- Emanuel Sonder
- Pierre-Luc Germain

## License

Artistic-2.0
