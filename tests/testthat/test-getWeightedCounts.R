test_that("getWeightedCounts returns a SummarizedExperiment with correct dims", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)

  peaks <- rowRanges(peakSE)
  genome <- BSgenome.Hsapiens.UCSC.hg38

  # Use fixed widths spanning multiple cut categories to ensure
  # consistent type columns across samples
  fixedWidths <- rep(c(80L, 200L, 400L), length.out=600)
  makeFrags <- function(seed) {
    set.seed(seed)
    idx <- sample(length(peaks), 600, replace=TRUE)
    s <- GenomicRanges::start(peaks)[idx]
    data.table::data.table(
      seqnames = as.character(GenomicRanges::seqnames(peaks))[idx],
      start    = pmax(1L, s),
      end      = pmax(1L, s + fixedWidths - 1L)
    )
  }
  atacFrag <- list(CTRL1=makeFrags(1), CTRL2=makeFrags(2))

  se <- getWeightedCounts(
    files      = NULL,
    atacFrag   = atacFrag,
    ranges     = peaks,
    genome     = genome,
    species    = "Homo sapiens",
    fragWeight = FALSE,
    peakWeight = FALSE
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(nrow(se), length(peaks))
  expect_equal(ncol(se), length(atacFrag))
  expect_true("counts" %in% assayNames(se))
  expect_true(all(assay(se, "counts") >= 0))
})

test_that("getWeightedCounts fragWeight=TRUE returns non-negative counts", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)

  peaks <- rowRanges(peakSE)[seq_len(200)]
  genome <- BSgenome.Hsapiens.UCSC.hg38

  fixedWidths <- rep(c(80L, 200L, 400L), length.out=300)
  makeFrags <- function(seed) {
    set.seed(seed)
    idx <- sample(length(peaks), 300, replace=TRUE)
    s <- GenomicRanges::start(peaks)[idx]
    data.table::data.table(
      seqnames = as.character(GenomicRanges::seqnames(peaks))[idx],
      start    = pmax(1L, s),
      end      = pmax(1L, s + fixedWidths - 1L)
    )
  }
  atacFrag <- list(S1=makeFrags(1), S2=makeFrags(2))

  se <- getWeightedCounts(
    files      = NULL,
    atacFrag   = atacFrag,
    ranges     = peaks,
    genome     = genome,
    species    = "Homo sapiens",
    fragWeight = TRUE,
    peakWeight = FALSE
  )

  expect_true(all(assay(se, "counts") >= 0))
})
