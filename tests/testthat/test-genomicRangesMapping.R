test_that("genomicRangesMapping with single byCols returns a sparse matrix", {
  peaks <- rowRanges(peakSE)[seq_len(100)]
  mm <- assay(motifMatches, "motifMatches")[seq_len(100), seq_len(3)]
  hits <- which(mm, arr.ind=TRUE)
  at <- data.table::data.table(
    seqnames = as.character(GenomicRanges::seqnames(peaks))[hits[,"row"]],
    start    = GenomicRanges::start(peaks)[hits[,"row"]],
    end      = GenomicRanges::end(peaks)[hits[,"row"]],
    sample   = colnames(mm)[hits[,"col"]],
    score    = 1L
  )
  res <- genomicRangesMapping(peaks, assayTable=at, byCols="sample",
                              scoreCol="score", type="equal")

  expect_true(inherits(res, "Matrix") || is.matrix(res))
  expect_equal(nrow(res), length(peaks))
})

test_that("genomicRangesMapping with two byCols returns a list of matrices", {
  peaks <- rowRanges(peakSE)[seq_len(100)]
  mm <- assay(motifMatches, "motifMatches")[seq_len(100), seq_len(3)]
  hits <- which(mm, arr.ind=TRUE)
  n <- nrow(hits)
  at <- data.table::data.table(
    seqnames = as.character(GenomicRanges::seqnames(peaks))[hits[,"row"]],
    start    = GenomicRanges::start(peaks)[hits[,"row"]],
    end      = GenomicRanges::end(peaks)[hits[,"row"]],
    motif    = colnames(mm)[hits[,"col"]],
    sample   = rep(c("A","B"), length.out=n),
    score    = 1L
  )
  res <- genomicRangesMapping(peaks, assayTable=at, byCols=c("motif","sample"),
                              scoreCol="score", type="equal")

  expect_type(res, "list")
  expect_true(all(sapply(res, function(m) inherits(m, "Matrix") || is.matrix(m))))
  expect_equal(nrow(res[[1]]), length(peaks))
  expect_equal(ncol(res[[1]]), 2L)
})

test_that("genomicRangesMapping values are non-negative", {
  peaks <- rowRanges(peakSE)[seq_len(50)]
  mm <- assay(motifMatches, "motifMatches")[seq_len(50), seq_len(2)]
  hits <- which(mm, arr.ind=TRUE)
  at <- data.table::data.table(
    seqnames = as.character(GenomicRanges::seqnames(peaks))[hits[,"row"]],
    start    = GenomicRanges::start(peaks)[hits[,"row"]],
    end      = GenomicRanges::end(peaks)[hits[,"row"]],
    motif    = colnames(mm)[hits[,"col"]],
    sample   = rep(c("S1","S2"), length.out=nrow(hits)),
    score    = runif(nrow(hits))
  )
  res <- genomicRangesMapping(peaks, assayTable=at, byCols=c("motif","sample"),
                              scoreCol="score", aggregationFun=sum, type="equal")
  expect_true(all(sapply(res, function(m) all(m >= 0))))
})
