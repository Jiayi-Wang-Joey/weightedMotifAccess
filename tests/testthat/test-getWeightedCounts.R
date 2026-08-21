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
    atac       = atacFrag,
    ranges     = peaks,
    genome     = genome,
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
    atac       = atacFrag,
    ranges     = peaks,
    genome     = genome,
    fragWeight = TRUE,
    peakWeight = FALSE
  )

  expect_true(all(assay(se, "counts") >= 0))
})

test_that("getWeightedCounts accepts file paths and GRanges for 'atac'", {
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
  ref <- assay(getWeightedCounts(atac=atacFrag, ranges=peaks, genome=genome))

  # the same fragments written out as plain .bed and as gzipped .tsv
  beds <- vapply(names(atacFrag), function(n) tempfile(fileext=".bed"), "")
  tsvs <- vapply(names(atacFrag), function(n) tempfile(fileext=".tsv.gz"), "")
  for (n in names(atacFrag)) {
    data.table::fwrite(atacFrag[[n]], beds[[n]], sep="\t", col.names=FALSE)
    data.table::fwrite(atacFrag[[n]], tsvs[[n]], sep="\t", col.names=FALSE)
  }

  expect_equal(assay(getWeightedCounts(atac=beds, ranges=peaks,
                                       genome=genome)), ref)
  expect_equal(assay(getWeightedCounts(atac=tsvs, ranges=peaks,
                                       genome=genome)), ref)

  # list of GRanges, and GRangesList
  grs <- lapply(atacFrag, GenomicRanges::makeGRangesFromDataFrame)
  expect_equal(assay(getWeightedCounts(atac=grs, ranges=peaks,
                                       genome=genome)), ref)
  expect_equal(assay(getWeightedCounts(atac=as(grs, "GRangesList"),
                                       ranges=peaks, genome=genome)), ref)

  # unnamed paths: sample names are derived from the file names
  se <- getWeightedCounts(atac=unname(tsvs), ranges=peaks, genome=genome)
  expect_equal(colnames(se), unname(sub("\\.tsv\\.gz$", "", basename(tsvs))))

  # unsupported extension
  expect_error(getWeightedCounts(atac=c(S1=tempfile(fileext=".foo")),
                                 ranges=peaks, genome=genome),
               "Cannot determine the format")
  # missing coordinate columns
  expect_error(getWeightedCounts(atac=list(S1=data.table::data.table(a=1)),
                                 ranges=peaks, genome=genome),
               "missing the column")
})

test_that("commented headers and extra columns are read without warnings", {
  # leading '#' comment lines plus columns beyond the coordinates: only the
  # first three are used, and no warning is emitted.
  lines <- c("# id=mysample", "# description=whatever",
             "chr1\t100\t200\textra\t2", "chr1\t300\t400\textra\t1")
  expected <- data.table::data.table(seqnames="chr1", start=c(100L, 300L),
                                     end=c(200L, 400L))

  plain <- tempfile(fileext=".tsv")
  writeLines(lines, plain)
  expect_silent(got <- weightedMotifAccess:::.importFragments(c(S1=plain)))
  expect_equal(got$S1, expected)

  gz <- paste0(plain, ".gz")
  con <- gzfile(gz, "wt"); writeLines(lines, con); close(con)
  expect_silent(got <- weightedMotifAccess:::.importFragments(c(S1=gz)))
  expect_equal(got$S1, expected)

  # a file with no comment header is still read correctly
  bare <- tempfile(fileext=".bed")
  writeLines(lines[3:4], bare)
  expect_equal(weightedMotifAccess:::.importFragments(c(S1=bare))$S1, expected)
})

test_that("a single sample can be passed to 'atac' without a list", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)

  peaks <- rowRanges(peakSE)[seq_len(200)]
  genome <- BSgenome.Hsapiens.UCSC.hg38

  set.seed(1)
  idx <- sample(length(peaks), 300, replace=TRUE)
  s <- GenomicRanges::start(peaks)[idx]
  dt <- data.table::data.table(
    seqnames = as.character(GenomicRanges::seqnames(peaks))[idx],
    start    = pmax(1L, s),
    end      = pmax(1L, s + rep(c(80L, 200L, 400L), length.out=300) - 1L))

  ref <- assay(getWeightedCounts(atac=list(dt), ranges=peaks, genome=genome))

  # a bare data.table is one sample, not one sample per column: a data.table
  # is also a list, so this guards the branch order in .importFragments().
  bare <- getWeightedCounts(atac=dt, ranges=peaks, genome=genome)
  expect_equal(ncol(bare), 1L)
  expect_equal(assay(bare), ref)

  # same for a bare data.frame, GRanges and file path
  expect_equal(assay(getWeightedCounts(atac=as.data.frame(dt), ranges=peaks,
                                       genome=genome)), ref)
  gr <- GenomicRanges::makeGRangesFromDataFrame(dt)
  expect_equal(assay(getWeightedCounts(atac=gr, ranges=peaks,
                                       genome=genome)), ref)
  # a bare path gives the same counts, but is named after the file rather
  # than falling back to "sample1"
  f <- tempfile(fileext=".bed")
  data.table::fwrite(dt, f, sep="\t", col.names=FALSE)
  fromFile <- getWeightedCounts(atac=f, ranges=peaks, genome=genome)
  expect_equal(colnames(fromFile), sub("\\.bed$", "", basename(f)))
  expect_equal(unname(assay(fromFile)), unname(ref))

  # 'chr' is accepted as an alias for 'seqnames'
  dt2 <- data.table::copy(dt)
  data.table::setnames(dt2, "seqnames", "chr")
  expect_equal(assay(getWeightedCounts(atac=dt2, ranges=peaks,
                                       genome=genome)), ref)
})
