local({
  dd <- betterChromVAR::getDummyData()
  dummy_counts <- dd$counts
  dummy_motifs <- dd$motifMatches
  dummy_counts$group <- rep(LETTERS[1:2], each=5)

  test_that("getInteractionsDeviations returns a SummarizedExperiment", {
    dev <- getInteractionsDeviations(dummy_counts, dummy_motifs, bait="motif1",
                                     minCount=5)
    expect_s4_class(dev, "SummarizedExperiment")
    expect_true("z" %in% assayNames(dev))
    expect_true("deviations" %in% assayNames(dev))
    expect_true("isBait" %in% colnames(rowData(dev)))
    expect_equal(sum(rowData(dev)$isBait), 1L)
  })

  test_that("getInteractionsDeviations adds overlap metadata columns", {
    dev <- getInteractionsDeviations(dummy_counts, dummy_motifs, bait="motif1",
                                     minCount=5)
    rd <- rowData(dev)
    expect_true("overlap.percentOfBait" %in% colnames(rd))
    expect_true("overlap.jaccard" %in% colnames(rd))
    expect_true("overlap.log2Enr" %in% colnames(rd))
  })

  test_that("getInteractionsDeviations errors on unknown bait", {
    expect_error(
      getInteractionsDeviations(dummy_counts, dummy_motifs, bait="notAMotif")
    )
  })

  test_that("discoverMotifInteractions returns a data.frame with expected columns", {
    dev <- getInteractionsDeviations(dummy_counts, dummy_motifs, bait="motif1",
                                     minCount=5)
    dev$group <- rep(LETTERS[1:2], each=5)
    res <- discoverMotifInteractions(dev, "group")

    expect_s3_class(res, "data.frame")
    expect_true("motif" %in% colnames(res))
    expect_true("adj.P.Val" %in% colnames(res))
    expect_true("difference" %in% colnames(res))
    expect_true("overlap.N" %in% colnames(res))
    expect_true(all(res$adj.P.Val >= 0 & res$adj.P.Val <= 1))
  })

  test_that("discoverMotifInteractions includes contrast column", {
    dev <- getInteractionsDeviations(dummy_counts, dummy_motifs, bait="motif1",
                                     minCount=5)
    dev$group <- rep(LETTERS[1:2], each=5)
    res <- discoverMotifInteractions(dev, "group")
    expect_true("contrast" %in% colnames(res))
  })
})
