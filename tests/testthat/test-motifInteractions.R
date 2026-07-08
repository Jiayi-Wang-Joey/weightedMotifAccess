set.seed(123)
data <- betterChromVAR::getDummyData()

dev <- getInteractionsDeviations(data$counts, data$motifMatches, bait="motif1",
                                 minCount=5)
dev$group <- rep(LETTERS[1:2], each=5)

test_that("motif interaction functions works", {
  expect_true(all(dim(dev)==c(ncol(motifMatches), ncol(counts))))
  res <- discoverMotifInteractions(dev, "group")
  
})



test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
