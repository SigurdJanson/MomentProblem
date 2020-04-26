setwd("..")
source("./R/FindPdf_Root.R")
setwd("./test")


test_that("Helper function: logseq", {
  # PRECONDITIONS
  expect_error(logseq(0, 1000, 1), "Values <= 0 not allowed in logarithmic scale")
  expect_error(logseq(1000, -1, by = -1), "Values <= 0 not allowed in logarithmic scale")
  
  # USING KNOWN RESULTS
  # Forwards using 'by'
  expect_equal(logseq(1, 1000, by = 1), c(1, 10, 100, 1000))
  # Backwards using 'by'
  expect_equal(logseq(10000, 1, by = -1), c(10000, 1000, 100, 10, 1))
  
  # USING EQUIVALENT EXPRESSION
  # Forwards using 'by'
  from <- 12
  to   <- 8634
  e <- 10^seq(log10(from), log10(to), 0.36)
  o <- logseq(from, to, 0.36)
  expect_equal(o, e)
  
  # Backwards using 'length.out'
  from <- 99
  to   <- 1
  e <- 10^seq(log10(from), log10(to), length.out = 27)
  o <- logseq(from, to, length.out = 27)
  expect_equal(o, e)
  
  # Backwards using 'along.with'
  from <- 12
  to   <- 5
  e <- 10^seq(log10(from), log10(to), along.with = 1:12)
  o <- logseq(from, to, along.with = 1:12)
  expect_equal(o, e)
})
