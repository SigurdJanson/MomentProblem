setwd("..")
source("./R/FindPdf_Root.R")
setwd("./test")


test_that("Helper function: logseq", {
  # exp(log(10)*seq(log10(from), log10(to), by = by))
  expect_equal(logseq(1, 1000, 1), c(1, 10, 100, 1000))
  
  from <- 12
  to   <- 8634
  e <- 10^seq(log10(from), log10(to), 0.36)
  o <- logseq(from, to, 0.36)
  expect_equal(o, e)
})
