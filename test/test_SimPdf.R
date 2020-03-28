library(testthat)
source("../R/ComparePdf.R")
source("./SimPdf.ori.R")

test_that("SimPdf - Compare with original version", {
  # Make sure that refactoring didn't have any side effects

  # Normal vs. Normal - Identical results
  for(R in c(2, 5, 10, 20)) {
    Range <- c(-R, R)
    o <- SimPdf(dnorm, dnorm, Range)
    e <- SimPdfOri(dnorm, dnorm, -R, R) / 100
    names(e) <- NULL
    expect_equal(o, e, info = paste("R =", R))
  }
  
  # Normal vs. Uniform
  for(R in c(1, 2, 5, 10, 20)) {
    Range <- c(-R, R)
    o <- SimPdf(dnorm, dunif, Range)
    e <- SimPdfOri(dnorm, dunif, -R, R) / 100
    names(e) <- NULL
    expect_equal(o, e, info = paste("R =", R))
  }
  
  # Normal vs. Log-Normal
  o <- SimPdf(dnorm, dlnorm)
  e <- SimPdfOri(dnorm, dlnorm) / 100
  names(e) <- NULL
  expect_equal(o, e)
})
  


test_that("SimPdf: ", {
  o <- SimPdf(dnorm, dnorm)
  e <- 1
  expect_identical(o, e)
  
  o <- SimPdf( dnorm, dt, Args2 = list(df = Inf))
  e <- 1
  expect_identical(o, e)
})




test_that("Moments", {

  r <- rnorm(1000)
  m <- MomentSet(r, nmax=8)$moments
  e <- momentsetOri(r, nmax=8)$moments
  expect_equal(m, e)

})
