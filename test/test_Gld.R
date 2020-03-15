library(testthat)
source("../MomentProblem.R")

test_that("GLD Helper Functions", {
  # .v1f
  # Must be zero for all lambda3 = lambda4
  e <- rep(0, 100)
  v <- runif(length(e), min = -0.25, max = 1000)
  expect_identical(.v1f(v, v), e)
  # Test the lower limit, lambda3 and 4 shall never be < -0.25
  expect_identical(.v1f(-0.25, -0.25), 0)
  
  # Must be positive for all lambda4 > lambda3
  e <- rep(TRUE, 100)
  v <- runif(length(e), min = -0.25, max = 1000)
  expect_identical(.v1f(0, 0+v) > 0, e)
  # ... and  negative for all lambda3 > lambda4
  e <- rep(FALSE, 100)
  expect_identical(.v1f(0+v, 0) > 0, e)
  
  #Must be Inf for any lambda = 0; NaN for both = 0
  expect_identical(.v1f(0, 1), Inf)
  expect_identical(.v1f(1, 0), -Inf)
  expect_identical(.v1f(0, 0), 0)
  #Must be Inf for any lambda = -1; NaN for both = -1
  expect_identical(.v1f(-1, 1), -Inf)
  expect_identical(.v1f(1, -1), Inf)
  expect_identical(.v1f(-1, -1), NaN)
  
  
  # .v2f
  #Must be NaN for any lambda = 0
  expect_identical(.v2f(0, 1), NaN)
  expect_identical(.v2f(1, 0), NaN)
  expect_identical(.v2f(0, 0), NaN)

  # Must be NaN for any lambda = -1 because of beta function
  expect_identical(.v2f(1, -1), Inf)
  expect_identical(.v2f(-1, 1), Inf)
  expect_identical(.v2f(-1, -1), -Inf)
  
  # .v3f
  # Must be zero for all lambda3 = lambda4
  e <- rep(0, 100)
  v <- runif(length(e), min = 1E-9, max = 1000)
  expect_equal(.v3f(v, v), e)
  # Must be NaN for any lambda = 0
  expect_identical(.v3f(0, 1), NaN)
  expect_identical(.v3f(1, 0), NaN)
  expect_identical(.v3f(0, 0), 0)
  
  #Must be ... for lambda3 = -0.5 or lambda4 = -1 because of beta function
  suppressWarnings({
    expect_identical(.v3f(1, -1), NaN)
  })
  suppressWarnings({
    expect_identical(.v3f(-1, 1), NaN)
  })
  expect_identical(.v3f(1, -0.5), Inf)
  expect_identical(.v3f(-0.5, 1), -Inf)
  suppressWarnings({
    expect_identical(.v3f(-0.5, -1), NaN)
  })

  # .v4f
  # Must be NaN for any lambda = 0
  expect_identical(.v4f(0, 1), NaN)
  expect_identical(.v4f(1, 0), NaN)
  expect_identical(.v4f(0, 0), NaN)
  # Must be NaN for any lambda = -0.5 (beta(3*-0.5+1, x) == NaN)
  suppressWarnings({
    expect_identical(.v4f(-0.5, 1), NaN)
  })
  suppressWarnings({
    expect_identical(.v4f(1, -0.5), NaN)
  })
  suppressWarnings({
    expect_identical(.v4f(-0.5, -0.5), NaN)
  })
  
})