setwd("..")
source("./R/FindPdf_GlD.R")
setwd("./test")

test_that("GLD: Helper Functions", {
  vk <- function(k, L3, L4) {
    j <- 0:k
    p1 <- (-1)^j / L3^(k-j) / L4^j
    p2 <- choose(k, j)
    p3 <- beta(L3*(k-j)+1, L4*j+1)
    return(sum(p1 * p2 * p3))
  }
  
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
  expect_identical(.v1f(0, 0), NaN)
  #Must be Inf for any lambda = -1; NaN for both = -1
  expect_identical(.v1f(-1, 1), -Inf)
  expect_identical(.v1f(1, -1), Inf)
  expect_identical(.v1f(-1, -1), NaN)

  # Specific values
  expect_identical(.v1f(2, 1), 1/6 - 1/2)
  expect_identical(.v1f(1, 2), 1/2 - 1/6)
  
  # Comparison to vk
  for(i in 1:40) {
    L3 <- runif(1, min = -0.25, 200)
    L4 <- runif(1, min = -0.25, 200)
    expect_equal(.v1f(L3, L4), vk(1, L3, L4))
  }
  
  # .v2f
  #Must be NaN for any lambda = 0
  expect_identical(.v2f(0, 1), NaN)
  expect_identical(.v2f(1, 0), NaN)
  expect_identical(.v2f(0, 0), NaN)

  # Must be NaN for any lambda = -1 because of beta function
  expect_identical(.v2f(1, -1), Inf)
  expect_identical(.v2f(-1, 1), Inf)
  expect_identical(.v2f(-1, -1), -Inf)
  
  # 
  expect_identical(.v2f(1, 1), 1/3 + 1/3 - 2/6)

  # Comparison to vk
  for(i in 1:40) {
    L3 <- runif(1, min = -0.25, 200)
    L4 <- runif(1, min = -0.25, 200)
    expect_equal(.v2f(L3, L4), vk(2, L3, L4))
  }
  
  # .v3f
  # Must be zero for all lambda3 = lambda4
  e <- rep(0, 100)
  v <- runif(length(e), min = 1E-9, max = 1000)
  expect_equal(.v3f(v, v), e)
  expect_identical(.v3f(1, 1), 1/4 - 1/4 - 3/12 + 3/12 )
  # Must be NaN for any lambda = 0
  expect_identical(.v3f(0, 1), NaN)
  expect_identical(.v3f(1, 0), NaN)
  expect_identical(.v3f(0, 0), NaN)
  
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
  
  # 
  #expect_identical(.v3f(2, 1), 1/56-1/4-3/4-1/40+1/20 )
  
  # Comparison to vk
  for(i in 1:40) {
    L3 <- runif(1, min = -0.25, 200)
    L4 <- runif(1, min = -0.25, 200)
    expect_equal(.v3f(L3, L4), vk(3, L3, L4))
  }
  

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
  
  # Comparison to vk
  for(i in 1:40) {
    L3 <- runif(1, min = -0.25, 200)
    L4 <- runif(1, min = -0.25, 200)
    expect_equal(.v4f(L3, L4), vk(4, L3, L4))
  }
  
})


test_that("GLD: Constructor", {
  # PRECONDITIONS
  expect_error(New_ByMomentPdf.gld(rep(1,5)),
               "Generalised lambda approximation can only handle up to 4 moments.")
  # Preconditions by parent class
  expect_error(New_ByMomentPdf.gld(rep(1,1)),
               "Mean and variance must be given at least.")
  expect_error(New_ByMomentPdf.gld(NULL),
               "At least one of 'TarMo' or 'TarFu' must be given.")
  expect_error(New_ByMomentPdf.gld(c(1,-1)),
               "Moments with even power must be >= 0.")
  
  # Created objects
  ObjectNames <- c("Function", "Moments", "ParamSolved", "TarFu", "TarMo",
                   "DistaFu", "DistaMo", "ParamSpace", "LaunchSpace", "Tolerance")
  Pdf <- New_ByMomentPdf.gld(rep(1,4))
  expect_named(Pdf, ObjectNames)
  expect_s3_class(Pdf, c("ByMomentPdf", "gld"))
  expect_identical(Pdf$Function, "gl")
})


test_that("GLD: SolutionMoments.gld", {
  # λ1 = -1.04895187e-9, λ2 = 1.4635577, λ3 = 0.13490936, λ4 = 0.13490936
  TargetMoments <- c(0, 1, 0, 3)
  Pdf <- New_ByMomentPdf.gld(TargetMoments)
  Pdf$ParamSolved <- list(1, list(c(-1.04895187e-9, 
                                    1.4635577,
                                    0.13490936,
                                    0.13490936)))
  # Check preconditions
  expect_error(SolutionMoments(Pdf, 9999), 
               "Solution 9999 not found")
  expect_error(SolutionMoments(Pdf, c(9999, 0)), 
               "Works only for a single solution")
  # Get moments
  TargetMoments <- c(0, 1, 0, 3)
  Pdf <- New_ByMomentPdf.gld(TargetMoments)
  Pdf$ParamSolved <- list(list(1, c(-1.04895187E-9, 1.4635577,
                                    0.13490936, 0.13490936)),
                          list(2, c(-1.04895187e-9, 1.4635577,
                                    0.13490936, 0.13490936)))
  o <- SolutionMoments(Pdf, 1)
  expect_length(o, 4)
  expect_equal(o, c(0, 1, 0, 3), tolerance = 0.005)
  o <- SolutionMoments(Pdf, 2)
  expect_length(o, 4)
  expect_equal(o, c(0, 1, 0, 3), tolerance = 0.005)
  # Test another index
  Pdf$ParamSolved <- list(list(99, c(-1.04895187e-9, 1.4635577,
                                     0.13490936, 0.13490936)))
  o <- SolutionMoments(Pdf, 99)
  expect_length(o, 4)
})