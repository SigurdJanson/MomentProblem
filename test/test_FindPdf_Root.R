library(testthat)
source("../R/FindPdf_Root.R")

test_that("New_ByMomentPdf.default", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  # 'New_ByMomentPdf' only calls 'New_ByMomentPdf.default'
  expect_identical(New_ByMomentPdf(TargetMoments), 
                   New_ByMomentPdf.default(TargetMoments))
  
  o <- New_ByMomentPdf.default(TargetMoments)
  expect_s3_class(o, c("list", "ByMomentPdf"))
  expect_named( o, c("Function", "Moments", "TarFu", "TarMo", 
                     "DistaFu", "DistaMo",
                     "ParamSpace", "LaunchSpace"))
  expect_identical(o$TarMo, TargetMoments)
  expect_type(o$ParamSpace,  "double")
  expect_equal(dim(o$ParamSpace),  c(1, 2))
  expect_type(o$LaunchSpace, "double")
  expect_equal(dim(o$LaunchSpace),  c(0, 1))
})



test_that("GetLaunchSpace.default: Preconditions", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default(TargetMoments)
  # Wrong method
  expect_error(GetLaunchSpace( Pdf, 5L, "Nonsense"),
               "")
  expect_error(GetLaunchSpace( Pdf, 5L, "R",
               "No parameter space specified"))
  # Expansion of arguments
  Pdf$ParamSpace <- matrix(c(-5, -5, 5, 5), ncol = 2, 
                           dimnames = list(NULL, c("from", "to")))
  expect_silent(GetLaunchSpace( Pdf, 5L, "R"))
  expect_silent(GetLaunchSpace( Pdf, 5L, "E"))
  expect_error(GetLaunchSpace( Pdf, 5L, "M"),
               "not yet implemented")
  
  #method = c("Random", "Even", "Manual")
})




test_that("GetLaunchSpace.default: Result", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default(TargetMoments)
  expect_equal( Pdf$ParamSpace, matrix(data = c(-Inf, Inf), ncol = 2L, 
                                       dimnames = list(NULL, c("from", "to"))) )

  
  N <- 2L # Number of launch points in param space
  Pdf$ParamSpace <- matrix(c(-5, -5, 5, 5), ncol = N, 
                           dimnames = list(NULL, c("from", "to")))
  Pdf <- GetLaunchSpace( Pdf, N, "Random")
  expect_type(Pdf$ParamSpace, "double") #
  expect_type(Pdf$LaunchSpace, "double") #
  expect_equal(sum(Pdf$LaunchSpace[1,] >= Pdf$ParamSpace[,"from"]), N)
  expect_equal(sum(Pdf$LaunchSpace[1,] <= Pdf$ParamSpace[,"to"]), N)
  expect_equal(sum(Pdf$LaunchSpace[2,] >= Pdf$ParamSpace[,"from"]), N)
  expect_equal(sum(Pdf$LaunchSpace[2,] <= Pdf$ParamSpace[,"to"]), N)
})