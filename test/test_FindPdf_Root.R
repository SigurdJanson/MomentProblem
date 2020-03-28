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
                     "ParamSpace", "SpaceLaunch"))
  expect_identical(o$TarMo, TargetMoments)
})


test_that("GetLaunchSpace.default: Preconditions", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default(TargetMoments)
  # Wrong method
  expect_error(GetLaunchSpace.default( Pdf, 5L, "Nonsense"),
               "")
  expect_silent(GetLaunchSpace.default( Pdf, 5L, "R"))
  expect_silent(GetLaunchSpace.default( Pdf, 5L, "E"))
  expect_error(GetLaunchSpace.default( Pdf, 5L, "M"),
               "Not yet implemented")
  
  #method = c("Random", "Even", "Manual") )
})
