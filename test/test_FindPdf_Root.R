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
  expect_equal(dim(o$ParamSpace),  c(2, 1))
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
  Pdf$ParamSpace <- matrix(c(-5, -5, 5, 5), ncol = 2, byrow = TRUE,
                           dimnames = list(c("from", "to"), NULL))
  expect_silent(GetLaunchSpace( Pdf, 5L, "R"))
  expect_silent(GetLaunchSpace( Pdf, 5L, "C"))
  expect_silent(GetLaunchSpace( Pdf, 5L, "H"))
  expect_error(GetLaunchSpace( Pdf, 5L, "M"),
               "not yet implemented")
  
  #method = c("Random", "Even", "Manual")
})




test_that("GetLaunchSpace.default: Result", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default(TargetMoments)
  expect_equal( Pdf$ParamSpace, matrix(data = c(-Inf, Inf), 
                                       nrow = 2L, byrow = TRUE,
                                       dimnames = list(c("from", "to"), NULL) ))

  for(m in c("Random", "Cubic", "Harmonic")) { # "Manual" not implemented
    N <- 2L   # Number of launch points in param space
    Dim <- 3L # Dimensions of parameter space = number of parameters
    Pdf$ParamSpace <- matrix(c(rep(-5, Dim), rep(5, Dim)), 
                             nrow = 2, byrow = TRUE, 
                             dimnames = list(c("from", "to"), NULL))
    Pdf <- GetLaunchSpace( Pdf, N, m )
    
    expect_type(Pdf$ParamSpace, "double") #
    expect_type(Pdf$LaunchSpace, "double") #
    
    # Check dimensions of result
    expect_equal(dim(Pdf$ParamSpace), c(2, Dim),
                 info = paste(m))

    if(m == "Cubic") 
      expect_equal(dim(Pdf$LaunchSpace), c(8, Dim),
                   info = paste(m))
    else
      expect_equal(dim(Pdf$LaunchSpace), c(2, Dim),
                   info = paste(m))
    
    # First launch point
    expect_equal(sum(Pdf$LaunchSpace[1,] >= Pdf$ParamSpace["from", ]), Dim,
                 info = paste(m))
    if(m == "Harmonic")
      expect_silent(sum(Pdf$LaunchSpace[1,] >= Pdf$ParamSpace["from", ]))
    expect_equal(sum(Pdf$LaunchSpace[1,] <= Pdf$ParamSpace["to", ]), Dim,
                 info = paste(m))
    # Second launch point
    expect_equal(sum(Pdf$LaunchSpace[2,] >= Pdf$ParamSpace["from", ]), Dim,
                 info = paste(m))
    expect_equal(sum(Pdf$LaunchSpace[2,] <= Pdf$ParamSpace["to", ]), Dim,
                 info = paste(m))
    
    if(m == "Cubic") {
      for(lp in 3:8) {
        expect_equal(sum(Pdf$LaunchSpace[lp,] >= Pdf$ParamSpace["from", ]), Dim,
                     info = paste(m))
        expect_equal(sum(Pdf$LaunchSpace[lp,] <= Pdf$ParamSpace["to", ]), Dim,
                     info = paste(m))
      }
    }
  }
})

#TODO: devise test with non-overlapping dimensions in param space
