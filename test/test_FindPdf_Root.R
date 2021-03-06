setwd("..")
source("./R/FindPdf_Root.R")
setwd("./test")


## CONSTRUCTOR ----
test_that("New_ByMomentPdf.default", {
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  # 'New_ByMomentPdf' only calls 'New_ByMomentPdf.default'
  expect_identical(New_ByMomentPdf(TargetMoments), 
                   New_ByMomentPdf.default(TargetMoments))
  
  o <- New_ByMomentPdf.default(TargetMoments)
  expect_s3_class(o, c("list", "ByMomentPdf"))
  expect_named( o, c("Function", "Moments", "ParamSolved", "TarFu", "TarMo",
                     "DistaFu", "DistaMo", 
                     "ParamSpace", "LaunchSpace", "Tolerance"))
  expect_identical(o$TarMo, TargetMoments)
  expect_type(o$ParamSpace,  "double")
  expect_equal(dim(o$ParamSpace),  c(2, 1))
  expect_type(o$LaunchSpace, "double")
  expect_equal(dim(o$LaunchSpace),  c(0, 1))
})





## PROBABILITY FUNCTION CALLS ----
test_that("Density, distribution, quantile & random function", {
  # Pdf <- New_ByMomentPdf(c(0, 1, 0, 3, 0, 15))
  # Pdf$Function <- "norm"
  # print(dPdf(Pdf, -5:5, ParamSet = NULL, mean = 2, sd = 0.1))
  # print(dPdf(Pdf, -5:5, ParamSet = list(mean = 2, sd = 0.1)))
  
  Pdf <- New_ByMomentPdf(c(0, 1, 0, 3, 0, 15))
  
  # First form, normal distribution
  # Not realistic but good enough for testing
  Pdf$Function <- "norm"
  expect_equal( dPdf(Pdf, -5:5, ParamSet = list(mean = 0, sd = 1)),
                dnorm(-5:5, mean = 0, sd = 1) )
  expect_equal( dPdf(Pdf, -5:5, ParamSet = list(mean = 2, sd = 0.1)),
                dnorm(-5:5, mean = 2, sd = 0.1) )
  # Second form, normal distribution
  expect_equal( dPdf(Pdf, -5:5, ParamSet = NULL, mean = 0, sd = 1),
                dnorm(-5:5, mean = 0, sd = 1) )
  expect_equal( dPdf(Pdf, -5:5, ParamSet = NULL, mean = 2, sd = 0.1),
                dnorm(-5:5, mean = 2, sd = 0.1) )
  # Third form, normal distribution
  Pdf$ParamSolved <- list(list(1, c(mean = 0, sd = 1)),
                          list(2, c(mean = 2, sd = 0.1)))
  expect_equal( dPdf(Pdf, -5:5, ParamSet = 1),
                dnorm(-5:5, mean = 0, sd = 1) )
  expect_equal( dPdf(Pdf, -5:5, ParamSet = 2),
                dnorm(-5:5, mean = 2, sd = 0.1) )
  
  # First form, Generalised Gamma distribution
  require(gld)
  Pdf$Function <- "gl"
  Lambda <- list(lambda1 = 0, lambda2 = 0.5, lambda3 = 0.5, lambda4 = 0.5)
  expect_equal( dPdf(Pdf, -5:5, ParamSet = Lambda),
                dgl(-5:5, unlist(Lambda)) )
  Lambda <- list(lambda1 = -0.1, lambda2 = 0.1, lambda3 = -0.1, lambda4 = 0.1)
  expect_equal( dPdf(Pdf, -5:5, ParamSet = Lambda),
                dgl(-5:5, unlist(Lambda)) )

  q <- seq(0, 1, by = 0.01)
  Lambda <- list(lambda1 = 0, lambda2 = 0.5, lambda3 = 0.5, lambda4 = 0.5)
  expect_equal( pPdf(Pdf, q, ParamSet = Lambda),
                pgl(q, unlist(Lambda)) )
  Lambda <- list(lambda1 = -0.1, lambda2 = 0.1, lambda3 = -0.1, lambda4 = 0.1)
  expect_equal( pPdf(Pdf, q, ParamSet = Lambda),
                pgl(q, unlist(Lambda)) )
  
  p <- seq(0, 1, by = 0.01)
  Lambda <- list(lambda1 = 0, lambda2 = 0.5, lambda3 = 0.5, lambda4 = 0.5)
  expect_equal( qPdf(Pdf, p, ParamSet = Lambda),
                qgl(p, unlist(Lambda)) )
  Lambda <- list(lambda1 = -0.1, lambda2 = 0.1, lambda3 = -0.1, lambda4 = 0.1)
  expect_equal( qPdf(Pdf, p, ParamSet = Lambda),
                qgl(p, unlist(Lambda)) )
  
  n <- 20
  Lambda <- list(lambda1 = 0, lambda2 = 0.5, lambda3 = 0.5, lambda4 = 0.5)
  set.seed(seed = 4321)
  o <- rPdf(Pdf, n, ParamSet = Lambda)
  set.seed(seed = 4321)
  e <- rgl(n, unlist(Lambda)) 
  expect_equal( o, e )
  
  Lambda <- list(lambda1 = -0.1, lambda2 = 0.1, lambda3 = -0.1, lambda4 = 0.1)
  set.seed(seed = 4321)
  o <- rPdf(Pdf, n, ParamSet = Lambda)
  set.seed(seed = 4321)
  e <- rgl(n, unlist(Lambda)) 
  expect_equal( o, e )
})



## LAUNCH SPACE ----
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




## SOLUTIONS ----
test_that("SolutionMoments", {
  succeed(message = "Tested via Gld class")
})


test_that("AddSolution", {
  # PRECONDITIONS
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default( TargetMoments )
  expect_error(AddSolution(Pdf, -1), 
               "Invalid launch point [(]must be > 0[)].") 
  expect_error(AddSolution(Pdf, 2), 
               "Index of launch point out of range.")
  
  # Add launching space to avoid out-of-range error
  Pdf$LaunchSpace <- matrix(c(0, 1, -1, 1, 1, 1, 0, 2, -1, 2, 1, 2),
                            nrow = 6, byrow = TRUE)
  Pdf$Function <- "norm"
  
  # Test
  Pdf <- AddSolution(Pdf, 1, c(0, 1))
  expect_equal(Pdf$ParamSolved, list(list(1, c(0, 1))))
  Pdf <- AddSolution(Pdf, 2, c(2, 3)) # since Append = FALSE
  expect_equal(Pdf$ParamSolved, list(list(2, c(2, 3))))
  # Append = TRUE, but launch point already exists: overwrite saved solution
  Pdf <- AddSolution(Pdf, 2, c(0, 1), Append = TRUE)
  expect_equal(Pdf$ParamSolved, list(list(2, c(0, 1))))
  Pdf <- AddSolution(Pdf, 4, c(2, 3), Append = TRUE)
  expect_equal(Pdf$ParamSolved, list(list(2, c(0, 1)),
                                     list(4, c(2, 3))))
  Pdf <- AddSolution(Pdf, 3, c(99, 1), Append = TRUE)
  expect_equal(Pdf$ParamSolved, list(list(2, c(0, 1)),
                                     list(3, c(99, 1)),
                                     list(4, c(2, 3))))
  Pdf <- AddSolution(Pdf, 1, c(7, 77), Append = TRUE)
  expect_equal(Pdf$ParamSolved, list(list(1, c(7, 77)),
                                     list(2, c(0, 1)),
                                     list(3, c(99, 1)),
                                     list(4, c(2, 3))))
})


test_that("FindPdf", {
  # PRECONDITIONS
  TargetMoments <- c(0, 1, 0, 3, 0, 15)
  Pdf <- New_ByMomentPdf.default( TargetMoments )
  expect_error(FindPdf(Pdf, -1), 
               "Invalid launch point [(]must be > 0[)].") 
  expect_error(FindPdf(Pdf, 1), 
               "Index of launch point out of range.")

  # RESULTS
  succeed(message = "Tested via Gld class")
})




test_that("EvaluatePdf / BestSolution", {
  CalcMoments <- function(a, b) c((b+a)/2, (b-a)^2/12, 0, -6/5)
  # Preconditions
  a <- 0
  b <- 1
  # Uniform distribution (https://t1p.de/q14s)
  TargetMoments <- CalcMoments(a, b)
  Pdf <- New_ByMomentPdf.default( TargetMoments )
  Pdf$Function <- "unif"
  
  # Evaluation using moments
  # Make sure that they are in the wrong order to test proper index search
  Pdf <- GetLaunchSpace( Pdf, Count = 4, Method = "Random" )
  Pdf$ParamSolved <- list(list(3, c(0.1, 1.1)))
  Pdf$ParamSolved <- c(Pdf$ParamSolved, list(list(1, c(-0.1, 1.1))))
  Pdf <- AddSolution( Pdf, 4, c(0.1, 0.9), 4, Append = TRUE)
  Pdf <- AddSolution( Pdf, 2, c(0, 1), 2, Append = TRUE)
  expect_length(Pdf$ParamSolved, 4L)# added this test because it was wrong once
  
  expect_null(Pdf$DistaMo) # test for side effect
  # Add moments because there is no built-in way to get them
  Pdf$Moments <- list(list(1, CalcMoments(-0.1, 1.1)),
                      list(2, CalcMoments(0, 1)),
                      list(3, CalcMoments(0.1, 1.1)),
                      list(4, CalcMoments(0.1, 0.9)))
  # Evaluate
  expect_silent(Pdf <- EvaluatePdf(Pdf, FALSE))
  expect_identical(nrow(Pdf$DistaMo), 4L)
  expect_identical(ncol(Pdf$DistaMo), 2L)
  expect_equal(Pdf$DistaMo[, "ID"], c(1:4))
  expect_equivalent(Pdf$DistaMo[2, "Delta"], 0)
  expect_gt(Pdf$DistaMo[1, "Delta"], Pdf$DistaMo[2, "Delta"])
  expect_gt(Pdf$DistaMo[3, "Delta"], Pdf$DistaMo[2, "Delta"])
  expect_gt(Pdf$DistaMo[4, "Delta"], Pdf$DistaMo[2, "Delta"])
  
  Best <- matrix(c(0, 1), nrow = 1, 
                 dimnames = list(NULL, paste0("Param", 1:2)))
  expect_equal(BestSolution(Pdf), Best)
  
  # Try list of moments in different order: result must be the same
  Pdf$Moments <- list(list(1, CalcMoments(-0.1, 1.1)),
                      list(3, CalcMoments(0.1, 1.1)),
                      list(2, CalcMoments(0, 1)),
                      list(4, CalcMoments(0.1, 0.9)))
  expect_equal(BestSolution(Pdf), Best)
  
  
  # Evaluation using distribution
  Pdf$Function <- "unif"
  Pdf$TarFu    <- list("dunif", c(0, 1))
  expect_silent(Pdf <- EvaluatePdf(Pdf, TRUE))

  Index <- which(Pdf$DistaFu[, "ID"] == 2)
  expect_equivalent(Pdf$DistaFu[Index, "Delta"], 0)
  expect_gt(Pdf$DistaMo[1, "Delta"], Pdf$DistaFu[Index, "Delta"])
  expect_gt(Pdf$DistaMo[3, "Delta"], Pdf$DistaFu[Index, "Delta"])
  expect_gt(Pdf$DistaMo[4, "Delta"], Pdf$DistaFu[Index, "Delta"])
  
  expect_equal(BestSolution(Pdf, TRUE), Best)
})