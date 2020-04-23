setwd("..")
source("./R/FindPdf_Root.R")
setwd("./test")

test_that("Summary", {
  # Same data as "AddSolution" test
  TargetMoments <- c(0, 1, 0, 3)
  Pdf <- New_ByMomentPdf.default( TargetMoments )
  Pdf$LaunchSpace <- matrix(c(0, 1, -1, 1, 0, 2, -1, 2),
                            nrow = 4, byrow = TRUE)
  Pdf$Function <- "norm"
  
  # Test without solutions
  expect_output(summary(Pdf), 
                "Solutions have not been evaluated")
  expect_output(summary(Pdf), 
                "Solutions[ ]+0")
  
  Pdf$DistaMo <- matrix(c(1:4, 0, sqrt(2^2+3^2), sqrt(2^2+3^2), 99), ncol = 2, 
                        dimnames = list(NULL, c("ID", "Delta")) )
  expect_output(summary(Pdf), 
                "Solutions[ ]+0")
  
  # New case: fill with solutions
  Pdf <- AddSolution(Pdf, 1, c(0, 1))
  Pdf <- AddSolution(Pdf, 2, c(2, 3), Append = TRUE) 
  Pdf <- AddSolution(Pdf, 4, c(2, 3), Append = TRUE)
  Pdf <- AddSolution(Pdf, 3, c(99, 1), Append = TRUE)
  
  expect_length(Pdf$DistaMo, 0) # side effect of adding solutions and solution count does not fit distance count
  Pdf$DistaMo <- matrix(c(1:4, 0, sqrt(2^2+3^2), sqrt(2^2+3^2), 99), ncol = 2, 
                        dimnames = list(NULL, c("ID", "Delta")) )
  
  # Tests
  expect_output(summary(Pdf), 
                "Solutions[ ]+4")
  expect_output(summary(Pdf), 
                "1[ ]+solutions with best characteristics")
  expect_output(summary(Pdf), 
                "Euclidean Distance[ ]+0")
  expect_output(summary(Pdf), 
                "Failed cycles[ ]+0")
  expect_output(summary(Pdf), 
                "Unique solutions[ ]+3")
  
  # New case
  Pdf <- AddSolution(Pdf, 2, c(0, 1), Append = TRUE)
  Pdf <- AddSolution(Pdf, 1, c(7, 77), Append = TRUE)
  Pdf$LaunchSpace <- rbind(Pdf$LaunchSpace, c(0, 0))
  Pdf <- AddSolution(Pdf, 5, c(7, 77), Append = TRUE)
  Pdf$LaunchSpace <- rbind(Pdf$LaunchSpace, c(0, 0))
  Pdf <- AddSolution(Pdf, 6, c(7, 77), Append = TRUE)
  Pdf$LaunchSpace <- rbind(Pdf$LaunchSpace, c(0, 0))
  Pdf <- AddSolution(Pdf, 7, c(NA, 77), Append = TRUE)
  expect_length(Pdf$ParamSolved, 7) # just to be on the safe side
  Pdf$DistaMo <- matrix(c(1:7, 0, 0, sqrt(4+9), sqrt(4+9), 76.32169, 76.32169, NA), 
                        ncol = 2, dimnames = list(NULL, c("ID", "Delta")) )

  # Tests
  expect_output(summary(Pdf), 
                "Solutions[ ]+7")
  expect_output(summary(Pdf), 
                "2[ ]+solutions with best characteristics")
  expect_output(summary(Pdf), 
                "Euclidean Distance[ ]+0")
  expect_output(summary(Pdf), 
                "Failed cycles[ ]+1")
  expect_output(summary(Pdf), 
                "Unique solutions[ ]+5")
  
})