setwd("..")
source("./R/RandomScatter.R")
setwd("./test")

test_that("GetRandomPoints: Preconditions", {
  expect_error(GetRandomPoints(0L, 2), 
               "Cannot handle cubes in less than 1 dimensions")
  expect_error(GetRandomPoints(1 - .Machine$double.eps, 2), 
               "Cannot handle cubes in less than 1 dimensions")
  
  expect_error(GetRandomPoints(2L, 0), 
               "At least 'N = 1' desired point must be set")
  expect_error(GetRandomPoints(2L, 1 - .Machine$double.eps), 
               "At least 'N = 1' desired point must be set")
})



test_that("GetRandomPoints: Results", {
  # Check boundaries of data points
  for(d in 1L:10L) {
    for(n in c(1L:10L, 25L, 50L, 100L)) {
      Result <- GetRandomPoints(d, n)

      # Check dimensions
      expect_equal(dim(Result), c(n, d), 
                   info = paste(d, n))
      
      # Check range of coordinates
      expect_true(all(Result <= 1))
      expect_true(all(Result >= 0))
    }
  }
  
  
  
})
