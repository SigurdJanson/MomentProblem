setwd("..")
source("./R/CubicScatter.R")
setwd("./test")

test_that("GetCubicPoints: Preconditions", {
  expect_error(GetCubicPoints(0L, 2), 
               "Cannot handle cubes in less than 1 dimensions")
  expect_error(GetCubicPoints(1 - .Machine$double.eps, 2), 
               "Cannot handle cubes in less than 1 dimensions")
  
  expect_error(GetCubicPoints(2L, 0), 
               "At least 'N = 1' desired point must be set")
  expect_error(GetCubicPoints(2L, 1 - .Machine$double.eps), 
               "At least 'N = 1' desired point must be set")
})



test_that("GetCubicPoints: Results", {
  # Check boundaries of data points
  for(d in 1L:10L) {
    for(n in c(1L:10L, 25L, 50L, 100L)) {
      Result <- GetCubicPoints(d, n)
      RealN <- ceiling(n^(1/d))^d
      
      # Check dimensions
      expect_equal(dim(Result), c(RealN, d), 
                   info = paste(d, n, RealN))
      
      # Check range of coordinates
      MinPos <- 1 / RealN^(1/d) / 2
      expect_equal( range(Result), c(MinPos, 1 - MinPos), 
                    info = paste(d, n, RealN) )
      
      # Check the range of distances between points
      if (RealN > 1) {
        Result <- dist(Result) # dist() expects points in rows
        MinDist <- 1 / RealN^(1/d)
        MaxDist <- sqrt(d) * (RealN^(1/d)-1) / RealN^(1/d) 
        expect_equal( range(Result), c(MinDist, MaxDist), 
                      info = paste(d, n, RealN) )
      }
    }
  }
  
  
  
})
