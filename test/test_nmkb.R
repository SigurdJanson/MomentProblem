setwd("..")
source("./R/nmkb.rcpp.R")
setwd("./test")


test_that("g", {
  # original function from dfoptim
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  
  for(x in 1:500) {
    c1 <- c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE)
    c3 <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)
    c4 <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE)
    lower <- c( 0.0,   2,   4,   8, -100, -200)
    input <- runif(length(lower), 0, 100) + lower
    upper <- runif(length(lower), 0, 100) + input
    lower[c4] <- rep(+Inf, length(lower[c4]))
    upper[c3] <- rep(-Inf, length(upper[c3]))
    
    e <- suppressWarnings( g(input) )
    o <- gcpp(input, c1, c3, c4, lower, upper)
    expect_equal(o, e)
  }
})


test_that("ginv", {
  # original function from dfoptim
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  
  for(x in 1:500) {
    c1 <- c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE)
    c3 <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)
    c4 <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE)
    lower <- c( 0.0,   2,   4,   8, -100, -200)
    input <- runif(length(lower), 0, 100) + lower
    upper <- runif(length(lower), 0, 100) + input
    lower[c4] <- rep(+Inf, length(lower[c4]))
    upper[c3] <- rep(-Inf, length(upper[c3]))
    
    e <- suppressWarnings( ginv(input) )
    o <- ginvcpp(input, c1, c3, c4, lower, upper)
    expect_equal(o, e)
  }
})


test_that("Preconditions of nmkbcpp", {
  rosbkext <- function(x) {
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
  
  np <- 10
  set.seed(123)
  p0 <- rnorm(np)
  p0[p0 >= +2] <- +2 - 1E-8
  p0[p0 <= -2] <- -2 + 1E-8
  
  badindex <- sample.int(np, 1)
  p0[badindex] <- -2
  expect_error(nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2), 
               "Starting vector 'par' must lie [*]strictly[*] between lower and upper bounds")
  p0[badindex] <- -2 + .Machine$double.eps
  expect_silent(nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2))
})



test_that("Rosenbrock: nmkb vs Rcpp function 'nmkb'", {
  rosbkext <- function(x) {
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
  
  # get original nmkb
  setwd("..")
  source("./R/nmkb.rcpp.R")
  nmkbcpp <- nmkb
  source("./R/nmkb.R")
  setwd("./test")
  
  np <- 10
  set.seed(123)
  p0 <- rnorm(np)
  p0[sample.int(np, 1)]
  p0[p0 > +2] <- +2 - 1E-8
  p0[p0 < -2] <- -2 + 1E-8
  e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2)
  o <- nmkbcpp(fn = rosbkext, par = p0, lower = -2, upper = 2)
  expect_identical(o, e)
  
  ctrl <- list(maxfeval = 50000)
  p0 <- rnorm(np)
  p0[p0 > +2] <- +2 - 1E-8
  p0[p0 < -2] <- -2 + 1E-8
  set.seed(123)
  e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
  set.seed(123)
  o <- nmkbcpp(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
  expect_identical(o, e)
})




test_that("non-smooth problem: nmkb vs nmkbcpp", {
  hald <- function(x) {
    #Hald J & Madsen K (1981), Combined LP and quasi-Newton methods 
    #for minimax optimization, Mathematical Programming, 20, p.42-62.
    i <- 1:21
    t <- -1 + (i - 1)/10
    f <- (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
    max(abs(f))
  }
  
  # get original nmkb
  setwd("..")
  source("./R/nmkb.rcpp.R")
  nmkbcpp <- nmkb
  source("./R/nmkb.R")
  setwd("./test")
  
  p0 <- runif(5, c(0,0,0,0,-2) - 1E-7, 4 - 1E-7)
  o <- nmkbcpp(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
  e <- nmkb(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
  expect_identical(o, e)
})

