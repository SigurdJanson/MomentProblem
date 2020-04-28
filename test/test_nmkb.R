setwd("..")
source("./R/nmkp.R")
source("./test/nmk.dfoptim.R")
source("./test/nmkb.dfoptim.R")
setwd("./test")

# Helper Functions ----
test_that("Helper .ReplaceVertex", {
  #.ReplaceVertex(V, Vnew, At)
  
  # Default 4 x 5 matrix
  V <- matrix(rep(1:4, 5) + rep(seq(0.1, 0.5, 0.1), each = 4), nrow = 4)
  
  # First column
  Vnew <- rep(1.0, nrow(V))
  e <- cbind( Vnew, V[, 1:(ncol(V)-1)], deparse.level = 0 )
  expect_equal(.ReplaceVertex(V, Vnew, 0), e) # At = 0
  expect_equal(.ReplaceVertex(V, Vnew, 1), e)
  
  # Last column
  Vnew <- rep(0.0, nrow(V))
  e <- cbind( V[, 1:(ncol(V)-1)], Vnew, deparse.level = 0 )
  expect_equal(.ReplaceVertex(V, Vnew, ncol(V)), e)
  expect_equal(.ReplaceVertex(V, Vnew, ncol(V)+1), e) # At > ncol

  # Middle columns
  Vnew <- rep(0.0, nrow(V))
  e <- cbind( V[, 1:2], Vnew, V[, 3:(ncol(V)-1)], deparse.level = 0 )
  expect_equal(.ReplaceVertex(V, Vnew, 3), e)
  e <- cbind( V[, 1:3], Vnew, V[, 4:(ncol(V)-1)], deparse.level = 0 )
  expect_equal(.ReplaceVertex(V, Vnew, 4), e)
  
  # Special case: min size of matrix has two columns
  V <- matrix(rep(1:4, 2), nrow = 4)
  for(c in 1:ncol(V)) {
    Vnew <- rep(c, nrow(V))
    e <- V
    e[, c] <- Vnew
    expect_equal(.ReplaceVertex(V, Vnew, c), e, 
                 info = paste0(c))
  }
  
})



# Testing nmkb/p ----
test_that("Maximise/Minimise using Rosenbrock: dfoptim::nmkb", {
  rosbkext <- function(x) {
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
  
  # Test maximisation by calling a negative 'rosbkext' and 
  # compare it with minimisation: results shall be equal
  np <- 10
  set.seed(123)
  p0 <- rnorm(np)
  p0[p0 > +2] <- +2 - 1E-8
  p0[p0 < -2] <- -2 + 1E-8
  
  # Start with dfoptim::nmkb
  ctrl <- list(maximize = FALSE)
  e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
  e$value <- e$value * -1
  ctrl <- list(maximize = TRUE)
  o <- nmkb(fn = function(x) -rosbkext(x), par = p0, lower = -2, upper = 2, control = ctrl)
  expect_equal(o, e) 

  # Verify with nmkp
  ctrl <- list(maximize = FALSE)
  e <- nmkp(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
  e$value <- e$value * -1
  ctrl <- list(maximize = TRUE)
  o <- nmkp(fn = function(x) -rosbkext(x), par = p0, lower = -2, upper = 2, control = ctrl)
  expect_equal(o, e) 
})


test_that("Minimise Common Test Functions: nmkp", {
  rosbkext <- function(x) {
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
  sphere <- function(x) sum(x^2)
  

  # Test maximisation by calling a negative 'rosbkext' and 
  # compare it with minimisation: results shall be equal
  tol <- 5e-5 # accepted tolerance
  for (np in c(2:6, 12, 24, 32)) {
    set.seed(123)
    p0 <- rnorm(np)
    p0[p0 > +2] <- +2 - 1E-8
    p0[p0 < -2] <- -2 + 1E-8
    
    ctrl <- list(maxfeval = 1E5, tol = 1E-8)
    o <- nmkp(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
    expect_equal(o$convergence, 0L, info = paste0(np))
    expect_equal(o$value, 0.0, tolerance = tol*np, info = paste0(np)) 
    expect_equal(o$par, rep(1.0, np), tolerance = tol*np, info = paste0(np)) 
  }
  
  # Sphere function
  tol <- 5e-5 # accepted tolerance
  box <- 5.2
  for (np in c(2:6, 12, 24, 32)) {
    set.seed(123)
    p0 <- rnorm(np)
    p0[p0 > +box] <- +box - 1E-8
    p0[p0 < -box] <- -box + 1E-8
    
    ctrl <- list(maxfeval = 1E5, tol = 1E-8)
    o <- nmkp(fn = sphere, par = p0, lower = -box, upper = box, control = ctrl)
    expect_equal(o$convergence, 0L, info = paste0(np)) 
    expect_equal(o$value, 0.0, tolerance = tol*np, info = paste0(np))
    expect_equal(o$par, rep(0.0, np), tolerance = tol*np, info = paste0(np))
  }
  
})



# Comparison with Original ----
test_that("Comparison 'dfoptim::nmkb' & 'nmk' vs 'nmkp'", {
  # Functions
  rosbk <- function(x) {
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
  sphere   <- function(x) sum(x^2)
  mccormick <- function(x) sin(sum(x)) + (x[1]-x[2])^2 - sum(c(1.5, 2.5)*x) +1
  rastrigin <- function(x) 10*length(x) + sum(x^2 - 10*cos(1*pi*x))
  schwefel <- function(x) 418.9829*length(x) - sum(x*sin(sqrt(abs(x))))
  sudipo <- function(x) sum((abs(x))^(2:(length(x)+1)))
  
  # Direct comparison only possible for two dimensions because of
  # improvement of transformation parameters rho, chi, ...
  np <- 2
  
  
  Functions <- c(rosbk, sphere, mccormick, rastrigin, schwefel, sudipo)
  loop <- 0
  for(f in Functions) {
    loop <- loop +1
    set.seed(123)
    p0 <- rnorm(np)
    # Unconstrained
    if (loop != 3) {
      e <- nmk(fn = f, par = p0)
      o <- nmkp(fn = f, par = p0)
      expect_equal(o, e, info = paste0(f, "(", p0, ")")) 
    }
    else { # for some reason it cannot handle the mccormick function
      # Error: "System ist für den Rechner singulär: reziproke Konditionszahl"
      expect_error(nmk(fn = f, par = p0))
      expect_error(nmkp(fn = f, par = p0))
    }
    # Box constrained
    p0[p0 > +2] <- +2 - 1E-8
    p0[p0 < -2] <- -2 + 1E-8
    e <- nmkb(fn = f, par = p0, lower = -2, upper = 2)
    o <- nmkp(fn = f, par = p0, lower = -2, upper = 2)
    expect_equal(o, e, info = paste0(f, "(", p0, ")")) 
  }

  
  # Rosenbrock terminates too soon with default 'control$maxfeval'.
  # Hence run it again. Just for diligence
  ctrl <- list(maxfeval = 50000)
  set.seed(123)
  p0 <- rnorm(np)
  p0[p0 > +2] <- +2 - 1E-8
  p0[p0 < -2] <- -2 + 1E-8
  e <- nmkb(fn = rosbk, par = p0, lower = -2, upper = 2, control = ctrl)
  o <- nmkp(fn = rosbk, par = p0, lower = -2, upper = 2, control = ctrl)
  expect_equal(o, e)
  
})




# setwd("..")
# source("./R/nmkb.rcpp.R")
# setwd("./test")
# 
# test_that("g", {
#   # original function from dfoptim
#   g <- function(x) {
#     gx <- x
#     gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)
#     gx[c3] <- log(x[c3] - lower[c3])
#     gx[c4] <- log(upper[c4] - x[c4])
#     gx
#   }
#   
#   # These are the cases that should happen
#   for(x in 1:1000) {
#     c1 <- c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE)
#     c3 <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)
#     c4 <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE)
#     lower <- c(0.0,   2,   4,   8, -100, -200)
#     input <- runif(length(lower), 0, 100) + lower
#     upper <- runif(length(lower), 0, 100) + input
#     lower[c4] <- rep(-Inf, length(lower[c4]))
#     upper[c3] <- rep(+Inf, length(upper[c3]))
#     
#     e <- suppressWarnings( g(input) )
#     o <- gcpp(input, c1, c3, c4, lower, upper)
#     expect_identical(o, e)
#   }
#   
#   # Edge cases
#   c1 <- rep(TRUE, 2)
#   c3 <- rep(FALSE, 2)
#   c4 <- rep(FALSE, 2)
#   lower <- rep(-2, 2)
#   upper <- rep( 5, 2)
#   input <- lower + .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
#   
#   input <- upper - .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
#   
#   # 
#   c1 <- rep(FALSE, 2)
#   c3 <- rep(TRUE, 2)
#   input <- lower + .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
#   
#   input <- upper - .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
#   
#   # 
#   c3 <- rep(FALSE, 2)
#   c4 <- rep(TRUE, 2)
#   input <- lower + .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
#   
#   input <- upper - .Machine$double.eps
#   e <- suppressWarnings( g(input) )
#   o <- gcpp(input, c1, c3, c4, lower, upper)
#   expect_identical(o, e)
# })
# 
# 
# test_that("ginv", {
#   # original function from dfoptim
#   ginv <- function(x) {
#     gix <- x
#     gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
#     gix[c3] <- lower[c3] + exp(x[c3])
#     gix[c4] <- upper[c4] - exp(x[c4])
#     gix
#   }
#   
#   for(x in 1:500) {
#     c1 <- c(TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE)
#     c3 <- c(FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)
#     c4 <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE)
#     lower <- c( 0.0,   2,   4,   8, -100, -200)
#     input <- runif(length(lower), 0, 100) + lower
#     upper <- runif(length(lower), 0, 100) + input
#     lower[c4] <- rep(+Inf, length(lower[c4]))
#     upper[c3] <- rep(-Inf, length(upper[c3]))
#     
#     e <- suppressWarnings( ginv(input) )
#     o <- ginvcpp(input, c1, c3, c4, lower, upper)
#     expect_identical(o, e)
#     if(any(is.nan(e))) expect_identical(o[is.nan(e)], NaN)
#   }
# })
# 
# 
# test_that("Preconditions of nmkbcpp", {
#   rosbkext <- function(x) {
#     n <- length(x)
#     sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
#   }
#   
#   np <- 10
#   set.seed(123)
#   p0 <- rnorm(np)
#   p0[p0 >= +2] <- +2 - 1E-8
#   p0[p0 <= -2] <- -2 + 1E-8
#   
#   badindex <- sample.int(np, 1)
#   p0[badindex] <- -2
#   expect_error(nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2), 
#                "Starting vector 'par' must lie [*]strictly[*] between lower and upper bounds")
#   p0[badindex] <- -2 + .Machine$double.eps
#   expect_silent(nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2))
# })
# 
# 
# 
# test_that("Rosenbrock: nmkb vs Rcpp function 'nmkb'", {
#   rosbkext <- function(x) {
#     n <- length(x)
#     sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
#   }
#   
#   # get original nmkb
#   setwd("..")
#   source("./R/nmkb.rcpp.R")
#   nmkbcpp <- nmkb
#   source("./R/nmkb.R")
#   setwd("./test")
#   
#   np <- 10
#   set.seed(123)
#   p0 <- rnorm(np)
#   p0[p0 > +2] <- +2 - 1E-8
#   p0[p0 < -2] <- -2 + 1E-8
#   e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2)
#   o <- nmkbcpp(fn = rosbkext, par = p0, lower = -2, upper = 2)
#   expect_identical(o, e)
#   
#   ctrl <- list(maxfeval = 50000)
#   p0 <- rnorm(np)
#   p0[p0 > +2] <- +2 - 1E-8
#   p0[p0 < -2] <- -2 + 1E-8
#   set.seed(123)
#   e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
#   set.seed(123)
#   o <- nmkbcpp(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
#   expect_identical(o, e)
# })
# 
# 
# 
# 
# test_that("non-smooth problem: nmkb vs nmkbcpp", {
#   hald <- function(x) {
#     #Hald J & Madsen K (1981), Combined LP and quasi-Newton methods 
#     #for minimax optimization, Mathematical Programming, 20, p.42-62.
#     i <- 1:21
#     t <- -1 + (i - 1)/10
#     f <- (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
#     max(abs(f))
#   }
#   
#   # get original nmkb
#   setwd("..")
#   source("./R/nmkb.rcpp.R")
#   nmkbcpp <- nmkb
#   source("./R/nmkb.R")
#   setwd("./test")
#   
#   p0 <- runif(5, c(0,0,0,0,-2) + 1E-7, 4 - 1E-7)
#   e <- nmkb(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
#   o <- nmkbcpp(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
#   expect_identical(o, e)
# })

