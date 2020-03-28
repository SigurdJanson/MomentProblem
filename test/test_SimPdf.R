library(testthat)
source("../R/SimPdf.R")

test_that("SimPdf - Compare with original version", {
  # Make sure that refactoring didn't have any side effects
  # Original function before refactoring
  # Almost original because the x <- ... expression had to be reordered
  # In it's original order there were deviations because of numeric precision
  SimPdfOri <- function( PDF1, PDF2, Xmin=-10, Xmax=10, Steps=1e5 ) {
    i <- 1:(Steps+1) #Step counter
    x <- Xmin + (Xmax-Xmin) / Steps *(i-1) #Location of steps

    f1 <- match.fun(PDF1) #Definition of pd function 1
    f2 <- match.fun(PDF2) #Definition of pd function 2
    rho1 <- f1(x) #Calculation of pd 1 at each step
    rho2 <- f2(x) #Calculation of pd 2 at each step
    rhomin <- pmin(rho1,rho2) #Minimum PD
    simil <- 200 * sum(rhomin) / (sum(rho1)+sum(rho2)) #Integration and standardization
    names(simil) = c("Similitude (%)")

    return(simil)
  }

  # Normal vs. Normal - Identical results
  for(R in c(2, 5, 10, 20)) {
    Range <- c(-R, R)
    o <- SimPdf(dnorm, dnorm, Range)
    e <- SimPdfOri(dnorm, dnorm, -R, R) / 100
    names(e) <- NULL
    expect_equal(o, e, info = paste("R =", R))
  }
  
  # Normal vs. Uniform
  for(R in c(1, 2, 5, 10, 20)) {
    Range <- c(-R, R)
    o <- SimPdf(dnorm, dunif, Range)
    e <- SimPdfOri(dnorm, dunif, -R, R) / 100
    names(e) <- NULL
    expect_equal(o, e, info = paste("R =", R))
  }
  
  # Normal vs. Log-Normal
  o <- SimPdf(dnorm, dlnorm)
  e <- SimPdfOri(dnorm, dlnorm) / 100
  names(e) <- NULL
  expect_equal(o, e)
})
  


test_that("SimPdf: ", {
  o <- SimPdf(dnorm, dnorm)
  e <- 1
  expect_identical(o, e)
  
  o <- SimPdf( dnorm, dt, Args2 = list(df = Inf))
  e <- 1
  expect_identical(o, e)
})




test_that("Moments", {
  # Original function before refactoring
  MnOri <- function(datasample, n=0) {
    N=length(datasample) #Sample size
    Mn=0 #Initialize moment
    for (i in 1:N){
      Mn=Mn+(datasample[i])^n #Accumulate moment
    }
    Mn=Mn/N #Averaging over sample size
    return(Mn)
  }


  # Original function before refactoring
  momentsetOri <- function(datasample, nmax=10){ #First 10 moments by default
    #Constructing a set of integer moments from 0 to "nmax" for a sample
    #Requires the function "Mn" defined in Appendix A.1.
    moments=c(1) #Zero-th moment value
    n=c(0) #Zero-th moment order
    for (i in 1:nmax){ #For each subsequent moment
      moments[i+1]=MnOri(datasample,i) #n-th moment value
      n[i+1]=i #n-th moment order
    }
    MS=data.frame(n, moments) #Constructing output data frame
    return(MS)
  }

  r <- rnorm(1000)
  m <- MomentSet(r, nmax=8)$moments
  e <- momentsetOri(r, nmax=8)$moments
  expect_equal(m, e)

})
