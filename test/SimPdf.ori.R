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
