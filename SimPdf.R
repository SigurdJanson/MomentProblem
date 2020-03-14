#' Mn
#' Determination of the n-th moment of a sample.
#' @param datasample A random sample
#' @param n The requested moment
#' @return A numeric value representing the requested moment.
#' @author Hugo Hernandez; refactoring by Jan Seifert
#' @references Hernandez, H. (2018). Comparison of Methods for the 
#' Reconstruction of Probability Density Functions from Data Samples.
#' Technical Report, DOI: 10.13140/RG.2.2.30177.35686
#' @examples
Mn <- function(datasample, n=0) {
  N <- length(datasample) #Sample size
  Mn <- 0 #Initialize moment
  for (i in 1:N){
    Mn <- Mn + (datasample[i])^n #Accumulate moment
  }
  Mn <- Mn/N #Averaging over sample size
  return(Mn)
}


#' MomentSet
#' Determination of the first integer moments from a sample
#' Constructing a set of integer moments from 0 to "nmax" for a sample
#'
#' @param datasample 
#' @param nmax Number of moments to be computed. First 10 moments by default.
#' @return A data frame with two columns. The first column denotes the 
#' order of the moment and the second column it's value.
#' @export
#' @author Hugo Hernandez; refactoring by Jan Seifert
#' @references Hernandez, H. (2018). Comparison of Methods for the 
#' Reconstruction of Probability Density Functions from Data Samples.
#' Technical Report, DOI: 10.13140/RG.2.2.30177.35686
#' @examples MomentSet(rnorm(1E4), nmax=8)
MomentSet <- function(datasample, nmax=10) { 
  moments <- 1 # Zero-th moment value
  n       <- 0 # Zero-th moment order
  for (i in 1:nmax) {
    moments[i+1] <- Mn(datasample, i) #n-th moment value
    n[i+1] <- i                       #n-th moment order
  }
  MS <- data.frame(n, moments) #Constructing output data frame
  return(MS)
}


#' SimPdf
#' PDF Similitude Assessment between a given "PDF1" function and a reference
#' "PDF2" function. "nsteps" between "xmin" and "xmax" will be used for the
#' assessment. By default the assessment range is from -10 to 10.
#' @param PDF1,PDF2 Two functions
#' @param xmin,xmax Lower and upper end of the range in which the 
#' distributions are compared
#' @param nsteps 
#' @param Args1,Args2 Extra arguments passed on to PDF1 / PDF2
#' @return Named numeric vector containing the percentage of similitude.
#' @export
#' @author Hugo Hernandez; refactoring by Jan Seifert
#' @references Hernandez, H. (2018). Comparison of Methods for the 
#' Reconstruction of Probability Density Functions from Data Samples.
#' Technical Report, DOI: 10.13140/RG.2.2.30177.35686
#' @examples
SimPdf <- function( PDF1, PDF2, Xmin=-10, Xmax=10, Steps=1e5, Args1 = NULL, Args2 = NULL ) {
  x <- seq(Xmin, Xmax, Steps)
  
  f1 <- match.fun(PDF1) # find PDF1 in environment
  f2 <- match.fun(PDF2) # find PDF2 in environment
  # Calculation of pd's
  rho1 <- do.call(f1, c(list(x), Args1))
  rho2 <- do.call(f2, c(list(x), Args2))
  
  rhomin <- pmin(rho1, rho2) # Minimum PD
  
  # Integration and standardization (200 = 100 (to get %) * 2 * ...)
  simil <- 200 * sum(rhomin) / (sum(rho1)+sum(rho2)) 
  names(simil) = c("Similitude (%)")
  
  return(simil)
}

