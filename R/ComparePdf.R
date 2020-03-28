#' Mn
#' Determines the n-th raw moment of a sample.
#' @param datasample A random sample
#' @param n The requested moment
#' @return A numeric value representing the requested moment.
#' @author Hugo Hernandez; refactoring by Jan Seifert
#' @references Hernandez, H. (2018). Comparison of Methods for the 
#' Reconstruction of Probability Density Functions from Data Samples.
#' Technical Report, DOI: 10.13140/RG.2.2.30177.35686
#' @examples
Mn <- function(datasample, n = 0) {
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



#' .DeltaPdf
#' Computes the Czekanowski distance between two empirical 
#' distribution functions.
#' @param Edf1 
#' @param Edf1 
#'
#' @return 
.DeltaPdf <- function(Edf1, Edf2) {
  # Integration and standardization (sum of minima / mean of values)
  rhomin <- pmin(Edf1, Edf2) # The smaller value of both
  Diff   <- 2 * sum(rhomin) / (sum(Edf1)+sum(Edf2)) 

  return(Diff)
}

#' SimPdf
#' PDF Similitude Assessment between a given "Pdf1" function and a reference
#' "Pdf2" function.
#' @param Pdf1,Pdf2 Two functions to be compared
#' @param Xmin,Xmax Lower and upper end of the range in which the 
#' distributions are compared (default is from -10 to 10)
#' @param Steps Defines the resolution by specifying the number of 
#' steps to be used between Xmin and Xmax.
#' @param Args1,Args2 list of extra arguments passed on to PDF1 / PDF2.
#' @details Simple symmetrical measure of the distance between 
#' two distributions.
#' @note This function implements the Czekanowski distance. It was 
#' designed as percentage but it can never become zero. Therefore, 
#' the normaing to 100% has been removed from the code. 
#' Ranges from 0 < result <= 1. Equal distributions result as 1.
#' @export
#' @author Hugo Hernandez; refactoring by Jan Seifert
#' @references Hernandez, H. (2018). Comparison of Methods for the 
#' Reconstruction of Probability Density Functions from Data Samples.
#' Technical Report, DOI: 10.13140/RG.2.2.30177.35686
#' @examples
SimPdf <- function( Pdf1, Pdf2, 
                    Xrange = c(-10, 10), Steps = 1e5, 
                    Args1 = NULL, Args2 = NULL ) {
  if(length(Xrange) == 1) Xrange <- c(-Xrange, Xrange) * sign(Xrange)
  x <- seq(Xrange[1], Xrange[2], length.out = Steps+1)
  
  f1 <- match.fun(Pdf1) # find PDF1 in environment
  f2 <- match.fun(Pdf2) # find PDF2 in environment
  # Calculation of pd's
  rho1 <- do.call(f1, c(list(x), Args1))
  rho2 <- do.call(f2, c(list(x), Args2))
  
  Diff <- .DeltaPdf(rho1, rho2)
  
  Result <- Diff
  return(Result)
}
#print(SimPdf(dnorm, dt, Xmin=-100, Xmax=100, Args1 = list(mean = 0), Args2 = list(df = Inf)))
#print(SimPdf(dnorm, dt, Xmin=-100, Xmax=100, Args1 = list(mean = 2), Args2 = list(df = Inf)))
