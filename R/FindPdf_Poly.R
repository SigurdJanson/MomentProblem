#' polyPDF
#' This function constructs a PD-function from the moments.
#' 
#' The output of this function is the PDF.
#' Moments can be generated using the "momentset" function given in Appendix A.2.
#' @param moments data frame, containing the moment order ("n") 
#' in one column and their values ("moments") in the other. 
#' The data frame must be labelled.
#' @param xmin,xmax Lower and upper limit of the PDF.
#' @param scale Scale factor can be used to transform data avoiding 
#' singularity (default = 1).
#' @param disp If TRUE, the PDF will be plotted and the polynomial 
#' coefficients are shown. 
#' @author Hugo Hernandez; small refactoring by Jan Seifert
#' @reference Hernandez, H. (2018). Comparison of Methods for the Reconstruction 
#' of Probability Density Functions from Data Samples. Technical Report 
#' @source https://t1p.de/wyhb
polyPDF <- function( moments, xmin = NULL, xmax = NULL, scale = 1, 
                     plot = FALSE) {
  if (is.null(xmin) || is.null(xmax))
    stop("Please input 'xmin' and 'xmax' estimated values")
  
  n <- length(moments$n) #Number of moments available
  degree <- n-1          #Degree of polynomial
  A <- matrix(0,n,n)     #Initialize matrix of coefficients
  B <- matrix(0,n,1)     #Initialize vector of independent terms

  for (i in 1:n) { #For each moment
    ni <- moments$n[i] #ni-th moment
    for (j in 1:n){ #For each power term
      A[i,j] <- (((xmax*scale)^(ni+j))-((xmin*scale)^(ni+j)))/(ni+j)
    }
    B[i] <- moments$moments[i] * scale^ni #Scale moments
  }
  a <- as.vector(solve(A,B)) #Find coefficients
  
  # Definition of PDF function
  rhof <- function(x){
    rho <- a[1] #Independent term
    for (i in 2:(degree+1)) {
      rho <- rho + a[i] * (x*scale)^(i-1) #Polynomial terms
    }
    
    #Density is zero when rho is negative or X is beyond boundaries
    rho <- rho * scale * as.integer(x >= xmin & x <= xmax) * as.integer(rho>0)
    return(rho)
  } #end of function
  
  if (isTRUE(plot)) {
    # Set coefficients names
    names(a)[1:2] <- c("(Intercept)", "x")
    if (length(a) > 2) {
      for (i in 3:length(a)) {
        names(a)[i] <- paste("x^", toString(i-1))
      }
    }
    cat("Polynomial model of the PDF:\n")
    print(a)
    
    #PDF plot
    i <- 1:1001
    y <- xmin + (xmax-xmin) * (i-1)/1000
    rho <- rhof(y) #Calculate density
    plot(y, rho, type="l", col="blue", xlim = c(xmin, xmax), ylim = c(0, max(rho)),
         xlab = "Measurement values", ylab = "Probability Density")
  }
  return(rhof)
}


#ListOfMoments <- list(data.frame(cbind(n = 1:8, moments = c(0, 1, 0, 3, 0, 15, 0, 105)))) #, 0, 945
#polyPDF(ListOfMoments[[1]], -16, +16, disp = TRUE)
#rm(ListOfMoments, polyPDF)