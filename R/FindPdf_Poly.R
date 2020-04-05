
New_ByMomentPdf.poly <- function( TarMo, TarFu = NULL ) {
  Pdf <- New_ByMomentPdf.default( TarMo, TarFu = NULL )

  Pdf$LaunchSpace <- matrix(c(xmin = TarMo[1]-TarMo[2]*5, 
                              xmax = TarMo[1]+TarMo[2]*5, 
                              scale = 1), ncol = 3)
  class(this) <- append(class(this), "poly")
}


dPdf.poly <- function(Pdf, x, ParamSet, ...) {
  if(is.list(ParamSet)) Param <- unlist(ParamSet)
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  xmin  <- Param["xmin"]
  xmax  <- Param["xmax"]
  scale <- Param["scale"]
  
  n <- length(Pdf$TarMo) # Number of moments available
  degree <- n-1          # Degree of polynomial
  rho <- Param[1]        # Independent term
  for (i in 2:(degree+1)) {
    rho <- rho + Param[i] * (x*scale)^(i-1) #Polynomial terms
  }
  
  # Density is zero when rho is negative or X is beyond boundaries
  rho <- rho * scale * as.integer(x >= xmin & x <= xmax) * as.integer(rho > 0)
  return(rho)
}

pPdf.poly <- function(Pdf, q, ...) {
  #TODO
}

qPdf.poly <- function(Pdf, p, ...) {
  #TODO
}

rPdf.poly <- function(Pdf, n, ...) {
  #TODO
}


SetLaunchSpace.poly <- function( Pdf, LaunchSpace ) {
  if (!is.matrix(LaunchSpace)) 
    LaunchSpace <- as.matrix(LaunchSpace, ncol = 3)
  if (ncol(LaunchSpace) != 3) 
    stop("Wrong launch space for polynomial estimate.")

  Pdf$LaunchSpace <- LaunchSpace
}


#' FindPdf.poly
#' This function constructs a PD-function from the moments.
#' 
#' The output of this function is the PDF.
#' @param Pdf `ByMomentPdf` object. 
#' The data frame must be labelled.
#' @param xmin,xmax Lower and upper limit of the PDF.
#' @param scale Scale factor can be used to transform data avoiding 
#' singularity (default = 1).
#' @param disp If `TRUE`, the PDF will be plotted and the polynomial 
#' coefficients are shown. 
#' @author Hugo Hernandez; small refactoring by Jan Seifert
#' @reference Hernandez, H. (2018). Comparison of Methods for the Reconstruction 
#' of Probability Density Functions from Data Samples. Technical Report 
#' @source https://t1p.de/wyhb
FindPdf.poly <- function( Pdf, xmin = NULL, xmax = NULL, scale = 1) {
  if(length(Pdf$LaunchSpace) == 0) {
    if (is.null(xmin) || is.null(xmax))
      stop("Pdf range 'xmin' and 'xmax' not fully specified.")
  } else {
    xmin  <- Pdf$LaunchSpace["xmin"]
    xmax  <- Pdf$LaunchSpace["xmax"]
    scale <- Pdf$LaunchSpace["scale"]
  }
  
  n <- length(Pdf$TarMo) #Number of moments available
  degree <- n-1          #Degree of polynomial
  A <- matrix(0, n, n)   #Initialize matrix of coefficients
  B <- matrix(0, n, 1)   #Initialize vector of independent terms

  for (m in 1:n) { #For each moment
    #_m <- m #moments$n[i] #ni-th moment
    for (pt in 1:n) { 
      #For each power term
      A[m,pt] <- (((xmax*scale)^(m+pt))-((xmin*scale)^(m+pt))) / (m+pt)
    }
    B[m] <- Pdf$TarMo[m] * scale^m #Scale moments
  }
  a <- as.vector(solve(A,B)) #Find coefficients
  
  # if (isTRUE(plot)) {
  #   # Set coefficients names
  #   names(a)[1:2] <- c("(Intercept)", "x")
  #   if (length(a) > 2) {
  #     for (i in 3:length(a)) {
  #       names(a)[i] <- paste("x^", toString(i-1))
  #     }
  #   }
  #   cat("Polynomial model of the PDF:\n")
  #   print(a)
  #   
  #   #PDF plot
  #   #-i <- 1:1001
  #   #-y <- xmin + (xmax-xmin) * (i-1)/1000
  #   y <- seq(xmin, xmax, length.out = 1000)
  #   
  #   rho <- rhof(y) #Calculate density
  #   plot(y, rho, type="l", col="blue", 
  #        xlim = c(xmin, xmax), ylim = c(0, max(rho)),
  #        xlab = "Measurement values", ylab = "Probability Density")
  # 
  #   invisible(rhof)
  # } else {
    return(rhof)
  # }
}


#ListOfMoments <- list(data.frame(cbind(n = 1:8, moments = c(0, 1, 0, 3, 0, 15, 0, 105)))) #, 0, 945
#polyPDF(ListOfMoments[[1]], -16, +16, disp = TRUE)
#rm(ListOfMoments, polyPDF)