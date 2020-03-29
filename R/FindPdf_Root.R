source("../R/LowDiscrepancyDistribution.R")

#' This is data to be included in my package
#'
#' @name ByMomentPdf
#' @docType data
#' @format List with classes "ByMomentPdf" and a second class to indicate
#' the generative function.
#' @author Jan Seifert 
#' @references 
#' @keywords data
NULL


New_ByMomentPdf <- function( TarMo ) {
  New_ByMomentPdf.default( TarMo )
}


New_ByMomentPdf.default <- function( TarMo ) {
  this <- list(
    # Pdf
    Function  = NULL,   # PDF: NULL if unknown or a list(Func, Args)
    Moments   = NULL,   # Moments
    # Target function
    TarFu     = NULL,   # PDF: NULL if unknown or a list(Func, Args)
    TarMo     = TarMo,  # Target moments
    # Distance to target
    DistaFu   = NULL,   # Distance between PDFs
    DistaMo   = NULL,   # Distance between moments
    # Dimensions of the parameter space of the PDF
    ParamSpace = matrix(data = c(-Inf, Inf), ncol = 2, 
                        dimnames = list(NULL, c("from", "to"))),
    # Starting points for approximation algorithms
    LaunchSpace = matrix(numeric(0)) # Matrix of coordinate vectors
  ) 
  
  class(this) <- append(class(this), "ByMomentPdf")
  return(this)
}


GetLaunchSpace <- function( Pdf, ... ) {
  UseMethod( "GetLaunchSpace" )
}

GetLaunchSpace.ByMomentPdf <- function( Pdf, Count,
                                    Method = c("Random", "Even", "Manual") ) {
  # PRECONDITIONS
  Method <- match.arg(Method)
  if (any(is.infinite(Pdf$ParamSpace))) 
    stop("No parameter space specified") # 
  
  # Dimensions of parameter space
  NDim <- ncol(Pdf$ParamSpace)
  
  # RUN
  # Random: 1 to n randomly set points
  if (Method == "Random") {
    LaunchPad <- numeric(0)
    for(c in 1:Count) { # rep(list(c(from = -0.25, to = 25)), NDim)
      Point <- runif(NDim, 
                         min = Pdf$ParamSpace[, "from"],
                         max = Pdf$ParamSpace[, "to"])
      LaunchPad <- cbind(LaunchPad, Point)
    }
    Pdf$LaunchSpace <- LaunchPad
  }
  
  # Even: distribute n starting points evenly across the space
  if (Method == "Even") {
    LaunchPad <- GetHarmonicPoints(NDim, Count)
    # Scale points to 'ParamSpace'
    Range <- Pdf$ParamSpace[, "to"] - Pdf$ParamSpace[, "from"]
    LaunchPad <- LaunchPad + Pdf$ParamSpace[, "from"]
    LaunchPad <- LaunchPad * Range
    Pdf$LaunchSpace <- LaunchPad
  }
  
  # Manual: select starting points in heat map
  if (Method == "Manual") {
    LaunchPad <- list()
    for(c in 1:Count) { # rep(list(c(from = -0.25, to = 25)), NDim)
      stop("not yet implemented")
    }
    Pdf$SpaceLaunch <- LaunchPad
  }
  
  return(Pdf)
}


### TESTING #####
x <- New_ByMomentPdf.default(c(0, 1))
x$ParamSpace <- matrix(c(-5, -5, 5, 5), ncol = 2,
                           dimnames = list(NULL, c("from", "to")))
GetLaunchSpace( x, 2L, "R")
