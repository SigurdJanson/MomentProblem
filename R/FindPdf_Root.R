source("../R/HarmonicScatter.R")
source("../R/CubicScatter.R")

#' Data structure contaning the information required to approximate
#' a PDF given specific moments.
#'
#' @name ByMomentPdf
#' @docType data
#' @format List with classes "ByMomentPdf" and a second class to indicate
#' the generative function.
#' \describe{
#'   \item{Function}{PDF}
#'   \item{Moments}{Observed moments (numeric vector).}
#'   \item{TarFu}{A function the new PDF shall mimic. Usually `NULL` 
#'   because unknown.}
#'   \item{TarMo}{The desired moments to be approximated 
#'   (numeric vector).}
#'   \item{DistaFu}{Distance measure between approximated function and 
#'   target function.}
#'   \item{DistaMo}{Distance measure between approximated moments and 
#'   target function.}
#'   \item{ParamSpace}{Range of definition that is allowed for 
#'   parameters of the PDF. Numeric matrix with two coordinate 
#'   vectors in the rows 'from' and 'to'. The number of columns 
#'   equals the number of parameters, i.e. the dimensions of the 
#'   paramter space.}
#'   \item{LaunchSpace}{Starting points in the parameter space 
#'   for approximative algorithms. Numeric matrix with as many
#'   rows as there are starting points. The number of columns 
#'   equals the number of parameters, i.e. the dimensions of the 
#'   paramter space.}
#' }
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
    ParamSpace = matrix(data = c(-Inf, Inf), nrow = 2,
                        dimnames = list(c("from", "to"), NULL)),
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
                                        Method = c("Random", "Cubic", 
                                                   "Harmonic", "Manual") ) {
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
    for(c in 1:Count) { 
      Point <- runif(NDim, 
                     min = Pdf$ParamSpace["from", ],
                     max = Pdf$ParamSpace["to", ])
      LaunchPad <- rbind(LaunchPad, Point)
    }
    dimnames(LaunchPad) <- NULL
    Pdf$LaunchSpace <- LaunchPad
  }
  
  # Cubic: distribute n starting points evenly across the space
  if (Method == "Cubic") {
    LaunchPad <- GetCubicPoints(NDim, Count)
    # Scale points to 'ParamSpace'
    Range <- Pdf$ParamSpace["to", ] - Pdf$ParamSpace["from", ]
    LaunchPad <- LaunchPad * Range
    LaunchPad <- LaunchPad + Pdf$ParamSpace["from", ]
    Pdf$LaunchSpace <- LaunchPad
  }
  
  # Harmonic: distribute n starting points evenly across the space
  if (Method == "Harmonic") {
    LaunchPad <- GetHarmonicPoints(NDim, Count)
    # Scale points to 'ParamSpace'
    Range <- Pdf$ParamSpace["to", ] - Pdf$ParamSpace["from", ]
    LaunchPad <- LaunchPad * Range
    LaunchPad <- LaunchPad + Pdf$ParamSpace["from", ]
    Pdf$LaunchSpace <- LaunchPad
  }
  
  # Manual: select starting points in heat map
  if (Method == "Manual") {
    LaunchPad <- list()
    for(c in 1:Count) { 
      stop("not yet implemented")
    }
    Pdf$SpaceLaunch <- LaunchPad
  }
  
  return(Pdf)
}


### TESTING #####
x <- New_ByMomentPdf.default(c(0, 1))
x$ParamSpace <- matrix( c(-5, -5, -5, 5, 5, 5), ncol = 3, byrow = TRUE,
                            dimnames = list(c("from", "to"), NULL) )
x <- GetLaunchSpace( x, Count = 2L, "Random")
