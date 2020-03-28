

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
    ParamSpace = NULL,
    # Starting points for approximation algorithms
    SpaceLaunch = NULL # List of coordinate vectors
    ) 
  
  class(this) <- append(class(this), "ByMomentPdf")
  return(this)
}




GetLaunchSpace.default <- function( PdfObj, Count,
                                    Method = c("Random", "Even", "Manual") ) {
  # PRECONDITIONS
  Method <- match.arg(Method)
  
  # Dimensions of parameter space
  NDim <- length(PdfObj$ParamSpace[[1]])
  
  # RUN
  # Random: 1 to n randomly set points
  if (Method == "Random") {
    LaunchPad <- list()
    for(c in 1:Count) { # rep(list(c(from = -0.25, to = 25)), NDim)
      LaunchPad[c] <- runif(NDim, 
                            min = PdfObj$ParamSpace[1:NDim]["from"],
                            max = PdfObj$ParamSpace[1:NDim]["to"])
    }
    PdfObj$SpaceLaunch <- LaunchPad
  }
  
  # Even: distribute n starting points evenly across the space
  if (Method == "Even") {
    LaunchPad <- list()
    for(c in 1:Count) { # rep(list(c(from = -0.25, to = 25)), NDim)
      LaunchPad[c] <- lapply(PdfObj$ParamSpace, 
                             function(x) seq(x["from"], x["to"], length.out = NDim)
      )
      #TODO: verify: has LunchPad the right format?
    }
    PdfObj$SpaceLaunch <- LaunchPad
  }
  
  # Manual: select starting points in heat map
  if (Method == "Manual") {
    LaunchPad <- list()
    for(c in 1:Count) { # rep(list(c(from = -0.25, to = 25)), NDim)
      stop("not yet implemented")
    }
    PdfObj$SpaceLaunch <- LaunchPad
  }
  
  return(PdfObj)
}

