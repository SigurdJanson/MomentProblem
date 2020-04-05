source("./R/HarmonicScatter.R")
source("./R/CubicScatter.R")

#' Data structure contaning the information required to approximate
#' a PDF given specific moments.
#'
#' @name ByMomentPdf
#' @docType data
#' @format List with classes `ByMomentPdf` and a second class to indicate
#' the generative function.
#' \describe{
#'   \item{Function}{String holding the base name of generating 
#'   distribution function. Must follow the R conventions of probability
#'   functions, i.e. for a uniform distribution `Function == "unif"`.
#'   \item{Moments}{{A list of moments for each solution. 
#'   The first list element contains a reference to the solution. The
#'   second holds a vectors of moments.}}
#'   \item{ParamSolved}{A list of parameters for each solution. 
#'   The first list element contains a reference to the solution. The
#'   second holds a vectors of parameters.}
#'   \item{TarFu}{A function the new PDF shall mimic. Usually `NULL` 
#'   because unknown. List with the format `list(Function, List-of-Args)`.}
#'   \item{TarMo}{The desired moments to be approximated 
#'   (numeric vector). Accordingly, `Pdf$TarMo[1]` is the mean, `Pdf$TarMo[2]`
#'   the variance, etc. The length (i.e. the highest allowed moment) is
#'   defined by each implementation of `ByMomentPdf`.}
#'   \item{DistaFu}{Distance measure between approximated solution 
#'   function and target function.}
#'   \item{DistaMo}{Distance measure between approximated solution 
#'   moments and target function.}
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
#'   \item{Tolerance}{Desired precision for converging algorithms.
#'   Default: `sqrt(.Machine$double.eps)`.}
#' }
#' @author Jan Seifert 
#' @references 
#' @keywords data



#' New_ByMomentPdf
#' Constructor for parent class `ByMomentPdf`.
#' @param TarMo Numeric vector of target moments.
#' @return Class `ByMomentPdf`
#' @export
#' @author Jan Seifert 
#' @examples
New_ByMomentPdf <- function( TarMo ) {
  New_ByMomentPdf.ByMomentPdf( TarMo )
}

#' New_ByMomentPdf.ByMomentPdf
#' @describeIn New_ByMomentPdf
New_ByMomentPdf.ByMomentPdf <- function( TarMo ) {
  this <- list(
    # Pdf
    Function    = NULL,   # PDF: NULL if unknown or a list(Func, Args)
    Moments     = NULL,   # Moments
    ParamSolved = NULL,
    # Target function
    TarFu       = NULL,   # PDF: NULL if unknown or a list(Func, Args)
    TarMo       = TarMo,  # Target moments
    # Distance to target
    DistaFu     = NULL,   # Distance between PDFs
    DistaMo     = NULL,   # Distance between moments
    # Dimensions of the parameter space of the PDF
    ParamSpace  = matrix(data = c(-Inf, Inf), nrow = 2,
                        dimnames = list(c("from", "to"), NULL)),
    # Starting points for approximation algorithms
    LaunchSpace = matrix(numeric(0)), # Matrix of coordinate vectors
    Tolerance   = sqrt(.Machine$double.eps)
  ) 
  
  class(this) <- append(class(this), "ByMomentPdf")
  return(this)
}


# Probability distribution functions ----
dPdf <- function(Pdf, x, ParamSet, ...) {
  UseMethod("dPdf")
}

dPdf.ByMomentPdf <- function(Pdf, x, ParamSet, ...) {
  FName <- paste0("d", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(x = x), c(Param, list(...))) )
}

pPdf <- function(Pdf, q, ParamSet, ...) {
  UseMethod("pPdf")
}

pPdf.ByMomentPdf <- function(Pdf, q, ParamSet, ...) {
  FName <- paste0("p", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(q = q), Param, ...) )
}

qPdf <- function(Pdf, p, ParamSet, ...) {
  UseMethod("qPdf")
}

qPdf.ByMomentPdf <- function(Pdf, p, ParamSet, ...) {
  FName <- paste0("q", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(p = p), Param, ...) )
}

rPdf <- function(Pdf, p, ParamSet, ...) {
  UseMethod("rPdf")
}

rPdf.ByMomentPdf <- function(Pdf, n, ParamSet, ...) {
  FName <- paste0("r", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(n = n), Param, ...) )
}



SolutionMoments <- function(Pdf) {
  UseMethod("SolutionMoments")
}



#' GetLaunchSpace
#'
#' @param Pdf An object of class `ByMomentPdf`.
#' @param Count Number of desired launch points.
#' @param Method A method to generate launch points.
#'
#' @return A class after adding the optimisation solution. 
#' The classes of `Pdf` are preserved.
#' @export
#'
#' @examples
GetLaunchSpace <- function( Pdf, ... ) {
  UseMethod( "GetLaunchSpace" )
}


#' GetLaunchSpace.ByMomentPdf
#' @describeIn GetLaunchSpace
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



#' AddSolution
#' Adds a new optimisation result to the solutions. 
#' @param Pdf An object of class `ByMomentPdf`.
#' @param LaunchPoint The index of a launching point referencing 
#' the point in `Pdf$LaunchSpace`.
#' @param SoluParam Numeric vector holding the parameters of
#' the pd-generating function.
#' @param Add If `FALSE` existing solutions will be overwritten (default).
#' Otherwise the new solution will be inserted in the list.
#' 
#' @note `AddSolution` makes sure that the list of solution is in the 
#' right order of launch points.
#' @return A class after adding the optimisation solution. 
#' The classes of `Pdf` are preserved.
#' @export
#'
#' @examples
AddSolution <- function( Pdf, LaunchPoint, SoluParam, Add = FALSE, ... ) {
  # RESULT
  UseMethod("AddSolution")
}


#' AddSolution.ByMomentPdf
#' @describeIn AddSolution
AddSolution.ByMomentPdf <-  function( Pdf, LaunchPoint, 
                                      SoluParam, Append = FALSE ) {
  # RESULT
  if (isFALSE(Append) ||
      is.null(Pdf$ParamSolved) || nrow(Pdf$Pdf$ParamSolved) <= 1) {
    Pdf$ParamSolved <- list(LaunchPoint, SoluParam)
  } else {
    # Check if already exists
    LPs <- lapply(x, `[[`, 1) # get existing launch points
    Pos <- which(LPs == LaunchPoint)
    if (length(Pos) != 0) {
      # Replace if it does
      Pdf$ParamSolved[[Pos]] <- list(LaunchPoint, SoluParam)
    } else {
      # Find position and insert it at the right position
      Pos <- which.max(LaunchPoint < LPs) - 1
      Pdf$ParamSolved <- append(Pdf$ParamSolved, 
                                list(LaunchPoint, SoluParam), 
                                after = Pos)
    }
  }
  #NOTE: At this point this class assumes that there is one distance
  # for each solution. If not, they cannot be referenced properly.
  Pdf$DistaFu <- NULL # not valid anymore
  Pdf$DistaMo <- NULL # not valid anymore
  return(Pdf)
}


#' FindPdf
#' Find solutions for probability-distribution generating 
#' function of `Pdf`.
#' Generic function for all sub-classes of `ByMomentPdf`.
#' @param Pdf An object of class `ByMomentPdf` or a sub-class.
#' @param LaunchPoint Index of the launch point is a reference to
#' the row of `Pdf$SpaceLaunch`.
#' @param ... Additional arguments to be passed to or from methods.
#' @note This method has no default and must be implemented by 
#' each child-class.
#' @return A class after adding the optimisation solution. 
#' The classes of `Pdf` are preserved.
#' @export
#'
#' @examples
FindPdf <- function( Pdf, LaunchPoint, Append = FALSE, ... ) {
  # PRECONDITIONS
  if(LaunchPoint < 1) stop("Invalid launch point (must be > 0).")
  if(LaunchPoint > nrow(Pdf$SpaceLaunch))
    stop("Index of launch point out of range.")
  
  # RESULT
  UseMethod("FindPdf")
}



#' EvaluatePdf
#' Assess the quality of solutions stored in `Pdf`. 
#' Generic function for all sub-classes of `ByMomentPdf`.
#' @param Pdf An object of class `ByMomentPdf` or a sub-class.
#' @param usePdf If `TRUE` the distance between found and desired solutions
#' is determined using the Pdf. If `FALSE` only the moment vectors will be 
#' used. Default is `FALSE`.
#' @param ... Additional arguments to be passed to or from methods.
#' @description By default, the Pdf itself is used to assess the quality
#' of any solutions. The range of evaluation is ± 5 stddev around the
#' mean. The distance function is the Jeffrey divergence which can be 
#' described as symmetric version of the Kullback-Leibler distance. The
#' `[jeffreys]()` function of the '[philentropy]' package is used to get it.
#' 
#' The distance between moment sets is an euclidean distance as given 
#' by Rs `[dist]()` function.
#' @note The code is not very efficient and may not be suited for highly 
#' iterative automated tasks.
#' @return A class after adding the optimisation solution. 
#' The classes of `Pdf` are preserved.
#' @references Deza, E., & Deza, M. M. (2009). Encyclopedia of Distances. 
#' Springer. https://doi.org/10.1007/978-3-642-00234-2
#' @export
#'
#' @examples
EvaluatePdf <- function( Pdf, ... ) {
  UseMethod("EvaluatePdf")
}

#' EvaluatePdf.ByMomentPdf
#' @describeIn EvaluatePdf
EvaluatePdf.ByMomentPdf <- function( Pdf, UsePdf = FALSE ) {
  if(is.null(Pdf$ParamSolved)) 
    stop("No solutions available in PDF")
  if(UsePdf && is.null(Pdf$TarFu)) 
    stop("'UsePdf' demands a target function in 'Pdf'")

  # Compute the distances
  if (UsePdf) {
    if(require(philentropy, quietly = TRUE) == FALSE)
      stop("Package 'philentropy' is not available")
    # 
    BaseRange <- c(-5, 5)
    EvalRange <- Pdf$TarMo[1] + BaseRange * Pdf$TarMo[2]
    EvalResolution <- diff(BaseRange) * 1E3
    x <- seq(EvalRange[1], EvalRange[2], length.out = EvalResolution)
    #
    Solutions <- lapply(Pdf$ParamSolved, `[`, 2)
    DResult <- list()
    #
    for(s in Solutions) {
      SoluResult  <- dPdf(x, s) # get PDF(x)
      TarFuResult <- do.call( Pdf$TarFu[[1]], append(list(x = x), Pdf$TarFu[2]) )
      D <- jeffreys(TarFuResult, SoluResult, TRUE, "log2")
      DResult <- append(DResult, list(D))
    }
    Pdf$DistaFu <- DResult
  } else { # Use moments
    # if moments are not available, generate
    if(is.null(Pdf$Moments)) {
      Pdf$Moments <- list()
      for(i in 1:length(Pdf$ParamSolved)) {
        Pdf$Moments <- append(Pdf$Moments, list(i, SolutionMoments(Pdf, i)))
      }
    }
    
    Moments <- lapply(Pdf$Moments, `[[`, 2)
    DResult <- list()
    Method <- 1 # i.e. "euclidean", see "?dist"
    Attrs <- list(Size = 2, Labels = "", Diag = FALSE, Upper = FALSE,
                  method = "euclidean", call = match.call(), class = "dist")
    for(m in Moments) {
      X <- rbind(Pdf$TarMo, m)
      D <- .Call(C_Cdist, X, Method, Attrs)
      DResult <- append(DResult, list(m[[1]], D))
    }
    Pdf$DistaMo <- DResult
  }
}

## TODO ####
# - documentation of pdf functions
# 

### TESTING #####
x <- New_ByMomentPdf.ByMomentPdf(c(0, 1))
x$ParamSpace <- matrix( c(-5, -5, -5, 5, 5, 5), ncol = 3, byrow = TRUE,
                           dimnames = list(c("from", "to"), NULL) )
x <- GetLaunchSpace( x, Count = 2L, "Random")


#' @param TarMo Target moments as vector. See details.
#' @param Method Function generating algorithm. Available are the 
#' generalised lambda distribution (`gld`), the Pearson distribution
#' (`pearson`) or a polynomial solution (`poly`).
#' @param Launch
#' @details `TarMo` expects the moments from 1 to m, zero-order moment
#' not included.
