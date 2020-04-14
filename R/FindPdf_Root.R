source("./R/HarmonicScatter.R")
source("./R/RandomScatter.R")
source("./R/CubicScatter.R")
source("./R/Helper_FindPdf_Root.R")

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
#'   second holds a vector of parameters.}
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
#' @param TarFu Name of target function as string.
#' @return Class `ByMomentPdf`
#' @export
#' @author Jan Seifert 
#' @examples
New_ByMomentPdf <- function( TarMo = NULL, TarFu = NULL ) {
  if (is.null(TarMo) && is.null(TarFu)) 
    stop("At least one of 'TarMo' or 'TarFu' must be given.")
  if (length(TarMo) < 2)
    stop("Mean and variance must be given at least.")
  if (any(TarMo[c(FALSE, TRUE)] < 0)) # check even positions
    stop("Moments with even power must be >= 0.")
  
  New_ByMomentPdf.default( TarMo, TarFu )
}

#' New_ByMomentPdf.ByMomentPdf
#' @describeIn New_ByMomentPdf
New_ByMomentPdf.default <- function( TarMo = NULL, TarFu = NULL ) {
  this <- list(
    # Pdf
    Function    = NULL,   # NULL if unknown or a list(Func, Args)
    Moments     = NULL,   # Moments
    ParamSolved = NULL,
    # Target FUnction and MOments
    TarFu       = TarFu,  # NULL if unknown or a list(Func, Args)
    TarMo       = TarMo,  # NULL if unknown or a numeric vector
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
#' Wrappers to call the underlying probability distribution, 
#' density, distribution function, quantile function and random generation. 
#' function. Calls the according function in `Pdf$Function`.
#' @param Pdf A `ByMomentPdf` object.
#' @param x,q Vector of quantiles.
#' @param p vector of probabilities.
#' @param n Number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param ParamSet The arguments passed on to the function. Either a list
#' or the index of the correct parameter in `Pdf$ParamSolved`.
#' @param ... Further arguments passed on to `Pdf$Function`.
#' @return Returns the result of the distribution function.
#' @export
#'
#' @author Jan Seifert
#' @examples
dPdf <- function(Pdf, x, ParamSet, ...) {
  UseMethod("dPdf")
}

dPdf.ByMomentPdf <- function(Pdf, x, ParamSet, ...) {
  FName <- paste0("d", Pdf$Function)
  if(is.list(ParamSet)) 
    Param <- ParamSet
  else {
    Indices   <- lapply(Pdf$ParamSolved, `[[`, 1)
    Solutions <- lapply(Pdf$ParamSolved, `[[`, 2)
    Param <- unlist(Solutions[Indices == ParamSet])
  }

  do.call( FName, c(list(x = x), c(Param, list(...))) )
}

#' pPdf
#' @describeIn dPdf Gives the distribution function
pPdf <- function(Pdf, q, ParamSet, ...) {
  UseMethod("pPdf")
}

pPdf.ByMomentPdf <- function(Pdf, q, ParamSet, ...) {
  FName <- paste0("p", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(q = q), c(Param, list(...))) )
}

#' qPdf
#' @describeIn dPdf Gives the quantile function
qPdf <- function(Pdf, p, ParamSet, ...) {
  UseMethod("qPdf")
}

qPdf.ByMomentPdf <- function(Pdf, p, ParamSet, ...) {
  FName <- paste0("q", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(p = p), c(Param, list(...))) )
}

#' rPdf
#' @describeIn dPdf Generates random deviates
rPdf <- function(Pdf, p, ParamSet, ...) {
  UseMethod("rPdf")
}

rPdf.ByMomentPdf <- function(Pdf, n, ParamSet, ...) {
  FName <- paste0("r", Pdf$Function)
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  
  do.call( FName, append(list(n = n), c(Param, list(...))) )
}



SolutionMoments <- function(Pdf, ...) {
  UseMethod("SolutionMoments")
}




# Launch Space ----

#' GetLaunchSpace
#' Generate a list of launch coordinates as starting points in the
#' search for solutions.
#' @param Pdf An object of class `ByMomentPdf`.
#' @param Count Number of desired launch points.
#' @param Method A method to generate launch points.
#'
#' @return A class after adding the optimisation solution. 
#' The classes of `Pdf` are preserved.
#' @export
#'
#' @author Jan Seifert
#' @examples
GetLaunchSpace <- function( Pdf, ... ) {
  if (is.null(Pdf$ParamSpace))
    stop("No parameter space specified") # 

  UseMethod( "GetLaunchSpace" )
}


#' GetLaunchSpace.ByMomentPdf
#' @describeIn GetLaunchSpace
GetLaunchSpace.ByMomentPdf <- function( Pdf, Count,
                                        Method = c("Random", "Cubic", 
                                                   "Harmonic", "Manual") ) {
  # PRECONDITIONS
  Method <- match.arg(Method)
  
  # Dimensions of parameter space
  NDim <- ncol(Pdf$ParamSpace)
  
  # RUN
  # Random: 1 to n randomly set points
  if (Method == "Random") {
    LaunchPad <- GetRandomPoints(NDim, Count)
    # Scale points to 'ParamSpace'
    Range <- Pdf$ParamSpace["to", ] - Pdf$ParamSpace["from", ]
    LaunchPad <- LaunchPad * Range
    LaunchPad <- LaunchPad + Pdf$ParamSpace["from", ]
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
    LaunchPad <- GetHarmonicPoints(NDim, Count, 
                                   Seed = 0.5 + runif(1, -0.1, 0.1))
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


#' SetLaunchSpace
#' Set a launch space instead of asking the object to generate
#' its own.
#' @param LaunchSpace A list of launch coordinates to set.
#' @return `ByMomentPdf`-object
#' @export
#' @author Jan Seifert
#' @examples
SetLaunchSpace <- function( LaunchSpace ) {
  UseMethod("SetLaunchSpace")
}



# Solutions ----

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
#' @author Jan Seifert
#' @examples
AddSolution <- function( Pdf, LaunchPoint, SoluParam, Add = FALSE, ... ) {
  # PRECONDITIONS
  if(LaunchPoint < 1L) stop("Invalid launch point (must be > 0).")
  if(LaunchPoint != 1 && LaunchPoint > nrow(Pdf$LaunchSpace))
    stop("Index of launch point out of range.")
  
  # RESULT
  UseMethod("AddSolution")
}


#' AddSolution.ByMomentPdf
#' @describeIn AddSolution
AddSolution.ByMomentPdf <-  function( Pdf, LaunchPoint, 
                                      SoluParam, Append = FALSE ) {
  # RESULT
  if (isFALSE(Append) ||
      is.null(Pdf$ParamSolved) || 
      length(Pdf$ParamSolved) < 1) {
    Pdf$ParamSolved <- list(list(LaunchPoint, SoluParam))
  } else {
    # Check if solution for this launch point already exists
    LPs <- unlist( lapply(Pdf$ParamSolved, `[[`, 1) )
    Pos <- which(LPs == LaunchPoint)
    
    if (length(Pos) != 0) {
      # Replace if it does
      Pdf$ParamSolved[[Pos]] <- list(LaunchPoint, SoluParam)
    } else {
      # Find position and insert it at the right position
      Len <- length(Pdf$ParamSolved)
      if (LaunchPoint < LPs[1]) {
        Pdf$ParamSolved <- c(list(list(LaunchPoint, SoluParam)), Pdf$ParamSolved)
      } else if (LaunchPoint > LPs[Len]) {
        Pdf$ParamSolved <- c(Pdf$ParamSolved, list(list(LaunchPoint, SoluParam)))
      } else {
        Pos <- which.min(LaunchPoint > LPs)
        Pdf$ParamSolved <- c(Pdf$ParamSolved[1L:(Pos-1L)], 
                             list(list(LaunchPoint, SoluParam)), 
                             Pdf$ParamSolved[Pos:Len])
      }
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
#' @author Jan Seifert
#' @examples
FindPdf <- function( Pdf, LaunchPoint, Append = FALSE, ... ) {
  # PRECONDITIONS
  if(LaunchPoint < 1L) stop("Invalid launch point (must be > 0).")
  if(LaunchPoint > nrow(Pdf$LaunchSpace))
    stop("Index of launch point out of range.")
  
  # RESULT
  UseMethod("FindPdf")
}



#' EvaluatePdf
#' Assess the quality of solutions stored in `Pdf` by assessing 
#' the distance to the target.
#' 
#' Generic function for all sub-classes of `[ByMomentPdf]`.
#' @param Pdf An object of class `[ByMomentPdf]` or a sub-class.
#' @param UsePdf If `TRUE` the distance between found and desired solutions
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
#' @return A class after adding the optimisation solution as element
#' `DistaMo` or `DistaFu`.
#' The classes of `Pdf` are preserved.
#' @references Deza, E., & Deza, M. M. (2009). Encyclopedia of Distances. 
#' Springer. https://doi.org/10.1007/978-3-642-00234-2
#' @export
#'
#' @author Jan Seifert
#' @examples
EvaluatePdf <- function(Pdf, ...) {
  UseMethod("EvaluatePdf")
}

#' EvaluatePdf.ByMomentPdf
#' @describeIn EvaluatePdf
EvaluatePdf.ByMomentPdf <- function(Pdf, UsePdf = FALSE) {
  if(is.null(Pdf$ParamSolved)) 
    stop("No solutions available in PDF")
  if(UsePdf && is.null(Pdf$TarFu)) 
    stop("'UsePdf' demands a target function in 'Pdf'")

  # Compute the distances
  if (UsePdf) {
    if(require(philentropy, quietly = TRUE) == FALSE)
      stop("Package 'philentropy' is not available")
    # 
    x <- .plotrange.ByMomentpdf(Pdf)
    #
    Solutions <- lapply(Pdf$ParamSolved, `[[`, 2)
    DResult <- list()
    Pos <- 1
    for(s in Solutions) {
      # Check if parameters are valid
      if (any(unlist(s) < Pdf$ParamSpace["from",]) || 
          any(unlist(s) > Pdf$ParamSpace["to",])) {
        DResult <- c(DResult, list(list(ID = Pdf$ParamSolved[[Pos]][[1]], 
                                   Delta = NA)))
      } else {
        # Compute distance
        SoluResult  <- dPdf(Pdf, x, as.list(s)) # get PDF(x)
        TarFuResult <- do.call( Pdf$TarFu[[1]], c(list(x = x), Pdf$TarFu[[2]]) )
        D <- jeffreys(TarFuResult, SoluResult, TRUE, "log2")
        DResult <- c(DResult, list(list(ID = Pdf$ParamSolved[[Pos]][[1]], 
                                   Delta = D)))
      }
      Pos <- Pos+1
    }
    Pdf$DistaFu <- DResult
  } else { # Use moments
    # if moments are not available, generate
    if(is.null(Pdf$Moments)) {
      Pdf$Moments <- list()
      for(i in 1:length(Pdf$ParamSolved)) {
        Pdf$Moments <- c(Pdf$Moments, list(list(i, SolutionMoments(Pdf, i))))
      }
    }
    
    Moments <- lapply(Pdf$Moments, `[[`, 2)
    DResult <- list()
    Method <- 1 # i.e. "euclidean", see "?dist"
    Attrs <- list(Size = 2, Labels = "", Diag = FALSE, Upper = FALSE,
                  method = "euclidean", call = match.call(), class = "dist")
    ID <- 1
    for(m in Moments) {
      X <- rbind(Pdf$TarMo, m)
      D <- .Call(C_Cdist, X, Method, Attrs)
      DResult <- c(DResult, list(ID = Pdf$Moments[ID][[1]], Delta = D))
      ID <- ID+1
    }
    Pdf$DistaMo <- DResult
  }
  return(Pdf)
}




#' BestSolution
#' Determines which solution stored in a `ByMomentPdf`-object is
#' closest to the target solution.
#' @param Pdf A `ByMomentPdf`-object
#' @param UsePdf If `TRUE` the distance between found and desired solutions
#' is determined using the Pdf. If `FALSE` only the moment vectors will be 
#' used. Default is `FALSE`.
#'
#' @return A single vector of solution parameters or a matrix 
#' if more than solution reaches the same optimal result.
#' @export
#' @author Jan Seifert
#' @examples
BestSolution <- function(Pdf, ...) {
  UseMethod("BestSolution")
}

#' BestSolution.ByMomentpdf
#' @describeIn BestSolution
BestSolution.ByMomentPdf <- function(Pdf, UsePdf = FALSE) {
  if(UsePdf) {
    if(is.null(Pdf$DistaFu)) Pdf <- EvaluatePdf(Pdf, UsePdf)
    DistaXx <- Pdf$DistaFu
  } else {
    if(is.null(Pdf$DistaMo)) Pdf <- EvaluatePdf(Pdf, UsePdf)
    DistaXx <- Pdf$DistaMo
  }
  
  # Get the index of the solutions with best results
  Distances <- unlist(lapply(DistaXx, `[`, 2))
  Best <- which(Distances == min(Distances, na.rm = TRUE))
  
  # Get all the parameter sets yielding best results
  BestParamSets <- lapply(Pdf$ParamSolved, `[[`, 2)
  BestParamSets <- BestParamSets[Best]
  UniBestParamSets <- unique(BestParamSets)
  NBest <- length(UniBestParamSets)

  if (length(UniBestParamSets) > 0) {
    UniBestParamSets <- matrix(unlist(UniBestParamSets), 
                               nrow = NBest, byrow = TRUE)
    colnames(UniBestParamSets) <- paste0("Param", 1:ncol(UniBestParamSets))
  }
  
  return(UniBestParamSets)
}




# Output ----

#' summary.ByMomentPdf
#' Result summaries of the results of optimisation function to identify 
#' probability distributions by moments.
#' @param object an object for which a summary is desired.
#' @param ... additional arguments affecting the summary produced.
#' @return Nothing
#' @export
#' @author Jan Seifert
#' @examples
summary.ByMomentPdf <- function( object, ...) {
  title  <- function(t) cat(paste0("\n", t, ":\n"))
  output <- function(x, y = "") {
    cat(paste0("   ", x, "  ", paste(y, collapse = ", "), "\n"))
  }
  
  title("Call")
  o <- rbind(Method = object$Function,
             TargetMoments = paste0(object$TarMo, collapse = ", "),
             StartingPoints = nrow(object$LaunchSpace))
  colnames(o) <- ""
  print(o, quote = FALSE)
  
  title("Results")
  output("Solutions", length(object$ParamSolved))
  
  title("Best solution")
  if (!is.null(object$DistaFu)) {
    Best <- BestSolution(object, UsePdf = TRUE)
    DistanceMethod <- "Jeffreys Distance"
    Distances <- unlist(lapply(object$DistaFu, `[`, 2))
    Distance <- min(Distances, na.rm = TRUE)
  }
  else if (!is.null(object$DistaFu)) {
    Best <- BestSolution(object, UsePdf = FALSE)
    DistanceMethod <- "Euclidean Distance"
    Distances <- unlist(lapply(object$DistaFu, `[`, 2))
    Distance <- min(Distances, na.rm = TRUE)
  }
  else
    output("No solutions available")
  
  if (!is.null(nrow(Best))) {
    cat("  ", nrow(Best), "solutions with best characteristics\n")
    output(DistanceMethod, round(Distance, digits = 4))
    output("")
    print(round(Best, digits = 4), rownames = FALSE)
  }
}



#' plot.ByMomentPdf
#' Plot a distributions function for a selected solution.
#' @param Pdf An object of class `ByMomentPdf`.
#' @param ParamSet Either a single integer as index for the selected 
#' `Pdf$ParamSolved` or a list with the parameters to be passed on to the
#' distribution function.
#' @param X The coordinates of points in the plot.
#' @param Type The type of distribution function is either "d" für probability 
#' density function, "p" probability distribution function, or "q" for 
#' the quantile function.
#' @param AddTarFu If `TRUE`the target function will be added to the plot 
#' (works only for dxxx-functions at the moment).
#' @param ... Additional arguments passed on to plot or the distribution 
#' function.
#' @return -
#' @export
#' @author Jan Seifert
#' @examples
plot.ByMomentPdf <- function( Pdf, ParamSet, X = NULL, Type = c(d, p, q), 
                              AddTarFu = FALSE, ...) {
  # PRECONDITIONS
  Type <- ifelse(missing(Type), "d", match.arg(Type))
  if(is.list(ParamSet)) Param <- ParamSet
  else Param <- Pdf$ParamSolved[ParamSet, ]
  X <- .plotrange.ByMomentpdf(Pdf, X)
  
  # RESULT
  if(Type == "d") {
    Y  <- dPdf(Pdf, X, Param) 
    XLab <- "x"
    YLab <- "Density"
  } else if (Type == "p") {
    Y  <- pPdf(Pdf, X, Param)
    XLab <- "q"
    YLab <- ""
  } else if (Type == "q") {
    Y  <- qPdf(Pdf, X, Param) 
    XLab <- "p"
    YLab <- ""
  }
  Main <- paste(Pdf$Function, format(Param, digits = 2), collapse = ", ")
  plot(X, Y, type = "l", main = Main, xlab = XLab, ylab = YLab, ...)
  
  if(AddTarFu) {
    #TODO: Type not supported here
    Y <- do.call( Pdf$TarFu[[1]], c(list(x = X), Pdf$TarFu[[2]]) )
    lines(X, Y, col = "cyan")
  }
}



#' hist.ByMomentPdf
#' A histogram of distances between the solutions and the target.
#' @param Pdf An object of class `ByMomentPdf`.
#' @param UsePdf  If `TRUE` the distance between found and desired solutions
#' is determined using the Pdf. If `FALSE` only the moment vectors will be 
#' used. Default is `FALSE`.
#' @return A histogram that shows the distribution of distance of all
#' solutions.
#' @export
#' @author Jan Seifert
#' @examples
hist.ByMomentPdf <- function(Pdf, UsePdf = FALSE) {
  if(UsePdf && is.null(Pdf$DistaFu)) 
    stop("Pdf has not been evaluated.")
  if(!UsePdf && is.null(Pdf$DistaMo)) 
    stop("Moments have not been evaluated.")
  
  if(UsePdf) {
    ID <- unlist(lapply(Pdf$DistaFu, `[`, "ID"))
    Distance <- unlist(lapply(Pdf$DistaFu, `[`, "Delta"))
    Data <- data.frame(ID, Distance)
  } else {
    Data <- data.frame(Distance = Pdf$DistaFu)
  }
  
  p <- ggplot(data = Data, aes(x = Distance, na.rm = TRUE)) + 
         geom_histogram()
  print(p)
  invisible(p)
}