# Helper functions for plots, printing and summaries

#' .plotrange.ByMomentpdf
#' Generates a sequence to plot a probability distribution.
#' @param Pdf A `ByMomentPdf` object
#' @param X A range around the expected mean. If `X` is a scalar the range will be set as 
#' `c(-X, X)` around the mean. If `X` is invalid `DefaultRange` will be used.
#' @param DefaultRange Use this range in case `X` is invalid.
#' @param DefaultResolution Desired length of the sequence.
#'
#' @return  Numeric vector.
#' @export
#'
#' @examples
.plotrange.ByMomentpdf <- function(Pdf, X = NULL, 
                                   DefaultRange = 5, DefaultResolution = 1E2) {
  # PRECONDITIONS
  if(is.null(X) || length(X) == 0)
    X <- DefaultRange
  
  # RESULTS
  if (length(X) == 1) {
    X <- c(-X, X)
    XRange <- Pdf$TarMo[1] + X * Pdf$TarMo[2]
    XDensity <- diff(X) * DefaultResolution
    X <- seq(XRange[1], XRange[2], length.out = XDensity)
  } else if (length(X) == 2) {
    XRange <- sort(X)
    XDensity <- diff(XRange) * DefaultResolution
    X <- seq(XRange[1], XRange[2], length.out = XDensity)
  }
  return(X)
}




#' logseq
#' Generate a sequence with logarithmically increasing intervals
#' @param from,to the starting and (maximal) end values of the sequence. 
#' Of length 1 unless just from is supplied as an unnamed argument.
#' @param by Increment of the sequence (numeric).
#' @return Numeric vector
#' @export
#' @examples
logseq <- function(from = 1, to = 1, by = ((log10(to) - log10(from))/(length.out - 1)),
                   length.out = NULL, along.with = NULL, ...) {
  if(from <= 0 || to <= 0) stop("Values <= 0 not allowed in logarithmic scale.")
  if(!is.null(along.with)) by <- ((log10(to) - log10(from))/(length(along.with) - 1))
  
  exp( log(10) * seq(log10(from), log10(to), by = by, ...) )
}

