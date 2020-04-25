# Helper functions for plots, printing and summaries

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


logseq <- function(from = 1, to = 1, by = 1) {
  exp(log(10)*seq(log10(from), log10(to), by = by))
}

