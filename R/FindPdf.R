library(ggplot2)
#source("FindPdf_Root.R")
source("./R/FindPdf_GlD.R")
source("./R/FindPdf_Poly.R")



# Question: are we talking about centralised or raw moments??? 
# (https://de.wikipedia.org/wiki/Moment_(Stochastik))

#' GetPdf
#' @param Pdf An initialised `ByMomentPdf` object with parameter
#' space and launch points.
#' @return A converged Pdf object of class `ByMomentPdf`.
#' @export
#'
#' @examples
GetPdf <- function( TarMo, Template = c("gld", "pearson", "poly"), 
                    Plot = FALSE ) {
  #Launch = c("Random", "Cubic", "Harmonic", "Manual")) {
  #Launch <- match.arg(Launch, several.ok = TRUE)
  Template <- match.arg(Template, several.ok = TRUE)
  
  # Initialise
  TargetFunction <- list("dnorm", list(mean = 0, sd = 1))
  Pdf <- New_ByMomentPdf.gld(TarMo, TarFu = TargetFunction)
  if (is.null(Pdf$TarFu)) stop("is.null(Pdf$TarFu) - 0")
  Pdf <- GetLaunchSpace(Pdf, 4, "Harmonic")
  if (is.null(Pdf$TarFu)) stop("is.null(Pdf$TarFu) - 1")
  
  # Converge for each launch point
  for(LP in 1:nrow(Pdf$LaunchSpace)) {
    Pdf <- FindPdf( Pdf, LP, Append = TRUE )
    if (is.null(Pdf$TarFu)) stop("is.null(Pdf$TarFu) - 2")
  }

  UsePdf <- TRUE
  Pdf  <- EvaluatePdf(Pdf, UsePdf = UsePdf)
  Best <- BestSolution(Pdf, Usepdf = UsePdf)
  
  # probably provide histogram of distances across solutions
  #if(Plot) {
  if(!is.null(Pdf$DistaFu)) {
    Data <- data.frame(Distance = Pdf$DistaFu)
    ggplot(data = Data, aes(x = Distance)) + 
      geom_histogram()
  }
  if(!is.null(Pdf$DistaMo)) {
    Data <- data.frame(Distance = Pdf$DistaFu)
    ggplot(data = Data, aes(x = Distance)) + 
      geom_histogram()
  }
  #}#Plot
  
  # Present best solution(s)
  #if(Plot) {
  #TODO: plot target function and Pdf
  #}#Plot
  return(Pdf)
}

Pdf <- GetPdf(c(0, 1, 0, 3), "gld")