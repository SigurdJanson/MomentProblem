source("FindPdf_Root.R")
source("FindPdf_GlD.R")
source("FindPdf_Poly.R")



# Question: are we talking about centralised or raw moments??? 
# (https://de.wikipedia.org/wiki/Moment_(Stochastik))




#' GetPdf
#' @param Pdf An initialised `ByMomentPdf` object with parameter
#' space and launch points.
#' @return A converged Pdf object of class `ByMomentPdf`.
#' @export
#'
#' @examples
GetPdf <- function( Pdf ) {
  #TarMo, Method = c("gld", "pearson", "poly"),
  #                  Launch = c("Random", "Cubic", "Harmonic", "Manual")) {
  
  Method <- match.arg(Method, several.ok = TRUE)
  Launch <- match.arg(Launch, several.ok = TRUE)
  
  Result <- list()
  
  # Converge for each launch point
  for(lp in 1:nrow(pdf$SpaceLaunch)) {
    #TODO: save result <- FindPdf( Pdf, lp )
  }
  
}