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
#' @author Jan Seifert
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
  Pdf <- GetLaunchSpace(Pdf, 1E4, "Harmonic")

  # Converge for each launch point
  for(LP in 1:nrow(Pdf$LaunchSpace)) {
    Pdf <- FindPdf( Pdf, LP, Append = TRUE )
    if(LP %% 1000 == 0) cat("o")
    else if(LP %% 100 == 0) cat(".")
  }

  UsePdf <- TRUE
  Pdf  <- EvaluatePdf(Pdf, UsePdf = UsePdf)
  Best <- BestSolution(Pdf, UsePdf = UsePdf)
  
  summary(Pdf)
  Champion <- as.list(Best[1,])
  names(Champion) <- paste0("lambda", 1:4)
  plot(Pdf, Champion, 5, AddTarFu = TRUE)
  hist(Pdf, UsePdf = TRUE)
  
  return(Pdf)
}

Pdf <- GetPdf(c(0, 1, 0, 3), "gld")

