library(gld)
source("./R/FindPdf_Root.R")
source("./R/nmkp.R")



#' .v1f
#' helper function to compute component v1
#' @param L3 GLD parameter lambda3
#' @param L4 GLD parameter lambda4
.v1f <- function( L3, L4 ) {
  s1 <- 1 / L3 / (L3+1)
  s2 <- 1 / L4 / (L4+1)
  return(s1-s2)
}


#' .v2f
#' helper function to compute component v2
#' @param L3 GLD parameter lambda3
#' @param L4 GLD parameter lambda4
.v2f <- function( L3, L4 ) {
  s1 <- 1 / L3^2 / (2*L3+1)
  s2 <- 1 / L4^2 / (2*L4+1)
  s3 <- 2 / L3 / L4 * suppressWarnings(beta(L3+1, L4+1))
  return(s1+s2-s3)
}


#' .v3f
#' helper function to compute component v3
#' @param L3 GLD parameter lambda3
#' @param L4 GLD parameter lambda4
.v3f <- function( L3, L4 ) {
  s1 <- 1 / L3^3 / (3*L3+1)
  s2 <- 1 / L4^3 / (3*L4+1)
  s3 <- 3 / L3^2 / L4 * suppressWarnings(beta(2*L3+1, L4+1))
  s4 <- 3 / L3 / L4^2 * suppressWarnings(beta(L3+1, 2*L4+1))
  return(s1-s2-s3+s4)
}


#' .v4f
#' helper function to compute component v4
#' @param L3 GLD parameter lambda3
#' @param L4 GLD parameter lambda4
.v4f <-  function( L3, L4 ) {
  s1 <- 1 / L3^4 / (4*L3+1)
  s2 <- 1 / L4^4 / (4*L4+1)
  s3 <- 6 / L3^2 / L4^2 * suppressWarnings(beta(2*L3+1, 2*L4+1))
  s4 <- 4 / L3^3 / L4 * suppressWarnings(beta(3*L3+1, L4+1))
  s5 <- 4 / L3 / L4^3 * suppressWarnings(beta(L3+1, 3*L4+1))
  return(s1+s2+s3-s4-s5)
}



#' .Alpha1
#' Helper function
.Alpha1 <- function( Lambda ) { 
  
  alpha1 <- Lambda[1] - 1/Lambda[2] * ( (1/(Lambda[3]+1)) - (1/(Lambda[4]+1)) )
  return(alpha1)
}

#' .Alpha2
#' Helper function
.Alpha2 <- function( v1, v2, Lambda ) { 
  # based on eqn 13: lambda2 <- sqrt(v2 - v1^2) / sqrt(Moments[2]) 
  alpha2 <- (v2 - v1^2) / Lambda[2]^2
  return(alpha2)
}


#' .Alpha3
#' helper function to compute component third moment
#' @param L vector containing GLD parameters 'lambda3' and 'lambda4'
#' @param A3 Expected value for 'alpha3' (3. moment)
#' @param my Expected mean of GLD
#' @note Solves equation 24 in Au-Yeung (2003)
.Alpha3 <- function( v ) { 
  
  alpha3 <- (v[3] - 3*v[1]*v[2] + 2*v[1]^3) / (v[2]-v[1]^2)^1.5
  
  return(alpha3)
}

#' .Alpha4
#' helper function to compute component third moment
#' @param L vector containing GLD parameters 'lambda3' and 'lambda4'
#' @param A4 Expected value for 'alpha4' (4. moment)
#' @param my Expected mean of GLD
#' @note Solves equation 25 in Au-Yeung (2003)
.Alpha4 <- function( v ) { 
  
  n <- (v[4] - 4*v[1]*v[3] + 6*v[1]^2*v[2] - 3*v[1]^4)
  d <- (v[2]-v[1]^2)^2
  alpha4 <- n / d
  
  return(alpha4)
}


# OBJECTIVE FUNCTION ----


#' .DeltaAllGLD
#' Helper function to compute the distance from the expected values
#' of all moments of the desired probability distribution.
#' The output of this function must be optimised.
#' @param L The four lambda parameters of the generalised lambda distribution
#' @param A The first four expected moments
#' @param MinMax Direction of distance, 0 (default) or 1.
#' @details This function can be used in classic optimisation algorithms 
#' like the 'optim' function. In that case the algorithm minimises the result and
#' 'MinMax' must be 0. Genetic algorithms maximise and 'MinMax' = 1 is correct.
#' @return The distance of the observed moments from the expected.
.DeltaAllGLD <- function(L, A, MinMax = 0) {
  L <- c(NA, NA, L)
  #if(A[3] != 0) L[4] <- L[4] 
  #else L[4] <- L[3] # When symmetry: L3 = L4
  v1 <- .v1f(L[3], L[4])
  v2 <- .v2f(L[3], L[4])
  v3 <- .v3f(L[3], L[4])
  v4 <- .v4f(L[3], L[4])

  L[2] <- sqrt(v2 - v1^2) / sqrt(A[2]) # eqn. 13
  L[1] <- A[1] + (1/(L[3]+1) - 1/(L[4]+1)) / L[2] # eqn 14
  
  Alpha2 <- .Alpha2(v1, v2, L) #(v2 - v1^2) / L[2]^2
  Alpha1 <- .Alpha1(L) # L[1] - 1/L[2] * ( (1/(L[3]+1)) - (1/(L[4]+1)) )

  Delta1 <- abs(A[1] - Alpha1)
  Delta2 <- abs(A[2] - Alpha2)
  Delta3 <- abs(A[3] - .Alpha3(c(v1, v2, v3, v4)))
  if(A[3] != 0) 
    Delta4 <- abs(A[4] - .Alpha4(c(v1, v2, v3, v4)))
  else
    Delta4 <- Delta3 # When symmetry: L3 = L4
  
  DeltaAll <- sqrt(Delta1^2 + Delta2^2 + Delta3^2 + Delta4^2)
  if(is.na(DeltaAll) || is.nan(DeltaAll) || is.infinite(DeltaAll))
    DeltaAll <- 999999
  if(MinMax == 1)
    DeltaAll <- 1 / DeltaAll
  return(DeltaAll)
}


# OPTIMISATION ----

#' .GLDTypeConstraints
#' The ideas is to select better init values for optimisation if 
#' the desired shape is known.
#' 
#' TODO: This is a STUB
.GLDTypeConstraints <- function( Shape = c("unimodal", "monotone", "u", "s"),
                                 Infinite = c("none", "lower", "upper", "both")) {
  L3 <- L4 <- NA # Init lambda3 and 4
  
  Shape <- match.arg(Shape)
  # Class I (lambda3 < 1, lambda4 < 1): Unimodal densities with continuous tails.
  if(Shape == "unimodal") {
    # Class Ia (lambda3, lambda4 ≤ ½)
    # Infinite slope in both directions
    
    # Class Ib (½ < lambda3 < 1, lambda4 ≤½), and
    
    # Class Ic (½ < lambda3 <1, ½ < lambda4 < 1).
    
    # Class V (lambda3 > 2, lambda4 > 2): 
    # Unimodal pdfs with both tails truncated.
    if(Infinite == "none") {
      L3 <- L4 <- 5
    }
  }
  if(Shape == "u") {
    # Class III (1 < lambda3 < 2 , 1< lambda4 < 2)
    # U-shaped densities with both tails truncated.
    L3 <- L4 <- 1 + ((2-1)/2)
    # symmetry is ignored but possible
    # is always truncated
  }
  
  if(Shape == "monotone") {
    # Class II (lambda3 > 1, lambda4 < 1)
    # Monotone pdfs similar to those of the exponential or chi²
    # distributions. Left tail is truncated.
    L3 <- 2
    L4 <- 1E-5
  }

  if(Shape == "s") {
    # Class IV (lambda3 >2, 1< lambda4 < 2)
    # Rarely occurring S-shaped pdfs with one mode and one
    # antimode. Both tails are truncated.
    L3 <- 3
    L4 <- 1
  }
}



#' .FindPDF_GLD_A3A4
#' 
#' @param Moments Vector holding the first four moments.
#' @param Tolerance numeric ≥ 0. Differences smaller than tolerance 
#' are not reported. The default value is close to 1.5e-8.
#' 
#' @details Moments are µ (mean), $\sigma$ (standard deviation), 
#' skewness and kurtosis.
#'
#' @return
#' @export
#'
#' @examples
.FindPDF_GLD_A3A4 <- function( Moments, Tolerance = sqrt(.Machine$double.eps) ) {
  Lambda <- numeric(4)
  names(Moments) <- c("Mean", "Var", "Skew", "Kurt")

  InitLambda <- runif(4, -0.25, 5) # initial values to start optimisation
  
  # Find Lambda3 and 4
  if(Moments["Skew"] == 999999) { # Symmetry: lambda3 == lambda4
    # Since (lambda3 == lambda4) we need only 1 dimension to optimise
    Tolerance <- 1E-4
    r <- optim(Start, method = c("Brent"), lower = -0.25, upper = 50, 
               fn = .DeltaA3A4GLD, A3 = Moments["Skew"], A4 = Moments["Kurt"], 
               control = list(reltol = Tolerance, maxit = 50000) )
    Lambda[3] <- r$par[1]
    Lambda[4] <- r$par[1]
  } else {
    # GA <- ga(type = "real-valued", fitness = function(x, ...) .DeltaA3A4GLD(x, ...),
    #          A3 = Moments["Skew"], A4 = Moments["Kurt"], MinMax = 1, 
    #          lower = c(-0.25, -0.25), upper = c(50, 50),
    #          popSize = 100, maxiter = 500, run = 100, pmutation = 0.15)
    #print(summary(GA))
    #plot(GA)
    #GA@solution[1,] - c(-0.15, -0.15)
    r <- optimx(GA@solution[1,], method = c("Nelder-Mead"), 
                fn = .DeltaA3A4GLD, A3 = Moments["Skew"], A4 = Moments["Kurt"], 
               control = list(reltol = Tolerance*100, maxit = 50000) )
    if (r$convcode == 0) {
      Lambda[3] <- r$x1 #r$par[1]
      Lambda[4] <- r$x2 #r$par[2]
      v1 <- .v1f(Lambda[3], Lambda[4])
      v2 <- .v2f(Lambda[3], Lambda[4])
      Lambda[2] <- sqrt(v2 - v1^2) / sqrt(Moments[2]) 
      Lambda[1] <- Moments[1] + 1/Lambda[2] * ( (1/(Lambda[3]+1)) + (1/(Lambda[4]+1)) )
    } else {
      Lambda <- rep(NA_integer_, 4)
    }
  }

  return(Lambda)
}


#' .FindPDF_GLD
#'
#' @param TarMo Target moments, vector with the first four moments.
#' @param InitLambda Starting point
#' @param Lower,Upper Lower and upper bounds on the parameters. A 
#' vector of the same length as the parameters.
#' @param Rigour How to enhance the 
#' @param Tolerance numeric ≥ 0. Differences smaller than tolerance 
#' are not reported. The default value is close to 1.5e-8.
#' 
#' @details Moments are µ (mean), $\sigma$ (standard deviation), 
#' skewness and kurtosis.
#' If `Rigour` is `"Genetic"` then the search starts with a genetic 
#' algorithm that is known to provide better coverage. After that a 
#' Nelder-Mead optimisation completes the optimisation using the result
#' of the genetic algorithm as starting point.
#' Other values of `Rigour` use a Nelder-Mead optimisation.
#' @return
#' @export
#'
#' @author Jan Seifert
#' @examples
.FindPDF_GLD <- function( TarMo, InitLambda, Lower, Upper,
                          Tolerance = 1E-6, Rigour = "None") { #sqrt(.Machine$double.eps) ) {
  Lambda <- numeric(4)
  #names(TarMo) <- c("Mean", "Var", "Skew", "Kurt")
  
  if (Rigour == "Genetics") {
    require(GA, quietly = TRUE)
    
    NelderMeadArgs <- list(method = "Nelder-Mead", 
                           poptim = 0.05,
                           pressel = 0.5,
                           control = list(reltol = Tolerance, maxit = 100))
    
    # Genetic algorithm first: locate the general vicinity of the minimum
    GA <- ga(type = "real-valued",
             fitness = function(x, ...) .DeltaAllGLD(x, ...), A = TarMo, MinMax = 1,
             suggestions = rbind(InitLambda),
             lower = Lower, upper = Upper,
             popSize = 100, maxiter = 1000, run = 100, pmutation = 0.15,
             monitor = NULL)
    ResetLambda <- GA@solution[1,]
  } else {
    ResetLambda <- InitLambda
  }

  # Run Nelder-Mead 
  if (length(ResetLambda) == length(InitLambda) ||
      all(ResetLambda >= Lower) || all(ResetLambda <= Upper)) {

    r <- nmkp(ResetLambda, fn = .DeltaAllGLD, A = TarMo, 
              lower = Lower, upper = Upper, 
              control = list(tol = 9e6, maxfeval = 5E5, trace = FALSE))
    
    if (r$convergence == 0) {
      Lambda[3] <- r$par[1]
      Lambda[4] <- r$par[2]
      Lambda[2] <- sqrt(.v2f(Lambda[3], Lambda[4]) - .v1f(Lambda[3], Lambda[4])^2) /
        sqrt(TarMo[2]) 
      Lambda[1] <- TarMo[1] + (1/(Lambda[3]+1) - 1/(Lambda[4]+1)) / Lambda[2] 
    } else {
      Lambda <- rep(NA_integer_, 4)
      r$value <- NA
    } 
  } else {
    Lambda <- rep(NA_integer_, 4)
    r$value <- NA
  }

  return(list(Lambda = Lambda, Distance = r$value))
}


# CLASS ----

#' New_ByMomentPdf.gld
#' Constructor of class `gld`.
#' @param TarMo Target moments, vector with the first four moments.
#' @param TarFu Target function is a list with the name of the 
#' desired density function as first entry and it's parameters as 
#' named list as the second.
#' @return An initialised object with S3 class `ByMomentPdf` and `gld`.
#' @export
#'
#' @author Jan Seifert
#' @examples
New_ByMomentPdf.gld <- function( TarMo, TarFu = NULL ) {
  if(length(TarMo) > 4) 
    stop("Generalised lambda approximation can only handle up to 4 moments.")
  this <- New_ByMomentPdf( TarMo, TarFu )
  
  this$Function   <- "gl"
  # Dimensions of the parameter space of the PDF
  # -+500 is arbitrary; -0.25 is the defined limit of the range of def.
  Dim34 <- c(-0.25, 500) # 
  this$ParamSpace <- matrix(data = c(Dim34, Dim34), nrow = 2, 
                            dimnames = list(c("from", "to"), NULL))
  
  class(this) <- append(class(this), "gld")
  return(this)
}



SolutionMoments.gld <- function(Pdf, Solution = 1) {
  if(length(Solution) == 1) {
    LPs <- lapply(Pdf$ParamSolved, `[[`, 1)
    Pos <- which(unlist(LPs) == Solution)
    if (length(Pos) == 1) {
      Solution <- Pdf$ParamSolved[[Pos]][[2]]
    } else {
      stop(paste("Solution", Solution, "not found."))
    }
  } else stop("Works only for a single solution.")
    
  Alpha <- numeric(4)
  
  v1 <- .v1f(Solution[3], Solution[4])
  v2 <- .v2f(Solution[3], Solution[4])
  v3 <- .v3f(Solution[3], Solution[4])
  v4 <- .v4f(Solution[3], Solution[4])
  
  Alpha[2] <- .Alpha2(v1, v2, Solution)
  Alpha[1] <- .Alpha1(Solution)
  Alpha[3] <- .Alpha3(c(v1, v2, v3, v4))
  Alpha[4] <- .Alpha4(c(v1, v2, v3, v4))

  return(Alpha)
}



FindPdf.gld <- function( Pdf, LaunchPoint, Append = FALSE ) {
  # PRECONDITIONS get tested by generic method 'FindPdf'
  # RESULTS
  LP <- Pdf$LaunchSpace[LaunchPoint, ]

  # Get a resulting PDF
  Result <- .FindPDF_GLD( Pdf$TarMo, LP, 
                          Lower = Pdf$ParamSpace["from",], 
                          Upper = Pdf$ParamSpace["to",], 
                          Pdf$Tolerance )
  
  # Determine distance to desired solution
  Pdf <- AddSolution( Pdf, LaunchPoint, SoluParam = Result$Lambda, 
                      DistaMo = Result$Distance, Append = Append )
  return(Pdf)
}

