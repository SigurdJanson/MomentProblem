library(gld)
library(GA)
library(optimx)
#library(plot3D)


#' .v1f
#' helper function to compute component v1
#' @param L3 GLD parameter lambda3
#' @param L4 GLD parameter lambda4
.v1f <- function( L3, L4 ) {
  if(is.numeric(c(L3, L4)))
    if(L3 == 0 && L4 == 0) 
      return(0) # TODO: verify this line
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
  if(L3 == 0 && L4 == 0) return(0) # TODO: verify this line
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


#' .Alpha3
#' helper function to compute component third moment
#' @param L vector containing GLD parameters 'lambda3' and 'lambda4'
#' @param A3 Expected value for 'alpha3' (3. moment)
#' @param my Expected mean of GLD
#' @note Solves equation 24 in Au-Yeung (2003)
.Alpha3 <- function( v, My = NA ) { 
  
  alpha3 <- (v[3] - 3*v[1]*v[2] + 2*v[1]^3) / (v[2]-v[1]^2)^1.5
  
  return(alpha3)
}

#' .Alpha4
#' helper function to compute component third moment
#' @param L vector containing GLD parameters 'lambda3' and 'lambda4'
#' @param A4 Expected value for 'alpha4' (4. moment)
#' @param my Expected mean of GLD
#' @note Solves equation 25 in Au-Yeung (2003)
.Alpha4 <- function( v, My = NA ) { 
  
  n <- (v[4] - 4*v[1]*v[3] + 6*v[1]^2*v[2] - 3*v[1]^4)
  d <- (v[2]-v[1]^2)^4
  alpha4 <- n / d
  
  return(alpha4)
}

#' .Delta_GLD
#' Helper function to compute the distance from the expected values.
#' The output of this function must be optimised.
#' @param L Vector containing lambda3 and 4
#' @param A3,A4 Desired third and fourth moment
#' @param MinMax 0 Indicates that the function shall minimise the result to
#' support classic optimisation algorithms; 1 indicates maximising for the 
#' use of genetic algorithms
.Delta_GLD <- function(L, A3, A4, MinMax = 0) {
  #if(min(L) < -0.25) return(Inf) # Required: min(lambda3,lambda4) > -1/4 
  if(length(L) == 1) L <- c(L, L)
  
  lambda3 <- L[1] # for readibility only
  lambda4 <- L[2] # for readibility only
  v1 <- .v1f(lambda3, lambda4) #ifelse(is.na(My), .v1f(lambda3, lambda4), My)
  v2 <- .v2f(lambda3, lambda4)
  v3 <- .v3f(lambda3, lambda4)
  v4 <- .v4f(lambda3, lambda4)
  
  Delta3 <- abs(A3 - .Alpha3(c(v1, v2, v3, v4)))
  Delta4 <- abs(A4 - .Alpha4(c(v1, v2, v3, v4)))

  if(MinMax == 0)
    return( Delta3 + Delta4 )
  else
    return( 1 / (Delta3 + Delta4) )
}


#' TODO: This is a STUB
.GLDTypeConstraints <- function( Type = c("")) {
  # Class I (lambda3 < 1, lambda4 < 1): Unimodal densities with continuous tails.
  
  # Class Ia (lambda3, lambda4 ≤ ½), 
  
  # Class Ib (½ < lambda3 < 1, lambda4 ≤½), and
  
  # Class Ic (½ < lambda3 <1, ½ < lambda4 < 1).
  
  # Class II (lambda3 > 1, lambda4 < 1): Monotone pdfs similar to those of the exponential or chi²
  # distributions. The left tail is truncated.
  
  # Class III (1 < lambda3 < 2 , 1< lambda4 < 2): U-shaped densities with both tails truncated.
  
  # Class IV (lambda3 >2, 1< lambda4 < 2): Rarely occurring S-shaped pdfs with one mode and one
  # antimode. Both tails are truncated.
  
  # Class V (lambda3 > 2, lambda4 > 2): Unimodal pdfs with both tails truncated.
}



#' FindPDF_GLD
#'
#' @param Moments Vector holding the first four moments.
#' @param Tolerance numeric ≥ 0. Differences smaller than tolerance 
#' are not reported. The default value is close to 1.5e-8.
#' 
#' @details Moments are µ (mean), $\sigma$ (standard deviation), skewness and kurtosis.
#'
#' @return
#' @export
#'
#' @examples
FindPDF_GLD <- function( Moments, Tolerance = sqrt(.Machine$double.eps) ) {
  Lambda <- numeric(4)
  names(Moments) <- c("Mean", "Var", "Skew", "Kurt")
  
  # Find Lambda3 and 4
  if(Moments["Skew"] == 999999) { # Symmetry: lambda3 == lambda4
    # Since (lambda3 == lambda4) we need only 1 dimension to optimise
    Start <- runif(1, -0.25, 50)
    Tolerance <- 1E-4
    r <- optim(Start, method = c("Brent"), lower = -0.25, upper = 50, 
               fn = .Delta_GLD, A3 = Moments["Skew"], A4 = Moments["Kurt"], 
               control = list(reltol = Tolerance, maxit = 50000) )
    Lambda[3] <- r$par[1]
    Lambda[4] <- r$par[1]
  } else {
    GA <- ga(type = "real-valued", fitness = function(x, ...) .Delta_GLD(x, ...),
             A3 = Moments["Skew"], A4 = Moments["Kurt"], MinMax = 1, 
             lower = c(-0.25, -0.25), upper = c(50, 50),
             popSize = 100, maxiter = 500, run = 100, pmutation = 0.15)
    #print(summary(GA))
    #plot(GA)
    #GA@solution[1,] - c(-0.15, -0.15)
    r <- optimx(GA@solution[1,], method = c("Nelder-Mead"), 
                fn = .Delta_GLD, A3 = Moments["Skew"], A4 = Moments["Kurt"], 
               control = list(reltol = Tolerance*100, maxit = 50000) )
    ###TODO: check if successful, first
    Lambda[3] <- r$x1 #r$par[1]
    Lambda[4] <- r$x2 #r$par[2]
  }
  
  v1 <- .v1f(Lambda[3], Lambda[4])
  v2 <- .v2f(Lambda[3], Lambda[4])
  
  # eqn 27
  Lambda[2] <- sqrt(v2 - v1^2) / sqrt(Moments[2]) 
  # eqn 28
  Lambda[1] <- Moments[1] + 1/Lambda[2] * ( (1/(Lambda[3]+1)) + (1/(Lambda[4]+1)) )
  
  print(Lambda)
  return(Lambda)
}


#  l <- FindPDF_GLD(c(0, 1, 0, 3), Tolerance = 1E-6)
# 
#  plotgld(lambda1 = l[1], lambda2 = l[2], lambda3 = l[3], lambda4 = l[4],  
#          param = "fmkl", lambda5 = NULL, add = NULL, truncate = 0,  
#          bnw = FALSE, col.or.type = 1, granularity = 10000, xlab = "x",  
#          ylab = NULL, quant.probs = seq(0,1,.25), new.plot = NULL)
# 
# source("./SimPdf.R")
# print(SimPdf( dnorm, dgl, Xmin=-10, Xmax=10, Steps=1e5, Args1 = NULL, Args2 = list(l) ))
#print(MomentSet( rgl(1e7, lambda1 = l[1], lambda2 = l[2], lambda3 = l[3], lambda4 = l[4]) ))
