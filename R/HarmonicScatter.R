
#' .HarmonyR
#' @describeIn .Harmony
.HarmonyR <- function( Dim, MaxDepth = 50L, .Depth = 0 ) {
  HarmonyConstants <- c(
    1 + sqrt(5/4)-0.5,   # Fibonacci 
    1.32471795724474602596090885447809, # Plastic number
    1.22074408460575947536168534910883  # 
  )
  if (Dim < length(HarmonyConstants))
    return(HarmonyConstants[Dim])
    
  if(.Depth == MaxDepth) 
    return(1)
  else {
    result <- .HarmonyR(Dim, MaxDepth = MaxDepth, .Depth = .Depth+1)
    return( (1 + result )^(1/(Dim+1)) )
  }
}


#' .Harmony
#' Compute harmonious numbers for a number of dimensions `Dim`.
#' @param Dim Number of dimensions of the space.
#' @param Tolerance Desired precision. The sequence stops if the 
#' following summand is smaller than tolerance. The default 
#' value is close to 1.5e-8 (numeric ≥ 0).
#' @param MaxDepth Number of recursive
#' @param .Depth Used by recursive algormithm
#' @note `.HarmonyR` function uses a recursive algorithm.
#' @source Roberts, M. (2018). The Unreasonable Effectiveness 
#' of Quasirandom Sequences. https://t1p.de/z30y
#' @author Jan Seifert
#' @return
#' @examples
.Harmony <- function( Dim, Tolerance = sqrt(.Machine$double.eps) ) {
  HarmonyConstants <- c(
    1 + sqrt(5/4)-0.5,   # Fibonacci 
    1.32471795724474602596090885447809, # Plastic number 
    1.22074408460575947536168534910883  # 
  )
  if (Dim < length(HarmonyConstants))
    return(HarmonyConstants[Dim])

  x <- 1.0 
  repeat {
    new <- (1 + x)^(1/(Dim+1))
    if(abs(new - x) < Tolerance) return(new)
    x <- new
  }
}



#' GetHarmonicPoints
#' Computes a `Dim`-dimensional low discrepancy 
#' quasirandom sequences using 'Harmonious numbers'.
#' @param Dim The number of dimensions.
#' @param N The number of samples to compute.
#' @param Seed Any number between 0 and 1 (default 0.5).
#' @details This function is uses quasi-random properties of 
#' infinite sequences. At the dimension of 1 it is equal to the
#' Finbonacci number. Because of it's properties it is recommended
#' to use this function with larger values of `N` (probably >= 50 
#' for 2 dimensions).
#' 
#' Seed can be any number. A common default setting 
#' is `Seed = 0`. but Roberts (2018) recommends `Seed = 0.5`. Change
#' `Seed` to get varying results.
#' @return A matrix with `Dim` columns and `N` rows (each row is 
#' a position vector). Values range from zero to one.
#' @export
#'
#' @references Roberts, M. (2018). The Unreasonable Effectiveness 
#' of Quasirandom Sequences. https://t1p.de/z30y
#' @examples
GetHarmonicPoints <- function( Dim, N = 20, Seed = 0.5 ) {
  # d=2   # dimensions  # n=50  # number of required points 
  G <- .Harmony(Dim) # g = phi(d)
  
  # for j in range(d): 
  #   alpha[j] = pow(1/g,j+1) %1 
  Alpha <- (1/G^(1:Dim)) %% 1
  
  Z <- matrix(nrow = N, ncol = Dim)
  for (d in 1:Dim) {
    Z[, d] <- (Seed + Alpha[d] * 1:N) %% 1
  }
  return(Z)
}

