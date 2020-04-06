

#' GetCubicPoints
#' Gives coordinates of point in a cubic space that are evenly distributed.
#' @param Dim Dimensions of the space
#' @param N Number of requested points
#' @details The space of a unit cubic is divided into many tiny 
#' (sub-) cubes, like a chessboard in 2d or a Rubic Cube in 3d.
#' @return A matrix with all point-vectors in rows and `Dim` columns.
#' Uses coordinates of a unit cube. The number of points may be greater 
#' than `N` because as number of points only a power of `Dim` is possible.
#' Number of points will be greater or equal than `N`.
#' @export
#' @examples 
#' GetCubicPoints(2, 2) # returns 4 points
#' GetCubicPoints(3, 2) # returns 8 points
#' GetCubicPoints(4, 2) # returns 16 points
GetCubicPoints <- function( Dim, N ) {
  # PRECONDITIONS
  if(Dim < 1) stop("Cannot handle cubes in less than 1 dimensions")
  Dim <- as.integer(Dim)
  if(N < 1) stop("At least 'N = 1' desired point must be set")
  N <- as.integer(N)

  # Special case: only 1 point, choose mid-point
  if(N == 1) {
    Result <- matrix(rep(0.5, Dim), nrow = 1)
    return( Result )
  }
  
  # Unit cube can be split into Count^Dim sub-cubes
  # Get 'CubeCount' that returns a number of points closest to 'N'
  EdgeLenCut <- ceiling( N^(1/Dim) ) #
  EdgeLen    <- 1 / EdgeLenCut       # edge length
  CubeCount  <- EdgeLenCut^N         # total number of sub-cubes
  
  # Get center coordinates of each sub-cube based on 1 dimension
  BaseVector1D <- EdgeLen * 0:(EdgeLenCut-1) + (EdgeLen/2)
  # Expand single dimension to all dimensions
  CoordinateCombns <- replicate(Dim, BaseVector1D, simplify = FALSE)
  Coordinates <- data.matrix(expand.grid(CoordinateCombns))
  dimnames(Coordinates) <- NULL
  
  return(Coordinates)
}


