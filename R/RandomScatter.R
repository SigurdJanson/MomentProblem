

#' GetRandomPoints
#' Gives coordinates of point in a cubic space that are evenly distributed.
#' @param Dim Dimensions of the space
#' @param N Number of requested points
#' @details Creates `N` random points in a `Dim`-dimensional cubi space.
#' @return A matrix with all point-vectors in rows and `Dim` columns.
#' Uses coordinates of a unit cube. The number of points may be greater 
#' than `N` because as number of points only a power of `Dim` is possible.
#' Number of points will be greater or equal than `N`.
#' @export
#' @examples 
#' GetRandomPoints(2, 2) # returns 2 points in 2-d space
#' GetRandomPoints(3, 2) # returns 2 points in 3-d space
#' GetRandomPoints(4, 2) # returns 2 points in 4-d space
GetRandomPoints <- function( Dim, N ) {
  # PRECONDITIONS
  if(Dim < 1) stop("Cannot handle cubes in less than 1 dimensions")
  Dim <- as.integer(Dim)
  if(N < 1) stop("At least 'N = 1' desired point must be set")
  N <- as.integer(N)
  
  Values <- runif(Dim * N, 0, 1)
  return(matrix(Values, nrow = N, ncol = Dim))
}