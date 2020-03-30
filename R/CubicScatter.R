

#' GetCubicPoints
#' Gives coordinates of point in a cubic space that are evenly distributed.
#' @param Dim Dimensions of the space
#' @param N Number of requested points
#' @details The space of a unit cubic is divided into many tiny 
#' (sub-) cubes, like a chessboard in 2d or a Rubic Cube in 3d.
#' @return A matrix with all point-vectors in rows and `Dim` columns.
#' Uses coordinates of a unit cube. The number of points may be greater 
#' than `N` because as number of points only a power of `Dim` is possible.
#' Number of points will be greater ir equal than `N`.
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



#' GetCubicPoints.Fragment
#' Gives point coordinates for a cubic space by splitting the space 
#' into concentric sub-cubes.
#' 
#' This function is unfinished. At this point this algorithm uses distances 
#' between vertices (0-face) along the edges (1-face) and the diagonal (n-face).
#' All diagonals for x-faces with 1 < x < n are ignored.
#' @param Dim 
#' @param N 
#' @details 
#' @return A matrix with all point-vectors in rows and `Dim` columns.
#' Uses coordinates of a unit cube.
#' @export
#' @source https://en.wikipedia.org/w/index.php?title=Hypercube&oldid=947339674
#' @examples
GetCubicPoints_Fragment <- function( Dim, N ) {
  Result <- matrix(rep(0.5, Dim), nrow = 1)
  
  # Special case: only 1 point, choose mid-point
  if(N == 1) return( Result )
  
  # Get number of edges & vertices of 'Dim'-cube
  Vertices <- 2L^Dim 
  Edges    <- 2L^(Dim - 2L) * choose(Dim, 2L)
  
  # Add first sub-cube (diagonal is always longer than un-sliced edge)
  Points   <- 1 + Vertices
  SubCubes <- 1
  Slices   <- c(0)
  repeat {
    if (Points > N) break
    
    # Determine longest distance (i.e. edge or diagonal)
    EdgeLen <- 2/(SubCubes+1) * seq(1, SubCubes, by=1) # Total len
    EdgeLen <- EdgeLen / (Slices+1) # Slice total length
    DiagLen <- sqrt(Dim) / SubCubes
    
    # Slice longest distance
    if(DiagLen >= max(EdgeLen)) {
      # Add sub-cube
      Points   <- Points + Vertices
      SubCubes <- SubCubes + 1
      Slices   <- c(Slices, 0)
    } else {
      FirstMax <- which.max(EdgeLen) # use first maximum and slice it
      Slices[FirstMax] <- Slices[FirstMax] + 1
      Points <- Points + Edges # We have 1 new point per edge
    }
  }#repeat
  
  # Compute point coordinates
  for (Cube in SubCubes) {
    # Vertex coordinates
    Delta <- 1/(SubCubes+1)
    CoordinateCombns <- replicate(2, c(Delta, 1-Delta), simplify = FALSE)
    VertexCoordinates <- data.matrix(expand.grid(CoordinateCombns))
    
    # Slicing points
    if(Slices[Cube] > 0) {
      for (Vertex1 in seq(1, nrow(VertexCoordinates)-1, by=1)) {
        for (Vertex2 in seq(2, nrow(VertexCoordinates), by=1)) {
          # How many coordinates differ between two vectors?
          Diff <- VertexCoordinates[Vertex1] == VertexCoordinates[Vertex2]
          # Use only neighbouring vertices
          if (sum(Diff) == 1) {
            #TODO: Compute slice coordinates
          }
        }
      }
    }
    
    # Add to result
    Result <- rbind(Result, VertexCoordinates)
    Result <- rbind(Result, SliceCoordinates)
  }
}

