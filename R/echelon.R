
#' Echelon
#' Compute the row echelon form of a matrix.
#' @param M A numeric matrix with at least two rows and columns
#' @param Reduced If `TRUE` the function will return the echelon in its 
#' reduced form.
#' @param Tolerance The smallest value that is considered different from
#' zero.
#' @return A numeric matrix in the requested form.
#' @author Jan Seifert
#' @source Instructions from https://stattrek.com/matrix-algebra/echelon-transform.aspx
#' @export
#' @examples
#' M <- matrix(c(1,2,3,2, 2,5,5,5, 2,6,4,6))
#' Echelon(M)
#' # Result: matrix(c(1,2,3,2, 0,1,-1,1, 0,0,0,0), nrow = 3, byrow = TRUE)
#' M <- matrix(c(1,2,3,2, 2,5,5,5, 2,6,4,6), nrow = 3, byrow = TRUE)
#' Echelon(M, Reduced = TRUE)
#' # Result: matrix(c(1,0,5,0, 0,1,-1,1, 0,0,0,0), nrow = 3, byrow = TRUE)
Echelon <- function(M, Reduced = FALSE, 
                    Tolerance = sqrt(.Machine$double.eps) ) {
  # PRECONDITIONS
  if (!is.numeric(M) && !is.matrix(M))
    stop("Input must be a numeric matrix.")
  if (length(dim(M)) != 2) stop("Matrix must be 2-dimensional")
  
  NR <- nrow(M)
  NC <- ncol(M)
  if(NR < 2 || NC < 2) stop("Matrix needs more rows and cols than 1")
  
  # RESULT
  # Find pivot (first non-zero entry in column of the matrix).
  # I.e. sort into right order first
  for(r in 1L:min(NR, NC)) {
    if(abs(M[r, r]) < Tolerance) {# sort only when necessary
      Pivot <- which.max(abs(M[r:NR, r]) > Tolerance)
      Pivot <- Pivot + (r-1) # make sure that Pivot in [r,NR]
      
      # Move pivot row to row 'r'
      if (abs(M[Pivot, r]) > Tolerance) { # only sort if Pivot is not zero
        if (r > 1)
          RowSort <- c( 1:(r-1), Pivot, (r:NR)[-(Pivot-r+1)] )
        else 
          RowSort <- c( Pivot, (r:NR)[-(Pivot-r+1)] )
        M <- M[RowSort,]
      }
    }
    # Multiply each element in pivot row by the inverse of the pivot, 
    # so the pivot equals 1.
    #\for(r in 1L:min(NR, NC))
    if(abs(M[r, r]) > Tolerance)
      M[r, ] <- M[r, ] / M[r, r] # from now on Pivot row is r
      
    # Add multiples of pivot row to each of the lower rows, so every 
    # element in pivot column of the lower rows equals 0.
    #for(r in 1:min(NR-1L, NC)) {
    if(r < NR) {
      Subtrahend <- M[(r+1):NR, r] #/ 1
      Subtrahend <- Subtrahend %*% t(M[r, r:NC])
      M[(r+1):NR, r:NC] <- M[(r+1):NR, r:NC] - Subtrahend
    }
  }

  # To get the matrix in reduced row echelon form, process 
  # non-zero entries above each pivot.
  if(Reduced) {
    for(r in NR:2L) {
      # Find pivot column
      PivCol <- match(1, M[r,]) #which.max(M[r,])
      if (is.na(PivCol)) next
      # Add multiples of the pivot row to each of the upper rows, 
      # until every element above the pivot equals 0.
      Subtrahend <- M[(r-1):1, PivCol]
      Subtrahend <- Subtrahend %*% t(M[r, PivCol:NC])
      M[(r-1):1, PivCol:NC] <- M[(r-1):1, PivCol:NC] - Subtrahend
    }
  }
  
  # RETURN
  M[abs(M) < Tolerance] <- 0
  return(M)
}



isEchelon <- function( M, Reduced = FALSE ) {
  NC <- ncol(M)
  NR <- nrow(M)
  # If there are more rows than cols, additional rows must be all 0
  if (NR > NC) {
    if (sum(M[(NC+1):NR,]) > 0) return(FALSE)
    # Drop zero rows to get a square matrix
    M <- M[1:NC,]
    NR <- NC
  }
  # At this point is guaranteed: NR <= NC
  
  # All values in lower triangle must be zero
  if (sum( M[lower.tri(M)] ) > 0) return(FALSE)
  # 
  R1 <- 1
  while (!is.na(R1)) {
    # Diagonal can only be 1 or zero
    D <- diag(M)
    if (!all(D %in% c(0, 1))) return(FALSE)
    # Ignore all rows where Pivot is in diagonal
    R1 <- match(0, D)
    if(is.na(R1)) return(TRUE)
    # ... and get new sub-matrix to test the diagonal
    # Get upper left corner of sub-atrix where M[R1,j] == 1
    C1 <- match(1, M[R1,])
    if (is.na(C1)) {  # We have left last pivot behind
      if (sum(M[R1:NR, ]) > 0)
        return(FALSE)
      else
        return(TRUE)
    }
    M <- M[R1:NR, C1:NC]
    NR <- NR-R1
    NC <- NC-C1
  }
}



isInconsistent <- function( M ) {
  lower.tri(mat)
}



