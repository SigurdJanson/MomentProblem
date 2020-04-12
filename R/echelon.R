#' Echelon
#' Compute the row echelon form of a matrix.
#' @param M A numeric matrix with at least two rows and columns.
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



#' isEchelon
#' Test if a matrix satisfies the echelon conditions.
#' @param M A numeric matrix.
#' @param Reduced If `TRUE`, this function tests for reduced echelon form.
#' `FALSE` (the default) tests for the plain form.
#'
#' @return `TRUE` if the echelon form is satsified.
#' @export
#' @author Jan Seifert
#' @examples
isEchelon <- function( M, Reduced = FALSE ) {
  # PRECONDITIONS
  if (!is.numeric(M) && !is.matrix(M))
    stop("Input must be a numeric matrix.")
  if (length(dim(M)) != 2) stop("Matrix must be 2-dimensional")
  NC <- ncol(M)
  NR <- nrow(M)
  if (NR == 1)
    return(M[1,1] == 1)
  if (NC == 1) {
    if(M[1,1] != 1)
      return(FALSE)
    else # [1,1] is the pivot, all values below must be 0
      return(all(M[2:nrow(M), 1] == 0))
  }
      
  # RESULT
  # If there are more rows than cols, additional rows must be all 0
  if (NR > NC) {
    if (any(M[(NC+1):NR,] > 0)) return(FALSE)
    # Drop zero rows to get a square matrix
    M <- M[1:NC,]
    NR <- NC
  }
  # At this point is guaranteed: NR <= NC
  
  # All values in lower triangle must be zero
  if (any( M[lower.tri(M)] > 0)) return(FALSE)
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
      return(all(M[R1:NR, ] == 0))
    }
    M <- M[R1:NR, C1:NC]
    NR <- NR-R1
    NC <- NC-C1
  }
  return(TRUE) # code should never reach this point
}



hasSolutions <- function( M ) {
  # PRECONDITIONS (other conditions are handled by called functions)
  if(!isEchelon(M)) M <- Echelon(M)
  # RESULT
  PivotCount <- 0L
  NC <- ncol(M)
  NR <- nrow(M)
  for(c in NC:1L) {
    PivotRow <- rev(M[,c]) == 1L
    PivotRow <- NR - which.max(PivotRow) +1
    # Verify
    if (M[PivotRow, c] == 1) {
      if (c == 1L)
        PivotCount <- PivotCount +1L
      else if (all(M[PivotRow, 1L:(c-1L)] == 0))
        PivotCount <- PivotCount +1L
    }
    
    # Pivot in last col, system is inconsistent
    if(PivotCount == 1L && c == NC)
      return(0L)
  }
  return(ifelse(PivotCount == NC-1L, 1L, Inf))
}


MatrixRank <- function( M ) {
  # PRECONDITIONS (other conditions are handled by called functions)
  if(!isEchelon(M)) M <- Echelon(M)
  # RESULT
  NotZero <- apply(M, 1, function(x) any(x != 0))
  return(sum(NotZero))
}


SolutionSpace <- function( M ) {
  
}

