
Echelon <- function(M, Reduced = FALSE, 
                    tolerance = sqrt(.Machine$double.eps) ) {
  # PRECONDITIONS
  if (!is.numeric(M) && !is.matrix(M))
    stop("Input must be a numeric matrix.")
  
  # RESULT
  NR <- nrow(M)
  NC <- ncol(M)
  
  # Find pivot (first non-zero entry in column of the matrix).
  # I.e. sort into right order first
  for(r in 1L:min(NR, NC)) {
    if(abs(M[r, r]) < tolerance) {# sort only when necessary
      Pivot <- which.max(abs(M[r:NR, r]) > tolerance)
      Pivot <- Pivot + (r-1) # make sure that Pivot in [r,NR]
      
      # Move pivot row to row 'r'
      if (abs(M[Pivot, r]) > tolerance) { # only sort if Pivot is not zero
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
    if(abs(M[r, r]) > tolerance)
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
      
      Subtrahend <- M[(r-1):1, PivCol]
      Subtrahend <- Subtrahend %*% t(M[r, PivCol:NC])
      M[(r-1):1, PivCol:NC] <- M[(r-1):1, PivCol:NC] - Subtrahend
    }
  }
  
  # RETURN
  M[abs(M) < tolerance] <- 0
  return(M)
}






