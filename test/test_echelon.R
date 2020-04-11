setwd("..")
source("./R/echelon.R")
setwd("./test")


test_that("Echelon (not reduced)", {
  # Echelon_Row_Reduction_in_R.pdf - self-computed solution
  M <- matrix(c(1,-1,1,1, 1,-1,-1, 2, 1,1,-1,3, 1,1,1,4), nrow = 4, byrow=TRUE)
  e <- matrix(c(1,-1,1,1, 0, 1, -1,1, 0,0,1,-0.5, 0,0,0,1), nrow = 4, byrow=TRUE)
  expect_equal(Echelon(M), e)
  # LinearDependence.pdf, example 3, S.2
  M <- matrix(c(1,2,3,2, 2,5,5,5, 2,6,4,6), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,2,3,2, 0,1,-1,1, 0,0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M), e)
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(0,1,2, 1,2,1, 2,7,8), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,2,1, 0,1,2, 0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M), e)
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(0,1, 1,2, 0,5), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,2, 0,1, 0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M), e)
  
  M <- matrix(c(1,2,3, 4,5,6), nrow = 2, byrow = TRUE)
  e <- matrix(c(1,2,3, 0,1,2), nrow = 2, byrow = TRUE)
  expect_equal(Echelon(M), e)

  # These are already in echelon form and should return the input without changes
  M <- matrix(c(1,2,3,4, 0,0,1,3, 0,0,0,1, 0,0,0,0), nrow = 4, byrow = TRUE)
  expect_identical(Echelon(M), M)
  
  M <- matrix(c(1,2,0,0, 0,0,1,0, 0,0,0,1, 0,0,0,0), nrow = 4, byrow = TRUE)
  expect_identical(Echelon(M), M)
  
  M <- matrix(c(1,1,3, 0,0,1), nrow = 2, byrow = TRUE)
  expect_identical(Echelon(M), M)

  M <- matrix(c(1,2,3, 0,1,2), nrow = 2, byrow = TRUE)
  expect_identical(Echelon(M), M)

  M <- matrix(c(1,0,1,2, 0,1,1,1, 0,0,0,1), nrow = 3, byrow = TRUE)
  expect_identical(Echelon(M), M)
  
  M <- matrix(c(1,0,2,-2, 0,1,0,2, 0,0,0,1), nrow = 3, byrow = TRUE)
  expect_identical(Echelon(M), M)
})


test_that("Reduced Echelon", {
  # Echelon_Row_Reduction_in_R.pdf
  M <- matrix(c(1,-1,1,1, 1,-1,-1,2, 1,1,-1,3, 1,1,1,4), nrow = 4, byrow=TRUE)
  e <- matrix(c(1,0,0,0,  0,1,0,0, 0,0,1,0, 0,0,0,1), nrow = 4, byrow=TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # LinearDependence.pdf, example 1, S.2 - self-computed solution
  M <- matrix(c(1,2,3,2, 2,5,5,5, 2,6,4,6), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,5,0, 0,1,-1,1, 0,0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # LinearDependance.pdf, example 2, S.3
  M <- matrix(c(1,2,4, 3,0,6, 2,2,6, 4,1,9), nrow = 4, byrow = TRUE)
  e <- matrix(c(1,0,2, 0,1,1, 0,0,0, 0,0,0), nrow = 4, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # LinearDependence.pdf, example 3, S.4
  M <- matrix(c(2,1,7, 5,6,7, 3,4,3), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,5, 0,1,-3, 0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # LinearDependence.pdf, example 4, S.5
  M <- matrix(c(1,2,2,   2,5,6, 3,5,4, 2,5,6), nrow = 4, byrow=TRUE)
  e <- matrix(c(1,0,-2,  0,1,2, 0,0,0, 0,0,0), nrow = 4, byrow=TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(0,1,2, 1,2,1, 2,7,8), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,-3, 0,1,2, 0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(0,1, 1,2, 0,5), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0, 0,1, 0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  # https://www.wolframalpha.com/widgets/view.jsp?id=7eab3cc8b79a0665f796eea7c14b2d90
  M <- matrix(c(1,2,3, 4,5,6), nrow = 2, byrow = TRUE)
  e <- matrix(c(1,0,-1, 0,1,2), nrow = 2, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  
  M <- matrix(c(1,0,1,2, 0,1,1,1, 0,0,0,1), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,1,0, 0,1,1,0, 0,0,0,1), nrow = 3, byrow = TRUE)
  expect_identical(Echelon(M, Reduced = TRUE), e)
  
  
  #
  # These are already in echelon form and will be reduced
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(1,2,3,4, 0,0,1,3, 0,0,0,1, 0,0,0,0), nrow = 4, byrow = TRUE)
  M <- matrix(c(1,2,0,0, 0,0,1,0, 0,0,0,1, 0,0,0,0), nrow = 4, byrow = TRUE)
  expect_identical(Echelon(M, Reduced = TRUE), M)
  # https://stattrek.com/matrix-algebra/echelon-transform.aspx
  M <- matrix(c(0,1,2, 1,2,1, 0,0,0), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,-3, 0,1,2, 0,0,0), nrow = 3, byrow = TRUE)
  expect_equal(Echelon(M, Reduced = TRUE), e)
  
  M <- matrix(c(1,2,3, 0,1,2), nrow = 2, byrow = TRUE)
  e <- matrix(c(1,0,-1, 0,1,2), nrow = 2, byrow = TRUE)
  expect_identical(Echelon(M, Reduced = TRUE), e)
  
  M <- matrix(c(1,0,2,-2, 0,1,0,2, 0,0,0,1), nrow = 3, byrow = TRUE)
  e <- matrix(c(1,0,2,0, 0,1,0,0, 0,0,0,1), nrow = 3, byrow = TRUE)
  expect_identical(Echelon(M, Reduced = TRUE), e)
  
  #
  # These are already in reuced echelon form and will be kept
  M <- matrix(c(1,1,3, 0,0,1), nrow = 2, byrow = TRUE)
  expect_identical(Echelon(M), M)

  M <- matrix(c(1,0,3, 0,1,1), nrow = 2, byrow = TRUE)
  expect_identical(Echelon(M), M)

  M <- matrix(c(1,2,0,0, 0,0,1,0, 0,0,0,1, 0,0,0,0), nrow = 4, byrow = TRUE)
  expect_identical(Echelon(M), M)
})


test_that("Echelon - Comparison with pracma::rref", {
  #' rref
  #' @param A 
  #' @source Taken from the pracma package (https://rdrr.io/cran/pracma/)
  #' @author Hans W. Borchers <hwborchers@googlemail.com>
  rref <- function(A) {
    stopifnot(is.numeric(A))
    if (!is.matrix(A))
      stop("Input parameter 'A' must be a matrix.")
    
    nr <- nrow(A); nc <- ncol(A)
    #tol <- eps() * max(nr, nc) * max(abs(A)) # original line
    tol <- .Machine$double.eps * max(nr, nc) * max(abs(A))
    
    r <- 1
    for (i in 1:nc) {
      pivot <- which.max(abs(A[r:nr, i]))
      pivot <- r + pivot - 1
      m <- abs(A[pivot, i])
      
      if (m <= tol) {
        A[r:nr, i] <- 0  # zeros(nr-r+1, 1)
      } else {
        A[c(pivot, r), i:nc] <- A[c(r, pivot), i:nc]
        A[r, i:nc] <- A[r, i:nc] / A[r, i]
        if (r == 1) {
          ridx <- c((r+1):nr)
        } else if (r == nr) {
          ridx <- c(1:(r-1))
        } else {
          ridx <- c(1:(r-1), (r+1):nr)
        }
        A[ridx, i:nc] <- A[ridx, i:nc] - 
          A[ridx, i, drop=FALSE] %*% A[r, i:nc, drop=FALSE]
        if (r == nr) break
        r <- r+1
      }
    }
    A[abs(A) < tol] <- 0
    return(A)
  }
  
  # Squared matrices
  for(Size in c(2, 5, 10, 11, 20, 50, 100)) {
    M <- matrix(
      sample.int(Size * Size, replace = TRUE),
      ncol = Size)
    Tolerance <- .Machine$double.eps * Size * max(abs(M))
    expect_equal(Echelon(M, Reduced = TRUE, tolerance = Tolerance), rref(M))
  }
  
})


test_that("Echelon - Comparison with rref by Fox", {
  #' @source https://markmail.org/thread/wjq2xro6sibzf5fm
  rref <- function(A, tol = sqrt(.Machine$double.eps),
                   verbose=FALSE, fractions=FALSE){
    ## A: coefficient matrix
    ## tol: tolerance for checking for 0 pivot
    ## verbose: if TRUE, print intermediate steps
    ## fractions: try to express nonintegers as rational numbers
    ## Written by John Fox
    if (fractions) {
      mass <- require(MASS)
      if (!mass) stop("fractions=TRUE needs MASS package")
    }
    if ((!is.matrix(A)) || (!is.numeric(A)))
      stop("argument must be a numeric matrix")
    n <- nrow(A)
    m <- ncol(A)
    for (i in 1:min(c(m, n))){
      col <- A[,i]
      col[1:n < i] <- 0
      # find maximum pivot in current column at or below current row
      which <- which.max(abs(col))
      pivot <- A[which, i]
      if (abs(pivot) <= tol) next     # check for 0 pivot
      if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows
      A[i,] <- A[i,]/pivot            # pivot
      row <- A[i,]
      A <- A - outer(A[,i], row)      # sweep
      A[i,] <- row                    # restore current row
      if (verbose)
        if (fractions) print(fractions(A))
      else print(round(A,round(abs(log(tol,10)))))
    }
    for (i in 1:n)
      if (max(abs(A[i,1:m])) <= tol)
        A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
    if (fractions) fractions (A)
    else round(A, round(abs(log(tol,10))))
  }
  
  succeed()
})

test_that("Echelon - test basic criteria", {
# A matrix that has undergone Gaussian elimination is said to be in row echelon form or, more properly, "reduced echelon form" or "row-reduced echelon form." Such a matrix has the following characteristics:
# 1. All zero rows are at the bottom of the matrix
# 2. The leading entry of each nonzero row after the first occurs to the right of the leading entry of the previous row.
# 3. The leading entry in any nonzero row is 1.
# 4. All entries in the column above and below a leading 1 are zero.
# Another common definition of echelon form only requires zeros below the leading ones, while the above definition also requires them above the leading ones.
  succeed()
})