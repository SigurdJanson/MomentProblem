#' nmkb
#' Nelder-Mead with Box constraints
#' @details see https://rdrr.io/cran/dfoptim/man/nmkb.html
#' @note The parameters controlling the behaviour of the simplex are chosen
#' according to Gao & Han (2010) instead of the original choice by Nelder & Mead.
#' @source This function is a copy of dfoptim::nmkb(). This algorithm is based 
#' on the [Matlab code of Prof. C.T. Kelley, given in his book "Iterative methods 
#' for optimization"](https://archive.siam.org/books/kelley/fr18/OPT_CODE/nelder.m). 
#' It was implemented for R by Ravi Varadhan with permission of 
#' Prof. Kelley and SIAM. However, there are some non-trivial modifications of 
#' the algorithm made by Ravi Varadhan. Some refactoring and further modifications
#' have been done by Jan Seifert.
#' @references 
#' Gao, F., & Han, L. (2010). Implementing the Nelder-Mead simplex algorithm with adaptive parameters. Computational Optimization and Applications, 51, 259–277. https://doi.org/10.1007/s10589-010-9329-3
#' Kelley, C. T. (1999). Iterative Methods for Optimization. Society for Industrial and Applied Mathematics. https://doi.org/10.1137/1.9781611970920
#' Nelder, J. A., & Mead, R. (1965). A simplex method for function minimization. Computer Journal, 7, 308–313.
# TODO: find speed optimisations
#
# Change re-ordering algorithm to enhance speed and to handle tie-breaks according to 
# Lagarias, J. et al (1998). "Convergence Properties of the Nelder–Mead Simplex Method in Low Dimensions"
# In a nonshrink condition only one point has changed and the ordering can be updated in linear time 
# (at most n comparisons) by one step of straight insertion sort. 
nmkb <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
  ctrl <- list(tol = 1e-06, maxfeval = min(5000, max(1500, 20 * length(par)^2)), 
               regsimp = TRUE, maximize = FALSE, restarts.max = 3, trace = FALSE)
  namc <- match.arg(names(control), choices = names(ctrl), several.ok = TRUE)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in 'control': ", namc[!(namc %in% names(ctrl))])
  if (!is.null(names(control))) ctrl[namc] <- control
  ftol <- ctrl$tol
  maxfeval <- ctrl$maxfeval
  regsimp <- ctrl$regsimp
  restarts.max <- ctrl$restarts.max
  maximize <- ifelse(ctrl$maximize, -1, 1)
  trace <- ctrl$trace

  # Dimensions of paramter space
  n <- length(par) 
  if (n == 1) 
    stop(call. = FALSE, "Use `optimize` for univariate optimization")
  if (n > 30) 
    warning("Nelder-Mead should not be used for high-dimensional optimization")
  
  # Spatial distortion functions to handle box constraints
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  
  # Set variables required to handle box constraints
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  
  if (any(c(par < lower, upper < par))) stop("Infeasible starting values!", call.=FALSE)
  low.finite <- is.finite(lower)
  upp.finite <- is.finite(upper)
  c1 <- low.finite & upp.finite  # both lower and upper bounds are finite 
  c2 <- !(low.finite | upp.finite) # both lower and upper bounds are infinite
  c3 <- !(c1 | c2) & low.finite # finite lower bound, infinite upper bound
  c4 <- !(c1 | c2) & upp.finite  # finite upper bound, infinite lower bound
  
  if (all(c2)) stop("Use `nmk()` for unconstrained optimization!", call.=FALSE)
  
  # Set objective function
  fnmb <- function(par, ...) { fn(ginv(par), ...) * maximize }

  x0 <- g(par) # x0 is the starting point
  
  # Initial values: V - vertices of the simplex. f - result of objective function.
  V <- cbind(rep(0, n), diag(n))
  f <- rep(0, n + 1)
  f[1] <- fnmb(x0, ...)
  V[, 1] <- x0
  scale <- max(1, sqrt(sum(x0^2)))

  if (regsimp) { # Create regular simplex
    alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, sqrt(n + 1) - 1)
    V[, -1] <- (x0 + alpha[2])
    diag(V[, -1]) <- x0[1:n] + alpha[1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
  } else {
    V[, -1] <- x0 + scale * V[, -1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
  }
  f[is.nan(f)] <- Inf
  nf <- n + 1     # number of calls of objective function
  # Sort simplex vertices according to f(v)
  ord <- order(f) #TODO: enhance for speed and tie-breaks
  f <- f[ord]
  V <- V[, ord]
  
  rho <- 1              # reflection  alpha
  chi <- 1 + 2/n        # expansion   beta
  gamma <- 0.75 - 1/2/n # contraction gamma
  sigma <- 1 - 1/n      # shrinkage   delta
  
  # Setting oshrink=0 gives vanilla Nelder-Mead, i.e. the simplex is merely
  # reflected without expansion, contraction, or shrinking
  oshrink <- 1L
  
  orth <- 0 #TODO: this variable has no explanation
  
  v <- V[, -1] - V[, 1]
  delf <- f[-1] - f[1]          # "delta f"
  diam <- sqrt(colSums(v^2))    # diameter of the simplex???
  sgrad <- c(solve(t(v), delf))
  alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
  
  # Termination criteria (`nf` has been defined above)
  simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 1])))
  dist <- f[n + 1] - f[1] # distance between f(highest) and f(lowest) vertex
  restarts <- 0 # number of restarts
  
  itc <- 0L # iteration counter
  while (nf < maxfeval & restarts < restarts.max & dist > ftol & simplex.size > 1e-06) {
    happy <- 0L  # `happy == 1` will indicate that new vertex is accepted
    itc <- itc + 1L
    
    # REFLECT (is always first try)
    fbc <- mean(f)
    xbar <- rowMeans(V[, 1:n])  # centroid of the simplex - TODO: enhance computation - see formula (7)
    xr <- (1 + rho) * xbar - rho * V[, n + 1]
    fr <- fnmb(xr, ...)
    nf <- nf + 1
    if (is.nan(fr)) fr <- Inf

    if (fr >= f[1] & fr < f[n]) {
      # Reflection was successful
      happy <- 1L
      xnew <- xr
      fnew <- fr
    } else if (fr < f[1]) { # if reflected point f(r) < f(old): EXPAND
      xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 1]
      fe <- fnmb(xe, ...)
      if (is.nan(fe)) fe <- Inf
      nf <- nf + 1
      if (fe < fr) {
        xnew <- xe
        fnew <- fe
        happy <- 1L
      } else {
        xnew <- xr
        fnew <- fr
        happy <- 1L
      }
    } else if (fr >= f[n] & fr < f[n + 1]) {
      # OUTSIDE CONTRACTION
      xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, n + 1]
      fc <- fnmb(xc, ...)
      if (is.nan(fc)) fc <- Inf
      nf <- nf + 1
      if (fc <= fr) {
        xnew <- xc
        fnew <- fc
        happy <- 1L
      }
    } else if (fr >= f[n + 1]) {
      # INSIDE CONTRACTION
      xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
      fc <- fnmb(xc, ...)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc < f[n + 1]) {
        xnew <- xc
        fnew <- fc
        happy <- 1L
      }
    }
    # Test for sufficient decrease; do oriented shrink if necessary
    if (happy == 1 & oshrink == 1L) {
      fbt <- mean(c(f[1:n], fnew))
      delfb <- fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) cat("Trouble - restarting: \n")
        restarts <- restarts + 1
        orth <- 1
        diams <- min(diam)
        sx <- sign(0.5 * sign(sgrad)) #TODO: isn't that just `sign(sgrad)` ?????
        happy <- 0L
        V[, -1] <- V[, 1]
        diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
      }
    }
    # New point accepted; remove old point and restart
    if (happy == 1) {
      V[, n + 1] <- xnew
      f[n + 1] <- fnew
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    } else if (happy == 0 & restarts < restarts.max) {
      if (orth == 0) orth <- 1
      V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
      for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
      nf <- nf + n
      # Sort
      ord <- order(f) #TODO: enhance for speed and tie-breaks
      V <- V[, ord]
      f <- f[ord]
    }
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
    f[is.nan(f)] <- Inf
    dist <- f[n + 1] - f[1]
    sgrad <- c(solve(t(v), delf))
    
    if (trace & !(itc%%2)) 
      cat("iter: ", itc, "\n", "value: ", f[1], "\n")
  } #while
  
  # Exit code and message
  if (dist <= ftol | simplex.size <= 1e-06) {
    conv <- 0
    message <- "Successful convergence"
  } else if (nf >= maxfeval) {
    conv <- 1
    message <- "Maximum number of fevals exceeded"
  } else if (restarts >= restarts.max) {
    conv <- 2
    message <- "Stagnation in Nelder-Mead"
  } else {
    conv <- 99
    message <- "Unexpected termination"
  }
  return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, feval = nf, 
              restarts = restarts, convergence = conv, message = message))
}




rosbkext <- function(x) {
  # Extended Rosenbrock function
  n <- length(x)
  sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}

np <- 10
#set.seed(123)
p0 <- rnorm(np)
p0[p0 > +2] <- +2 - 1E-8
p0[p0 < -2] <- -2 + 1E-8
e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2)
# o <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2)
# 
 ctrl <- list(maxfeval = 50000)
 #set.seed(123)
 p0 <- rnorm(np)
 p0[p0 > +2] <- +2 - 1E-8
 p0[p0 < -2] <- -2 + 1E-8
 e <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)
# o <- nmkb(fn = rosbkext, par = p0, lower = -2, upper = 2, control = ctrl)



