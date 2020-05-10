
.ReplaceVertex <- function(V, Vnew, At) {
  len <- ncol(V)
  if (At <= 1)
    Result <- cbind(Vnew, V[, 1L:(len-1L)], deparse.level = 0)
  else if (At >= len) {
    V[, len] <- Vnew
    Result <- V
  } else {
    Result <- cbind(V[, 1L:(At-1L)], Vnew, V[, At:(len-1)], deparse.level = 0)
  }
  
  return(Result)
}


#' nmkp
#' An implementation of the Nelder-Mead algorithm for derivative-free 
#' optimization. This version refactors, optimises, corrects and extends 
#' [dfoptim::nmk()][dfoptim::nmk()] and [dfoptim::nmkb()][dfoptim::nmkb()].
#' @usage nmkp(par, fn, lower=-Inf, upper=Inf, control = list(), ...)
#' @param par A starting vector of parameter values. Must be feasible, i.e. 
#' lie strictly between lower and upper bounds.
#' @param fn Nonlinear objective function that is to be optimized. A scalar 
#' function that takes a real vector as argument and returns a scalar that 
#' is the value of the function at that point (see details).
#' @param lower,upper Lower/upper bounds on the parameters. A vector of 
#' the same length as the parameters. If a single value is specified, it is 
#' assumed that the same lower/upper bound applies to all parameters.
#' @param control A list of control parameters. List entries must be named.
#' See *Details* for more information.
#' @param ... Additional arguments passed to `fn`.
#' @details 
#' If the arguments `lower` and `upper` are `+-Inf` this function runs an unconstrained 
#' algorithm. With bounds it inflates the space and runs a traditional simplex search 
#' in this inflated space. In other words, bounds are enforced by means of a parameter 
#' transformation.
#' 
#' ## Control Arguments
#' Argument control is a list specifing any changes to default values of 
#' algorithm control parameters for the outer loop. Note that the names of 
#' these must be specified completely. Partial matching will not work. 
#' The list items are as follows:
#' \describe{
#'   \item{tol}{Convergence tolerance. Iteration is terminated when the 
#'   absolute difference in function value between successive iteration 
#'   is below tol. Default is 1.e-06.}
#'   \item{maxfeval}{Maximum number of objective function evaluations 
#'   allowed (in other words every call of `fn`). Default is 
#'   `min(5000, max(1500, 20*length(par)^2))`.}
#'   \item{regsimp}{A logical variable indicating whether the starting 
#'   parameter configuration is a regular simplex. Default is TRUE.}
#'   \item{maximize}{A logical variable indicating whether the objective 
#'   function should be maximized. Default is FALSE.}
#'   \item{restarts.max}{Maximum number of times the algorithm should be 
#'   restarted before declaring failure. Default is 3.}
#'   \item{trace}{A logical variable indicating whether the starting 
#'   parameter configuration is a regular simplex. Default is FALSE.}
#'   \item{adapt.behavior}{The parameters controlling the behaviour of the
#'   simplex are chosen according to Gao & Han (2010) instead of the original
#'   choice by Nelder & Mead.}
#' }
#' @return A list with the following components:
#' \describe{
#'   \item{par}{Best estimate of the parameter vector found by the algorithm.}
#'   \item{value}{The value of the objective function at termination.}
#'   \item{feval}{The number of times the objective `fn` was evaluated.}
#'   \item{restarts}{The number of times the algorithm had to be restarted 
#'   when it stagnated.}
#'   \item{convergence}{An integer code indicating type of convergence. 
#'   0 indicates successful convergence. Positive integer codes indicate 
#'   failure to converge. The `message` will explain.}
#'   \item{message}{Text message indicating the type of convergence or failure.}
#' }
#' @source This function is a heavily refactored copy of dfoptim::nmkb(). This 
#' algorithm is based on the [Matlab code of Prof. C.T. Kelley, given in his book 
#' "Iterative methods for optimization"](https://archive.siam.org/books/kelley/fr18/OPT_CODE/nelder.m). 
#' It was implemented for R by Ravi Varadhan with permission of 
#' Prof. Kelley and SIAM. However, there are some non-trivial modifications of 
#' the algorithm made by Ravi Varadhan. Some refactoring and further modifications
#' have been done by Jan Seifert.
#' @references 
#' Gao, F., & Han, L. (2010). Implementing the Nelder-Mead simplex algorithm with adaptive parameters. Computational Optimization and Applications, 51, 259–277. https://doi.org/10.1007/s10589-010-9329-3
#' Kelley, C. T. (1999). Iterative Methods for Optimization. Society for Industrial and Applied Mathematics. https://doi.org/10.1137/1.9781611970920
#' Nelder, J. A., & Mead, R. (1965). A simplex method for function minimization. Computer Journal, 7, 308–313.
#' Lagarias, J., Reeds, J., Wright, M., & Wright, P. (1998). Convergence Properties of the Nelder–Mead Simplex Method in Low Dimensions. SIAM Journal on Optimization, 9, 112–147. https://doi.org/10.1137/S1052623496303470
nmkp <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
  # BASIC CONSTANTS THAT DETERMINE THE BEHAVIOUR OF THE ALGORITHM
  ctrl <- list(tol = 1e-06, maxfeval = min(5000, max(1500, 20 * length(par)^2)), 
               regsimp = TRUE, maximize = FALSE, restarts.max = 3, 
               adapt.behavior = TRUE, trace = FALSE)
  namc <- match.arg(names(control), choices = names(ctrl), several.ok = TRUE)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in 'control': ", namc[!(namc %in% names(ctrl))])
  if (!is.null(names(control))) ctrl[namc] <- control
  ftol <- ctrl$tol
  stol <- 1E-6
  maxfeval <- ctrl$maxfeval
  restarts.max <- ctrl$restarts.max
  maximize <- ifelse(ctrl$maximize, -1, 1)
  trace <- ctrl$trace
  # Setting oshrink = FALSE gives vanilla Nelder-Mead, i.e. the simplex is merely
  # reflected expanded, or contracted without shrinking
  oshrink <- TRUE
  # Dimensions of paramter space
  n <- length(par)
  if (n == 1) 
    stop(call. = FALSE, "Use `optimize` for univariate optimization")
  if (n > 32) 
    warning("Nelder-Mead should not be used for high-dimensional optimization")
  
  # Basic Nelder-Mead transformation parameters
  if (ctrl$adapt.behavior) {
    rho <- 1              # reflection  alpha
    chi <- 1 + 2/n        # expansion   beta
    gamma <- 0.75 - 1/2/n # contraction gamma
    sigma <- 1 - 1/n      # shrinkage   delta
  } else {
    rho <- 1
    chi <- 2
    gamma <- 0.5
    sigma <- 0.5
  }
  
  # SETUP
  # Spatial distortion functions to handle box constraints
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1]) / rangex[c1] - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + rangex[c1]/2 * (1 + tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  
  getsimplexsize <- function(V, DiffTo1 = NULL) {
    #Vold <- V # FOR TESTING ONLY
    #if(!all(c2)) V <- apply(V, MARGIN = 2, ginv)
    if (is.null(DiffTo1)) DiffTo1 <- V[, -1] - V[, 1]
    size <- sum(abs(DiffTo1)) / max(1, sum(abs(V[, 1])))
    return(size)
  }
  
  # Set variables required to handle box constraints
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  
  if (any(c(par < lower, upper < par))) stop("Infeasible starting values!", call.=FALSE)
  low.finite <- is.finite(lower)
  upp.finite <- is.finite(upper)
  c1 <- low.finite & upp.finite   # both lower and upper bounds are finite 
  c2 <- !(low.finite | upp.finite) # both lower and upper bounds are infinite
  c3 <- !(c1 | c2) & low.finite   # finite lower bound, infinite upper bound
  c4 <- !(c1 | c2) & upp.finite   # finite upper bound, infinite lower bound
  rangex <- (upper - lower)       # range needed for g and ginv
  
  # Set objective function
  if(all(c2)) {
    fnmb <- function(par, ...) { fn(par, ...) * maximize }
    x0 <- par # x0 is starting point without transformation
  }
  else {
    fnmb <- function(par, ...) { fn(ginv(par), ...) * maximize }
    x0 <- g(par) # x0 is the starting point
  }

  
  # RUN
  # Initial values: V - vertices of the simplex. f - result of objective function.
  V <- cbind(rep(0, n), diag(n))
  f <- rep(0, n + 1)
  f[1] <- fnmb(x0, ...)
  V[, 1] <- x0
  scale <- max(1, sqrt(sum(x0^2)))

  if (ctrl$regsimp) { # Create regular simplex
    alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, sqrt(n + 1) - 1)
    V[, -1] <- (x0 + alpha[2])
    diag(V[, -1]) <- x0[1:n] + alpha[1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
  } else {
    V[, -1] <- x0 + scale * V[, -1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
  }
  f[is.nan(f)] <- Inf
  nf <- n + 1L     # number of calls of objective function
  # Initial sorting of simplex vertices according to f(v)
  ord <- order(f)
  f <- f[ord]
  V <- V[, ord]
  # Start values
  v <- V[, -1] - V[, 1]
  delf <- f[-1] - f[1]          # "delta f"
  diam <- sqrt(colSums(v^2))    # diameter of the simplex???
  sgrad <- c(solve(t(v), delf))
  alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
  xbar <- rowMeans(V[, 1:n]) # initial centroid of the simplex
  
  # Termination criteria (`nf` has been defined above)
  ##/simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 1])))
  simplex.size <- getsimplexsize(V)
  dist <- f[n + 1] - f[1] # distance between f(highest) and f(lowest) vertex
  restarts <- 0L # number of restarts
  
  itc <- 0L # iteration counter
  while (nf < maxfeval && restarts < restarts.max && dist > ftol && simplex.size > stol) {
    happy <- FALSE  # `happy == 1` will indicate that new vertex is accepted
    itc <- itc + 1L
    nonshrinkit <- TRUE
    
    # REFLECT (is always first try)
    fbc <- mean(f)
    xr <- (1 + rho) * xbar - rho * V[, n + 1]
    fr <- fnmb(xr, ...)
    nf <- nf + 1L
    if (is.nan(fr)) fr <- Inf

    if (fr >= f[1] && fr < f[n]) {
      # Reflection was successful
      happy <- TRUE
      xnew <- xr
      fnew <- fr
    } else if (fr < f[1]) { # if reflected point f(r) < f(old): EXPAND
      xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 1]
      fe <- fnmb(xe, ...)
      if (is.nan(fe)) fe <- Inf
      nf <- nf + 1L
      if (fe < fr) {
        xnew <- xe
        fnew <- fe
        happy <- TRUE
      } else {
        xnew <- xr
        fnew <- fr
        happy <- TRUE
      }
    } else if (fr >= f[n] && fr < f[n + 1]) {
      # OUTSIDE CONTRACTION
      xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, n + 1]
      fc <- fnmb(xc, ...)
      if (is.nan(fc)) fc <- Inf
      nf <- nf + 1L
      if (fc <= fr) {
        xnew <- xc
        fnew <- fc
        happy <- TRUE
      }
    } else if (fr >= f[n + 1]) {
      # INSIDE CONTRACTION
      xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
      fc <- fnmb(xc, ...)
      if (is.nan(fc))
        fc <- Inf
      nf <- nf + 1L
      if (fc < f[n + 1]) {
        xnew <- xc
        fnew <- fc
        happy <- TRUE
      }
    }
    # Test for sufficient decrease; do oriented shrink if necessary
    if (happy && oshrink) {
      # Check if average function value has improved
      #/fbt <- mean(c(f[1:n], fnew))
      delfb <- mean(c(f[1:n], fnew)) - fbc #/fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) cat("Trouble - restarting: \n")
        restarts <- restarts + 1L
        diams <- min(diam)
        # According to the original matlab code the instruction below must be:
        # sx=.5+sign(sgrad); sx=sign(sx);
        sx <- sign(0.5 + sign(sgrad))
        happy <- FALSE
        V[, -1L] <- V[, 1L]
        diag(V[, -1L]) <- diag(V[, -1L]) - diams * sx[1L:n]
        nonshrinkit <- FALSE
      }
    }
    # New point accepted; remove old point and restart
    if (happy) {
      if(nonshrinkit) {
        # Insert new vertex and reorder
        InsertAt <- which.max(fnew < f)
        V <- .ReplaceVertex(V, xnew, At = InsertAt)
        f <- append(f[1:n], fnew, InsertAt-1)
        # Update centroid of the simplex
        if(InsertAt < n+1)
          xbar <- xbar + (xnew - V[, n+1])/n 
      } else {
        # Insert new vertex and reorder
        V[, n + 1] <- xnew
        f[n + 1] <- fnew
        ord <- order(f)
        V <- V[, ord]
        f <- f[ord]
        # Update centroid of the simplex
        xbar <- rowMeans(V[, 1:n])
      }
    } else if (!happy && restarts < restarts.max) {
      V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
      for (j in 2:ncol(V)) f[j] <- fnmb(V[, j], ...)
      nf <- nf + n
      ord <- order(f) # Order the vertices
      V <- V[, ord]
      f <- f[ord]
      xbar <- rowMeans(V[, 1:n]) # Centroid
    }
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    ##/simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1]))) #TODO: change simplex.size to rationale to Singer & Singer formula 16
    simplex.size <- getsimplexsize(V, v)
    f[is.nan(f)] <- Inf
    dist <- f[n + 1] - f[1]
    sgrad <- c(solve(t(v), delf))
    
    if (trace & !(itc%%2)) 
      cat("iter:", itc, "\n", "value:", f[1], "\n")
  } #while
  
  # Exit code and message
  if (dist <= ftol || simplex.size <= stol) {
    conv <- 0L
    message <- "Successful convergence"
    #if (simplex.size <= stol) print("simple") #TESTING
    #else print(simplex.size)
  } else if (nf >= maxfeval) {
    conv <- 1L
    message <- "Maximum number of fevals exceeded"
  } else if (restarts >= restarts.max) {
    conv <- 2L
    message <- "Stagnation in Nelder-Mead"
  } else {
    conv <- 99L
    message <- "Unexpected termination"
  }
  return(list(par = ifelse(c2, V[, 1], ginv(V[, 1])), value = f[1] * maximize, 
              feval = nf, restarts = restarts, 
              convergence = conv, message = message))
}

#### Testing code --------------
# rosbkext <- function(x) {
#   # Extended Rosenbrock function
#   n <- length(x)
#   sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
# }
# 
# np <- 6 #12
# for (box in c(2, 4, 12, 24, 32, 64, 128)) {
#   #set.seed(123)
#   p0 <- rnorm(np)
#   p0[p0 > +2] <- +2 - 1E-8
#   p0[p0 < -2] <- -2 + 1E-8
#   
#   #box <- c(rep(box, np/2), rep(Inf, np/2))
#   
#   ctrl <- list(maxfeval = 5E4, tol = 1E-12)
#   #o <- nmkp(fn = rosbkext, par = p0, lower = -2, upper = +box, control = ctrl)
#   o <- nmkp(fn = rosbkext, par = p0, control = ctrl)
#   print(o$message)
#   cat("f(", format(o$par, digits = 2), ") =", format(o$value, digits=3), "\n")
# }
