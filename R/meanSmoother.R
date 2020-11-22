#' Windowed mean smoother
#'
#' TODO: farm out to C++ and test, at least when there are no NAs
#' 
#' @name meanSmoother
#'
#' @param x     Input data matrix: samples are columns, regions/loci are rows
#' @param k     Number of windows to use (default k=1, i.e., 3 windows)
#' @param iter  Number of iterations to smooth (default is 2)
#' @param na.rm Whether to remove NAs prior to smoothing (TRUE)
#' @param delta Convergence threshhold (overrides iter if > 0; default is 0) 
#'
#' @return      Smoothed data matrix
#'
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(10000), ncol=25)
#' smooth.dummy <- meanSmoother(dummy)
#' smooth.dummy <- meanSmoother(dummy, iter=3) 
#' smooth.dummy <- meanSmoother(dummy, delta=1e-3) 
#' 
meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE, delta=0) {

  if (k == 0) {
    message("Returning unsmoothed x. This is probably an error.") 
    return(x) 
  } 

  stopifnot(length(x) >= k)

  i <- 0 
  eps <- delta + 1 
  while (i < iter & eps > delta) {

    x0 <- x 
    i <- i + 1 
    if (!anyNA(x0)) {
      x <- .meanSmoother.rcpp(x0, k=k)
    } else {
      x <- .meanSmoother.internal(x0, k=k, na.rm=na.rm)
    }
    eps <- median(abs(x - x0)) # R builtin is fastish

  }
  
  return(x)

}


# helper fn
.meanSmoother.internal <- function(x, k, na.rm=TRUE){

  n <- length(x) 
  y <- rep(NA, n)

  first <- k + 1                # first eligible position to smooth
  last <- n - k                 # last eligible position to smooth
  excess <- seq((last + 1), n)  # excess bins beyond eligible

  # why, it even looks like C++ now. note that na.rm can create issues
  for (i in first:last) y[i] <- .window.mean(x, i=i, k=k, na.rm=na.rm)
  
  # it is possible for k to be 0 in this loop, it appears
  for (i in 1:k) y[i] <- .window.mean(x, i=i, k=i-1, na.rm=na.rm)

  # it is definitely possible for k to be 0 in this loop
  for (i in excess) y[i] <- .window.mean(x, i=i, k=n-i, na.rm=na.rm)

  return(y)

}


# helper fn
.window.mean <- function(x, j, k, na.rm=na.rm) {

  ifelse(k > 0, mean(x[(j-k-1):(j+k)], na.rm=na.rm), x[j])
  # startpos = j - k - 1
  # endpos = j + k
  # span = k + 2 

}


# helper fn; farm out to Rcpp
.meanSmoother.rcpp <- function(v, k) { 

  stop(".meanSmoother.rcpp is not yet tested")
  
  n <- length(v) 
  y <- rep(NA, n)               # has to be a better way to do this 
  first <- k + 1                # first eligible position to smooth
  last <- n - k                 # last eligible position to smooth
  excess <- seq((last + 1), n)  # excess bins beyond eligible

  # k should never be 0 in the following loop:
  for (i in first:last) y[i] <- .window.mean.rcpp(v, i=i, k=k)
  
  # it is possible for k to be 0 in this loop:
  for (i in 1:k) y[i] <- .window.mean.rcpp(v, i=i, k=i-1)

  # it is definitely possible for k to be 0 in this loop
  for (i in excess) y[i] <- .window.mean.rcpp(v, i=i, k=n-i)

  return(y)

}


# helper fn; farm out to Rcpp
.window.mean.rcpp <- function(v, j, k, na.rm=na.rm) {

  stop(".window.mean.rcpp is not yet tested")
  ifelse(k > 0, mean(v[(j-k-1):(j+k)]), v[j])

}
