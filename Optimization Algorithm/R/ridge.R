#' Ridge Regression
#' 
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values.
#' 
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
ridge_regression <- function(y, X, lambda) {
  n = dim(X)
  if (length(y) != n[1]) stop("X and y are not conformable")
  if (lambda < 0) stop("The tuning parameter is negative")
  return(solve(t(X)%*%X+lambda*diag(n[2]))%*%t(X)%*%y)
}
