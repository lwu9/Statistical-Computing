#' Sweep k
#'
#' \code{sweep_k} applies the sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
sweep_k <- function(A, k) {
  if (A[k,k]==0 || sum((t(A)-A)^2) > 10^(-10) )
    return("error")
  A_swp <- A
  len <- dim(A)[1]
  A_swp <- A-matrix(A[,k],len,1)%*%matrix(A[k,],1,len)/A[k,k]
  A_swp[,k] <- A[,k]/A[k,k]
  A_swp[k,] <- A[k,]/A[k,k]
  A_swp[k,k] <- -1/A[k,k]
  return(A_swp)
}