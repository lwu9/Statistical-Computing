#' Check Matrix is k-banded
#' 
#' \code{check_k_banded} returns a Boolean variable indicating whether or
#' not a matrix is k-banded.
#' 
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
check_k_banded <- function(A, k) {
  if (k < 0 )
    stop("k is negative")
  if (round(k)-k != 0)
    stop("k is not an integer")
  n = dim(A)[1]
  zero_row <- sapply((k+2):n, function(i) return(sum(abs(A[i,1:(i-k-1)]))))
  nonzero <- sum(sapply((k+2):n, function(i) return(abs(A[i,(i-k)]))))
  if (isSymmetric(A) && nonzero!=0 && sum(zero_row) == 0)
    return(TRUE)
  else return(FALSE)
}

#' Banded Cholesky
#' 
#' \code{chol_banded} computes the Cholesky decomposition of a banded 
#' positive definite matrix.
#' 
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
chol_banded <- function(A, k, check=TRUE) {
  
  p = 1
  if (check==TRUE) {
    #A <- Matrix(A,sparse = T)
    if (check_k_banded(A,k)==TRUE) {
      m = dim(A)[1]
      R <- matrix(0,m,m)
      for (i in 1:(m-k)) {
        n = dim(A)[1]
        r = matrix(0,n,n)
        #p = 1
        r[p,p] <- sqrt(A[p,p])
        r[p,(p+1):(p+k)] <- A[p,(p+1):(p+k)]/r[p,p]
        #length(r[p,(p+1):n])
        A <- A[(p+1):n,(p+1):n]-matrix(r[p,(p+1):n],n-p,1)%*%r[p,(p+1):n]
        R[i,i:(n+i-1)] <- r[p,]
      }
      for (i in (m-k+1):(m-1)) {
        n = dim(A)[1]
        r = matrix(0,n,n)
        #p = 1
        r[p,p] <- sqrt(A[p,p])
        r[p,(p+1):n] <- A[p,(p+1):n]/r[p,p]
        A <- A[(p+1):n,(p+1):n]-matrix(r[p,(p+1):n],n-p,1)%*%r[p,(p+1):n]
        R[i,i:(n+i-1)] <- r[p,]
      }
      R[m,m] <- sqrt(A)
      R <- Matrix(R, sparse=TRUE)
      return(R)
    }
    else if (!isSymmetric(A))
      stop(paste0("A is not symmetric, A is not ", k,"-banded"))
    else stop(paste0("A is not ", k,"-banded"))
  }
  else {
    A <- Matrix(A,sparse = T)
    m = dim(A)[1]
    R <- matrix(0,m,m)
    for (i in 1:(m-k)) {
      n = dim(A)[1]
      r = matrix(0,1,n)
      #p = 1
      r[p,p] <- sqrt(A[p,p])
      r[p,(p+1):(p+k)] <- A[p,(p+1):(p+k)]/r[p,p]
      #length(r[p,(p+1):n])
      A <- A[(p+1):n,(p+1):n]-
        Matrix(matrix(r[p,(p+1):n],n-p,1)%*%r[p,(p+1):n],sparse=TRUE)
      R[i,i:(n+i-1)] <- r[p,]
    }
    for (i in (m-k+1):(m-1)) {
      n = dim(A)[1]
      r = matrix(0,n,n)
      #p = 1
      r[p,p] <- sqrt(A[p,p])
      r[p,(p+1):n] <- A[p,(p+1):n]/r[p,p]
      A <- A[(p+1):n,(p+1):n]-matrix(r[p,(p+1):n],n-p,1)%*%r[p,(p+1):n]
      R[i,i:(n+i-1)] <- r[p,]
    }
    R[m,m] <- sqrt(A)
    R <- Matrix(R, sparse=TRUE)
    return(R)
  }
}


