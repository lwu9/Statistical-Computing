#' Inverse Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
isweep <- function(A, k=NULL) {
  if (typeof(A)=="character")
    return("error")
  if(length(k)==0)
    k=1:dim(A)[1]
  for(i in k) {
    iswp <- isweep_k(A,i)
    if(typeof(iswp)=="character")
      return("error")
    else
      A <- iswp
  }
  return(A)
}