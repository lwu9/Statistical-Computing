#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
sweep <- function(A, k=NULL) {
  if (typeof(A)=="character")
    return("error")
  else{
    if(length(k)==0)
      k=1:dim(A)[1]
    for(i in k) {
      swp <- sweep_k(A,i)
      if(typeof(swp)=="character")
        return("error")
      else
        A <- swp
    }
    return(A)
  }
}
