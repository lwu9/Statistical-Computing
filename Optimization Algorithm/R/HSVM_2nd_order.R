#' Compute quasi-Newton Step (Naive) for HSVM
#'
#' @param X Design matrix
#' @param y Binary response vector
#' @param b Regression vector
#' @param g Gradient vector
#' @param delta HSVM parameter
#' @param lambda Regularization parameter
newton_step_naive <- function(X, y, b, g, delta, lambda=1e-4) {
  n = length(y); p = dim(X)[2]
  tt = X%*%b*y
  diagline <- array(0,n)
  ind = ((tt>(1-delta)) & (tt<=1))
  diagline[ind] <- 1/delta
  W <- diag(diagline)
  # Hessian = lambda*diag(p)
  # for (i in 1:n) {
  #   if (((y[i]*X[i,]%*%b)<=1) && ((y[i]*X[i,]%*%b)>(1-delta)))
  #     Hessian = Hessian + as.matrix(X[i,])%*%X[i,]/delta
  # }
  # t(chol(X%*%t(X)))%*%chol(X%*%t(X)) (lower %*% uppper triangular) is X%*%t(X)
  R = chol( lambda*diag(p) + t(X)%*%W%*%X)
  yy = forwardsolve(t(R), g)
  return(backsolve(R, yy))
}

#' Compute quasi-Newton Step (Sherman-Morrison-Woodbury) for HSVM
#'
#' @param X Design matrix
#' @param y Binary response vector
#' @param b Regression vector
#' @param g Gradient vector
#' @param delta HSVM parameter
#' @param lambda Regularization parameter
newton_step_smw <- function(X, y, b, g, delta, lambda=1e-4) {
  p = dim(X)[2]
  tt = X%*%b*y
  ind = ((tt>(1-delta)) & (tt<=1))
  X_tilde <- X[ind,]
  cma <- chol(delta*diag(sum(ind))+X_tilde%*%t(X_tilde)/lambda)
  cma_inv <- chol2inv(cma)
  Hessian_inv <- (diag(p)/lambda-t(X_tilde)%*%cma_inv%*%X_tilde/lambda^2)
  return(Hessian_inv%*%g)
}

#' Backtracking for steepest descent
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param g Gradient vector
#' @param d descent direction vector
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
backtrack_descent <- function(fx, x, g, d, tt=1, alpha=0.5, beta=0.9) {
  while (fx(x-tt*d) > fx(x)-alpha*tt*t(g)%*%d) {
    tt <- beta*tt
  }
  return(tt)
}

#' Damped Newton's Method for Fitting HSVM
#'
#' @param y Binary response
#' @param X Design matrix
#' @param b Initial regression coefficient vector
#' @param delta HSVM parameter
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param naive Whether to use Newton Naive or Newton Sherman-Morrison-Woodbury method
#' @export
hsvm_newton <- function(X, y, b, delta, lambda=1e-4, max_iter=1e2, tol=1e-3, naive = 1) {
  # relchg_ite <- c()
  # relchg_fun <- c()
  # obj_fun <- c()
  # norm_grad <- c()
  fx <- function(b) return(fx_hsvm(y, X, b, delta, lambda))
  for (i in 1:max_iter) {
    g = gradf_hsvm(y, X, b, delta, lambda)
    if(naive == 0) d = newton_step_smw(X, y, b, g, delta, lambda=1e-4)
    else d = newton_step_naive(X, y, b, g, delta, lambda=1e-4)
    tt = backtrack_descent(fx, b, g, d, tt=1, alpha=0.5, beta=0.9)
    d_ite =  tt*d
    # relchg_ite <- c(relchg_ite, norm(d_ite,"2")/norm(b,"2"))
    beta = b - d_ite
    # f1 <- fx_hsvm(y, X, beta, delta, lambda)
    # f0 <- fx_hsvm(y, X, b, delta, lambda)
    # d_fun <- f1-f0
    decrement_square <- t(g)%*%d
    # obj_fun[i] <- fx_hsvm(y, X, beta, delta, lambda)
    # norm_grad[i] <- norm(gradf_hsvm(y, X, beta, delta, lambda),"2")
    # relchg_fun <- c(relchg_fun,norm(d_fun,"2")/norm(f0,"2"))
    # result = list(beta,obj_fun,norm_grad,
    #               relchg_fun,relchg_ite, decrement_square/2,i)
    #print(decrement_square)
    #print(paste0(i,":",decrement_square))
    if (decrement_square/2 < tol) return(beta)
    else b = beta
  }
}

