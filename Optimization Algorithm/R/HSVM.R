#' Objective Function for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
fx_hsvm <- function(y, X, beta, delta=0.4, lambda=1e-4) {
  n = dim(X)[1]
  tt = X%*%beta*y
  ind2 = ((tt>(1-delta)) & (tt<=1))
  phi2 <- sum((1-tt[ind2])^2)/2/delta
  phi3 <- sum(1-tt[tt<(1-delta)]-delta/2)
  phi <- phi2+phi3
  # phi = c()
  # for(i in 1:n) {
  #   tt <- y[i]*X[i,]%*%beta
  #   if(tt > 1) phi[i] = 0
  #   else if (tt<=1 & tt>(1-delta))
  #     phi[i] = (1-tt)^2/2/delta
  #   else phi[i] = 1-tt-delta/2
  # }
  penal  <- lambda/2*norm(beta,"2")^2
  return(sum(phi)+penal)
}

#' Gradient for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
gradf_hsvm <- function(y, X, beta, delta=0.4, lambda=1e-4) {
  n = dim(X)[1]
  tt = X%*%beta*y
  ind <- (tt<=1 & tt>(1-delta))
  if (sum(ind)==1) grad <- as.matrix(X[ind,])%*%(X[ind,]%*%beta-y[ind])/delta
  else grad <- t(X[ind,])%*%(X[ind,]%*%beta-y[ind])/delta
  ind2 <- (tt<=(1-delta))
  if (sum(ind2)==1) grad = grad - as.matrix(X[ind2,])%*%y[ind2]
  else grad = grad - t(X[ind2,])%*%y[ind2]

  # grad = beta*0
  # for(i in 1:n) {
  #   tt <- y[i]*X[i,]%*%beta
  #   if (tt<=1 & tt>(1-delta))
  #     grad <- grad+X[i,]*(as.numeric(y[i]^2*(X[i,]%*%beta)-y[i])/delta)
  #   else if(tt<=1-delta) grad = grad-X[i,]*y[i]
  # }
  penal  <- lambda*beta
  return(grad+penal)
}

#' Gradient Descent (Fixed Step-Size)
#'
#' @param y binary response
#' @param X design matrix
#' @param b0 initial parameter estimate
#' @param delta shape parameter
#' @param t step-size
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(y, X, b0, delta, tt, lambda=1e-4, max_iter=1e2, tol=1e-3) {
  relchg_ite <- c()
  relchg_fun <- c()
  obj_fun <- c()
  norm_grad <- c()
  for (i in 1:max_iter) {
    # relative change in the iterate values
    d_ite =  tt*gradf_hsvm(y, X, b0, delta, lambda)
    relchg_ite <- c(relchg_ite,norm(d_ite,"2")/norm(b0,"2"))
    beta = b0 - d_ite
    f1 <- fx_hsvm(y, X, beta, delta, lambda)
    f0 <- fx_hsvm(y, X, b0, delta, lambda)
    d_fun <- f1-f0
    obj_fun[i] <- fx_hsvm(y, X, beta, delta, lambda)
    norm_grad[i] <- norm(gradf_hsvm(y, X, beta, delta, lambda),"2")
    relchg_fun <- c(relchg_fun,norm(d_fun,"2")/norm(f0,"2"))
    result = list(beta,obj_fun,norm_grad,
                  relchg_fun,relchg_ite)
    if (abs(d_fun)<tol) {
      return(result)
    }
    else {b0 = beta}
  }
  return(result)
}

#' Gradient Descent (Backtracking Step-Size)
#'
#' @param y binary response
#' @param X design matrix
#' @param b0 initial parameter estimate
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_backtrack <- function(y, X, b0, delta, lambda=1e-4, max_iter=1e4, tol=1e-3) {
  relchg_ite <- c()
  relchg_fun <- c()
  obj_fun <- c()
  norm_grad <- c()
  t_size = 0.5; beta_size= 0.5
  for (i in 1:max_iter) {
    # relative change in the iterate values
    d_ite = t_size*gradf_hsvm(y, X, b0, delta, lambda)
    beta = b0 - d_ite
    f1 <- fx_hsvm(y, X, beta, delta, lambda)
    f0 <- fx_hsvm(y, X, b0, delta, lambda)
    d_fun <- f1-f0
    obj_fun[i] <- fx_hsvm(y, X, beta, delta, lambda)
    norm_grad[i] <- norm(gradf_hsvm(y, X, beta, delta, lambda),"2")
    # If loss fun is increasing, then use backtracking
    if (d_fun > 0) {
      t_size = beta_size*t_size
      d_ite = - t_size*gradf_hsvm(y, X, b0, delta, lambda)
      beta = b0 + d_ite
      f1 <- fx_hsvm(y, X, beta, delta, lambda)
      f0 <- fx_hsvm(y, X, b0, delta, lambda)
      while(d_fun > -0.5*t_size*norm(gradf_hsvm(y, X, b0, delta, lambda),"2")) {
        t_size = beta_size*t_size
        d_ite = - t_size*gradf_hsvm(y, X, b0, delta, lambda)
        beta = b0 + d_ite
        f1 <- fx_hsvm(y, X, beta, delta, lambda)
        f0 <- fx_hsvm(y, X, b0, delta, lambda)
        d_fun <- f1-f0
      }
    }
    relchg_ite <- c(relchg_ite,norm(d_ite,"2")/norm(b0,"2"))
    relchg_fun <- c(relchg_fun,norm(d_fun,"2")/norm(f0,"2"))
    result = list(beta,obj_fun,norm_grad,
                  relchg_fun,relchg_ite)
    if (abs(d_fun)<tol) {

      return(result)
    }
    else b0 = beta
  }
  return(result)
}

#' Gradient Descent
#'
#' @param y binary response
#' @param X design matrix
#' @param b0 initial parameter estimate
#' @param delta shape parameter
#' @param t step-size
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent <- function(y, X, b0, delta, t=NULL, lambda=1e-4, max_iter=1e2, tol=1e-3) {
  if (!is.null(t)) {
    L=sapply(1:length(y), function(j) {y[j]^2*norm(X[j,],"2")**2})
    L=sum(L)/delta+abs(lambda)
    t = 1/2/L
    return(gradient_descent_fixed(y, X, b0, delta, t, lambda, max_iter, tol))
  }
  else
    return(gradient_descent_backtrack(y, X, b0, delta, lambda, max_iter, tol))
}
