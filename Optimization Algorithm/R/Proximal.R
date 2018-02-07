#' Softthreshold operator
#'
#' @param x input vector
#' @param lambda threshold
#' @export
#' @examples
#' softthreshold(seq(-10, 10, length.out=10), 1)
softthreshold <- function(x, lambda) {
  st <- numeric(length=length(x))
  st[which(x > lambda)] <- x[which(x > lambda)] - lambda
  st[which(x < -lambda)] <- x[which(x < -lambda)] + lambda
  return(st)
}

#' Proximal mapping of the scaled nuclear norm
#'
#' @param X input matrix
#' @param gamma scale parameter
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' Y <- prox_nuc(X, 1)
prox_nuc <- function(X, gamma) {
  SVD <- svd(X)
  return(SVD$u%*%diag(softthreshold(SVD$d, gamma))%*%t(SVD$v))
}

#' Project matrix onto observed index set
#'
#' @param X input matrix
#' @param Omega index set
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' PX <- project_omega(X, Omega)
project_omega <- function(X, Omega) {
  X[-Omega] = 0
  return(X)
}

#' Proximal-Gradient Step
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param t step size
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' t <- 1
#' U1 <- prox_step(X, U, Omega, lambda, tt)
prox_step <- function(X, U, Omega, lambda, tt) {
  x_half <- U-tt*(project_omega(U, Omega) - project_omega(X, Omega))
  return(proj_nuc(x_half, lambda))
}

#' Matrix Completion Loss
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' matrix_completion_loss(X, U, Omega, lambda)
matrix_completion_loss <- function(X, U, Omega, lambda) {
 nuclear <- sum(svd(U)$d)
 return(sum((project_omega(X, Omega)-project_omega(U, Omega))^2)/2+lambda*nuclear)
}

#' Proximal Gradient Descent
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' pg_sol <- prox_grad(X, U, Omega, lambda)
prox_grad <- function(X, U, Omega, lambda, t, max_iter=1e2, tol=1e-3) {
  U_new <- prox_step(X, U, Omega, lambda, t)
  for (i in 1:max_iter) {
    loss <- matrix_completion_loss(X, U, Omega, lambda)
    loss_new = matrix_completion_loss(X, U_new, Omega, lambda)
    obj_chg <- abs((loss_new - loss)/loss)
    diff_iter <- norm(U_new-U,"F")/norm(U,"F")
    if (obj_chg < tol) {
      return(list(iter_val=U_new, obj_val=loss_new,
                  obj_chg = obj_chg, iter_chg=diff_iter))
    }
    U <- U_new
    U_new <- prox_step(X, U, Omega, lambda, t)
  }
  return(list(iter_val=U_new, obj_val=loss_new,
              obj_chg = obj_chg, iter_chg=diff_iter))
}

#' Soft-Impute Algorithm
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' si_sol <- soft_impute(X, U, Omega, lambda)
soft_impute <- function(X, U, Omega, lambda, max_iter=1e2, tol=1e-3) {
  for (i in 1:max_iter) {
    U_old <- U
    Y <- U
    obj_old <- matrix_completion_loss(X, U, Omega, lambda)
    Y[Omega] <- X[Omega]
    SVD <- svd(Y)
    D_tilde <- softthreshold(diag(SVD$d), lambda)
    U <- SVD$u%*%diag(D_tilde)%*%t(SVD$v)
    obj_new <- matrix_completion_loss(X, U, Omega, lambda)
    obj_chg <- abs((obj_old - obj_new)/obj_old)
    rel_chg_ier <- norm(U-U_old, "F")/norm(U_old, "F")
    if (rel_chg < tol) return(list(iter_val=U, obj = obj_new, obj_chg =obj_chg, iter_chg=rel_chg_ier))
  }
  return(list(iter_val=U, obj = obj_new, obj_chg =obj_chg, iter_chg=rel_chg_ier))
}

imageA <- as.matrix(read.table('/Users/lili/Desktop/NCSU/ST758/HW/imageA.csv', sep=",", header=FALSE))
image(imageA, col=gray(0:3/3))
n <- dim(imageA)[1]; p <- dim(imageA)[2]
Omega <- 1:(n*p)
Omega <- Omega[-Omega[is.na(imageA)] ]
lambda_max <- max(svd(project_omega(imageA, Omega))$d)
lambdas <- exp(seq(-6, log(lambda_max,10), length.out = 5))
U <- matrix(rnorm(n*p), n, p)
U[Omega] <- imageA[Omega]
si_sol <- list()
for (lambda  in lambdas) {
  si_sol[[1]] <- soft_impute(imageA, U, Omega, lambda)
  print(lambda)
}