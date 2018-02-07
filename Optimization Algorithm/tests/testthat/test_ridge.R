#library(testhat)
context("ridge regression")
test_that("ridge_regression throws an error if k is negative or X, y non-conformable",{
  X <- matrix(1:6,2,3)
  y <- 1:4
  expect_error( ridge_regression(y,X,lambda=0.3), "X and y are not conformable" )
  y <- 1:2
  expect_error( ridge_regression(y,X,lambda=-0.3),"The tuning parameter is negative" )
})
test_that("check the correctness of the estimated regression coefficients",{
  X <- matrix(1:6,2,3)
  n <- dim(X)
  y <- 1:2
  lambda <- 0.3
  bhat <- ridge_regression(y,X,lambda)
  expect_equal(norm((t(X)%*%X+lambda*diag(n[2]))%*%bhat-t(X)%*%y),
               expected=0.0,tolerance=1e-5)
})
