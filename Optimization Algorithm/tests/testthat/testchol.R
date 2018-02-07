context("K_banded Cholesky Decomposition")

test_that("check_k_banded throws an error if k is negative or is not an interger", {
  
  i <- c(1:10,1:9,1:8); j <- c(1:10,2:10,3:10); x <- c(1:10,5 * (2:10),3*(3:10))
  A = sparseMatrix(i, j, x = x)
  A <- as.matrix(round(10*A)/10)
  A = t(A)%*%A+0.001*diag(10)
  #A <- as.matrix(round(10*A)/10)
  expect_error( check_k_banded(A, -1),"k is negative" )
  expect_error( check_k_banded(A, 0.1),"k is not an integer" )
})

test_that("check_k_banded returns T or F to check whether A is k-banded", {

  i <- c(1:10,1:9,1:8); j <- c(1:10,2:10,3:10); x <- c(1:10,5 * (2:10),3*(3:10))
  A = sparseMatrix(i, j, x = x)
  A <- as.matrix(round(10*A)/10)
  A = t(A)%*%A+0.001*diag(10)
  expect_equal( check_k_banded(A, 2),TRUE )
  expect_equal( check_k_banded(A, 3),FALSE)
  A[1,2] <- A[1,2]-1
  expect_equal( check_k_banded(A, 2),FALSE)
})

test_that("chol_banded throws an error if k is negative or is not an interger", {
  # create a positive definite A
  i <- c(1:10,1:9,1:8); j <- c(1:10,2:10,3:10); x <- c(1:10,5 * (2:10),3*(3:10))
  A = sparseMatrix(i, j, x = x)
  A <- as.matrix(round(10*A)/10)
  A = t(A)%*%A+0.001*diag(10)
  expect_error( chol_banded(A, -1),"k is negative" )
  expect_error( chol_banded(A, 0.1),"k is not an integer" )
  expect_equal(sum(abs(chol_banded(A,2)-chol(A))),expected=0.0,tolerance=1e-7)
})


