test_that("kernel_se returns symmetric positive definite matrix", {
  X <- matrix(rnorm(10 * 2), 10, 2)
  K <- kernel_se(X)
  expect_equal(dim(K), c(10, 10))
  expect_true(isSymmetric(K))
  expect_true(all(eigen(K)$values > -1e-8))
})

test_that("kernel_se_ard returns correct dimensions", {
  X1 <- matrix(rnorm(5 * 3), 5, 3)
  X2 <- matrix(rnorm(7 * 3), 7, 3)
  K  <- kernel_se_ard(X1, X2, lengthscales = c(1, 1, 1))
  expect_equal(dim(K), c(5, 7))
})
