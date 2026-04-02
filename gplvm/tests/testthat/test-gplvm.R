test_that("bgplvm runs and returns correct structure", {
  set.seed(42)
  Y <- matrix(rnorm(30 * 5), 30, 5)
  fit <- bgplvm(Y, Q = 2, n_iter = 10, verbose = 0)

  expect_s3_class(fit, "bgplvm")
  expect_equal(dim(fit$mu_X), c(30, 2))
  expect_equal(length(fit$elbo_trace), 10)
  expect_true(is.finite(tail(fit$elbo_trace, 1)))
})
