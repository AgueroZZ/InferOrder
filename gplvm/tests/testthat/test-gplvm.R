test_that("bgplvm returns correct class and structure", {
  set.seed(42)
  Y   <- matrix(rnorm(30 * 5), 30, 5)
  fit <- bgplvm(Y, Q = 2, M = 8, n_iter = 5, verbose = 0)

  expect_s3_class(fit, "bgplvm")
  expect_equal(dim(fit$mu_X), c(30, 2))
  expect_equal(dim(fit$S_X),  c(30, 2))
  expect_equal(dim(fit$Z),    c(8, 2))
  expect_length(fit$lengthscales, 2)
  expect_true(is.finite(fit$variance))
  expect_true(is.finite(fit$sigma2))
  expect_true(all(fit$S_X > 0))
})

test_that("bgplvm ELBO trace has correct length", {
  set.seed(43)
  Y   <- matrix(rnorm(20 * 4), 20, 4)
  fit <- bgplvm(Y, Q = 1, M = 5, n_iter = 10, verbose = 0, patience = 999L)
  expect_length(fit$elbo_trace, 10)
  expect_true(all(is.finite(fit$elbo_trace)))
})

test_that("extract_pseudotime works for Q=1 fit", {
  set.seed(44)
  Y   <- matrix(rnorm(40 * 6), 40, 6)
  fit <- bgplvm(Y, Q = 1, M = 8, n_iter = 5, verbose = 0)
  pt  <- extract_pseudotime(fit, scale01 = TRUE)

  expect_s3_class(pt, "data.frame")
  expect_equal(nrow(pt), 40)
  expect_named(pt, c("cell", "pseudotime", "lower95", "upper95", "s"))
  expect_true(all(is.finite(pt$pseudotime)))
  expect_true(all(pt$s >= 0))
  expect_true(min(pt$pseudotime) >= 0 - 1e-10)
  expect_true(max(pt$pseudotime) <= 1 + 1e-10)
})

test_that("extract_pseudotime works for Q=2 fit (auto-selects dim)", {
  set.seed(45)
  Y   <- matrix(rnorm(25 * 5), 25, 5)
  fit <- bgplvm(Y, Q = 2, M = 6, n_iter = 5, verbose = 0)
  pt  <- extract_pseudotime(fit)
  expect_equal(nrow(pt), 25)
})

test_that("print.bgplvm runs without error", {
  set.seed(46)
  Y   <- matrix(rnorm(20 * 4), 20, 4)
  fit <- bgplvm(Y, Q = 2, M = 5, n_iter = 3, verbose = 0)
  expect_output(print(fit), "Bayesian GPLVM")
})

test_that("plot.bgplvm runs without error for Q>=2", {
  set.seed(47)
  Y   <- matrix(rnorm(20 * 4), 20, 4)
  fit <- bgplvm(Y, Q = 2, M = 5, n_iter = 3, verbose = 0)
  expect_no_error(plot(fit))
})

test_that("plot.bgplvm runs without error for Q=1", {
  set.seed(48)
  Y   <- matrix(rnorm(20 * 4), 20, 4)
  fit <- bgplvm(Y, Q = 1, M = 5, n_iter = 3, verbose = 0)
  expect_no_error(plot(fit))
})

test_that("bgplvm respects M parameter", {
  set.seed(49)
  Y   <- matrix(rnorm(50 * 5), 50, 5)
  fit <- bgplvm(Y, Q = 2, M = 12, n_iter = 3, verbose = 0)
  expect_equal(fit$M, 12)
  expect_equal(nrow(fit$Z), 12)
})

test_that("check_convergence runs without error", {
  set.seed(50)
  Y   <- matrix(rnorm(20 * 4), 20, 4)
  fit <- bgplvm(Y, Q = 2, M = 5, n_iter = 10, verbose = 0)
  expect_output(check_convergence(fit))
})
