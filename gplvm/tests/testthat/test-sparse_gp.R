test_that("psi1 reduces to kernel_se_ard when posterior variance -> 0", {
  set.seed(1)
  N <- 5; M <- 3; Q <- 2
  mu_X <- matrix(rnorm(N * Q), N, Q)
  Z    <- matrix(rnorm(M * Q), M, Q)
  ls   <- c(1.0, 0.5)
  var  <- 2.0
  # Near-zero posterior variance
  s2_X <- matrix(1e-10, N, Q)

  psi1   <- .expected_K_NM(mu_X, s2_X, Z, ls, var)
  K_exact <- kernel_se_ard(mu_X, Z, lengthscales = ls, variance = var)

  expect_equal(psi1, K_exact, tolerance = 1e-6)
})

test_that("psi1 has correct dimensions", {
  N <- 10; M <- 4; Q <- 3
  mu_X <- matrix(0, N, Q)
  s2_X <- matrix(0.1, N, Q)
  Z    <- matrix(0, M, Q)
  res  <- .expected_K_NM(mu_X, s2_X, Z, rep(1, Q), 1.0)
  expect_equal(dim(res), c(N, M))
})

test_that("psi1 entries are positive and bounded by variance", {
  set.seed(2)
  N <- 8; M <- 3; Q <- 2
  mu_X <- matrix(rnorm(N * Q), N, Q)
  s2_X <- matrix(runif(N * Q, 0.01, 1), N, Q)
  Z    <- matrix(rnorm(M * Q), M, Q)
  ls   <- c(1.0, 1.0); var <- 1.5
  psi1 <- .expected_K_NM(mu_X, s2_X, Z, ls, var)
  expect_true(all(psi1 > 0))
  expect_true(all(psi1 <= var + 1e-10))
})

test_that("psi2 is symmetric and positive semidefinite", {
  set.seed(3)
  N <- 6; M <- 3; Q <- 2
  mu_X <- matrix(rnorm(N * Q), N, Q)
  s2_X <- matrix(0.1, N, Q)
  Z    <- matrix(rnorm(M * Q), M, Q)
  psi2 <- .psi2_matrix(mu_X, s2_X, Z, c(1, 1), 1.0)
  expect_equal(dim(psi2), c(M, M))
  expect_true(isSymmetric(psi2, tol = 1e-8))
  expect_true(all(eigen(psi2, only.values = TRUE)$values >= -1e-8))
})

test_that("expected_trace_K_NN returns N * variance", {
  expect_equal(.expected_trace_K_NN(100L, 2.5), 250)
  expect_equal(.expected_trace_K_NN(1L, 1.0), 1.0)
})

test_that(".compute_K_MM returns valid Cholesky", {
  set.seed(4)
  Z  <- matrix(rnorm(10 * 2), 10, 2)
  res <- .compute_K_MM(Z, c(1.0, 1.0), 1.0)
  expect_equal(dim(res$K), c(10, 10))
  expect_true(isSymmetric(res$K))
  expect_true(all(eigen(res$K, only.values = TRUE)$values > 0))
  # Verify L L^T ≈ K
  expect_equal(res$L %*% t(res$L), res$K, tolerance = 1e-10)
})

test_that(".chol_with_jitter handles near-singular matrix", {
  K <- matrix(1, 5, 5)  # rank-1, needs jitter
  expect_no_error({
    res <- .chol_with_jitter(K)
    expect_true(!is.null(res$L))
  })
})
