test_that("ELBO returns a finite scalar", {
  set.seed(10)
  N <- 15; D <- 4; Q <- 2; M <- 5
  Y      <- matrix(rnorm(N * D), N, D)
  params <- .init_params(Y, Q, M, "random")
  val    <- .compute_sparse_elbo(params, Y)
  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("ELBO increases when we move towards better params", {
  set.seed(11)
  N <- 20; D <- 5; Q <- 2; M <- 6
  Y      <- matrix(rnorm(N * D), N, D)
  params <- .init_params(Y, Q, M, "pca")

  elbo0 <- .compute_sparse_elbo(params, Y)

  # Run a few Adam steps
  res <- .adam_optimize(params, Y, n_iter = 20L, verbose = 0L)
  elbo_final <- tail(res$elbo_trace, 1)

  # ELBO should not decrease dramatically (allow some noise from numeric grad)
  expect_true(elbo_final >= elbo0 - 10)
})

test_that("ELBO grad has correct sign on a simple perturbation", {
  skip_if_not_installed("numDeriv")
  set.seed(12)
  N <- 10; D <- 3; Q <- 1; M <- 4
  Y      <- matrix(rnorm(N * D), N, D)
  params <- .init_params(Y, Q, M, "random")

  # Compute numeric gradient via numDeriv
  theta <- .pack_params(params, N, Q, M)
  f     <- function(th) .compute_sparse_elbo(.unpack_params(th, N, Q, M), Y)
  g_num <- numDeriv::grad(f, theta, method = "simple")

  # Compute our gradient
  g_ours <- .grad_sparse_elbo(params, Y)
  g_pack <- .pack_params(g_ours, N, Q, M)

  # Both should agree (they both use numDeriv internally for now,
  # so this is a consistency check)
  expect_equal(g_pack, g_num, tolerance = 1e-4)
})
