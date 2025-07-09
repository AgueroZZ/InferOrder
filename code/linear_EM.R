# internal functions for EM algorithm with linear model
solve_beta <- function(t_grid, Nk, xbar_k, betaprec) {
  K <- length(t_grid)
  d <- ncol(xbar_k)

  # scalar
  sum_N_t2 <- sum(Nk * t_grid^2)
  rhs <- colSums(matrix(Nk * t_grid, nrow=K, ncol=d) * xbar_k)

  # Solve
  lhs <- betaprec * diag(d) + sum_N_t2 * diag(d)
  beta_hat <- solve(lhs, rhs)
  return(beta_hat)
}

# EM algorithm for linear model
EM_algorithm_linear <- function(data, K, betaprec = 1,
                                t_grid   = NULL, seed = 1,
                                sigma_init = NULL,
                                max_iter = 100, tol = 1e-4,
                                verbose  = TRUE) {

  data_center <- colMeans(data)
  data <- sweep(data, 2, data_center)

  n <- nrow(data); d <- ncol(data)
  if (is.null(t_grid)) t_grid <- seq(-1, 1, length.out = K)

  set.seed(seed)
  beta       <- rnorm(d)
  pi_k       <- rep(1/K, K)
  sigma_list <- if (is.null(sigma_init)) lapply(1:K, function(k) diag(d)) else sigma_init

  loglik_trace <- numeric(max_iter)
  mu_list      <- vector("list", K)
  gamma        <- matrix(0, n, K)

  for (iter in 1:max_iter) {
    # 1) component means
    for (k in 1:K) mu_list[[k]] <- t_grid[k] * beta

    # 2) fast E-step must have invSigma/logdet in params
    params <- list(pi = pi_k, mu = mu_list, sigma = sigma_list)
    params <- init_cov_cache_fast(params)   # <-- cache inverses & logdets
    gamma  <- ESTEP(data, params)      # or just ESTEP(data, params)

    # 3) update pi
    Nk   <- colSums(gamma)
    pi_k <- Nk / n

    # 4) weighted means
    xbar_k <- t(sapply(1:K, function(k) colSums(gamma[,k] * data) / Nk[k]))

    # 5) solve beta
    beta <- solve_beta(t_grid, Nk, xbar_k, betaprec)

    # 6) update means
    for (k in 1:K) mu_list[[k]] <- t_grid[k] * beta

    # 7) update sigma
    for (k in 1:K) {
      diff <- sweep(data, 2, mu_list[[k]])
      W    <- sqrt(gamma[,k])
      sigma_list[[k]] <- crossprod(diff * W, diff * W)/Nk[k] + 1e-6*diag(d)
    }

    # 8) compute loglik + penalty
    ll <- sum(log(rowSums(sapply(1:K, function(k) {
      pi_k[k] * mvtnorm::dmvnorm(data, mean=mu_list[[k]], sigma=sigma_list[[k]])
    })) + 1e-300))
    ll <- ll - 0.5 * betaprec * sum(beta^2)

    loglik_trace[iter] <- ll
    if (verbose) {
      cat(sprintf("Iteration %3d: objective = %.6f\n", iter, ll))
    }
    if (iter>1 && (ll - loglik_trace[iter-1]) < tol) {
      if (verbose) cat(sprintf("Converged at iteration %d with objective %.6f\n", iter, ll))
      loglik_trace <- loglik_trace[1:iter]
      break
    }
  }

  # restore centering
  for (k in 1:K) mu_list[[k]] <- mu_list[[k]] + data_center

  params <- list(pi = pi_k,
                 mu = mu_list,
                 sigma = sigma_list,
                 beta = beta,
                 t_grid = t_grid)

  list(params = params,
       gamma = gamma,
       loglik_trace = loglik_trace)
}
