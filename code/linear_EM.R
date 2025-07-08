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

# EM algorithm using a linear prior
EM_algorithm_linear <- function(data, K, betaprec=1,
                                t_grid=NULL, seed=1,
                                sigma_init=NULL,
                                max_iter=100, tol=1e-4, verbose=TRUE) {

  # center the data
  data_center <- colMeans(data)
  data <- sweep(data, 2, data_center)

  n <- nrow(data)
  d <- ncol(data)
  if (is.null(t_grid)) {
    t_grid <- seq(-1, 1, length.out=K)
  }

  # Initialize
  set.seed(seed)
  beta <- rnorm(d)
  pi_k <- rep(1/K, K)
  if (is.null(sigma_init)) {
    sigma_list <- lapply(1:K, function(k) diag(d))
  } else {
    sigma_list <- sigma_init
  }

  loglik_trace <- numeric(max_iter)

  for (iter in 1:max_iter) {
    # ===================
    # Compute current u_k
    # ===================
    mu_list <- lapply(1:K, function(k) t_grid[k] * beta)

    # ===================
    # E-step
    # ===================
    gamma <- ESTEP(data, list(pi=pi_k, mu=mu_list, sigma=sigma_list))

    # ===================
    # M-step
    # ===================
    Nk <- colSums(gamma)
    pi_k <- Nk / n

    # Weighted means
    xbar_k <- matrix(NA, nrow=K, ncol=d)
    for (k in 1:K) {
      weighted_sum <- colSums(gamma[,k] * data)
      xbar_k[k,] <- weighted_sum / Nk[k]
    }

    # Solve for beta
    beta <- solve_beta(t_grid, Nk, xbar_k, betaprec)

    # Update mu_list with new beta
    mu_list <- lapply(1:K, function(k) t_grid[k] * beta)

    # Update sigma
    for (k in 1:K) {
      diff <- sweep(data, 2, mu_list[[k]])
      weighted_outer <- matrix(0, d, d)
      for (i in 1:n) {
        weighted_outer <- weighted_outer + gamma[i,k] * (diff[i,] %o% diff[i,])
      }
      sigma_list[[k]] <- weighted_outer / Nk[k] + 1e-6 * diag(d)
    }

    # ===================
    # Compute objective
    # ===================
    ll <- 0
    for (i in 1:n) {
      temp <- 0
      for (k in 1:K) {
        dens <- mvtnorm::dmvnorm(data[i,], mean=mu_list[[k]], sigma=sigma_list[[k]])
        temp <- temp + pi_k[k] * dens
      }
      ll <- ll + log(temp + 1e-300)
    }

    # Add prior penalty on beta
    penalty <- 0.5 * betaprec * sum(beta^2)
    ll <- ll - penalty

    loglik_trace[iter] <- ll

    if (verbose) {
      cat(sprintf("Iteration %d: objective = %.4f\n", iter, ll))
    }

    if (iter > 1 && (loglik_trace[iter] - loglik_trace[iter-1]) < tol) {
      loglik_trace <- loglik_trace[1:iter]
      break
    }
  }

  # add centering back
  mu_list <- lapply(mu_list, function(mu) mu + data_center)

  # ===================
  # Return in standard format
  # ===================
  params <- list(
    pi = pi_k,
    mu = mu_list,
    sigma = sigma_list,
    beta = beta,
    t_grid = t_grid
  )

  return(list(
    params = params,
    gamma = gamma,
    loglik_trace = loglik_trace
  ))
}
