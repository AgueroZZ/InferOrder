# internal functions for EM algorithm with linear model
# solve_beta <- function(t_grid, Nk, xbar_k, betaprec) {
#   K <- length(t_grid)
#   d <- ncol(xbar_k)
#
#   # scalar
#   sum_N_t2 <- sum(Nk * t_grid^2)
#   rhs <- colSums(matrix(Nk * t_grid, nrow=K, ncol=d) * xbar_k)
#
#   # Solve
#   lhs <- betaprec * diag(d) + sum_N_t2 * diag(d)
#   beta_hat <- solve(lhs, rhs)
#   return(beta_hat)
# }

# EM algorithm for linear model
# EM_algorithm_linear <- function(data, K, betaprec = 1,
#                                 t_grid   = NULL, seed = 1,
#                                 sigma_init = NULL,
#                                 max_iter = 100, tol = 1e-4,
#                                 verbose  = TRUE) {
#
#   data_center <- colMeans(data)
#   data <- sweep(data, 2, data_center)
#
#   n <- nrow(data); d <- ncol(data)
#   if (is.null(t_grid)) t_grid <- seq(-1, 1, length.out = K)
#
#   set.seed(seed)
#   beta       <- rnorm(d)
#   pi_k       <- rep(1/K, K)
#   sigma_list <- if (is.null(sigma_init)) lapply(1:K, function(k) diag(d)) else sigma_init
#
#   loglik_trace <- numeric(max_iter)
#   mu_list      <- vector("list", K)
#   gamma        <- matrix(0, n, K)
#
#   for (iter in 1:max_iter) {
#     # 1) component means
#     for (k in 1:K) mu_list[[k]] <- t_grid[k] * beta
#
#     # 2) fast E-step must have invSigma/logdet in params
#     params <- list(pi = pi_k, mu = mu_list, sigma = sigma_list)
#     params <- init_cov_cache_fast(params)   # <-- cache inverses & logdets
#     gamma  <- ESTEP(data, params)      # or just ESTEP(data, params)
#
#     # 3) update pi
#     Nk   <- colSums(gamma)
#     pi_k <- Nk / n
#
#     # 4) weighted means
#     xbar_k <- t(sapply(1:K, function(k) colSums(gamma[,k] * data) / Nk[k]))
#
#     # 5) solve beta
#     beta <- solve_beta(t_grid, Nk, xbar_k, betaprec)
#
#     # 6) update means
#     for (k in 1:K) mu_list[[k]] <- t_grid[k] * beta
#
#     # 7) update sigma
#     for (k in 1:K) {
#       diff <- sweep(data, 2, mu_list[[k]])
#       W    <- sqrt(gamma[,k])
#       sigma_list[[k]] <- crossprod(diff * W, diff * W)/Nk[k] + 1e-6*diag(d)
#     }
#
#     # 8) compute loglik + penalty
#     ll <- sum(log(rowSums(sapply(1:K, function(k) {
#       pi_k[k] * mvtnorm::dmvnorm(data, mean=mu_list[[k]], sigma=sigma_list[[k]])
#     })) + 1e-300))
#     ll <- ll - 0.5 * betaprec * sum(beta^2)
#
#     loglik_trace[iter] <- ll
#     if (verbose) {
#       cat(sprintf("Iteration %3d: objective = %.6f\n", iter, ll))
#     }
#     if (iter>1 && (ll - loglik_trace[iter-1]) < tol) {
#       if (verbose) cat(sprintf("Converged at iteration %d with objective %.6f\n", iter, ll))
#       loglik_trace <- loglik_trace[1:iter]
#       break
#     }
#   }
#
#   # restore centering
#   for (k in 1:K) mu_list[[k]] <- mu_list[[k]] + data_center
#
#   params <- list(pi = pi_k,
#                  mu = mu_list,
#                  sigma = sigma_list,
#                  beta = beta,
#                  t_grid = t_grid)
#
#   list(params = params,
#        gamma = gamma,
#        loglik_trace = loglik_trace)
# }




#’ Solve for beta with full Σ_k^{-1} weighting
#’
#’ @param t_grid         Numeric vector of length K with the grid points t_k
#’ @param Nk             Numeric vector of length K with responsibility sums N_k
#’ @param xbar_k         K×d matrix of weighted sample means \bar x_k
#’ @param betaprec       Scalar prior precision on beta
#’ @param invSigma_list  List of length K, each a d×d matrix Σ_k^{-1}
#’ @return                Numeric vector of length d (the updated beta)
solve_beta <- function(t_grid, Nk, xbar_k, betaprec, invSigma_list) {
  K <- length(t_grid)
  d <- ncol(xbar_k)

  # Build LHS = betaprec*I + ∑_k N_k t_k^2 Σ_k^{-1}
  LHS <- betaprec * diag(d)
  for (k in seq_len(K)) {
    LHS <- LHS + Nk[k] * (t_grid[k]^2) * invSigma_list[[k]]
  }

  # Build RHS = ∑_k t_k Σ_k^{-1} (N_k * xbar_k[k,])
  RHS <- numeric(d)
  for (k in seq_len(K)) {
    RHS <- RHS + t_grid[k] * (invSigma_list[[k]] %*% (Nk[k] * xbar_k[k, ]))
  }

  # Solve for beta
  beta_hat <- solve(LHS, RHS)
  beta_hat <- as.numeric(beta_hat)
  return(beta_hat)
}


#’ EM algorithm for the linear‐prior model using fast E‐step
#’
#’ @param data         n×d data matrix
#’ @param K            Number of grid points / components
#’ @param betaprec     Scalar prior precision on beta
#’ @param t_grid       Optional K‐vector of grid points (defaults to seq(-1,1))
#’ @param seed         RNG seed for initialization
#’ @param sigma_init   Optional list of K initial covariance matrices
#’ @param max_iter     Maximum EM iterations
#’ @param tol          Convergence tolerance on log‐likelihood change
#’ @param verbose      Whether to print iteration logs
#’ @return             List with params (pi, mu, sigma, beta, t_grid), gamma, loglik_trace
EM_algorithm_linear <- function(data, K, betaprec = 1,
                                t_grid     = NULL, seed = 1,
                                sigma_init = NULL,
                                max_iter   = 100, tol = 1e-4,
                                verbose    = TRUE) {

  # 1) Center the data
  data_center <- colMeans(data)
  data        <- sweep(data, 2, data_center)

  n <- nrow(data); d <- ncol(data)
  if (is.null(t_grid)) t_grid <- seq(-1, 1, length.out = K)

  # 2) Initialize
  set.seed(seed)
  beta       <- rnorm(d)
  pi_k       <- rep(1/K, K)
  sigma_list <- if (is.null(sigma_init))
    replicate(K, diag(d), simplify = FALSE)
  else
    sigma_init

  loglik_trace <- numeric(max_iter)
  mu_list      <- vector("list", K)
  gamma        <- matrix(0, n, K)

  for (iter in seq_len(max_iter)) {
    # 3) Means = t_k * beta
    for (k in seq_len(K)) {
      mu_list[[k]] <- as.numeric(t_grid[k] * beta)
    }

    # 4) Fast E-step
    params <- list(pi = pi_k, mu = mu_list, sigma = sigma_list)
    params <- init_cov_cache_fast(params)
    gamma  <- ESTEP(data, params)

    # 5) Update weights
    Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
    pi_k <- Nk / n

    # 6) Weighted means
    xbar_k <- t(sapply(seq_len(K), function(k) {
      colSums(gamma[,k] * data) / Nk[k]
    }))

    # 7) Solve beta
    beta <- solve_beta(t_grid, Nk, xbar_k, betaprec, params$invSigma)

    # 8) Update means again
    for (k in seq_len(K)) {
      mu_list[[k]] <- as.numeric(t_grid[k] * beta)
    }

    # 9) EEI covariances (shared diagonal)
    V <- numeric(d)
    for (k in seq_len(K)) {
      diff <- sweep(data, 2, mu_list[[k]])
      V   <- V + colSums(gamma[,k] * (diff^2))
    }
    V_shared     <- V / n + 1e-6 + 1e-8
    Sigma_shared <- diag(V_shared, d)
    sigma_list   <- replicate(K, Sigma_shared, simplify = FALSE)

    # 10) Compute log-likelihood safely
    mixdens <- rowSums(sapply(seq_len(K), function(k) {
      pi_k[k] * mvtnorm::dmvnorm(data,
                                 mean  = mu_list[[k]],
                                 sigma = sigma_list[[k]])
    }))
    mixdens <- pmax(mixdens, 1e-300)         # floor to avoid log(0)
    ll      <- sum(log(mixdens))
    ll      <- ll - 0.5 * betaprec * sum(beta^2)

    # NaN check
    if (!is.finite(ll)) {
      warning(sprintf("Iteration %d produced NaN log-lik; stopping early.", iter))
      break
    }

    loglik_trace[iter] <- ll

    # 11) Verbose logging
    if (verbose) {
      cat(sprintf("Iteration %3d: objective = %.6f\n", iter, ll))
    }
    if (iter > 1 && (ll - loglik_trace[iter-1]) < tol) {
      if (verbose) cat(sprintf("Converged at iteration %d with objective %.6f\n",
                               iter, ll))
      break
    }
  }

  # 12) Restore centering on means
  for (k in seq_len(K)) {
    mu_list[[k]] <- mu_list[[k]] + data_center
  }

  # 13) Return (crop trace if early break)
  loglik_trace <- loglik_trace[seq_len(iter)]
  params_out   <- list(pi     = pi_k,
                       mu     = mu_list,
                       sigma  = sigma_list,
                       beta   = beta,
                       t_grid = t_grid)

  list(params       = params_out,
       gamma        = gamma,
       loglik_trace = loglik_trace)
}


