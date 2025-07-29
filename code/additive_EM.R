library(matrixStats)
library(mvtnorm)
library(Matrix)

# Compute generalized log-determinant (sum of log of nonzero eigenvalues)
generalized_logdet <- function(Q, tol = 1e-10) {
  # Q: symmetric matrix
  # tol: eigenvalue threshold; if NULL, use a heuristic
  #
  # Returns sum(log(λ_i)) over all eigenvalues λ_i > tol.

  # 1) compute eigenvalues
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)
  vals <- ev$values

  # 2) set a default tolerance if none given
  if (is.null(tol)) {
    # anything smaller than n*max(vals)*eps is treated as zero
    tol <- length(vals) * max(vals) * .Machine$double.eps
  }

  # 3) keep only those above tol
  positive <- vals[vals > tol]
  if (length(positive) == 0) {
    warning("All eigenvalues are below tolerance; returning 0")
    return(0)
  }

  # 4) sum logs
  sum(log(positive))
}

# Generate default initialization parameters for GMM
make_default_init_add <- function(X, Q_list) {
  d <- ncol(X)

  K_vec <- sapply(Q_list, function(Qj) nrow(Qj) / d)  # Qj is K_j * d x K_j * d
  K <- prod(K_vec)

  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  output <- list(
    pi    = rep(1/K, K),
    mu    = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = lapply(1:K, function(k) diag(d))
  )

  # sort by the pivot mu (i.e. the mu with largest range)
  pivot_index <- which.max(apply(do.call(rbind, output$mu), 2, function(mu) max(mu) - min(mu)))

  mu_first <- sapply(output$mu, function(mu) mu[pivot_index])
  order_idx <- order(mu_first)
  output$pi    <- output$pi[order_idx]
  output$mu    <- output$mu[order_idx]
  output$sigma <- output$sigma[order_idx]
  return(output)
}

# 1. Cache covariance inverses and log-determinants
init_cov_cache_fast <- function(params) {
  K <- length(params$pi)
  invSigma <- vector("list", K)
  logdet   <- numeric(K)
  for (k in seq_len(K)) {
    S <- params$sigma[[k]]
    L <- chol(S)
    invSigma[[k]] <- chol2inv(L)
    logdet[k]     <- 2 * sum(log(diag(L)))
  }
  params$invSigma <- invSigma
  params$logdet   <- logdet
  params
}

# 2. Vectorized E-step
ESTEP <- function(data, params) {
  n <- nrow(data); K <- length(params$pi); d <- ncol(data)
  log_gamma <- matrix(0, n, K)
  Mu_mat <- do.call(cbind, params$mu)
  for (k in seq_len(K)) {
    Xc <- sweep(data, 2, Mu_mat[,k], "-")
    Qf <- rowSums((Xc %*% params$invSigma[[k]]) * Xc)
    log_gamma[,k] <- log(params$pi[k]) - 0.5*(Qf + params$logdet[k] + d*log(2*pi))
  }
  lse <- rowLogSumExps(log_gamma)
  exp(log_gamma - lse)
}

# 3. M-step for additive model
MSTEP_additive <- function(data, gamma, params, A, Q_theta, modelName = "EEI") {
  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)
  DK <- d * K

  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Step 1: Weighted average for each cluster
  Wx <- t(data) %*% gamma   # d × K
  Wx_vec <- as.numeric(Wx) # dK vector

  if (modelName != "EEI") stop("Currently only EEI supported for additive M-step")

  # Shared diagonal covariance
  sigma_vec <- numeric(d)
  mu_mat <- A  # A: (dK × Dtheta), e.g., Dtheta = K1 + K2

  # Step 2: Estimate θ by solving penalized least squares
  blocks <- vector("list", K)
  rhs <- numeric(d * K)

  # Construct weight matrix W = diag(N_k / sigma)
  sigma_inv_diag <- rep(1 / diag(params$sigma[[1]]), times = K)
  W_diag <- rep(Nk, each = d) * sigma_inv_diag
  W_mat <- Diagonal(x = W_diag)

  lhs <- t(mu_mat) %*% W_mat %*% mu_mat + Q_theta
  rhs <- t(mu_mat) %*% (sigma_inv_diag * Wx_vec)

  fac <- Cholesky(lhs, LDL = FALSE)
  theta_hat <- solve(fac, rhs)

  # Step 3: Generate mu_list = A %*% theta
  mu_vec <- as.vector(mu_mat %*% theta_hat)
  mu_list <- split(mu_vec, rep(seq_len(K), each = d))

  # Step 4: Shared diagonal covariance (same as before)
  V <- numeric(d)
  for (k in seq_len(K)) {
    diff <- sweep(data, 2, mu_list[[k]])
    V <- V + colSums(gamma[,k] * (diff^2))
  }
  V_shared <- V / n + 1e-6
  Sigma_shared <- diag(V_shared, d)
  sigma_list <- replicate(K, Sigma_shared, simplify = FALSE)

  list(pi = pi_new, mu = mu_list, sigma = sigma_list, theta = theta_hat)
}

# 4. Create the observation matrix A
make_additive_A <- function(K_vec) {
  K <- prod(K_vec)  # Total number of clusters
  total_theta_dim <- sum(K_vec)

  # Resulting A matrix is K x sum(K_vec)
  A <- Matrix::Matrix(0, nrow = K, ncol = total_theta_dim, sparse = TRUE)

  # Each row corresponds to an index combination (i1, i2, ..., ip)
  index_mat <- as.matrix(expand.grid(lapply(K_vec, seq_len)))  # K x p

  for (row in seq_len(K)) {
    idxs <- index_mat[row, ]  # vector of length p
    col_pos <- mapply(function(i, offset) i + offset,
                      i = idxs,
                      offset = cumsum(c(0, K_vec[-length(K_vec)])))
    A[row, col_pos] <- 1
  }

  return(A)
}

# hlep function to split theta by Q_list
split_theta_by_Qlist <- function(theta, Q_list, d = NULL) {
  # theta: stacked vector of all additive components
  # Q_list: list of Q matrices, each of size (Kj·d × Kj·d)
  # d: number of coordinates

  if (is.null(d)) {
    # infer d from Q_list structure
    d_vec <- sapply(Q_list, function(Qj) {
      nrow(Qj) / round(nrow(Qj) / ncol(Qj))  # extract d
    })
    if (length(unique(d_vec)) != 1) stop("Inconsistent d across Q_list")
    d <- d_vec[1]
  }

  theta_list <- list()
  cursor <- 1
  for (i in seq_along(Q_list)) {
    Qj <- Q_list[[i]]
    Kj <- nrow(Qj) / d
    n_elem <- Kj * d
    theta_block <- theta[cursor:(cursor + n_elem - 1)]
    theta_mat <- matrix(theta_block, nrow = Kj, ncol = d, byrow = TRUE)
    theta_list[[i]] <- theta_mat
    cursor <- cursor + n_elem
  }

  theta_list
}


# 5. additive smooth EM algorithm
additive_smoothEM <- function(data, Q_list,
                              init_params     = NULL,
                              modelName       = "EEI",
                              max_iter        = 100,
                              tol             = 1e-4,
                              relative_lambda = FALSE,
                              verbose         = TRUE) {
  n <- nrow(data)
  d <- ncol(data)

  # Each Qj is (Kj·d) × (Kj·d), so Kj = nrow(Qj) / d
  K_vec   <- sapply(Q_list, function(Qj) nrow(Qj) / d)
  K       <- prod(K_vec)
  D_theta <- sum(sapply(Q_list, nrow))

  # Build A and original Q_theta
  A_scalar    <- make_additive_A(K_vec)             # K × sum Kj
  A_block     <- kronecker(A_scalar, Diagonal(d))   # (K·d) × (sum Kj·d)
  Q_theta_orig <- bdiag(Q_list)                     # prior precision

  # Initialize parameters
  if (is.null(init_params)) {
    theta0    <- rnorm(D_theta)
    mu_mat    <- matrix(A_block %*% theta0, nrow = K, ncol = d, byrow = TRUE)
    mu_list   <- split(mu_mat, rep(seq_len(K), each = d))
    pi0       <- rep(1/K, K)
    sigma0    <- diag(apply(data, 2, var) + 1e-6)
    sigma_list0 <- replicate(K, sigma0, simplify = FALSE)
    params    <- init_cov_cache_fast(list(mu = mu_list, pi = pi0, sigma = sigma_list0))
  } else {
    params <- init_cov_cache_fast(init_params)
  }

  # Storage for log-likelihood components
  loglik_data    <- numeric(max_iter)
  loglik_penalty <- numeric(max_iter)
  loglik_trace   <- numeric(max_iter)

  # start iterations
  for (iter in seq_len(max_iter)) {
    # keep previous full state in case we need to roll back
    prev_params  <- params
    prev_gamma   <- if (iter > 1) gamma else NULL
    prev_theta   <- if (iter > 1) theta_hat else NULL
    prev_obj     <- if (iter > 1) loglik_trace[iter - 1] else -Inf

    # 1) optionally rescale prior
    Q_theta <- Q_theta_orig
    if (relative_lambda) {
      current_sigma <- diag(params$sigma[[1]])
      scale_vec     <- rep(1 / sqrt(current_sigma), times = sum(K_vec))
      Sscale        <- Diagonal(x = scale_vec)
      Q_theta       <- Sscale %*% Q_theta_orig %*% Sscale
    }

    # 2) E-step
    gamma <- ESTEP(data, params)

    # 3) M-step
    new_params <- MSTEP_additive(
      data     = data,
      gamma    = gamma,
      params   = params,
      A        = A_block,
      Q_theta  = Q_theta,
      modelName = modelName
    )
    theta_hat <- new_params$theta
    params     <- init_cov_cache_fast(new_params)

    # 4) compute penalized log-likelihood
    dens_mat <- sapply(seq_len(K), function(k) {
      params$pi[k] * dmvnorm(data,
                             mean  = params$mu[[k]],
                             sigma = params$sigma[[k]])
    })
    loglik_data[iter]    <- sum(log(rowSums(dens_mat) + 1e-300))
    loglik_penalty[iter] <- -0.5 * crossprod(theta_hat, Q_theta %*% theta_hat) + 0.5 * generalized_logdet(Q_theta)
    loglik_trace[iter]   <- loglik_data[iter] + loglik_penalty[iter]

    if (verbose) {
      cat(sprintf("Iter %3d: data_ll = %.4f, penalty = %.4f, total = %.4f\n",
                  iter, loglik_data[iter], loglik_penalty[iter], loglik_trace[iter]))
    }

    # 5) stop if objective decreased
    if (iter > 1 && loglik_trace[iter] < prev_obj) {
      if (verbose) {
        cat(sprintf("Objective decreased at iter %d (%.4f → %.4f), rolling back and stopping.\n",
                    iter - 1, prev_obj, loglik_trace[iter]))
      }
      # roll back to previous parameters
      params       <- prev_params
      gamma        <- prev_gamma
      theta_hat    <- prev_theta
      # truncate logs
      loglik_data    <- loglik_data[1:(iter - 1)]
      loglik_penalty <- loglik_penalty[1:(iter - 1)]
      loglik_trace   <- loglik_trace[1:(iter - 1)]
      break
    }

    # 6) check convergence
    if (iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol) {
      if (verbose) cat(sprintf("Converged at iter %d\n", iter))
      loglik_data    <- loglik_data[1:iter]
      loglik_penalty <- loglik_penalty[1:iter]
      loglik_trace   <- loglik_trace[1:iter]
      break
    }
  }

  # 7) split θ into list of matrices
  theta_list <- split_theta_by_Qlist(theta_hat, Q_list, d = d)

  list(
    params         = params,
    theta          = theta_list,
    A              = A_block,
    gamma          = gamma,
    loglik_data    = loglik_data,
    loglik_penalty = loglik_penalty,
    loglik_trace   = loglik_trace
  )
}
