# compute log-sum-exp in a numerically stable way
logsumexp <- function(x) {
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}

# E-step function for the EM algorithm
ESTEP <- function(data, params) {
  n <- nrow(data)
  K <- length(params$pi)
  log_gamma <- matrix(NA, n, K)

  for (k in 1:K) {
    log_density <- mvtnorm::dmvnorm(data, mean=params$mu[[k]], sigma=params$sigma[[k]], log=TRUE)
    log_gamma[, k] <- log(params$pi[k] + 1e-300) + log_density
  }

  # log-sum-exp across clusters for each observation
  logsumexp_vec <- apply(log_gamma, 1, logsumexp)

  # Normalize responsibilities
  gamma <- exp(log_gamma - logsumexp_vec)

  return(gamma)
}

# M-step function for the EM algorithm
MSTEP_general <- function(data, gamma, Q_prior = NULL, modelName = "VVV") {
  n <- nrow(data)
  d <- ncol(data)
  K <- ncol(gamma)

  Nk <- colSums(gamma)
  Nk[Nk < 1e-8] <- 1e-8

  # ======== 1. Compute mu ========
  # Solve (H + Q_prior) * U = b
  b_vec <- numeric(d * K)
  H <- matrix(0, d * K, d * K)
  for (k in 1:K) {
    weighted_sum <- colSums(gamma[, k] * data)
    b_vec[((k - 1) * d + 1):(k * d)] <- weighted_sum
    H[((k - 1) * d + 1):(k * d), ((k - 1) * d + 1):(k * d)] <- Nk[k] * diag(d)
  }
  if (!is.null(Q_prior)) {
    H <- H + Q_prior
  }
  mu_vec <- solve(H, b_vec)
  mu_list <- lapply(1:K, function(k) mu_vec[((k - 1) * d + 1):(k * d)])

  # ======== 2. Compute Sigma ========
  sigma_list <- vector("list", K)
  if (modelName == "VVV") {
    # Full covariance for each
    for (k in 1:K) {
      diff <- sweep(data, 2, mu_list[[k]])
      weighted_outer <- matrix(0, d, d)
      for (i in 1:n) {
        weighted_outer <- weighted_outer + gamma[i,k] * (diff[i,] %o% diff[i,])
      }
      sigma_list[[k]] <- weighted_outer / Nk[k] + 1e-6 * diag(d)
    }
  } else if (modelName == "VII") {
    # Spherical with separate variance
    for (k in 1:K) {
      diff <- sweep(data, 2, mu_list[[k]])
      weighted_sqsum <- 0
      for (i in 1:n) {
        weighted_sqsum <- weighted_sqsum + gamma[i,k] * sum(diff[i,]^2)
      }
      lambda_k <- weighted_sqsum / (Nk[k]*d)
      sigma_list[[k]] <- diag(lambda_k, d)
    }
  } else if (modelName == "EII") {
    # Shared spherical covariance
    total_sqsum <- 0
    total_N <- sum(Nk)
    for (k in 1:K) {
      diff <- sweep(data, 2, mu_list[[k]])
      for (i in 1:n) {
        total_sqsum <- total_sqsum + gamma[i,k] * sum(diff[i,]^2)
      }
    }
    lambda_shared <- total_sqsum / (total_N * d)
    for (k in 1:K) {
      sigma_list[[k]] <- diag(lambda_shared, d)
    }
  } else {
    stop(sprintf("Unsupported modelName: %s", modelName))
  }

  # ======== 3. Mixing proportions ========
  pi_new <- Nk / n

  return(list(pi = pi_new, mu = mu_list, sigma = sigma_list))
}

# EM algorithm function
EM_algorithm <- function(data, Q_prior = NULL, init_params,
                         modelName = "VVV",
                         max_iter=100, tol=1e-4, verbose=TRUE) {
  params <- init_params
  K <- length(params$pi)
  d <- ncol(data)
  n <- nrow(data)

  loglik_trace <- numeric(max_iter)

  for (iter in 1:max_iter) {
    # E-STEP
    gamma <- ESTEP(data, params)

    # M-STEP
    new_params <- MSTEP_general(data, gamma, Q_prior, modelName=modelName)

    # Compute log-likelihood + prior penalty
    ll <- 0
    for (i in 1:n) {
      temp <- 0
      for (k in 1:K) {
        dens <- mvtnorm::dmvnorm(data[i,], mean=new_params$mu[[k]], sigma=new_params$sigma[[k]])
        temp <- temp + new_params$pi[k] * dens
      }
      ll <- ll + log(temp + 1e-300)
    }

    # Add prior penalty if needed
    if (!is.null(Q_prior)) {
      U_vec <- unlist(lapply(new_params$mu, as.numeric))
      penalty <- 0.5 * t(U_vec) %*% Q_prior %*% U_vec
      logdet_Q <- 0.5 * as.numeric(determinant(Q_prior, logarithm=TRUE)$modulus)
      ll <- ll - as.numeric(penalty) + logdet_Q
    }

    loglik_trace[iter] <- ll

    # Convergence
    if (iter > 1 && (loglik_trace[iter] - loglik_trace[iter - 1]) < tol) {
      if (verbose) cat(sprintf("Converged at iteration %d with objective %.4f\n", iter, ll))
      loglik_trace <- loglik_trace[1:iter]
      break
    }
    if (verbose) cat(sprintf("Iteration %d: objective = %.4f\n", iter, ll))

    params <- new_params
  }

  return(list(params = params, gamma = gamma, loglik_trace = loglik_trace))
}


# Generate default initialization parameters
make_default_init <- function(X, K) {
  d <- ncol(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  list(
    pi = rep(1/K, K),
    mu = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = lapply(1:K, function(k) diag(d))
  )
}
