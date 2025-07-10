library(matrixStats)
library(mvtnorm)
library(Matrix)

# Numerically stable log-sum-exp
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Generate default initialization parameters for GMM
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

# 1. Cache covariance inverses and log-determinants
#
# This function computes and stores Sigma_k^{-1} and log|Sigma_k| for each component.
init_cov_cache_fast <- function(params) {
  K <- length(params$pi)
  invSigma <- vector("list", K)
  logdet   <- numeric(K)
  for (k in seq_len(K)) {
    S <- params$sigma[[k]]
    L <- chol(S)
    invSigma[[k]] <- chol2inv(L)
    logdet[k] <- 2 * sum(log(diag(L)))
  }
  params$invSigma <- invSigma
  params$logdet   <- logdet
  params
}



# 2. Vectorized E-step
#
# Computes posterior responsibilities gamma_{ik} in a single loop over K using
# matrix operations and cached inverses.
ESTEP <- function(data, params) {
  # data: n x d matrix, params: list with pi, mu, invSigma, logdet
  n <- nrow(data)
  K <- length(params$pi)
  d <- ncol(data)
  log_gamma <- matrix(0, n, K)

  # Combine component means into d x K matrix
  Mu_mat <- do.call(cbind, params$mu)
  for (k in seq_len(K)) {
    Xc <- sweep(data, 2, Mu_mat[,k], FUN = "-")  # center data by mu_k
    Qf <- rowSums((Xc %*% params$invSigma[[k]]) * Xc)  # quadratic form
    # log responsibility up to normalization
    log_gamma[,k] <- log(params$pi[k]) -
      0.5 * (Qf + params$logdet[k] + d * log(2*pi))
  }

  # Normalize with log-sum-exp for numerical stability
  lse <- rowLogSumExps(log_gamma)
  gamma <- exp(log_gamma - lse)
  gamma
}


# 3. Vectorized M-step
MSTEP <- function(data, gamma, Q_prior = NULL, modelName = "VVV") {
  n <- nrow(data)
  d <- ncol(data)
  K <- ncol(gamma)
  DK <- d * K

  Nk <- colSums(gamma)
  Nk[Nk < 1e-8] <- 1e-8

  # 1) Update mixture weights
  pi_new <- Nk / n

  # 2) Prepare right-hand side: concatenated weighted sums
  b_mat <- t(data) %*% gamma      # d × K
  b_vec <- as.numeric(b_mat)      # length d*K

  # 3) Choose solve strategy
  if (DK > 200) {
    # --- Sparse block-diagonal solve ---
    H_sp <- Diagonal(x = rep(Nk, each = d))
    if (!is.null(Q_prior)) {
      # Q_prior should be a dK×dK sparse block-diagonal dgCMatrix
      H_sp <- H_sp + Q_prior
    }
    fac    <- Cholesky(H_sp, LDL = FALSE)
    mu_vec <- solve(fac, b_vec)
  } else {
    # --- Dense block-diagonal solve via LAPACK (fast for small DK) ---
    H <- diag(rep(Nk, each = d))
    if (!is.null(Q_prior)) H <- H + as.matrix(Q_prior)
    mu_vec <- solve(H, b_vec)
  }
  # split back into list of length K
  mu_list <- split(mu_vec, rep(seq_len(K), each = d))

  # 4) Update covariances exactly as before
  sigma_list <- vector("list", K)
  for (j in seq_len(K)) {
    Xc <- sweep(data, 2, mu_list[[j]])
    if (modelName == "VVV") {
      W   <- sqrt(gamma[,j])
      S_j <- crossprod(Xc * W, Xc * W) / Nk[j] + 1e-6 * diag(d)
    } else if (modelName == "VII") {
      lam <- sum((Xc^2) * gamma[,j])/(Nk[j] * d)
      S_j <- diag(lam, d)
    } else if (modelName == "EII") {
      sigma_list[[j]] <- sum((Xc^2) * gamma[,j])
      next
    } else stop("Unsupported modelName")
    sigma_list[[j]] <- S_j
  }
  if (modelName == "EII") {
    total_sq <- sum(unlist(sigma_list))
    lam_sh   <- total_sq / (sum(Nk) * d)
    sigma_list <- lapply(seq_len(K), function(i) diag(lam_sh, d))
  }

  list(pi = pi_new, mu = mu_list, sigma = sigma_list)
}




# 4. Fast EM algorithm
EM_algorithm <- function(data, init_params,
                              Q_prior   = NULL,
                              modelName = "VVV",
                              max_iter  = 100,
                              tol       = 1e-4,
                              verbose   = TRUE) {
  params <- init_params
  # Cache inverses and log-dets
  params <- init_cov_cache_fast(params)
  loglik_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    # E-step
    gamma <- ESTEP(data, params)

    # M-step
    new_params <- MSTEP(data, gamma, Q_prior, modelName)
    # Update cache
    new_params <- init_cov_cache_fast(new_params)

    # Compute log-likelihood
    dens_mat <- sapply(seq_along(new_params$pi), function(k) {
      dmvnorm(data, mean = new_params$mu[[k]], sigma = new_params$sigma[[k]])
    })
    ll <- sum(log(dens_mat %*% new_params$pi + 1e-300))

    # Add prior penalty if provided
    if (!is.null(Q_prior)) {
      U_vec <- unlist(new_params$mu)
      penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
      logdetQ <- 0.5 * as.numeric(determinant(Q_prior, TRUE)$modulus)
      ll <- ll - as.numeric(penalty) + logdetQ
    }

    loglik_trace[iter] <- ll
    if (iter > 1 && (loglik_trace[iter] - loglik_trace[iter - 1]) < tol) {
      if (verbose) cat(sprintf("Converged at iteration %d with objective %.4f\n", iter, ll))
      loglik_trace <- loglik_trace[1:iter]
      break
    }
    if (verbose) cat(sprintf("Iteration %d: objective = %.4f\n", iter, ll))
    params <- new_params
  }

  list(params = params, gamma = gamma, loglik_trace = loglik_trace[1:iter])
}
