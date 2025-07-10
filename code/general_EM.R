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
    pi    = rep(1/K, K),
    mu    = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = lapply(1:K, function(k) diag(d))
  )
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

# 3. M‐step with “EEI” (shared diagonal) and support for VVV/VII/EII/EEI prior‐weighted means
MSTEP <- function(data, gamma, params, Q_prior = NULL, modelName = "VVV") {
  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)
  DK <- d * K

  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Weighted sums ∑ γ_ik x_i
  Wx <- t(data) %*% gamma  # d × K

  # 1) Means (unregularized)
  mu_list <- lapply(seq_len(K), function(k) Wx[,k] / Nk[k])

  # 2) Covariances
  sigma_list <- vector("list", K)
  if (modelName == "EEI") {
    V <- numeric(d)
    for (k in seq_len(K)) {
      diff <- sweep(data, 2, mu_list[[k]])
      V   <- V + colSums(gamma[,k] * (diff^2))
    }
    V_shared    <- V / n + 1e-6
    Sigma_shared <- diag(V_shared, d)
    sigma_list  <- replicate(K, Sigma_shared, simplify = FALSE)
  } else {
    for (k in seq_len(K)) {
      diff <- sweep(data, 2, mu_list[[k]])
      W    <- sqrt(gamma[,k])
      if (modelName == "VVV") {
        sigma_list[[k]] <- crossprod(diff * W, diff * W) / Nk[k] + 1e-6 * diag(d)
      } else if (modelName == "VII") {
        lam <- sum((diff^2) * gamma[,k]) / (Nk[k] * d)
        sigma_list[[k]] <- diag(lam, d)
      } else if (modelName == "EII") {
        sigma_list[[k]] <- sum((diff^2) * gamma[,k])
      } else {
        stop("Unsupported modelName: ", modelName)
      }
    }
    if (modelName == "EII") {
      total <- sum(unlist(sigma_list))
      lam_sh <- total / (sum(Nk) * d)
      sigma_list <- replicate(K, diag(lam_sh, d), simplify = FALSE)
    }
  }

  # 3) Precision‐weighted solve for VVV/VII/EII/EEI when prior present
  if (modelName %in% c("VVV", "VII", "EII", "EEI") && !is.null(Q_prior)) {
    blocks <- vector("list", K)
    rhs    <- numeric(DK)
    for (k in seq_len(K)) {
      Sinv        <- params$invSigma[[k]]
      blocks[[k]] <- Nk[k] * Sinv
      rhs[((k-1)*d + 1):(k*d)] <- Sinv %*% Wx[,k]
    }
    if (DK > 200) {
      Hsp   <- bdiag(blocks) + Q_prior
      fac   <- Cholesky(Hsp, LDL = FALSE)
      mu_vec <- solve(fac, rhs)
    } else {
      H      <- as.matrix(bdiag(blocks)) + as.matrix(Q_prior)
      mu_vec <- solve(H, rhs)
    }
    mu_list <- split(mu_vec, rep(seq_len(K), each = d))
  }

  list(pi    = pi_new,
       mu    = mu_list,
       sigma = sigma_list)
}

# 4. EM wrapper
EM_algorithm <- function(data, init_params,
                         Q_prior   = NULL,
                         modelName = "VVV",
                         max_iter  = 100,
                         tol       = 1e-4,
                         verbose   = TRUE) {
  params       <- init_cov_cache_fast(init_params)
  loglik_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    gamma      <- ESTEP(data, params)
    new_params <- MSTEP(data, gamma, params, Q_prior, modelName)
    params     <- init_cov_cache_fast(new_params)

    dens_mat <- sapply(seq_along(params$pi), function(k) {
      params$pi[k] * dmvnorm(data,
                             mean  = params$mu[[k]],
                             sigma = params$sigma[[k]])
    })
    ll <- sum(log(rowSums(dens_mat) + 1e-300))
    if (!is.null(Q_prior)) {
      U_vec   <- unlist(params$mu)
      penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
      logdetQ <- 0.5 * as.numeric(determinant(Q_prior, TRUE)$modulus)
      ll <- ll - as.numeric(penalty) + logdetQ
    }

    loglik_trace[iter] <- ll
    if (verbose) cat(sprintf("Iteration %3d: objective = %.6f\n", iter, ll))
    if (iter > 1 && (ll - loglik_trace[iter-1]) < tol) {
      if (verbose) cat(sprintf("Converged at iteration %d with objective %.6f\n", iter, ll))
      loglik_trace <- loglik_trace[1:iter]
      break
    }
  }

  list(params       = params,
       gamma        = gamma,
       loglik_trace = loglik_trace)
}



