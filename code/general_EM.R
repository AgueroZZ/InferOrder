library(matrixStats)
library(mvtnorm)
library(Matrix)

compute_permuted_energy <- function(U, Q, P = NULL) {
  # U: K × D numeric matrix of row-means
  # Q: (K*D) × (K*D) precision matrix (from kronecker)
  # P: optional integer vector of length K; a permutation of 1:K
  #
  # Returns: (P U)ᵀ Q (P U) with U row-stacked into length-K*D

  # Input checks omitted for brevity…

  # 1) Permute rows if needed
  Up <- if (is.null(P)) U else U[P, , drop = FALSE]

  # 2) Row-major stack: c(Up[1,], Up[2,], …, Up[K,])
  u_vec <- as.numeric(t(Up))

  # 3) Quadratic form
  as.numeric(crossprod(u_vec, Q %*% u_vec))
}


optimal_permutation <- function(mu_list, Q, verbose = FALSE, random_times = 10) {
  # mu_list: list of length K of length-D numeric vectors
  # Q: (K*D) × (K*D) precision matrix
  # verbose: report only if energy improves
  # random_times: number of random permutations to try in addition to column-based

  # Stack mu_list into a K x D matrix
  U <- do.call(rbind, mu_list)
  K <- nrow(U)
  D <- ncol(U)

  # Compute base energy (no permutation)
  u_base      <- as.numeric(t(U))  # row-major stack
  energy_base <- as.numeric(crossprod(u_base, Q %*% u_base))

  best_energy <- energy_base
  best_perm   <- seq_len(K)

  # 1) Try each column-based ordering
  for (j in seq_len(D)) {
    Pj  <- order(U[, j])
    u_j <- as.numeric(t(U[Pj, , drop = FALSE]))
    Ej  <- as.numeric(crossprod(u_j, Q %*% u_j))
    if (Ej < best_energy) {
      best_energy <- Ej
      best_perm   <- Pj
    }
  }

  # 2) Try random permutations if requested
  if (random_times > 0) {
    for (t in seq_len(random_times)) {
      Pr  <- sample.int(K)  # random permutation
      u_r <- as.numeric(t(U[Pr, , drop = FALSE]))
      Er  <- as.numeric(crossprod(u_r, Q %*% u_r))
      if (Er < best_energy) {
        best_energy <- Er
        best_perm   <- Pr
      }
    }
  }

  # Verbose report
  if (verbose) {
    if (best_energy < energy_base) {
      cat(sprintf(
        "Energy improved: %.6g -> %.6g\n",
        energy_base, best_energy
      ))
    } else {
      cat("No improvement in energy (remains ",
          sprintf("%.6g", energy_base), ")\n", sep = "")
    }
  }

  return(best_perm)
}


plot_sequence <- function(U, P = NULL, tol = 1e-6) {
  # U: K × 2 numeric matrix of 2D points
  # P: optional integer vector of length K representing a permutation of 1:K
  #    If NULL, no permutation is applied
  #
  # This function draws two side-by-side plots:
  #   Left:  original sequence of points connected by arrows
  #   Right: permuted sequence (if P is provided) connected by arrows

  if (!is.matrix(U) || ncol(U) != 2) {
    stop("`U` must be a K x 2 matrix")
  }
  K <- nrow(U)

  if (!is.null(P)) {
    if (length(P) != K || !setequal(P, seq_len(K))) {
      stop("`P` must be a permutation vector of 1:K")
    }
    Up <- U[P, , drop = FALSE]
    right_title <- "Permuted sequence"
  } else {
    Up <- U
    right_title <- "Original sequence (no P)"
  }

  # Set up side-by-side plotting
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par), add = TRUE)

  # Left panel: original sequence
  plot(U, pch = 19, asp = 1,
       xlab = "U[,1]", ylab = "U[,2]",
       main = "Original sequence")
  for (k in seq_len(K - 1)) {
    if (sqrt(sum((U[k + 1, ] - U[k, ])^2)) > tol) {
      arrows(U[k, 1], U[k, 2],
             U[k + 1, 1], U[k + 1, 2],
             col = "orange", lwd = 2, length = 0.1)
    }
  }

  # Right panel: permuted sequence
  plot(Up, pch = 19, asp = 1,
       xlab = "Up[,1]", ylab = "Up[,2]",
       main = right_title)
  for (k in seq_len(K - 1)) {
    if (sqrt(sum((Up[k + 1, ] - Up[k, ])^2)) > tol) {
      arrows(Up[k, 1], Up[k, 2],
             Up[k + 1, 1], Up[k + 1, 2],
             col = "steelblue", lwd = 2, length = 0.1)
    }
  }
}

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

# initialization using k-means clustering on ordering_vec
make_init <- function(X, ordering_vec, K, assume_EEI = TRUE) {
  d <- ncol(X)

  # K-means on ordering vector
  cl <- kmeans(ordering_vec, centers = K)
  cluster_ids <- cl$cluster

  # Reorder clusters from left to right
  cluster_order <- order(cl$centers)
  cluster_rank  <- match(cluster_ids, cluster_order)

  # Means
  mu_list <- lapply(1:K, function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0) rep(0, d) else colMeans(X[idx, , drop = FALSE])
  })

  # EEI-style sigma
  if (assume_EEI) {
    # Compute global diagonal variance
    diffs <- do.call(rbind, lapply(1:K, function(k) {
      idx <- which(cluster_rank == k)
      sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
    }))
    diag_var <- apply(diffs^2, 2, mean) + 1e-6
    Sigma_shared <- diag(diag_var, d)
    sigma_list   <- replicate(K, Sigma_shared, simplify = FALSE)
  } else {
    # Full covariance per cluster
    sigma_list <- lapply(1:K, function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) <= 1) {
        diag(rep(1, d))
      } else {
        cov(X[idx, , drop = FALSE]) + 1e-6 * diag(d)
      }
    })
  }

  Nk <- table(factor(cluster_rank, levels = 1:K))
  pi_vec <- as.numeric(Nk) / length(ordering_vec)

  list(pi = pi_vec, mu = mu_list, sigma = sigma_list)
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
MSTEP <- function(data, gamma, params, Q_prior = NULL,
                  modelName = "VVV") {
  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)
  DK <- d * K

  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Step 1: weighted sum ∑ γ_ik x_i
  Wx <- t(data) %*% gamma  # d × K

  # Step 2: Regularized mean estimation
  mu_list <- lapply(seq_len(K), function(k) Wx[,k] / Nk[k])

  if (modelName %in% c("VVV", "VII", "EII", "EEI") && !is.null(Q_prior)) {
    blocks <- vector("list", K)
    rhs    <- numeric(DK)
    for (k in seq_len(K)) {
      Sinv <- params$invSigma[[k]]
      blocks[[k]] <- Nk[k] * Sinv
      rhs[((k-1)*d + 1):(k*d)] <- Sinv %*% Wx[,k]
    }

    if (DK > 200) {
      Hsp     <- bdiag(blocks) + Q_prior
      fac     <- Cholesky(Hsp, LDL = FALSE)
      mu_vec  <- solve(fac, rhs)
    } else {
      H       <- as.matrix(bdiag(blocks)) + as.matrix(Q_prior)
      mu_vec  <- solve(H, rhs)
    }
    mu_list <- split(mu_vec, rep(seq_len(K), each = d))
  }

  # Step 3: Estimate covariances using regularized mu_list
  sigma_list <- vector("list", K)
  if (modelName == "EEI") {
    V <- numeric(d)
    for (k in seq_len(K)) {
      diff <- sweep(data, 2, mu_list[[k]])
      V <- V + colSums(gamma[,k] * (diff^2))
    }
    V_shared <- V / n + 1e-6
    Sigma_shared <- diag(V_shared, d)
    sigma_list <- replicate(K, Sigma_shared, simplify = FALSE)

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

  return(list(pi = pi_new, mu = mu_list, sigma = sigma_list))
}

# 4. EM wrapper
EM_algorithm <- function(data, init_params,
                         Q_prior   = NULL,
                         modelName = "VVV",
                         max_iter  = 100,
                         tol       = 1e-4,
                         relative_lambda = FALSE,
                         permutation = FALSE,
                         random_times = 10,
                         verbose   = TRUE) {

  params       <- init_cov_cache_fast(init_params)
  K <- length(params$pi)
  loglik_trace <- numeric(max_iter)

  if (relative_lambda) {
    # store the original Q_prior for rescaling
    Q_prior_orig <- Q_prior
  }

  for (iter in seq_len(max_iter)) {

    # Check if the permutation is requested
    if (permutation) {
      # Compute the optimal permutation
      perm_idx <- optimal_permutation(params$mu, Q_prior, verbose = verbose, random_times = random_times)
      # Apply the permutation to the parameters
      params$pi    <- params$pi[perm_idx]
      params$mu    <- params$mu[perm_idx]
      params$sigma <- params$sigma[perm_idx]
      # # Apply the same permutation to gamma
      # gamma <- gamma[, perm_idx, drop = FALSE]
    }

    # Optional rescaling of Q_prior if relative_lambda is TRUE
    if (relative_lambda) {
      sigma_vec <- diag(params$sigma[[1]])  # EEI: shared diagonal
      if (length(unique(sapply(params$sigma, function(S) all.equal(S, params$sigma[[1]])))) > 1) {
        stop("Relative lambda requires shared diagonal covariance (EEI) or identical diagonal covariances (EII) across clusters.")
      }
      scale_vec <- rep(1 / sqrt(sigma_vec), times = K)  # correct cluster-major order
      Sscale <- Diagonal(x = scale_vec)
      # rescale the original Q_prior
      Q_prior <- Sscale %*% Q_prior_orig %*% Sscale
    }

    gamma      <- ESTEP(data, params)
    new_params <- MSTEP(data, gamma, params, Q_prior, modelName)

    last_params <- params
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
      logdetQ <- generalized_logdet(Q_prior)
      ll <- ll - as.numeric(penalty) + 0.5 * logdetQ
    }

    loglik_trace[iter] <- ll
    if (verbose) cat(sprintf("Iteration %3d: objective = %.6f\n", iter, ll))
    if (iter > 1 && (ll - loglik_trace[iter-1]) < tol) {
      # Check if the last log likelihood is larger than the previous one
      if (ll < loglik_trace[iter - 1]) {
        if (verbose) cat(sprintf("Objective decreased at iteration %d, final objective = %.6f\n", iter, loglik_trace[iter - 1]))
        loglik_trace <- loglik_trace[1:(iter - 1)]
        params <- last_params
        gamma  <- ESTEP(data, params)
        break
      }else{
        if (verbose) cat(sprintf("Converged at iteration %d with objective %.6f\n", iter, ll))
        loglik_trace <- loglik_trace[1:iter]
        break
      }
    }
  }

  list(params       = params,
       gamma        = gamma,
       loglik_trace = loglik_trace)
}



