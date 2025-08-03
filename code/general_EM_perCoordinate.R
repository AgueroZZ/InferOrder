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
make_default_init <- function(X, K, ordering = TRUE) {
  d <- ncol(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  output <- list(
    pi    = rep(1/K, K),
    mu    = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = lapply(1:K, function(k) diag(d))
  )

  if(ordering){
    # sort by the pivot mu (i.e. the mu with largest range)
    pivot_index <- which.max(apply(do.call(rbind, output$mu), 2, function(mu) max(mu) - min(mu)))

    mu_first <- sapply(output$mu, function(mu) mu[pivot_index])
    order_idx <- order(mu_first)
    output$pi    <- output$pi[order_idx]
    output$mu    <- output$mu[order_idx]
    output$sigma <- output$sigma[order_idx]
  }

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

# 3. M‐step coordinate-wise
MSTEP_coordinatewise <- function(data, gamma, Q_prior, eps = 1e-6, parallel = FALSE, ncores = 2) {
  n <- nrow(data)
  d <- ncol(data)
  K <- ncol(gamma)

  Nk <- colSums(gamma)
  Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Step 1: update mu per coordinate
  update_fn <- function(j) {
    compute_mu_posterior(xj = data[, j], Gamma = gamma, Q_1D = Q_prior)
  }

  if (parallel) {
    mu_list_j <- parallel::mclapply(1:d, update_fn, mc.cores = ncores)
  } else {
    mu_list_j <- lapply(1:d, update_fn)
  }

  mu_mat <- do.call(rbind, mu_list_j)  # d x K
  mu_list <- lapply(1:K, function(k) mu_mat[, k])

  # Step 2: update shared diagonal variance (EEI)
  V <- numeric(d)
  for (k in 1:K) {
    diff <- sweep(data, 2, mu_list[[k]])
    V <- V + colSums(gamma[, k] * diff^2)
  }
  V_shared <- V / n + eps
  Sigma_shared <- diag(V_shared, d)
  sigma_list <- replicate(K, Sigma_shared, simplify = FALSE)

  return(list(
    pi    = pi_new,
    mu    = mu_list,
    sigma = sigma_list
  ))
}



