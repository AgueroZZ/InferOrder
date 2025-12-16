library(matrixStats)
library(mvtnorm)
library(Matrix)

# Compute generalized log-determinant (sum of log of nonzero eigenvalues)
generalized_logdet <- function(Q, eigen_tol = 0, rank_deficiency = 0) {
  # Q: symmetric matrix
  # eigen_tol: eigenvalue threshold; if NULL, use a heuristic
  #
  # Returns sum(log(λ_i)) over all eigenvalues λ_i > eigen_tol

  # 1) compute eigenvalues
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)
  vals <- ev$values

  # if rank_deficiency > 0, take out the last `rank_deficiency` eigenvalues
  if (rank_deficiency > 0) {
    vals <- vals[-seq(length(vals) - rank_deficiency + 1, length(vals))]
  }

  # 2) set a default tolerance if none given
  if (is.null(eigen_tol)) {
    # anything smaller than n*max(vals)*eps is treated as zero
    eigen_tol <- length(vals) * max(vals) * .Machine$double.eps
  }

  # 3) keep only those above tol
  positive <- vals[vals > eigen_tol]
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

# initialization using ordering_vec
make_init <- function(X, ordering_vec, K, assume_EEI = TRUE, nugget = 0, discretization = "equal") {
  d <- ncol(X)

  if(discretization == "kmeans"){
    # K-means on ordering vector
    cl <- kmeans(ordering_vec, centers = K)
  }
  else if(discretization == "quantile"){
    # based on cutting ordering_vec into K quantiles
    quantiles <- quantile(ordering_vec, p = seq(0,1, length.out = (K+1)))
    cl <- list()
    cl$centers <- (quantiles[-1] + quantiles[-K])/2
    assign_cluster <- function(value) which(quantiles >= value)[1]
    cl$cluster <- Vectorize(assign_cluster)(ordering_vec)
  }else{
    cutoff <- seq(min(ordering_vec), max(ordering_vec), length.out = (K + 1))
    cl <- list()
    cl$centers <- (cutoff[-1] + cutoff[-K])/2
    assign_cluster <- function(value) which(cutoff >= value)[1]
    cl$cluster <- Vectorize(assign_cluster)(ordering_vec)
  }

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
    diag_var <- apply(diffs^2, 2, mean) + nugget
    Sigma_shared <- diag(diag_var, d)
    sigma_list   <- replicate(K, Sigma_shared, simplify = FALSE)
  } else {
    # Full covariance per cluster
    sigma_list <- lapply(1:K, function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) <= 1) {
        diag(rep(1, d))
      } else {
        cov(X[idx, , drop = FALSE]) + nugget * diag(d)
      }
    })
  }

  Nk <- table(factor(cluster_rank, levels = 1:K))
  pi_vec <- as.numeric(Nk) / length(ordering_vec)

  list(pi = pi_vec, mu = mu_list, sigma = sigma_list)
}

# 0. Record convergence
compute_log_joint_observed <- function(X, params, Q_prior = NULL, eigen_tol = 0, rank_deficiency = 0) {
  K <- length(params$pi)
  log_dens_mat <- sapply(seq_len(K), function(k) {
    dmvnorm(X,
            mean  = params$mu[[k]],
            sigma = params$sigma[[k]],
            log   = TRUE)
  })
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K

  loglik <- sum(rowLogSumExps(log_gamma))

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu)
    penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    loglik  <- loglik - as.numeric(penalty) + 0.5 * logdetQ
  }

  return(loglik)
}

# 2. Penalized ELBO = Q-function + entropy + log-prior
compute_penalized_ELBO <- function(X, Gamma, params, Q_prior = NULL, eigen_tol = 0, rank_deficiency = 0 ) {
  K <- length(params$pi)
  log_dens_mat <- sapply(seq_len(K), function(k) {
    dmvnorm(X,
            mean  = params$mu[[k]],
            sigma = params$sigma[[k]],
            log   = TRUE)
  })
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K

  # 2b. Q-function: E_q[log p(X,Z)]
  Qval <- sum(Gamma * log_gamma)

  # 2c. entropy: −E_q[log q]
  entropy <- - sum(Gamma * log(Gamma + 1e-300))

  elbo <- Qval + entropy

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu)
    penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    elbo    <- elbo - as.numeric(penalty) + 0.5 * logdetQ
  }

  return(elbo)
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
# (still need to confirm the computation of optimal Sigma when EEI with relative lambda...)
MSTEP <- function(data, gamma, params, Q_prior = NULL,
                  relative_lambda = FALSE,
                  modelName    = "VVV",
                  iterate_once = TRUE,
                  nugget      = 0,
                  rank_deficiency = 0,
                  tol_inner    = 1e-6,
                  max_inner    = 20,
                  verbose      = FALSE) {

  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)
  DK <- d * K

  # 1. Update mixing proportions
  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # 2. Precompute weighted sums
  Wx <- t(data) %*% gamma   # d × K

  # 3. Initialize mu_list and sigma_list
  mu_list    <- lapply(seq_len(K), function(k) Wx[,k] / Nk[k])
  sigma_list <- params$sigma

  # 4. Handle relative_lambda scaling if requested
  if (relative_lambda && !is.null(Q_prior)) {
    if (length(unique(sapply(sigma_list, function(S) all.equal(S, sigma_list[[1]])))) > 1) {
      stop("Relative lambda requires shared diagonal covariance (EEI/EII) across clusters.")
    }
    sigma_vec   <- diag(sigma_list[[1]])
    scale_vec   <- rep(1 / sqrt(sigma_vec), times = K)
    Sscale      <- Diagonal(x = scale_vec)
    Q_prior_orig <- Q_prior
    Q_prior      <- Sscale %*% Q_prior_orig %*% Sscale
  }

  # 5. Single‐shot if iterate_once or no penalization/model unsupported
  if (iterate_once || is.null(Q_prior) ||
      !(modelName %in% c("VVV","VII","EII","EEI"))) {
    # 5a. μ‐update if penalized
    if (!is.null(Q_prior) && modelName %in% c("VVV","VII","EII","EEI")) {
      blocks <- vector("list", K); rhs <- numeric(DK)
      for (k in seq_len(K)) {
        Sinv         <- params$invSigma[[k]]
        blocks[[k]]  <- Nk[k] * Sinv
        rhs[((k-1)*d+1):(k*d)] <- Sinv %*% Wx[,k]
      }
      Hmat <- bdiag(blocks) + Q_prior
      if (DK > 200) {
        fac    <- Cholesky(Hmat, LDL = FALSE)
        mu_vec <- solve(fac, rhs)
      } else {
        mu_vec <- solve(as.matrix(Hmat), rhs)
      }
      mu_list <- split(mu_vec, rep(seq_len(K), each = d))
    }
    # 5b. covariance update once
    sigma_list <- lapply(seq_len(K), function(k) {
      diff <- sweep(data, 2, mu_list[[k]])
      W    <- sqrt(gamma[,k])
      switch(modelName,
             VVV = crossprod(diff * W, diff * W) / Nk[k] + nugget * diag(d),
             VII = { lam <- sum((diff^2)*gamma[,k])/(Nk[k]*d); diag(lam,d) },
             EII = sum((diff^2)*gamma[,k]),
             sigma_list[[k]]  # for EEI handled below
      )
    })
    if (modelName == "EEI") {
      if(relative_lambda){
        GX  <- t(data) %*% gamma
        GX2 <- t(data^2) %*% gamma
        mu_mat <- do.call(cbind, mu_list)
        V <- rowSums( GX2 - 2 * (mu_mat * GX) + sweep(mu_mat^2, 2, Nk, "*") )

        U_vec <- unlist(mu_list)
        U_vec_d <- c()
        for (dd in seq_len(d)) {
          u_index <- seq(from = dd, to = length(U_vec), by = d)
          U_vec_d <- c(U_vec_d, as.numeric(crossprod(U_vec[u_index], Q_prior_orig[u_index,u_index] %*% U_vec[u_index])))
        }
        sharedSig <- (diag(V, d) + diag(U_vec_d, d)) / (n + K - rank_deficiency) + nugget
      }else{
        GX  <- t(data) %*% gamma
        GX2 <- t(data^2) %*% gamma
        mu_mat <- do.call(cbind, mu_list)
        V <- rowSums( GX2 - 2 * (mu_mat * GX) + sweep(mu_mat^2, 2, Nk, "*") )

        Vs <- V / n + nugget
        sharedSig <- diag(Vs, d)
      }
      sigma_list <- replicate(K, sharedSig, simplify = FALSE)
    }
    return(list(pi = pi_new, mu = mu_list, sigma = sigma_list))
  }

  # 6. Otherwise, run inner block‐coordinate until convergence
  for (inner in seq_len(max_inner)) {
    mu_old    <- do.call(c, mu_list)
    sigma_old <- unlist(sigma_list)

    # 6a. Recompute invSigma and logdet from current sigma_list
    invSigma_list <- vector("list", K)
    for (k in seq_len(K)) {
      S <- sigma_list[[k]]
      L <- chol(S)
      invSigma_list[[k]] <- chol2inv(L)
    }

    if (relative_lambda) {
      # Rescale Q_prior if relative_lambda is TRUE
      sigma_vec   <- diag(sigma_list[[1]])
      scale_vec   <- rep(1 / sqrt(sigma_vec), times = K)
      Sscale      <- Diagonal(x = scale_vec)
      Q_prior      <- Sscale %*% Q_prior_orig %*% Sscale
    }

    # 6b. μ‐update using freshly computed invSigma_list
    blocks <- vector("list", K); rhs <- numeric(DK)
    for (k in seq_len(K)) {
      Sinv        <- invSigma_list[[k]]
      blocks[[k]] <- Nk[k] * Sinv
      rhs[((k-1)*d+1):(k*d)] <- Sinv %*% Wx[,k]
    }
    Hmat <- bdiag(blocks) + Q_prior
    if (DK > 200) {
      fac    <- Cholesky(Hmat, LDL = FALSE)
      mu_vec <- solve(fac, rhs)
    } else {
      mu_vec <- solve(as.matrix(Hmat), rhs)
    }
    mu_list <- split(mu_vec, rep(seq_len(K), each = d))

    # 6c. Σ‐update given new mu_list
    if (modelName == "EEI") {
      if(relative_lambda){
        GX  <- t(data) %*% gamma
        GX2 <- t(data^2) %*% gamma
        mu_mat <- do.call(cbind, mu_list)
        V <- rowSums( GX2 - 2 * (mu_mat * GX) + sweep(mu_mat^2, 2, Nk, "*") )

        U_vec <- unlist(mu_list)
        U_vec_d <- c()
        for (dd in seq_len(d)) {
          u_index <- seq(from = dd, to = length(U_vec), by = d)
          U_vec_d <- c(U_vec_d, as.numeric(crossprod(U_vec[u_index], Q_prior_orig[u_index,u_index] %*% U_vec[u_index])))
        }
        sharedSig <- (diag(V, d) + diag(U_vec_d, d)) / (n + K - rank_deficiency) + nugget
      }else{
        GX  <- t(data) %*% gamma
        GX2 <- t(data^2) %*% gamma
        mu_mat <- do.call(cbind, mu_list)
        V <- rowSums( GX2 - 2 * (mu_mat * GX) + sweep(mu_mat^2, 2, Nk, "*") )

        Vs  <- V / n + nugget
        sharedSig  <- diag(Vs, d)
      }
      sigma_list <- replicate(K, sharedSig, simplify = FALSE)
    } else {
      for (k in seq_len(K)) {
        diff <- sweep(data, 2, mu_list[[k]])
        W    <- sqrt(gamma[,k])
        sigma_list[[k]] <- switch(modelName,
                                  VVV = crossprod(diff * W, diff * W) / Nk[k] + nugget * diag(d),
                                  VII = { lam <- sum((diff^2)*gamma[,k])/(Nk[k]*d); diag(lam,d) },
                                  EII = sum((diff^2)*gamma[,k]),
                                  sigma_list[[k]]
        )
      }
      if (modelName == "EII") {
        total <- sum(unlist(sigma_list))
        lam_sh <- total/(sum(Nk)*d)
        sigma_list <- replicate(K, diag(lam_sh,d), simplify = FALSE)
      }
    }

    # 6d. Check inner convergence
    mu_diff    <- max(abs(mu_vec - mu_old))
    sigma_diff <- max(abs(unlist(sigma_list) - sigma_old))
    if (verbose) {
      cat(sprintf("  MSTEP inner iter %d: Δμ=%.3e, ΔΣ=%.3e\n",
                  inner, mu_diff, sigma_diff))
    }
    if (max(mu_diff, sigma_diff) < tol_inner) break
  }

  list(pi    = pi_new,
       mu    = mu_list,
       sigma = sigma_list)
}


# 4. EM wrapper
EM_algorithm <- function(data, init_params,
                         Q_prior   = NULL,
                         iterate_once = TRUE,
                         max_inner = 10,
                         modelName = "VVV",
                         max_iter  = 100,
                         tol       = 1e-4,
                         inner_tol = 1e-6,
                         eigen_tol = 0,
                         rank_deficiency = 0, # should be the rank deficiency of entire (joint) Q_prior
                         nugget    = 0,
                         relative_lambda = FALSE,
                         random_times = 10,
                         verbose   = TRUE) {

  params       <- init_cov_cache_fast(init_params)
  K <- length(params$pi)
  loglik_trace <- numeric(max_iter)
  elbo_trace   <- numeric(max_iter)

  if (relative_lambda) {
    # store the original Q_prior for rescaling
    Q_prior_orig <- Q_prior
  }

  for (iter in seq_len(max_iter)) {

    # E-Step
    gamma      <- ESTEP(data, params)

    # M-Step
    # check if relative_lambda is TRUE and Q_prior is not NULL
    if (relative_lambda && !is.null(Q_prior)) {
      new_params <- MSTEP(data = data, gamma = gamma, params = params,
                          Q_prior = Q_prior_orig,
                          relative_lambda = relative_lambda,
                          modelName = modelName,
                          iterate_once = iterate_once,
                          max_inner = max_inner,
                          tol_inner = inner_tol,
                          rank_deficiency = rank_deficiency,
                          nugget = nugget,
                          verbose = verbose)
    }
    else{
      new_params <- MSTEP(data = data, gamma = gamma, params = params,
                          Q_prior = Q_prior,
                          modelName = modelName,
                          iterate_once = iterate_once,
                          max_inner = max_inner,
                          tol_inner = tol,
                          verbose = verbose)
    }
    last_params <- params
    params     <- init_cov_cache_fast(new_params)

    # Optional rescaling of Q_prior if relative_lambda is TRUE
    if (relative_lambda && !is.null(Q_prior)) {
      sigma_vec <- diag(params$sigma[[1]])  # EEI: shared diagonal
      if (length(unique(sapply(params$sigma, function(S) all.equal(S, params$sigma[[1]])))) > 1) {
        stop("Relative lambda requires shared diagonal covariance (EEI) or identical diagonal covariances (EII) across clusters.")
      }
      scale_vec <- rep(1 / sqrt(sigma_vec), times = K)  # correct cluster-major order
      Sscale <- Diagonal(x = scale_vec)
      # rescale the original Q_prior
      Q_prior <- Sscale %*% Q_prior_orig %*% Sscale
    }

    # Record log-likelihood and ELBO
    ll <- compute_log_joint_observed(data, params, Q_prior, eigen_tol, rank_deficiency)
    elbo <- compute_penalized_ELBO(data, gamma, params, Q_prior, eigen_tol, rank_deficiency)

    loglik_trace[iter] <- ll
    elbo_trace[iter]   <- elbo

    if (verbose) {
      cat(sprintf("Iteration %3d: penLogLik = %.6f, ELBO = %.6f\n", iter, ll, elbo))
    }
    if (iter > 1 && abs(elbo - elbo_trace[iter-1]) < tol) {
      # Check if elbo decreased
      if (elbo_trace[iter - 1] - elbo > 1e-10) {
        if (verbose) {
          cat(sprintf("ELBO decreased at iteration %d: %.6f -> %.6f\n", iter, elbo_trace[iter - 1], elbo))
        }
        loglik_trace <- loglik_trace[1:(iter - 1)]
        elbo_trace   <- elbo_trace[1:(iter - 1)]
        params <- last_params
        gamma  <- ESTEP(data, params)
        break
      }else{
        if (verbose) {
          cat(sprintf("Converged at iteration %d with ELBO %.6f\n", iter, elbo))
          cat(sprintf("Final pLoglik: %.6f\n", ll))
        }
        loglik_trace <- loglik_trace[1:iter]
        elbo_trace   <- elbo_trace[1:iter]
        break
      }
    }
  }

  list(params       = params,
       gamma        = gamma,
       elbo_trace   = elbo_trace,
       loglik_trace = loglik_trace)
}



