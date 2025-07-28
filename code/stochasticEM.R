make_default_init <- function(X, K) {
  d <- ncol(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  list(
    pi    = rep(1/K, K),
    mu    = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = rep(1, d)
  )
}


ESTEP <- function(X, U_list, Sigma_diag, pi_vec = NULL) {
  n <- nrow(X)
  d <- ncol(X)
  K <- length(U_list)
  log_gamma <- matrix(NA, n, K)

  if (is.null(pi_vec)) pi_vec <- rep(1/K, K)
  Sinv <- diag(1 / Sigma_diag)
  logdet <- sum(log(Sigma_diag))

  for (k in seq_len(K)) {
    Xc <- sweep(X, 2, U_list[[k]])
    Qf <- rowSums((Xc %*% Sinv) * Xc)
    log_gamma[, k] <- log(pi_vec[k]) - 0.5 * (Qf + logdet + d * log(2 * pi))
  }

  lse <- matrixStats::rowLogSumExps(log_gamma)
  Gamma <- exp(log_gamma - lse)
  Gamma
}


SSTEP <- function(Gamma, threshold = c("random", "hard")) {
  threshold <- match.arg(threshold)
  n <- nrow(Gamma)
  K <- ncol(Gamma)
  Z <- matrix(0, n, K)

  if (threshold == "random") {
    for (i in 1:n) {
      Z[i, sample(1:K, 1, prob = Gamma[i, ])] <- 1
    }
  } else if (threshold == "hard") {
    idx <- max.col(Gamma)
    Z[cbind(1:n, idx)] <- 1
  }
  Z
}


smooth_MSTEP <- function(data, gamma, params, Q_prior = NULL) {
  n <- nrow(data)
  d <- ncol(data)
  K <- ncol(gamma)
  DK <- d * K

  # Mixture weights
  Nk <- colSums(gamma)
  Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Weighted sums
  Wx <- t(data) %*% gamma

  # Shared Σ from previous step
  Sinv <- diag(1 / params$sigma)

  # Build linear system
  blocks <- vector("list", K)
  rhs <- numeric(DK)
  for (k in seq_len(K)) {
    blocks[[k]] <- Nk[k] * Sinv
    rhs[((k - 1) * d + 1):(k * d)] <- Sinv %*% Wx[, k]
  }

  if (!is.null(Q_prior)) {
    H <- as.matrix(bdiag(blocks)) + as.matrix(Q_prior)
  } else {
    H <- as.matrix(bdiag(blocks))
  }

  mu_vec <- solve(H, rhs)
  mu_list <- split(mu_vec, rep(seq_len(K), each = d))

  # Update shared diagonal covariance
  V <- numeric(d)
  for (k in seq_len(K)) {
    diff <- sweep(data, 2, mu_list[[k]])
    V <- V + colSums(gamma[, k] * (diff^2))
  }
  V_shared <- V / n + 1e-6

  list(
    pi = pi_new,
    mu = mu_list,
    sigma = V_shared
  )
}


MSTEP <- function(X, params = NULL, Gamma = NULL, Z = NULL,
                  smoother = NULL, Q_prior = NULL) {
  if (!is.null(Gamma)) {
    if (is.null(params)) stop("Need previous params for smooth_MSTEP!")
    return(smooth_MSTEP(data = X, gamma = Gamma, params = params, Q_prior = Q_prior))
  }
  if (!is.null(Z)) {
    if (is.null(smoother)){
      warning("No smoother provided, using the Q-matrix provided for smoothing.")
      if (is.null(Q_prior)){
        stop("Either provide a smoother or a Q_prior for smoothing.")
      }
      return(smooth_MSTEP(data = X, gamma = Z, params = params, Q_prior = Q_prior))
    }
    else{
      return(smoother(Z, X, params))
    }
  }
  stop("Must provide either Gamma or Z.")
}

stochastic_EM <- function(X,
                          init_params,
                          Q_prior = NULL,
                          smoother = NULL,
                          threshold = NULL,
                          max_iter = 100,
                          tol = 1e-4,
                          verbose = TRUE) {
  # Initialize
  params <- init_params
  mu_current <- params$mu
  sigma_current <- diag(params$sigma)

  # Trace of relative changes
  rel_change_mu_trace <- numeric()
  rel_change_sigma_trace <- numeric()
  rel_change_gamma_trace <- numeric()

  for (iter in 1:max_iter) {
    ## 1️⃣ E-step
    Sigma_diag <- params$sigma
    Gamma <- ESTEP(X, params$mu, Sigma_diag, pi_vec = params$pi)

    ## 2️⃣ S-step
    if (!is.null(threshold)) {
      Z <- SSTEP(Gamma, threshold)
      if (is.null(smoother)) {
        warning("No smoother provided, using the Q-matrix provided for smoothing.")
        res <- MSTEP(X, params = params, Gamma = Z, Q_prior = Q_prior)
        if (is.null(Q_prior)) {
          stop("Either provide a smoother or a Q_prior for smoothing.")
        }
      }else{
        res <- MSTEP(X, Z = Z, smoother = smoother)
      }
    } else {
      res <- MSTEP(X, params = params, Gamma = Gamma, Q_prior = Q_prior)
    }

    ## 3️⃣ Compute relative changes
    if (iter > 1) {
      # mu
      delta_U_num <- sqrt(sum((unlist(res$mu) - unlist(mu_current))^2))
      delta_U_den <- sqrt(sum(unlist(mu_current)^2)) + 1e-8
      rel_delta_U <- delta_U_num / delta_U_den

      # sigma
      sigma_new <- diag(res$sigma)
      delta_Sigma_num <- sqrt(sum((sigma_new - sigma_current)^2))
      delta_Sigma_den <- sqrt(sum(sigma_current^2)) + 1e-8
      rel_delta_Sigma <- delta_Sigma_num / delta_Sigma_den

      # gamma
      delta_Gamma_num <- sqrt(sum((Gamma - Gamma_prev)^2))
      delta_Gamma_den <- sqrt(sum(Gamma_prev^2)) + 1e-8
      rel_delta_Gamma <- delta_Gamma_num / delta_Gamma_den

      # Store
      rel_change_mu_trace <- c(rel_change_mu_trace, rel_delta_U)
      rel_change_sigma_trace <- c(rel_change_sigma_trace, rel_delta_Sigma)
      rel_change_gamma_trace <- c(rel_change_gamma_trace, rel_delta_Gamma)

      if (verbose) {
        cat(sprintf("Iteration %d:\n", iter))
        cat(sprintf("  relΔU     = %.6e\n", rel_delta_U))
        cat(sprintf("  relΔSigma = %.6e\n", rel_delta_Sigma))
        cat(sprintf("  relΔGamma = %.6e\n", rel_delta_Gamma))
      }

      # Convergence check
      if (rel_delta_U < tol && rel_delta_Sigma < tol && rel_delta_Gamma < tol) {
        cat("Converged.\n")
        break
      }
    } else if (verbose) {
      cat(sprintf("Iteration %d (initialization) \n", iter))
    }

    ## 4️⃣ Update for next iteration
    params <- res
    mu_current <- params$mu
    sigma_current <- diag(params$sigma)
    Gamma_prev <- Gamma
  }

  ## Return
  list(
    params = params,
    Gamma = Gamma,
    rel_change_mu_trace = rel_change_mu_trace,
    rel_change_sigma_trace = rel_change_sigma_trace,
    rel_change_gamma_trace = rel_change_gamma_trace
  )
}




