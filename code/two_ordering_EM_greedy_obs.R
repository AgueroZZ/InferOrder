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

compute_mu_posterior <- function(xj, Gamma, Q_1D) {
  K <- ncol(Gamma)
  Nk <- colSums(Gamma)
  rhs <- colSums(Gamma * xj)
  H <- diag(Nk) + Q_1D
  solve(H, rhs)
}

compute_posterior_energy <- function(xj, Gamma, mu_jk, Q_1D, scaled = TRUE) {

  # xj: length-N vector
  # Gamma: N x K matrix
  # mu_jk: length-K vector
  # Q_1D: K x K prior precision matrix

  pred  <- Gamma %*% mu_jk
  resid <- xj - pred
  sigma2_hat <- mean(resid^2)  # MLE estimate of residual variance

  # Data fit term
  pred_mat <- matrix(mu_jk, nrow = nrow(Gamma), ncol = length(mu_jk), byrow = TRUE)
  resid2   <- (xj - pred_mat)^2
  data_fit <- sum(Gamma * resid2) / (2 * sigma2_hat)

  # Normalization constant of observation model
  norm_obs <- sum(Gamma) * log(2 * pi * sigma2_hat) / 2

  if (scaled) {
    # Scale the prior precision matrix
    Q_1D <- Q_1D / sigma2_hat
  }

  # Prior penalty
  penalty_prior <- as.numeric(0.5 * crossprod(mu_jk, Q_1D %*% mu_jk))

  # Prior normalization constant
  prior_norm <- (length(mu_jk) / 2) * log(2 * pi) - 0.5 * as.numeric(generalized_logdet(Q_1D))

  # Total posterior negative log-density
  total <- data_fit + norm_obs + penalty_prior + prior_norm
  return(total)
}

compute_residual <- function(xj, mu_list, Gamma, j_local = NULL) {
  K <- length(mu_list)
  n <- length(xj)
  mu_jk <- if (!is.null(j_local)) {
    sapply(mu_list, function(mu_k) mu_k[j_local])
  } else {
    approx_mu_for_j(xj, Gamma)
  }
  pred <- Gamma %*% mu_jk
  mean((xj - pred)^2)
}

approx_mu_for_j <- function(xj, Gamma) {
  Nk <- colSums(Gamma)
  colSums(xj * Gamma) / Nk
}

two_ordering_smoothEM <- function(
    X,
    K = 5,
    max_iter_backfitting = 20,
    max_iter_EM = 100,
    init_method = c("random", "pca", "tsne"),
    greedy = FALSE,
    reassign_by = c("likelihood", "posterior"),
    make_Q_prior = NULL,
    lambda_vec = c(10, 10),
    rw_order_vec = c(2, 2),
    relative_lambda = TRUE,
    parallel = FALSE,
    ncores = 2,
    verbose = TRUE,
    init_coord_assign = NULL,
    em_args = list()
) {
  # Select initialization method and reassignment criterion
  init_method <- match.arg(init_method)
  reassign_by  <- match.arg(reassign_by)
  n <- nrow(X); d <- ncol(X)

  # Helper function for generating initial EM parameters
  get_init <- function(Xj) {
    if (init_method == "random") {
      make_default_init(Xj, K)
    } else if (init_method == "pca") {
      pc <- prcomp(Xj, center = TRUE, scale. = TRUE)$x[, 1]
      make_init(Xj, ordering_vec = pc, K = K)
    } else if (init_method == "tsne") {
      tsne <- Rtsne::Rtsne(
        Xj, dims = 1, perplexity = 50,
        verbose = FALSE, check_duplicates = FALSE,
        theta = 0.5
      )$Y[, 1]
      make_init(Xj, ordering_vec = tsne, K = K)
    } else {
      make_default_init(Xj, K)
    }
  }

  # Initialize coordinate assignment if not using greedy
  coord_assign <- if (is.null(init_coord_assign)) {
    sample(1:2, d, replace = TRUE)
  } else {
    init_coord_assign
  }

  ######
  # GREEDY INITIALIZATION
  ######
  if (greedy) {
    # 1) Choose the most variable coordinate as seed for ordering 1
    var_vec <- apply(X, 2, var)
    j1 <- which.max(var_vec)
    Q1 <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1])
    init1 <- get_init(X[, j1, drop = FALSE])
    Gamma1 <- ESTEP(
      X[, j1, drop = FALSE],
      init_cov_cache_fast(init1)
    )

    # 2) Choose the coordinate with highest posterior energy under ordering 1 as seed for ordering 2
    mu1_seed <- compute_mu_posterior(X[, j1], Gamma1, Q1)
    energies <- sapply(setdiff(1:d, j1), function(j) {
      compute_posterior_energy(
        X[, j], Gamma1, mu1_seed, Q1,
        scaled = relative_lambda
      )
    })
    j2 <- setdiff(1:d, j1)[which.max(energies)]
    Q2 <- make_Q_prior(
      K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2]
    )
    init2 <- get_init(X[, j2, drop = FALSE])
    Gamma2 <- ESTEP(
      X[, j2, drop = FALSE],
      init_cov_cache_fast(init2)
    )

    if (verbose) {
      cat(sprintf(
        "[Greedy Init] seed coordinates: j1 = %d, j2 = %d\n",
        j1, j2
      ))
    }

    # Initialize assignments for seeds
    coord_assign <- rep(NA_integer_, d)
    coord_assign[j1] <- 1
    coord_assign[j2] <- 2

    # 3) Greedily add remaining coordinates and update Gamma1/Gamma2 after each addition
    remaining <- setdiff(1:d, c(j1, j2))
    for (iter_g in seq_along(remaining)) {
      j <- remaining[iter_g]
      xj <- X[, j]

      # Compute posterior energy for both orderings
      mu1j <- compute_mu_posterior(xj, Gamma1, Q1)
      E1 <- compute_posterior_energy(
        xj, Gamma1, mu1j, Q1,
        scaled = relative_lambda
      )
      mu2j <- compute_mu_posterior(xj, Gamma2, Q2)
      E2 <- compute_posterior_energy(
        xj, Gamma2, mu2j, Q2,
        scaled = relative_lambda
      )

      if (E1 < E2) {
        coord_assign[j] <- 1
        J1_cur <- which(coord_assign == 1)
        init1_cur <- get_init(X[, J1_cur, drop = FALSE])
        Gamma1 <- ESTEP(
          X[, J1_cur, drop = FALSE],
          init_cov_cache_fast(init1_cur)
        )
        if (verbose) cat(sprintf(
          "[Greedy Iter %d] coord %d assigned to ordering 1 (E1=%.3f, E2=%.3f)\n",
          iter_g, j, E1, E2
        ))
      } else {
        coord_assign[j] <- 2
        J2_cur <- which(coord_assign == 2)
        init2_cur <- get_init(X[, J2_cur, drop = FALSE])
        Gamma2 <- ESTEP(
          X[, J2_cur, drop = FALSE],
          init_cov_cache_fast(init2_cur)
        )
        if (verbose) cat(sprintf(
          "[Greedy Iter %d] coord %d assigned to ordering 2 (E1=%.3f, E2=%.3f)\n",
          iter_g, j, E1, E2
        ))
      }
    }

    if (verbose) {
      cat(sprintf(
        "[Greedy Init Completed] ordering1 size = %d, ordering2 size = %d\n",
        sum(coord_assign == 1), sum(coord_assign == 2)
      ))
    }
  }

  # If only greedy init with no backfitting, return directly
  if (greedy && max_iter_backfitting == 0) {
    return(list(
      coord_assign = coord_assign,
      J1 = which(coord_assign == 1),
      J2 = which(coord_assign == 2)
    ))
  }

  ######
  # BACKFITTING LOOP
  ######
  res1 <- res2 <- NULL
  for (iter in seq_len(max_iter_backfitting)) {
    J1 <- which(coord_assign == 1)
    J2 <- which(coord_assign == 2)
    X1 <- X[, J1, drop = FALSE]
    X2 <- X[, J2, drop = FALSE]

    # Refit ordering 1
    if (length(J1) > 0) {
      Q1 <- make_Q_prior(
        K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1]
      )
      res1 <- list()
      res1$gamma <- ESTEP(
        X1, init_cov_cache_fast(get_init(X1))
      )
      res1$params <- MSTEP_coordinatewise(
        X1, res1$gamma, Q1, parallel = parallel, ncores = ncores
      )
    }
    # Refit ordering 2
    if (length(J2) > 0) {
      Q2 <- make_Q_prior(
        K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2]
      )
      res2 <- list()
      res2$gamma <- ESTEP(
        X2, init_cov_cache_fast(get_init(X2))
      )
      res2$params <- MSTEP_coordinatewise(
        X2, res2$gamma, Q2, parallel = parallel, ncores = ncores
      )
    }

    # Reassign all coordinates
    new_assign <- coord_assign
    for (j in seq_len(d)) {
      xj <- X[, j]
      if (reassign_by == "likelihood") {
        r1 <- if (!is.null(res1))
          compute_residual(xj, res1$params$mu, res1$gamma,
                           j_local = which(J1 == j)) else Inf
        r2 <- if (!is.null(res2))
          compute_residual(xj, res2$params$mu, res2$gamma,
                           j_local = which(J2 == j)) else Inf
        new_assign[j] <- if (r1 < r2) 1 else 2
      } else {
        Q1_j <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1])
        mu1_j <- compute_mu_posterior(xj, res1$gamma, Q1_j)
        E1 <- compute_posterior_energy(
          xj, res1$gamma, mu1_j, Q1_j,
          scaled = relative_lambda
        )
        Q2_j <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2])
        mu2_j <- compute_mu_posterior(xj, res2$gamma, Q2_j)
        E2 <- compute_posterior_energy(
          xj, res2$gamma, mu2_j, Q2_j,
          scaled = relative_lambda
        )
        new_assign[j] <- if (E1 < E2) 1 else 2
      }
    }

    reassign_count <- sum(new_assign != coord_assign)
    if (verbose) {
      cat(sprintf(
        "[Backfitting %2d] reassigned %d coords (J1=%d, J2=%d)\n",
        iter, reassign_count,
        sum(coord_assign == 1), sum(coord_assign == 2)
      ))
    }
    if (reassign_count == 0) break
    coord_assign <- new_assign
  }

  ######
  # FINAL FIT
  ######
  J1 <- which(coord_assign == 1)
  J2 <- which(coord_assign == 2)
  X1 <- X[, J1, drop = FALSE]
  X2 <- X[, J2, drop = FALSE]

  res1 <- res2 <- NULL
  if (length(J1) > 0) {
    Q1 <- make_Q_prior(
      K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1]
    )
    res1 <- list()
    res1$gamma <- ESTEP(
      X1, init_cov_cache_fast(get_init(X1))
    )
    res1$params <- MSTEP_coordinatewise(
      X1, res1$gamma, Q1, parallel = parallel, ncores = ncores
    )
  }
  if (length(J2) > 0) {
    Q2 <- make_Q_prior(
      K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2]
    )
    res2 <- list()
    res2$gamma <- ESTEP(
      X2, init_cov_cache_fast(get_init(X2))
    )
    res2$params <- MSTEP_coordinatewise(
      X2, res2$gamma, Q2, parallel = parallel, ncores = ncores
    )
  }

  return(list(
    coord_assign = coord_assign,
    J1           = J1,
    J2           = J2,
    result1      = res1,
    result2      = res2,
    Gamma1       = if (!is.null(res1)) res1$gamma else NULL,
    Gamma2       = if (!is.null(res2)) res2$gamma else NULL,
    mu1_list     = if (!is.null(res1)) res1$params$mu else NULL,
    mu2_list     = if (!is.null(res2)) res2$params$mu else NULL
  ))
}
