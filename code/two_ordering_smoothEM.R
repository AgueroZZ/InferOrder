two_ordering_smoothEM <- function(X, K = 5,
                                  max_iter_EM = 100,
                                  max_iter = 20,
                                  modelName = "EEI",
                                  init_method = c("random", "pca", "warm_start", "tsne"),
                                  make_Q_prior = NULL,
                                  lambda_vec = c(10, 10),
                                  rw_order_vec = c(2, 2),
                                  reassign_by = c("likelihood", "posterior"),
                                  init_coord_assign = NULL,
                                  relative_lambda = TRUE,
                                  em_args = list(), verbose = TRUE) {

  init_method <- match.arg(init_method)

  reassign_by <- match.arg(reassign_by)
  n <- nrow(X)
  d <- ncol(X)

  coord_assign <- if (is.null(init_coord_assign)) sample(1:2, d, replace = TRUE) else init_coord_assign

  res1 <- res2 <- NULL  # initialize to allow reuse if skipped

  for (iter in 1:max_iter) {
    J1 <- which(coord_assign == 1)
    J2 <- which(coord_assign == 2)
    X1 <- X[, J1, drop = FALSE]
    X2 <- X[, J2, drop = FALSE]

    if(init_method == "random"){
      initial_params1 <- make_default_init(X1, K)
      initial_params2 <- make_default_init(X2, K)
    }else if(init_method == "pca"){
      pc1_metric1 <- prcomp(X1,center = TRUE, scale. = TRUE)$x[,1]
      pc1_metric2 <- prcomp(X2,center = TRUE, scale. = TRUE)$x[,1]
      initial_params1 <- make_init(X = X1, ordering_vec = pc1_metric1, K = K)
      initial_params2 <- make_init(X = X2, ordering_vec = pc1_metric2, K = K)
    }else if(init_method == "warm_start"){
      # not yet implemented
      restop("Warm start initialization not implemented yet.")
    }else if(init_method == "tsne"){
      tsne <- Rtsne(X1, dims = 1, perplexity = 50, verbose = F, check_duplicates = FALSE, theta = 0.5, pca = TRUE, ncol = 1)
      tsne_metric1 <- tsne$Y[,1]
      tsne <- Rtsne(X2, dims = 1, perplexity = 50, verbose = F, check_duplicates = FALSE, theta = 0.5, pca = TRUE, ncol = 1)
      tsne_metric2 <- tsne$Y[,1]
      initial_params1 <- make_init(X = X1, ordering_vec = tsne_metric1, K = K)
      initial_params2 <- make_init(X = X2, ordering_vec = tsne_metric2, K = K)
    }else{
      stop("Unknown initialization method: ", init_method)
    }

    if (length(J1) > 0) {
      Q1 <- make_Q_prior(K = K, d = length(J1), lambda = lambda_vec[1], q = rw_order_vec[1])
      res1 <- do.call(EM_algorithm, c(list(
        data = X1,
        max_iter = max_iter_EM,
        init_params = initial_params1,
        Q_prior = Q1,
        relative_lambda = relative_lambda,
        modelName = modelName,
        verbose = FALSE
      ), em_args))
    } else {
      if (verbose) cat("Ordering 1 skipped (empty)\n")
    }

    if (length(J2) > 0) {
      Q2 <- make_Q_prior(K = K, d = length(J2), lambda = lambda_vec[2], q = rw_order_vec[2])
      res2 <- do.call(EM_algorithm, c(list(
        data = X2,
        max_iter = max_iter_EM,
        init_params = initial_params2,
        Q_prior = Q2,
        relative_lambda = relative_lambda,
        modelName = modelName,
        verbose = FALSE
      ), em_args))
    } else {
      if (verbose) cat("Ordering 2 skipped (empty)\n")
    }

    new_assign <- coord_assign
    for (j in 1:d) {
      xj <- X[, j]

      if (reassign_by == "likelihood") {
        j1 <- which(J1 == j)
        j2 <- which(J2 == j)
        r1 <- if (length(J1) > 0) compute_residual(xj, res1$params$mu, res1$gamma, j_local = if (length(j1) > 0) j1 else NULL) else Inf
        r2 <- if (length(J2) > 0) compute_residual(xj, res2$params$mu, res2$gamma, j_local = if (length(j2) > 0) j2 else NULL) else Inf
        new_assign[j] <- if (r1 < r2) 1 else 2
      } else {
        Q1_j <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1])
        mu1_j <- compute_mu_posterior(xj, res1$gamma, Q1_j)
        E1 <- compute_posterior_energy(xj, res1$gamma, mu1_j, Q1_j, scaled = relative_lambda)

        Q2_j <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2])
        mu2_j <- compute_mu_posterior(xj, res2$gamma, Q2_j)
        E2 <- compute_posterior_energy(xj, res2$gamma, mu2_j, Q2_j, scaled = relative_lambda)
        new_assign[j] <- if (E1 < E2) 1 else 2
      }
    }

    reassign_count <- sum(new_assign != coord_assign)
    if (verbose) cat(sprintf("Iteration %2d: reassigned %d coordinates (J1 = %d, J2 = %d)\n",
                             iter, reassign_count, length(J1), length(J2)))
    if (reassign_count == 0) break
    coord_assign <- new_assign
  }

  list(
    coord_assign = coord_assign,
    J1 = J1, J2 = J2,
    result1 = res1,
    result2 = res2,
    Gamma1 = if (!is.null(res1)) res1$gamma else NULL,
    Gamma2 = if (!is.null(res2)) res2$gamma else NULL,
    mu1_list = if (!is.null(res1)) res1$params$mu else NULL,
    mu2_list = if (!is.null(res2)) res2$params$mu else NULL
  )
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


