# Compute state-level (occupancy) entropy from a responsibility matrix Gamma
# Gamma: n x K matrix, rows sum ~ 1 (responsibilities p(z_i = k | x_i))
# base:  log base ("e", "2", or "10")
# normalize: if TRUE, return H/log(K) in addition to H (so it's in [0,1])
# eps: small positive number for numerical stability in log
compute_state_entropy <- function(Gamma,
                                  base = c("e", "2", "10"),
                                  normalize = TRUE,
                                  eps = 1e-12) {
  base <- match.arg(base)
  if (!is.matrix(Gamma)) stop("Gamma must be a matrix.")
  if (ncol(Gamma) < 1L || nrow(Gamma) < 1L) stop("Gamma must be n x K with n, K >= 1.")

  # Row-normalize defensively (in case of tiny row sums)
  rs <- rowSums(Gamma)
  rs[rs < eps] <- eps
  G <- sweep(Gamma, 1L, rs, "/")

  # State occupancy p_k = average responsibility per state
  pk <- colMeans(G)
  pk <- pmax(pk, 0)                     # clip negatives due to numeric issues
  pk <- pk / sum(pk + eps)              # normalize to sum 1

  # Choose log base
  logfun <- switch(base,
                   "e"  = log,
                   "2"  = function(x) log(x) / log(2),
                   "10" = function(x) log(x, base = 10))

  # State-level entropy H_state = - sum_k p_k log p_k
  H <- -sum(pk * logfun(pk + eps))

  # Effective number of states
  Neff <- switch(base,
                 "e"  = exp(H),
                 "2"  = 2^H,
                 "10" = 10^H)

  out <- list(
    H_state = H,
    H_state_normalized = if (normalize) H / logfun(ncol(G)) else NULL,
    Neff = Neff,
    pk = pk
  )
  return(out)
}





compute_entropy_summary <- function(X,
                                    K = 50,
                                    lambda = 500,
                                    q = 2,
                                    modelName = "EEI",
                                    relative_lambda = TRUE,
                                    rank_deficiency = 2,
                                    iterate_once = TRUE,
                                    tol = 1e-3,
                                    max_iter = 100,
                                    verbose = FALSE) {
  # Construct RW prior
  Q_prior_RW <- make_random_walk_precision(K = K, d = 1, lambda = lambda, q = q)

  max_vec <- numeric(0)
  mean_vec <- numeric(0)

  # Loop over columns of X
  for (i in seq_len(ncol(X))) {
    selected_X <- X[, i, drop = FALSE]
    init_params <- make_default_init(selected_X, K = K)

    result_RW <- EM_algorithm(
      data = selected_X,
      Q_prior = Q_prior_RW,
      init_params = init_params,
      max_iter = max_iter,
      modelName = modelName,
      relative_lambda = relative_lambda,
      rank_deficiency = rank_deficiency,
      iterate_once = iterate_once,
      tol = tol,
      verbose = verbose
    )

    gamma <- result_RW$gamma
    entropy_vec <- -rowSums(gamma * log(gamma + 1e-12))

    max_vec <- c(max_vec, max(entropy_vec))
    mean_vec <- c(mean_vec, mean(entropy_vec))
  }

  return(list(max_vec = max_vec, mean_vec = mean_vec))
}




compute_min_slope_smooth <- function(x, sort_x = TRUE, n_grid = 200L, q = 0.2) {
  stopifnot(is.numeric(x))
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L) return(NA_real_)

  # scale sorted data to [0,1]
  x <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) + 1e-12)

  # sort data
  if (sort_x) {
    x <- sort(x, na.last = NA)
    t <- seq_len(n)                # treat order index as "latent z"
  } else {
    t <- seq_len(n)
  }

  # spline interpolator (natural cubic spline)
  spline_fit <- splinefun(t, x, method = "natural")

  # evaluate derivative on a regular grid
  grid <- seq(min(t), max(t), length.out = n_grid)
  slope_est <- spline_fit(grid, deriv = 1)

  # return the 10% quantile of absolute slope estimates
  return(quantile(abs(slope_est), probs = q, na.rm = TRUE))
  # # return minimum absolute slope (avoid NaN)
  # min(abs(slope_est), na.rm = TRUE)
}
