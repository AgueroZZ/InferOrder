simulate_2d_processes <- function(f = NULL,
                                  nrow         = 50,
                                  ncol         = 50,
                                  grid         = NULL,
                                  d            = 1,
                                  epsilon_sd   = 1,
                                  seed         = 123,
                                  matern_nu    = 2,
                                  matern_scale = 0.3,
                                  matern_var   = 1) {
  if (!is.null(seed)) set.seed(seed)

  # 1. build or accept the grid ----------------------------------
  if (is.null(grid)) {
    xs   <- seq(0, 1, length.out = nrow)
    ys   <- seq(0, 1, length.out = ncol)
    grid <- expand.grid(x = xs, y = ys)
  } else {
    # accept data.frame or matrix with x,y columns
    if (!all(c("x","y") %in% colnames(grid))) {
      stop("`grid` must have columns named 'x' and 'y'.")
    }
    grid <- as.data.frame(grid)[,c("x","y")]
  }
  n <- nrow(grid)

  # 2. generate true surfaces Fmat (n × d) -----------------------
  if (is.null(f)) {
    # simulate d independent Matérn fields via fields + MASS
    if (!requireNamespace("fields", quietly = TRUE) ||
        !requireNamespace("MASS",   quietly = TRUE)) {
      stop("Please install.packages(c('fields','MASS')) to simulate Matérn fields.")
    }
    locs    <- as.matrix(grid)
    distmat <- fields::rdist(locs)
    C <- matern_var * fields::Matern(distmat,
                                     range      = matern_scale,
                                     smoothness = matern_nu)
    samples <- MASS::mvrnorm(n = d,
                             mu = rep(0, n),
                             Sigma = C)
    Fmat <- t(samples)  # now n × d
  } else {
    # user-specified function or list of functions
    if (is.function(f)) {
      f_list <- rep(list(f), d)
    } else if (is.list(f) && length(f) == d && all(vapply(f, is.function, FALSE))) {
      f_list <- f
    } else {
      stop("`f` must be NULL, a single function, or a list of `d` functions.")
    }
    Fmat <- matrix(NA_real_, nrow = n, ncol = d)
    for (k in seq_len(d)) {
      vals <- with(grid, f_list[[k]](x, y))
      if (length(vals) == 1) {
        Fmat[, k] <- rep(vals, n)
      } else if (length(vals) == n) {
        Fmat[, k] <- vals
      } else {
        stop(sprintf("Function %d returned %d values, but grid has %d locations.",
                     k, length(vals), n))
      }
    }
  }

  # 3. add independent Gaussian noise -----------------------------
  Y <- matrix(NA_real_, nrow = n, ncol = d)
  for (k in seq_len(d)) {
    Y[, k] <- Fmat[, k] + rnorm(n, mean = 0, sd = epsilon_sd)
  }
  colnames(Y) <- paste0("proc", seq_len(d))

  # 4. return results ---------------------------------------------
  list(
    grid = grid,
    Y    = Y
  )
}
