#' Initialize parameters for bgplvm
#'
#' @param Y    N x D data matrix
#' @param Q    Latent dimensionality
#' @param M    Number of inducing points
#' @param init "pca" (default), "random", or a named list of parameters
#' @return Named list with fields: mu_X, log_s_X, Z, log_ls, log_var, log_sigma2
#' @keywords internal
.init_params <- function(Y, Q, M, init = "pca") {
  N <- nrow(Y)

  if (is.list(init)) {
    # User-supplied warm start — validate fields
    required <- c("mu_X", "log_s_X", "Z", "log_ls", "log_var", "log_sigma2")
    missing_fields <- setdiff(required, names(init))
    if (length(missing_fields) > 0L) {
      stop("init list missing fields: ", paste(missing_fields, collapse = ", "))
    }
    return(init)
  }

  mu_X <- if (init == "pca") {
    pca <- prcomp(Y, center = TRUE, scale. = FALSE)
    pca$x[, seq_len(Q), drop = FALSE]
  } else if (init == "random") {
    matrix(rnorm(N * Q, sd = 0.1), N, Q)
  } else {
    stop("init must be 'pca', 'random', or a named list")
  }

  # Inducing points: subsample rows of mu_X
  idx_z <- sample.int(N, M)
  Z <- mu_X[idx_z, , drop = FALSE]

  list(
    mu_X       = mu_X,
    log_s_X    = matrix(log(0.1), N, Q),  # initial posterior std dev = 0.1
    Z          = Z,
    log_ls     = rep(0, Q),               # lengthscale = 1 per dim
    log_var    = 0,                        # signal variance = 1
    log_sigma2 = log(0.1)                 # noise variance = 0.1
  )
}

#' Pack all parameters into a flat numeric vector (for L-BFGS-B / grad checks)
#'
#' @param params Named list from .init_params
#' @param N,Q,M Integer dimensions
#' @return Numeric vector of length 2*N*Q + M*Q + Q + 2
#' @keywords internal
.pack_params <- function(params, N, Q, M) {
  c(
    as.vector(params$mu_X),     # N*Q
    as.vector(params$log_s_X),  # N*Q
    as.vector(params$Z),        # M*Q
    params$log_ls,              # Q
    params$log_var,             # 1
    params$log_sigma2           # 1
  )
}

#' Unpack a flat numeric vector back into a named parameter list
#'
#' @param theta Numeric vector (output of .pack_params)
#' @param N,Q,M Integer dimensions
#' @return Named list matching .init_params structure
#' @keywords internal
.unpack_params <- function(theta, N, Q, M) {
  idx <- 0L
  take <- function(n) {
    out <- theta[idx + seq_len(n)]
    idx <<- idx + n
    out
  }
  list(
    mu_X       = matrix(take(N * Q), N, Q),
    log_s_X    = matrix(take(N * Q), N, Q),
    Z          = matrix(take(M * Q), M, Q),
    log_ls     = take(Q),
    log_var    = take(1L),
    log_sigma2 = take(1L)
  )
}

#' Add diagonal jitter to a matrix in-place
#'
#' @param K      Square numeric matrix
#' @param jitter Scalar to add to diagonal
#' @return Modified matrix
#' @keywords internal
.add_jitter <- function(K, jitter) {
  diag(K) <- diag(K) + jitter
  K
}

#' Cholesky with automatic jitter fallback
#'
#' Attempts chol(K + jitter * I), doubling jitter on failure.
#'
#' @param K             Square positive (semi-)definite matrix
#' @param jitter_start  Initial jitter (default 1e-6)
#' @param max_attempts  Maximum doubling attempts (default 6)
#' @return list(L = lower triangular Cholesky, jitter = final jitter used)
#' @keywords internal
.chol_with_jitter <- function(K, jitter_start = 1e-6, max_attempts = 6L) {
  jitter <- jitter_start
  for (i in seq_len(max_attempts)) {
    Kj <- .add_jitter(K, jitter)
    L <- tryCatch(
      t(chol(Kj)),  # lower triangular
      error = function(e) NULL
    )
    if (!is.null(L)) return(list(L = L, jitter = jitter))
    jitter <- jitter * 10
  }
  stop(
    "Cholesky failed after ", max_attempts, " attempts with final jitter ",
    jitter, ". Matrix may not be positive definite."
  )
}

#' Clamp parameters to numerically safe ranges (in-place update)
#'
#' Prevents variance collapse or explosion during Adam optimization.
#'
#' @param params Named list of parameters
#' @return Same list with values clamped
#' @keywords internal
.clamp_params <- function(params) {
  params$log_s_X    <- pmax(pmin(params$log_s_X,    2),  -10)
  params$log_var    <- pmax(pmin(params$log_var,     5),   -5)
  params$log_sigma2 <- pmax(pmin(params$log_sigma2,  2),  -10)
  params$log_ls     <- pmax(pmin(params$log_ls,      3),   -3)
  params
}
