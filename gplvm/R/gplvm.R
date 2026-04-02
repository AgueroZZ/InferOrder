#' Fit a Sparse Bayesian GPLVM
#'
#' Fits a sparse Bayesian Gaussian Process Latent Variable Model (BGPLVM) using
#' variational inference with inducing points (Titsias 2009 / Titsias & Lawrence 2010).
#'
#' The model:
#'   - Prior: p(X) = prod_i N(x_i; 0, I_Q)
#'   - Likelihood: p(Y | X) = prod_d GP(y_d; 0, K(X,X) + sigma^2 I)
#'   - Variational posterior: q(X) = prod_i N(x_i; mu_i, diag(s_i^2))
#'   - Sparse approximation via M inducing inputs Z (M << N)
#'
#' Scales to N > 5000 cells at O(NM^2) per ELBO evaluation.
#'
#' @param Y         N x D data matrix (e.g., NMF or PCA loadings). Rows = cells.
#' @param Q         Latent dimensionality. Use Q=2 for visualization, Q=1 for pseudo-time.
#' @param M         Number of inducing points. Default: min(100, ceiling(sqrt(N))).
#' @param n_iter    Maximum optimization iterations (default 1000).
#' @param lr        Adam learning rate (default 0.01).
#' @param optimizer "adam" (default) or "lbfgsb".
#' @param init      "pca" (default), "random", or a named list of initial parameters.
#' @param verbose   Print ELBO every this many iterations; 0 = silent (default 100).
#' @param tol_elbo  Early stopping tolerance (default 1e-4).
#' @param patience  Early stopping patience in iterations (default 20).
#' @param ...       Additional arguments passed to the optimizer.
#'
#' @return An object of class \code{"bgplvm"} with components:
#'   \item{mu_X}{N x Q matrix of posterior latent means}
#'   \item{S_X}{N x Q matrix of posterior latent variances}
#'   \item{Z}{M x Q matrix of learned inducing inputs}
#'   \item{lengthscales}{Length-Q vector of ARD kernel lengthscales}
#'   \item{variance}{Signal variance (scalar)}
#'   \item{sigma2}{Noise variance (scalar)}
#'   \item{K_MM}{M x M inducing kernel matrix at convergence}
#'   \item{K_MM_inv}{M x M inverse of K_MM}
#'   \item{elbo_trace}{ELBO values per iteration}
#'   \item{N, D, Q, M}{Dimensions}
#'   \item{converged}{Logical; TRUE if early stopping criterion was met}
#'   \item{n_iter_run}{Number of iterations actually run}
#'
#' @references
#'   Titsias, M. (2009). Variational learning of inducing variables in sparse
#'   Gaussian processes. AISTATS.
#'
#'   Titsias, M. & Lawrence, N. (2010). Bayesian Gaussian process latent
#'   variable model. AISTATS.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' Y <- matrix(rnorm(200 * 10), 200, 10)
#'
#' # 2D embedding
#' fit <- bgplvm(Y, Q = 2, M = 20, n_iter = 100, verbose = 0)
#' plot(fit)
#'
#' # Pseudo-time (1D)
#' fit1d <- bgplvm(Y, Q = 1, M = 20, n_iter = 100, verbose = 0)
#' pt <- extract_pseudotime(fit1d)
#' }
#' @export
bgplvm <- function(Y,
                   Q         = 2L,
                   M         = NULL,
                   n_iter    = 1000L,
                   lr        = 0.01,
                   optimizer = c("adam", "lbfgsb"),
                   init      = "pca",
                   verbose   = 100L,
                   tol_elbo  = 1e-4,
                   patience  = 20L,
                   ...) {
  Y <- as.matrix(Y)
  N <- nrow(Y)
  D <- ncol(Y)

  if (N < 2L) stop("Y must have at least 2 rows")
  if (D < Q)  stop("Q (", Q, ") cannot exceed D (", D, ")")

  if (is.null(M)) M <- min(100L, ceiling(sqrt(N)))
  M <- as.integer(M)
  if (M > N)  stop("M (", M, ") cannot exceed N (", N, ")")
  if (M < 2L) stop("M must be at least 2")

  optimizer <- match.arg(optimizer)
  cl        <- match.call()

  if (verbose > 0L) {
    message(sprintf(
      "Fitting sparse BGPLVM: N=%d, D=%d, Q=%d, M=%d, optimizer=%s",
      N, D, Q, M, optimizer
    ))
  }

  params <- .init_params(Y, Q, M, init)

  opt_result <- if (optimizer == "adam") {
    .adam_optimize(
      params_init = params,
      Y           = Y,
      n_iter      = n_iter,
      lr          = lr,
      tol_elbo    = tol_elbo,
      patience    = patience,
      verbose     = verbose,
      ...
    )
  } else {
    .lbfgsb_optimize(
      params_init = params,
      Y           = Y,
      n_iter      = n_iter,
      verbose     = verbose > 0L,
      ...
    )
  }

  p <- opt_result$params

  # Cache K_MM and its inverse at converged parameters
  ls_final  <- exp(p$log_ls)
  var_final <- exp(p$log_var)
  kmm_res   <- .compute_K_MM(p$Z, ls_final, var_final)

  structure(
    list(
      mu_X         = p$mu_X,
      log_s_X      = p$log_s_X,
      S_X          = exp(2 * p$log_s_X),
      Z            = p$Z,
      lengthscales = ls_final,
      variance     = var_final,
      sigma2       = exp(p$log_sigma2),
      K_MM         = kmm_res$K,
      K_MM_inv     = chol2inv(t(kmm_res$L)),
      elbo_trace   = opt_result$elbo_trace,
      N            = N,
      D            = D,
      Q            = Q,
      M            = M,
      Y            = Y,
      call         = cl,
      converged    = opt_result$converged,
      n_iter_run   = length(opt_result$elbo_trace)
    ),
    class = "bgplvm"
  )
}

#' Extract pseudo-time from a bgplvm fit
#'
#' Returns a per-cell pseudo-time coordinate with posterior uncertainty.
#' If the model was fit with Q=1, the single latent dimension IS pseudo-time.
#' If Q>1, the dimension with the smallest ARD lengthscale (most "active") is used.
#'
#' @param object   A \code{bgplvm} object
#' @param dim      Which latent dimension to use (default NULL = auto-select)
#' @param scale01  Rescale pseudo-time to [0, 1] (default TRUE)
#' @param flip     Reverse the time direction (default FALSE)
#'
#' @return A data.frame with columns:
#'   \item{cell}{Integer cell index (1:N)}
#'   \item{pseudotime}{Posterior mean pseudo-time coordinate}
#'   \item{lower95}{2.5th percentile: pseudotime - 1.96 * s}
#'   \item{upper95}{97.5th percentile: pseudotime + 1.96 * s}
#'   \item{s}{Posterior standard deviation}
#'
#' @export
extract_pseudotime <- function(object, dim = NULL, scale01 = TRUE, flip = FALSE) {
  if (!inherits(object, "bgplvm")) {
    stop("object must be of class 'bgplvm'")
  }

  # Auto-select: dimension with smallest lengthscale = most "on" in ARD
  if (is.null(dim)) {
    dim <- which.min(object$lengthscales)
  }
  if (dim < 1L || dim > object$Q) {
    stop("dim must be between 1 and Q (", object$Q, ")")
  }

  pt_mean <- object$mu_X[, dim]
  pt_sd   <- sqrt(object$S_X[, dim])

  if (flip) pt_mean <- -pt_mean

  if (scale01) {
    rng <- range(pt_mean)
    if (diff(rng) < .Machine$double.eps) {
      warning("Pseudo-time range is essentially zero; scale01 skipped.")
    } else {
      scale_factor <- diff(rng)
      pt_mean <- (pt_mean - rng[1]) / scale_factor
      pt_sd   <- pt_sd / scale_factor
    }
  }

  data.frame(
    cell       = seq_len(object$N),
    pseudotime = pt_mean,
    lower95    = pt_mean - 1.96 * pt_sd,
    upper95    = pt_mean + 1.96 * pt_sd,
    s          = pt_sd
  )
}

#' Predict latent embedding for new observations
#'
#' Fits q(X_new) for new data Ynew while fixing all kernel hyperparameters
#' at their trained values.
#'
#' @param object  A \code{bgplvm} object
#' @param newdata N_new x D matrix of new observations (same D as training Y)
#' @param n_iter  Optimization iterations for q(X_new) (default 200)
#' @param verbose Print interval (default 0)
#' @param ...     Additional arguments to .adam_optimize
#'
#' @return A data.frame with columns mu_1, ..., mu_Q, s_1, ..., s_Q
#' @export
predict.bgplvm <- function(object, newdata, n_iter = 200L, verbose = 0L, ...) {
  Ynew <- as.matrix(newdata)
  N_new <- nrow(Ynew)

  if (ncol(Ynew) != object$D) {
    stop("newdata has ", ncol(Ynew), " columns but model expects ", object$D)
  }

  Q <- object$Q
  M <- object$M

  # Fix hyperparameters; only optimize mu_X and log_s_X for new cells
  fixed_log_ls     <- log(object$lengthscales)
  fixed_log_var    <- log(object$variance)
  fixed_log_sigma2 <- log(object$sigma2)
  fixed_Z          <- object$Z

  # Initialize new latent coords via regression onto training latent space
  # Simple: project Ynew onto training principal components
  pca_new <- tryCatch({
    pc <- prcomp(object$Y, center = TRUE, scale. = FALSE)
    scale(Ynew, center = pc$center, scale = FALSE) %*% pc$rotation[, seq_len(Q), drop = FALSE]
  }, error = function(e) {
    matrix(0, N_new, Q)
  })

  init_new <- list(
    mu_X       = pca_new,
    log_s_X    = matrix(log(0.1), N_new, Q),
    Z          = fixed_Z,
    log_ls     = fixed_log_ls,
    log_var    = fixed_log_var,
    log_sigma2 = fixed_log_sigma2
  )

  # Optimize only mu_X and log_s_X; freeze other params by zeroing their gradients
  frozen_grad <- function(params, Y, YtY_trace = NULL) {
    g <- .grad_sparse_elbo(params, Y, YtY_trace)
    # Zero out gradients for frozen hyperparameters
    g$Z          <- g$Z * 0
    g$log_ls     <- g$log_ls * 0
    g$log_var    <- g$log_var * 0
    g$log_sigma2 <- g$log_sigma2 * 0
    g
  }

  res <- .adam_optimize(
    params_init = init_new,
    Y           = Ynew,
    n_iter      = n_iter,
    verbose     = verbose,
    grad_fn     = frozen_grad,
    ...
  )

  p     <- res$params
  s_mat <- sqrt(exp(2 * p$log_s_X))

  df <- as.data.frame(p$mu_X)
  names(df) <- paste0("mu_", seq_len(Q))
  for (q in seq_len(Q)) df[[paste0("s_", q)]] <- s_mat[, q]

  df
}

# ---- S3 Methods ---------------------------------------------------------------

#' @export
print.bgplvm <- function(x, ...) {
  cat("Bayesian GPLVM (sparse, variational)\n")
  cat(sprintf("  N = %d cells, D = %d dims, Q = %d latent, M = %d inducing\n",
              x$N, x$D, x$Q, x$M))
  cat(sprintf("  Noise sigma^2       : %.4f\n", x$sigma2))
  cat(sprintf("  Signal variance     : %.4f\n", x$variance))
  cat(sprintf("  ARD lengthscales    : %s\n",
              paste(sprintf("%.3f", x$lengthscales), collapse = ", ")))
  cat(sprintf("  Iterations run      : %d  (converged: %s)\n",
              x$n_iter_run, x$converged))
  if (!is.null(x$elbo_trace) && length(x$elbo_trace) > 0L) {
    cat(sprintf("  Final ELBO          : %.4f\n", tail(x$elbo_trace, 1)))
  }
  invisible(x)
}

#' @export
summary.bgplvm <- function(object, ...) {
  print(object)
  if (object$Q >= 1L) {
    cat("\nPseudo-time summary (dim with smallest lengthscale):\n")
    pt <- extract_pseudotime(object, scale01 = TRUE)
    cat(sprintf("  Mean uncertainty (s): %.4f\n", mean(pt$s)))
    cat(sprintf("  Pseudo-time range   : [%.4f, %.4f]\n",
                min(pt$pseudotime), max(pt$pseudotime)))
  }
  invisible(object)
}

#' Plot the latent space of a bgplvm fit
#'
#' @param x           A \code{bgplvm} object
#' @param dims        Integer vector of length 2: which latent dims to plot (default c(1,2))
#' @param col         Optional color vector of length N
#' @param uncertainty If TRUE, draw uncertainty ellipses scaled by posterior std dev
#' @param ...         Additional arguments passed to \code{plot()}
#' @export
plot.bgplvm <- function(x, dims = c(1L, 2L), col = NULL,
                        uncertainty = FALSE, ...) {
  if (x$Q < 2L) {
    # 1D: plot pseudo-time with uncertainty ribbon
    pt <- extract_pseudotime(x, scale01 = TRUE)
    ord <- order(pt$pseudotime)
    plot(pt$pseudotime[ord], pt$s[ord],
         type = "l",
         xlab = "Pseudo-time",
         ylab = "Posterior std dev",
         main = "BGPLVM Pseudo-time Uncertainty",
         ...)
    return(invisible(x))
  }

  mu <- x$mu_X[, dims, drop = FALSE]

  plot_args <- list(
    x    = mu[, 1],
    y    = mu[, 2],
    xlab = paste0("Latent dim ", dims[1]),
    ylab = paste0("Latent dim ", dims[2]),
    main = "BGPLVM Latent Space",
    col  = if (is.null(col)) "black" else col,
    pch  = 16L,
    cex  = 0.7
  )
  do.call(plot, c(plot_args, list(...)))

  if (uncertainty) {
    s <- sqrt(x$S_X[, dims, drop = FALSE])
    for (i in seq_len(x$N)) {
      graphics::symbols(
        x      = mu[i, 1],
        y      = mu[i, 2],
        circles = sqrt(s[i, 1]^2 + s[i, 2]^2),
        add    = TRUE,
        inches = FALSE,
        fg     = grDevices::adjustcolor(
          if (is.null(col)) "black" else col[i], alpha.f = 0.2
        )
      )
    }
  }

  invisible(x)
}

#' Check convergence diagnostics for a bgplvm fit
#'
#' Prints a summary of convergence and suggests remedies if not converged.
#'
#' @param object A \code{bgplvm} object
#' @param ... Ignored
#' @export
check_convergence <- function(object, ...) UseMethod("check_convergence")

#' @export
check_convergence.bgplvm <- function(object, ...) {
  cat("=== Convergence Diagnostics ===\n")
  cat(sprintf("  Iterations run : %d\n", object$n_iter_run))
  cat(sprintf("  Converged      : %s\n", object$converged))

  if (!is.null(object$elbo_trace) && length(object$elbo_trace) >= 10L) {
    n  <- length(object$elbo_trace)
    tail_elbo <- object$elbo_trace[max(1L, n - 9L):n]
    cat(sprintf("  Last 10 ELBOs  : %s\n",
                paste(sprintf("%.2f", tail_elbo), collapse = ", ")))

    # Check if ELBO is still increasing
    delta <- diff(tail(object$elbo_trace, 20L))
    if (mean(delta) > 1e-2) {
      cat("  WARNING: ELBO still increasing — consider more iterations.\n")
    } else {
      cat("  ELBO appears stable.\n")
    }
  }

  if (!object$converged) {
    cat("\nSuggestions:\n")
    cat("  - Increase n_iter\n")
    cat("  - Try a lower learning rate (lr = 0.005)\n")
    cat("  - Increase M (more inducing points)\n")
  }

  invisible(object)
}
