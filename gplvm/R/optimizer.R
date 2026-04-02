#' Adam optimizer for sparse BGPLVM
#'
#' Maximizes the ELBO using the Adam adaptive gradient algorithm.
#' Parameters are stored as named lists (no pack/unpack per step),
#' so Adam moment states mirror the parameter structure.
#'
#' @param params_init  Named list of initial parameters (from .init_params)
#' @param Y            N x D data matrix
#' @param n_iter       Maximum iterations (default 1000)
#' @param lr           Learning rate (default 0.01)
#' @param beta1        First moment decay (default 0.9)
#' @param beta2        Second moment decay (default 0.999)
#' @param eps          Numerical epsilon (default 1e-8)
#' @param tol_elbo     Early stopping: converge if |delta ELBO| < tol_elbo
#'                     for `patience` consecutive iterations (default 1e-4)
#' @param patience     Patience for early stopping (default 20)
#' @param verbose      Print every this many iterations; 0 = silent (default 100)
#' @param grad_fn      Gradient function (default .grad_sparse_elbo)
#' @return list(params = final params, elbo_trace = numeric vector, converged = logical)
#' @keywords internal
.adam_optimize <- function(params_init, Y,
                           n_iter   = 1000L,
                           lr       = 0.01,
                           beta1    = 0.9,
                           beta2    = 0.999,
                           eps      = 1e-8,
                           tol_elbo = 1e-4,
                           patience = 20L,
                           verbose  = 100L,
                           grad_fn  = .grad_sparse_elbo) {
  params      <- params_init
  elbo_trace  <- numeric(n_iter)
  YtY_trace   <- sum(Y^2)

  # Initialize Adam moment states (same structure as params)
  zero_like <- function(p) lapply(p, function(x) x * 0)
  m <- zero_like(params)  # first moments
  v <- zero_like(params)  # second moments

  no_improve_count <- 0L
  best_elbo        <- -Inf
  converged        <- FALSE

  for (iter in seq_len(n_iter)) {
    # Gradient (ascent direction)
    g <- grad_fn(params, Y, YtY_trace)

    # Update moments for each parameter field
    for (nm in names(params)) {
      m[[nm]] <- beta1 * m[[nm]] + (1 - beta1) * g[[nm]]
      v[[nm]] <- beta2 * v[[nm]] + (1 - beta2) * g[[nm]]^2
    }

    # Bias-corrected estimates
    m_hat <- lapply(m, function(x) x / (1 - beta1^iter))
    v_hat <- lapply(v, function(x) x / (1 - beta2^iter))

    # Gradient ascent step
    for (nm in names(params)) {
      params[[nm]] <- params[[nm]] + lr * m_hat[[nm]] / (sqrt(v_hat[[nm]]) + eps)
    }

    # Clamp to numerically safe ranges
    params <- .clamp_params(params)

    # Evaluate ELBO
    elbo_val       <- .compute_sparse_elbo(params, Y, YtY_trace)
    elbo_trace[iter] <- elbo_val

    if (verbose > 0L && iter %% verbose == 0L) {
      message(sprintf("Iter %4d / %d  |  ELBO: %10.4f", iter, n_iter, elbo_val))
    }

    # Early stopping
    if (elbo_val > best_elbo + tol_elbo) {
      best_elbo        <- elbo_val
      no_improve_count <- 0L
    } else {
      no_improve_count <- no_improve_count + 1L
      if (no_improve_count >= patience) {
        if (verbose > 0L) {
          message(sprintf("Converged at iter %d (no improvement for %d iters)",
                          iter, patience))
        }
        converged   <- TRUE
        elbo_trace  <- elbo_trace[seq_len(iter)]
        break
      }
    }
  }

  if (!converged) elbo_trace <- elbo_trace[seq_len(n_iter)]

  list(params = params, elbo_trace = elbo_trace, converged = converged)
}

#' L-BFGS-B optimizer wrapper for sparse BGPLVM
#'
#' Alternative to Adam; better for small N where the full parameter vector
#' fits comfortably in memory and second-order curvature helps.
#'
#' @param params_init  Named list of initial parameters
#' @param Y            N x D data matrix
#' @param n_iter       Maximum iterations passed to stats::optim (default 500)
#' @param verbose      Logical; print optim output (default FALSE)
#' @return list(params, elbo_trace = NULL, converged = logical)
#' @keywords internal
.lbfgsb_optimize <- function(params_init, Y, n_iter = 500L, verbose = FALSE) {
  N <- nrow(Y)
  Q <- ncol(params_init$mu_X)
  M <- nrow(params_init$Z)
  YtY_trace <- sum(Y^2)

  theta0 <- .pack_params(params_init, N, Q, M)

  neg_elbo <- function(th) {
    p <- .unpack_params(th, N, Q, M)
    -.compute_sparse_elbo(p, Y, YtY_trace)
  }

  neg_grad <- function(th) {
    p  <- .unpack_params(th, N, Q, M)
    gp <- .grad_sparse_elbo(p, Y, YtY_trace)
    -.pack_params(gp, N, Q, M)
  }

  ctrl <- list(maxit = n_iter, trace = as.integer(verbose))
  res  <- stats::optim(
    par     = theta0,
    fn      = neg_elbo,
    gr      = neg_grad,
    method  = "L-BFGS-B",
    control = ctrl
  )

  final_params <- .unpack_params(res$par, N, Q, M)
  converged    <- res$convergence == 0L

  list(params = final_params, elbo_trace = NULL, converged = converged)
}
