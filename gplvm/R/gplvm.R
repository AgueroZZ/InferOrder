#' Fit a Bayesian GPLVM
#'
#' Fits a Bayesian Gaussian Process Latent Variable Model using variational
#' inference. The model assumes:
#'   Y ~ GP(0, K(X, X) + sigma^2 * I)
#' where X are the latent coordinates with prior X ~ N(0, I).
#'
#' @param Y Data matrix of size N x D (N observations, D dimensions)
#' @param Q Dimensionality of the latent space
#' @param kernel Kernel function to use (default: kernel_se_ard)
#' @param n_iter Number of optimization iterations
#' @param lr Learning rate for gradient-based optimization
#' @param verbose Print progress every this many iterations (0 = silent)
#'
#' @return A list of class "bgplvm" with components:
#'   \item{mu_X}{Posterior mean of latent coordinates (N x Q)}
#'   \item{S_X}{Posterior variance of latent coordinates (N x Q)}
#'   \item{kernel_params}{Fitted kernel hyperparameters}
#'   \item{sigma2}{Fitted noise variance}
#'   \item{elbo_trace}{ELBO values per iteration}
#' @export
bgplvm <- function(Y, Q = 2, kernel = kernel_se_ard,
                   n_iter = 200, lr = 0.01, verbose = 50) {
  Y <- as.matrix(Y)
  N <- nrow(Y)
  D <- ncol(Y)

  # Initialize latent means via PCA
  pca <- prcomp(Y, center = TRUE, scale. = FALSE)
  mu_X <- pca$x[, seq_len(Q), drop = FALSE]
  S_X  <- matrix(0.1, nrow = N, ncol = Q)  # log-variance init

  # Hyperparameters (log-parameterized for unconstrained optimization)
  log_lengthscales <- rep(0, Q)   # lengthscale = 1 per dim
  log_variance     <- 0           # signal variance = 1
  log_sigma2       <- log(0.1)    # noise variance = 0.1

  elbo_trace <- numeric(n_iter)

  for (i in seq_len(n_iter)) {
    params <- list(
      mu_X = mu_X, S_X = S_X,
      log_ls = log_lengthscales,
      log_var = log_variance,
      log_sigma2 = log_sigma2
    )

    elbo_val <- .compute_elbo(Y, params, kernel)
    elbo_trace[i] <- elbo_val

    if (verbose > 0 && (i %% verbose == 0)) {
      message(sprintf("Iter %d / %d | ELBO: %.4f", i, n_iter, elbo_val))
    }

    # Gradient step (finite differences placeholder — replace with analytic grads)
    grad <- .finite_diff_elbo(Y, params, kernel, eps = 1e-4)
    mu_X             <- mu_X             + lr * grad$mu_X
    S_X              <- S_X              + lr * grad$S_X
    log_lengthscales <- log_lengthscales + lr * grad$log_ls
    log_variance     <- log_variance     + lr * grad$log_var
    log_sigma2       <- log_sigma2       + lr * grad$log_sigma2
  }

  structure(
    list(
      mu_X = mu_X,
      S_X  = exp(S_X),
      kernel_params = list(
        lengthscales = exp(log_lengthscales),
        variance     = exp(log_variance)
      ),
      sigma2      = exp(log_sigma2),
      elbo_trace  = elbo_trace,
      Y = Y, Q = Q
    ),
    class = "bgplvm"
  )
}

#' Compute the Variational ELBO
#'
#' @param Y Data matrix N x D
#' @param params List of variational and kernel parameters
#' @param kernel Kernel function
#' @keywords internal
.compute_elbo <- function(Y, params, kernel) {
  N <- nrow(Y)
  D <- ncol(Y)

  ls     <- exp(params$log_ls)
  var_k  <- exp(params$log_var)
  sigma2 <- exp(params$log_sigma2)
  mu_X   <- params$mu_X
  S_X    <- exp(params$S_X)   # variances (positive)

  # Expected kernel using only posterior means (0th-order approximation)
  K <- kernel(mu_X, lengthscales = ls, variance = var_k)
  K <- K + (sigma2 + 1e-6) * diag(N)

  # Log-likelihood term: -0.5 * D * [log|K| + tr(K^{-1} YY^T/D)]
  chol_K <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(chol_K)) return(-Inf)

  log_det_K <- 2 * sum(log(diag(chol_K)))
  K_inv_Y   <- chol2inv(chol_K) %*% Y
  trace_term <- sum(Y * K_inv_Y)

  ll <- -0.5 * (D * log_det_K + trace_term + N * D * log(2 * pi))

  # KL divergence KL[q(X) || p(X)], p(X) = N(0, I)
  # KL = 0.5 * sum(mu^2 + S - 1 - log(S))
  kl <- 0.5 * sum(mu_X^2 + S_X - 1 - log(S_X + 1e-8))

  ll - kl
}

#' Finite-difference gradient of ELBO (placeholder)
#' @keywords internal
.finite_diff_elbo <- function(Y, params, kernel, eps = 1e-4) {
  f0 <- .compute_elbo(Y, params, kernel)

  grad_mu_X <- matrix(0, nrow = nrow(params$mu_X), ncol = ncol(params$mu_X))
  for (i in seq_len(nrow(params$mu_X))) {
    for (j in seq_len(ncol(params$mu_X))) {
      p2 <- params; p2$mu_X[i, j] <- p2$mu_X[i, j] + eps
      grad_mu_X[i, j] <- (.compute_elbo(Y, p2, kernel) - f0) / eps
    }
  }

  # Scalar gradients
  scalar_grad <- function(name) {
    p2 <- params; p2[[name]] <- p2[[name]] + eps
    (.compute_elbo(Y, p2, kernel) - f0) / eps
  }

  list(
    mu_X       = grad_mu_X,
    S_X        = matrix(0, nrow = nrow(params$S_X), ncol = ncol(params$S_X)),
    log_ls     = sapply(seq_along(params$log_ls), function(k) {
      p2 <- params; p2$log_ls[k] <- p2$log_ls[k] + eps
      (.compute_elbo(Y, p2, kernel) - f0) / eps
    }),
    log_var    = scalar_grad("log_var"),
    log_sigma2 = scalar_grad("log_sigma2")
  )
}

#' Print method for bgplvm
#' @export
print.bgplvm <- function(x, ...) {
  cat("Bayesian GPLVM\n")
  cat(sprintf("  N = %d observations, D = %d dimensions, Q = %d latent dims\n",
              nrow(x$Y), ncol(x$Y), x$Q))
  cat(sprintf("  Noise sigma^2 = %.4f\n", x$sigma2))
  cat(sprintf("  Kernel lengthscales: %s\n",
              paste(round(x$kernel_params$lengthscales, 3), collapse = ", ")))
  cat(sprintf("  Final ELBO: %.4f\n", tail(x$elbo_trace, 1)))
  invisible(x)
}

#' Plot latent space of a bgplvm fit
#'
#' @param x A bgplvm object
#' @param dims Which two latent dimensions to plot (default c(1, 2))
#' @param ... Additional arguments passed to plot()
#' @export
plot.bgplvm <- function(x, dims = c(1, 2), ...) {
  mu <- x$mu_X[, dims]
  plot(mu[, 1], mu[, 2],
       xlab = paste0("Latent dim ", dims[1]),
       ylab = paste0("Latent dim ", dims[2]),
       main = "BGPLVM Latent Space",
       ...)
}
