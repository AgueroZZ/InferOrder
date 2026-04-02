#' Compute the sparse Bayesian GPLVM ELBO
#'
#' Implements the collapsed variational lower bound from Titsias (2009) and
#' Titsias & Lawrence (2010), with expectations over q(X) computed analytically
#' via psi-statistics for the ARD squared-exponential kernel.
#'
#' ELBO = ELL - KL[q(X) || p(X)]
#'
#' After marginalizing the inducing outputs U analytically, the expected
#' log-likelihood (summed over D output dimensions) is:
#'
#'   ELL = - N*D/2 * log(2*pi*sigma2)
#'         - D/2 * (log|A| - log|K_MM|)
#'         - 1/(2*sigma2) * ||Y||_F^2
#'         + 1/(2*sigma2^2) * ||L_A^{-T} psi1^T Y||_F^2
#'         - D/(2*sigma2) * (psi0 - tr(K_MM^{-1} psi2))
#'
#' where A = K_MM + sigma^{-2} * psi2, and L_A is its lower Cholesky factor.
#'
#' @param params       Named list: mu_X (N x Q), log_s_X (N x Q), Z (M x Q),
#'                     log_ls (length Q), log_var (scalar), log_sigma2 (scalar)
#' @param Y            N x D data matrix
#' @param YtY_trace    Precomputed sum(Y^2); if NULL, computed internally
#' @return Scalar ELBO value (to maximize)
#' @keywords internal
.compute_sparse_elbo <- function(params, Y, YtY_trace = NULL) {
  N <- nrow(Y)
  D <- ncol(Y)

  # --- Natural-scale parameters ---
  ls     <- exp(params$log_ls)
  var_k  <- exp(params$log_var)
  sigma2 <- exp(params$log_sigma2)
  s2_X   <- exp(2 * params$log_s_X)  # N x Q posterior variances
  mu_X   <- params$mu_X
  Z      <- params$Z
  M      <- nrow(Z)

  # --- Psi statistics ---
  psi0 <- .expected_trace_K_NN(N, var_k)                    # scalar = N * var_k
  psi1 <- .expected_K_NM(mu_X, s2_X, Z, ls, var_k)         # N x M
  psi2 <- .psi2_matrix(mu_X, s2_X, Z, ls, var_k)           # M x M

  # --- K_MM Cholesky: K_MM = L L^T, L lower triangular ---
  kmm_res <- .compute_K_MM(Z, ls, var_k)
  L       <- kmm_res$L

  # --- A = K_MM + sigma^{-2} psi2, Cholesky: A = L_A L_A^T ---
  A      <- kmm_res$K + psi2 / sigma2
  La_res <- .chol_with_jitter(A)
  L_A    <- La_res$L

  # --- Log-determinants (via Cholesky diagonals, never log(det(.))) ---
  log_det_K_MM <- 2 * sum(log(diag(L)))
  log_det_A    <- 2 * sum(log(diag(L_A)))

  # --- Data fit: 1/(2*sigma2^2) * ||L_A^{-T} psi1^T Y||_F^2 ---
  # backsolve(t(L_A), X) solves t(L_A) * result = X  =>  result = L_A^{-T} X
  psi1tY <- crossprod(psi1, Y)          # M x D: psi1^T Y
  c_mat  <- backsolve(t(L_A), psi1tY)  # M x D: L_A^{-T} psi1^T Y

  # --- tr(K_MM^{-1} psi2) ---
  # tr(K_MM^{-1} psi2) = tr(L^{-T} L^{-1} psi2)
  # Step 1: W = L^{-1} psi2   (forwardsolve: solves L * W = psi2)
  # Step 2: tr(L^{-T} W) = sum(diag(backsolve(t(L), W)))  (backsolve: solves t(L) * V = W)
  W_psi2            <- forwardsolve(L, psi2)
  trace_KMMInv_psi2 <- sum(diag(backsolve(t(L), W_psi2)))

  # --- ||Y||_F^2 ---
  if (is.null(YtY_trace)) YtY_trace <- sum(Y^2)

  # --- Expected log-likelihood (Titsias & Lawrence 2010) ---
  ell <- (
    - N * D / 2 * log(2 * pi * sigma2)         # -ND/2 * log(2pi*sigma^2)
    - D / 2 * (log_det_A - log_det_K_MM)        # -D/2 * log|A / K_MM|
    - YtY_trace / (2 * sigma2)                  # -1/(2sigma^2) * ||Y||^2
    + sum(c_mat^2) / (2 * sigma2^2)             # +1/(2sigma^4) * ||L_A^{-T} psi1^T Y||^2
    - D * (psi0 - trace_KMMInv_psi2) / (2 * sigma2)  # trace correction
  )

  # --- KL[q(X) || p(X)], p(X) = N(0, I) ---
  # KL = 0.5 * sum(mu^2 + s2 - 1 - log(s2))
  #    = 0.5 * sum(mu^2 + exp(2*log_s) - 1 - 2*log_s)
  kl <- 0.5 * sum(mu_X^2 + s2_X - 1 - 2 * params$log_s_X)

  ell - kl
}

#' Gradient of the sparse BGPLVM ELBO w.r.t. all parameters
#'
#' Uses numerical differentiation (numDeriv) for all parameters. This is
#' correct and numerically stable. Analytic gradients for mu_X and log_s_X
#' are a planned optimization for future versions.
#'
#' @param params       Named list (same as .compute_sparse_elbo)
#' @param Y            N x D data matrix
#' @param YtY_trace    Precomputed sum(Y^2); if NULL, computed internally
#' @return Named list with same structure as params, each field is the gradient
#' @keywords internal
.grad_sparse_elbo <- function(params, Y, YtY_trace = NULL) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop(
      "Package 'numDeriv' is required for gradient computation. ",
      "Install with: install.packages('numDeriv')"
    )
  }

  N <- nrow(Y)
  Q <- ncol(params$mu_X)
  M <- nrow(params$Z)

  if (is.null(YtY_trace)) YtY_trace <- sum(Y^2)

  theta <- .pack_params(params, N, Q, M)

  f <- function(th) {
    p <- .unpack_params(th, N, Q, M)
    .compute_sparse_elbo(p, Y, YtY_trace)
  }

  g <- numDeriv::grad(f, theta, method = "simple")
  .unpack_params(g, N, Q, M)
}
