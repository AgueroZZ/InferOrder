#' Compute the inducing kernel matrix K_MM with Cholesky
#'
#' @param Z           M x Q matrix of inducing inputs
#' @param lengthscales Length-Q vector of ARD lengthscales
#' @param variance    Signal variance (scalar)
#' @return list(K = M x M jittered matrix, L = lower Cholesky, jitter)
#' @keywords internal
.compute_K_MM <- function(Z, lengthscales, variance) {
  K <- kernel_se_ard(Z, Z, lengthscales = lengthscales, variance = variance)
  chol_res <- .chol_with_jitter(K)
  list(K = .add_jitter(K, chol_res$jitter), L = chol_res$L, jitter = chol_res$jitter)
}

#' Analytical E_q(X)[K_NM] for ARD squared-exponential kernel
#'
#' For mean-field q(X) = prod_i N(mu_i, diag(s2_i)), computes the N x M matrix
#' whose (i, m) entry is E_{x_i ~ q}[k(x_i, z_m)].
#'
#' For the SE-ARD kernel k(x, z) = var * exp(-0.5 * sum_q (x_q - z_q)^2 / ls_q^2),
#' with x_i ~ N(mu_i, diag(s2_i)), each 1D integral is a Gaussian-times-Gaussian:
#'
#'   E[exp(-0.5*(x - z)^2 / ls^2)] = exp(-0.5*(mu - z)^2 / (ls^2 + s2)) * sqrt(ls^2 / (ls^2 + s2))
#'
#' The full result factorizes over dimensions.
#'
#' @param mu_X  N x Q posterior means
#' @param s2_X  N x Q posterior variances (positive)
#' @param Z     M x Q inducing inputs
#' @param lengthscales Length-Q vector
#' @param variance     Scalar signal variance
#' @return N x M matrix (psi1 statistic)
#' @keywords internal
.expected_K_NM <- function(mu_X, s2_X, Z, lengthscales, variance) {
  N <- nrow(mu_X)
  M <- nrow(Z)
  Q <- ncol(mu_X)

  log_psi1 <- matrix(0, N, M)

  for (q in seq_len(Q)) {
    ls2_q <- lengthscales[q]^2
    # Outer differences: (N x 1) - (1 x M) => N x M
    diff2  <- outer(mu_X[, q], Z[, q], "-")^2
    denom  <- outer(s2_X[, q], rep(ls2_q, M), "+")  # N x M, ls^2 + s2_i
    log_psi1 <- log_psi1 -
      0.5 * diff2 / denom +
      0.5 * log(ls2_q) -
      0.5 * log(denom)
  }

  variance * exp(log_psi1)
}

#' Analytical E_q(X)[K_MN K_NM] for ARD squared-exponential kernel
#'
#' Computes the M x M "psi2" matrix:
#'   psi2_{m,m'} = sum_i E_{x_i}[k(z_m, x_i) * k(x_i, z_m')]
#'
#' For SE-ARD, the closed-form per cell i, per pair (m, m') involves:
#'   - A term from the product of two SE kernels integrated against N(x_i; mu_i, S_i)
#'   - Factorizes over Q dimensions
#'
#' Per dimension q, the integrand is:
#'   exp(-0.5*(z_mq - x)^2/ls_q^2) * exp(-0.5*(x - z_m'q)^2/ls_q^2) * N(x; mu_iq, s2_iq)
#' which equals:
#'   exp(-0.25*(z_mq - z_m'q)^2/ls_q^2) *
#'   N(0.5*(z_mq + z_m'q); mu_iq, ls_q^2/2 + s2_iq) * sqrt(pi * (ls_q^2 + 2*s2_iq) / ls_q^2... )
#'
#' Full derivation: Titsias & Lawrence (2010), Appendix A, equations A.4--A.6.
#'
#' @param mu_X  N x Q
#' @param s2_X  N x Q
#' @param Z     M x Q
#' @param lengthscales Length-Q
#' @param variance     Scalar
#' @return M x M symmetric matrix (psi2 statistic)
#' @keywords internal
.psi2_matrix <- function(mu_X, s2_X, Z, lengthscales, variance) {
  N <- nrow(mu_X)
  M <- nrow(Z)
  Q <- ncol(mu_X)

  psi2 <- matrix(0, M, M)

  for (i in seq_len(N)) {
    # Contribution of cell i: M x M matrix
    log_contrib <- matrix(0, M, M)
    for (q in seq_len(Q)) {
      ls2_q  <- lengthscales[q]^2
      s2_iq  <- s2_X[i, q]
      mu_iq  <- mu_X[i, q]

      # (z_mq - z_m'q)^2 cross term
      zdiff2 <- outer(Z[, q], Z[, q], "-")^2  # M x M
      # midpoint of z_mq and z_m'q
      zmid   <- outer(Z[, q], Z[, q], "+") / 2  # M x M
      # (midpoint - mu_iq)^2
      mdiff2 <- (zmid - mu_iq)^2

      denom_q <- ls2_q / 2 + s2_iq

      log_contrib <- log_contrib -
        0.25 * zdiff2 / ls2_q -
        0.5  * mdiff2 / denom_q +
        0.5  * log(ls2_q / 2) -
        0.5  * log(denom_q)
    }
    psi2 <- psi2 + variance^2 * exp(log_contrib)
  }

  # Symmetrize for numerical cleanliness
  0.5 * (psi2 + t(psi2))
}

#' Expected trace of K_NN under q(X) for SE-ARD kernel
#'
#' For SE-ARD, k(x_i, x_i) = variance for all x_i (diagonal of kernel = constant).
#' Therefore E[tr(K_NN)] = N * variance.
#'
#' @param N        Integer number of observations
#' @param variance Scalar signal variance
#' @return Scalar N * variance
#' @keywords internal
.expected_trace_K_NN <- function(N, variance) {
  N * variance
}
