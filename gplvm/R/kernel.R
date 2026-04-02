#' Squared Exponential (RBF) Kernel
#'
#' Computes the squared exponential covariance matrix between two sets of points.
#'
#' @param X1 Matrix of size n1 x D
#' @param X2 Matrix of size n2 x D (defaults to X1)
#' @param lengthscale Kernel lengthscale parameter
#' @param variance Signal variance (amplitude squared)
#'
#' @return Covariance matrix of size n1 x n2
#' @export
kernel_se <- function(X1, X2 = NULL, lengthscale = 1.0, variance = 1.0) {
  if (is.null(X2)) X2 <- X1
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  # Squared Euclidean distances
  D <- outer(rowSums(X1^2), rowSums(X2^2), "+") - 2 * X1 %*% t(X2)
  variance * exp(-0.5 * D / lengthscale^2)
}

#' ARD Squared Exponential Kernel
#'
#' Squared exponential kernel with Automatic Relevance Determination (ARD),
#' using a separate lengthscale per input dimension.
#'
#' @param X1 Matrix of size n1 x D
#' @param X2 Matrix of size n2 x D (defaults to X1)
#' @param lengthscales Vector of lengthscales, length D
#' @param variance Signal variance
#'
#' @return Covariance matrix of size n1 x n2
#' @export
kernel_se_ard <- function(X1, X2 = NULL, lengthscales, variance = 1.0) {
  if (is.null(X2)) X2 <- X1
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  # Scale each dimension by its lengthscale
  X1s <- sweep(X1, 2, lengthscales, "/")
  X2s <- sweep(X2, 2, lengthscales, "/")

  D <- outer(rowSums(X1s^2), rowSums(X2s^2), "+") - 2 * X1s %*% t(X2s)
  variance * exp(-0.5 * D)
}
