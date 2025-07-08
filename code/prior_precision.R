# generate prior precision matrices for VAR1
make_VAR1_precision <- function(K, d, A=NULL, Q=NULL) {
  if (is.null(A)) A <- diag(d) * 0.8
  if (is.null(Q)) Q <- diag(d) * 0.5

  Qinv <- solve(Q)
  At_Qinv <- t(A) %*% Qinv

  H <- matrix(0, d*K, d*K)

  # Prior for U1 ~ N(0, I)
  H[1:d, 1:d] <- diag(d)

  # VAR(1) structure
  for (k in 2:K) {
    idx_k   <- ((k-1)*d + 1):(k*d)
    idx_km1 <- ((k-2)*d + 1):( (k-1)*d )

    # Diagonal blocks
    H[idx_k, idx_k]     <- H[idx_k, idx_k] + Qinv
    H[idx_km1, idx_km1] <- H[idx_km1, idx_km1] + At_Qinv %*% A

    # Off-diagonal blocks
    H[idx_k, idx_km1]   <- H[idx_k, idx_km1] - Qinv %*% A
    H[idx_km1, idx_k]   <- H[idx_km1, idx_k] - At_Qinv
  }
  return(H)
}


# generate prior precision matrices for p-th order random walk
make_random_walk_precision <- function(K, d, q=1, lambda=1) {
  # Construct q-th order difference matrix
  D_q <- diag(K)
  for (i in 1:q) {
    D_q <- diff(D_q)
  }

  # Precision scalar part
  Q_scalar <- t(D_q) %*% D_q

  # Add tiny ridge for numerical stability
  Q_scalar <- Q_scalar + 1e-6 * diag(K)

  # Final precision matrix
  Q_U <- lambda * kronecker(Q_scalar, diag(d))

  return(Q_U)
}
