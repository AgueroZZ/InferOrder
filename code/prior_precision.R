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
make_random_walk_precision <- function(K, d, q=1, lambda=1, ridge=0) {
  # Construct q-th order difference matrix
  D_q <- diag(K)
  for (i in 1:q) {
    D_q <- diff(D_q)
  }

  # Precision scalar part
  Q_scalar <- t(D_q) %*% D_q

  # Add tiny ridge for numerical stability
  Q_scalar <- Q_scalar + ridge * diag(K)

  # Final precision matrix
  Q_U <- lambda * kronecker(Q_scalar, diag(d))

  return(Q_U)
}


# generate prior precision matrices for q-th order random walk on a KxK lattice
make_lattice_rwq_precision <- function(K, d, q = 1, lambda = 1, ridge = 0) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install.packages('Matrix')")
  }
  # Validate inputs
  if (K <= 0 || K != as.integer(K))    stop("K must be a positive integer.")
  if (d <= 0 || d != as.integer(d))    stop("d must be a positive integer.")
  if (q <= 0 || q != as.integer(q) || q >= K) {
    stop("q must be a positive integer less than K.")
  }

  # 1. Build the q-th order difference matrix in 1D
  D <- diag(K)
  for (i in seq_len(q)) {
    D <- diff(D)
  }

  # 2. 1D precision: Q1D = t(D) %*% D
  Q1D <- crossprod(D)
  Q1D <- Matrix::Matrix(Q1D, sparse = TRUE)

  # 3. Kronecker sum to get 2D precision on K×K lattice
  I_K <- Matrix::Diagonal(K)
  Q2D_base <- Matrix::kronecker(Q1D, I_K) +
    Matrix::kronecker(I_K, Q1D)

  # 4. Add ridge term if requested
  if (ridge > 0) {
    Q2D_base <- Q2D_base + ridge * Matrix::Diagonal(K * K)
  }

  # 5. Replicate across d independent processes
  Q_U <- lambda * Matrix::kronecker(Q2D_base, Matrix::Diagonal(d))

  return(Q_U)
}
