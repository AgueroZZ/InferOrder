compute_permuted_energy <- function(U, Q, P = NULL) {
  # U: K × D numeric matrix of row‐means
  # Q: (K*D) × (K*D) precision matrix (constructed via kronecker)
  # P: optional integer vector of length K; a permutation of 1:K
  #
  # Returns: the scalar (P U)ᵀ Q (P U), where U is first permuted by P (if given),
  # then column-stacked into a length-K*D vector.

  # --- Input checks ---
  if (!is.matrix(U)) {
    stop("`U` must be a numeric matrix")
  }
  K <- nrow(U)
  D <- ncol(U)
  if (!is.null(P)) {
    if (length(P) != K || !setequal(P, seq_len(K))) {
      stop("`P` must be a permutation vector of 1:K")
    }
  }
  if (!is.matrix(Q) || any(dim(Q) != c(K * D, K * D))) {
    stop("`Q` must be a (K*D) x (K*D) matrix, where K = nrow(U), D = ncol(U)")
  }

  # --- Apply permutation to rows of U, if provided ---
  Up <- if (is.null(P)) U else U[P, , drop = FALSE]

  # --- Column-major stack into a vector of length K*D ---
  # as.numeric(Up) produces c(Up[1,1], Up[2,1], …, Up[K,1],
  #                          Up[1,2], Up[2,2], …, Up[K,2], …)
  u_vec <- as.numeric(Up)

  # --- Compute quadratic form ---
  # crossprod(u_vec, Q %*% u_vec) gives a 1×1 matrix
  energy <- as.numeric(crossprod(u_vec, Q %*% u_vec))

  return(energy)
}

optimal_permutation <- function(mu_list, Q, verbose = FALSE) {
  # mu_list: list of length K, each element is a numeric vector of length D
  # Q: (K*D) x (K*D) precision matrix
  # verbose: if TRUE, report only whether energy improved and from what to what
  #
  # Returns the permutation (vector of length K) that minimizes (P U)ᵀ Q (P U).

  # Stack mu_list into U (K x D)
  U <- do.call(rbind, mu_list)
  K <- nrow(U)
  D <- ncol(U)

  # Compute energy with no permutation
  u_base <- as.numeric(U)
  energy_base <- as.numeric(crossprod(u_base, Q %*% u_base))

  # Initialize best
  best_energy <- energy_base
  best_perm   <- seq_len(K)

  # Try each column-based ordering
  for (j in seq_len(D)) {
    Pj <- order(U[, j])
    u_vec_j <- as.numeric(U[Pj, , drop = FALSE])
    Ej <- as.numeric(crossprod(u_vec_j, Q %*% u_vec_j))
    if (Ej < best_energy) {
      best_energy <- Ej
      best_perm   <- Pj
    }
  }

  # Verbose report
  if (verbose) {
    if (best_energy < energy_base) {
      cat(sprintf(
        "Energy improved: %.6g -> %.6g\n",
        energy_base, best_energy
      ))
    } else {
      cat("No improvement in energy (remains ", sprintf("%.6g", energy_base), ")\n", sep = "")
    }
  }

  return(best_perm)
}

plot_sequence <- function(U, P = NULL, tol = 1e-6) {
  # U: K × 2 numeric matrix of 2D points
  # P: optional integer vector of length K representing a permutation of 1:K
  #    If NULL, no permutation is applied
  #
  # This function draws two side-by-side plots:
  #   Left:  original sequence of points connected by arrows
  #   Right: permuted sequence (if P is provided) connected by arrows

  if (!is.matrix(U) || ncol(U) != 2) {
    stop("`U` must be a K x 2 matrix")
  }
  K <- nrow(U)

  if (!is.null(P)) {
    if (length(P) != K || !setequal(P, seq_len(K))) {
      stop("`P` must be a permutation vector of 1:K")
    }
    Up <- U[P, , drop = FALSE]
    right_title <- "Permuted sequence"
  } else {
    Up <- U
    right_title <- "Original sequence (no P)"
  }

  # Set up side-by-side plotting
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par), add = TRUE)

  # Left panel: original sequence
  plot(U, pch = 19, asp = 1,
       xlab = "U[,1]", ylab = "U[,2]",
       main = "Original sequence")
  for (k in seq_len(K - 1)) {
    if (sqrt(sum((U[k + 1, ] - U[k, ])^2)) > tol) {
      arrows(U[k, 1], U[k, 2],
             U[k + 1, 1], U[k + 1, 2],
             col = "orange", lwd = 2, length = 0.1)
    }
  }

  # Right panel: permuted sequence
  plot(Up, pch = 19, asp = 1,
       xlab = "Up[,1]", ylab = "Up[,2]",
       main = right_title)
  for (k in seq_len(K - 1)) {
    if (sqrt(sum((Up[k + 1, ] - Up[k, ])^2)) > tol) {
      arrows(Up[k, 1], Up[k, 2],
             Up[k + 1, 1], Up[k + 1, 2],
             col = "steelblue", lwd = 2, length = 0.1)
    }
  }
}


# set.seed(42)
# U <- matrix(rnorm(40), ncol = 2)
# Q <- make_random_walk_precision(K = 20, d = 2, q = 1, lambda = 1)
# # order based on the first column
# P <- order(U[,2])
# plot_sequence(U, P)
#
#
# energy_original <- compute_permuted_energy(U, Q)
# energy_permuted <- compute_permuted_energy(U, Q, P)
# cat("Original energy:", energy_original, "\n")
# cat("Permuted energy:", energy_permuted, "\n")
