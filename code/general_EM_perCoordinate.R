library(matrixStats)
library(mvtnorm)
library(Matrix)

# Compute generalized log-determinant (sum of log of nonzero eigenvalues)
generalized_logdet <- function(Q, eigen_tol = 0, rank_deficiency = 0) {
  # Q: symmetric matrix
  # eigen_tol: eigenvalue threshold; if NULL, use a heuristic
  #
  # Returns sum(log(λ_i)) over all eigenvalues λ_i > eigen_tol

  # 1) compute eigenvalues
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)
  vals <- ev$values

  # if rank_deficiency > 0, take out the last `rank_deficiency` eigenvalues
  if (rank_deficiency > 0) {
    vals <- vals[-seq(length(vals) - rank_deficiency + 1, length(vals))]
  }

  # 2) set a default tolerance if none given
  if (is.null(eigen_tol)) {
    # anything smaller than n*max(vals)*eps is treated as zero
    eigen_tol <- length(vals) * max(vals) * .Machine$double.eps
  }

  # 3) keep only those above tol
  positive <- vals[vals > eigen_tol]
  if (length(positive) == 0) {
    warning("All eigenvalues are below tolerance; returning 0")
    return(0)
  }

  # 4) sum logs
  sum(log(positive))
}

# Numerically stable log-sum-exp
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Generate default initialization parameters for GMM
make_default_init <- function(X, K, ordering = TRUE) {
  d <- ncol(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  output <- list(
    pi    = rep(1/K, K),
    mu    = lapply(1:K, function(k) runif(d, mins, maxs)),
    sigma = lapply(1:K, function(k) diag(d))
  )

  if(ordering){
    # sort by the pivot mu (i.e. the mu with largest range)
    pivot_index <- which.max(apply(do.call(rbind, output$mu), 2, function(mu) max(mu) - min(mu)))

    mu_first <- sapply(output$mu, function(mu) mu[pivot_index])
    order_idx <- order(mu_first)
    output$pi    <- output$pi[order_idx]
    output$mu    <- output$mu[order_idx]
    output$sigma <- output$sigma[order_idx]
  }

  return(output)
}

# initialization using k-means clustering on ordering_vec
make_init <- function(X, ordering_vec, K, assume_EEI = TRUE, nugget = 0) {
  d <- ncol(X)

  # K-means on ordering vector
  cl <- kmeans(ordering_vec, centers = K)
  cluster_ids <- cl$cluster

  # Reorder clusters from left to right
  cluster_order <- order(cl$centers)
  cluster_rank  <- match(cluster_ids, cluster_order)

  # Means
  mu_list <- lapply(1:K, function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0) rep(0, d) else colMeans(X[idx, , drop = FALSE])
  })

  # EEI-style sigma
  if (assume_EEI) {
    # Compute global diagonal variance
    diffs <- do.call(rbind, lapply(1:K, function(k) {
      idx <- which(cluster_rank == k)
      sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
    }))
    diag_var <- apply(diffs^2, 2, mean) + nugget
    Sigma_shared <- diag(diag_var, d)
    sigma_list   <- replicate(K, Sigma_shared, simplify = FALSE)
  } else {
    # Full covariance per cluster
    sigma_list <- lapply(1:K, function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) <= 1) {
        diag(rep(1, d))
      } else {
        cov(X[idx, , drop = FALSE]) + nugget * diag(d)
      }
    })
  }

  Nk <- table(factor(cluster_rank, levels = 1:K))
  pi_vec <- as.numeric(Nk) / length(ordering_vec)

  list(pi = pi_vec, mu = mu_list, sigma = sigma_list)
}

# 0. Record convergence
compute_log_joint_observed <- function(X, params, Q_prior = NULL, eigen_tol = 0, rank_deficiency = 0) {
  K <- length(params$pi)
  log_dens_mat <- sapply(seq_len(K), function(k) {
    dmvnorm(X,
            mean  = params$mu[[k]],
            sigma = params$sigma[[k]],
            log   = TRUE)
  })
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K

  loglik <- sum(rowLogSumExps(log_gamma))

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu)
    penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    loglik  <- loglik - as.numeric(penalty) + 0.5 * logdetQ
  }

  return(loglik)
}

# 2. Penalized ELBO = Q-function + entropy + log-prior
compute_penalized_ELBO <- function(X, Gamma, params, Q_prior = NULL, eigen_tol = 0, rank_deficiency = 0 ) {
  K <- length(params$pi)
  log_dens_mat <- sapply(seq_len(K), function(k) {
    dmvnorm(X,
            mean  = params$mu[[k]],
            sigma = params$sigma[[k]],
            log   = TRUE)
  })
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K

  # 2b. Q-function: E_q[log p(X,Z)]
  Qval <- sum(Gamma * log_gamma)

  # 2c. entropy: −E_q[log q]
  entropy <- - sum(Gamma * log(Gamma + 1e-300))

  elbo <- Qval + entropy

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu)
    penalty <- 0.5 * crossprod(U_vec, Q_prior %*% U_vec)
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    elbo    <- elbo - as.numeric(penalty) + 0.5 * logdetQ
  }

  return(elbo)
}

# 1. Cache covariance inverses and log-determinants
init_cov_cache_fast <- function(params) {
  K <- length(params$pi)
  invSigma <- vector("list", K)
  logdet   <- numeric(K)
  for (k in seq_len(K)) {
    S <- params$sigma[[k]]
    L <- chol(S)
    invSigma[[k]] <- chol2inv(L)
    logdet[k]     <- 2 * sum(log(diag(L)))
  }
  params$invSigma <- invSigma
  params$logdet   <- logdet
  params
}

# 2. Vectorized E-step
ESTEP <- function(data, params) {
  n <- nrow(data); K <- length(params$pi); d <- ncol(data)
  log_gamma <- matrix(0, n, K)
  Mu_mat <- do.call(cbind, params$mu)
  for (k in seq_len(K)) {
    Xc <- sweep(data, 2, Mu_mat[,k], "-")
    Qf <- rowSums((Xc %*% params$invSigma[[k]]) * Xc)
    log_gamma[,k] <- log(params$pi[k]) - 0.5*(Qf + params$logdet[k] + d*log(2*pi))
  }
  lse <- rowLogSumExps(log_gamma)
  exp(log_gamma - lse)
}

#' Coordinate-wise M-step for EEI with separable prior
#'
#' @description
#' Per-coordinate M-step that is algebraically equivalent to your current
#' block M-step (`MSTEP`) when:
#'   (i) covariance model is EEI (shared diagonal across components),
#'   (ii) the smoothness prior is separable: Q_full = I_d ⊗ Q_K.
#'
#' It updates mu by solving d independent K×K linear systems and then
#' updates the shared diagonal covariance exactly as in your EEI branch,
#' including handling of `relative_lambda`, `rank_deficiency`, and `nugget`.
#'
#' @param data  (n × d) matrix of observations.
#' @param gamma (n × K) responsibility matrix; each row sums to 1.
#' @param params list with at least `sigma` (list of length K, each a d×d diag matrix).
#'               If `iterate_once = TRUE`, it should also be consistent with
#'               the current Sigma used to weight the mu-update, exactly as in `MSTEP`.
#' @param Q_K   (K × K) precision matrix acting along components for a single coordinate.
#'              Use NULL to indicate no penalization.
#' @param relative_lambda logical; if TRUE, scale the prior by 1/sigma_j^2 for each coordinate j,
#'                        matching your `MSTEP`'s `relative_lambda` behavior.
#' @param iterate_once    logical; if TRUE, perform a single μ update (weighted by current Sigma in params)
#'                        and one Σ update, then return. If FALSE and Q_K is not NULL,
#'                        run inner coordinate-wise μ/Σ updates until convergence.
#' @param nugget          nonnegative scalar added to diagonal variances (EEI).
#' @param rank_deficiency integer; used in the EEI variance denominator as (n + K - rank_deficiency)
#'                        when `relative_lambda = TRUE`, exactly as in your `MSTEP`.
#' @param tol_inner       tolerance for inner-loop convergence on (μ, Σ).
#' @param max_inner       maximum number of inner iterations.
#' @param verbose         logical; print inner-loop deltas if TRUE.
#'
#' @return list with elements:
#'   - pi    : length-K vector of mixing proportions
#'   - mu    : list of length K; each element is a length-d mean vector
#'   - sigma : list of length K; each element is the shared diagonal covariance (d × d)
#'
#' @details
#' Algebraic equivalence to the block update holds for EEI and Q_full = I_d ⊗ Q_K because
#' the μ update decouples across coordinates:
#'   [diag(Nk/sigma_j^2) + Q_use] μ_{j·} = (t(Gamma) x_{·j}) / sigma_j^2,
#' where Q_use = Q_K / sigma_j^2 if `relative_lambda=TRUE`, otherwise Q_use = Q_K.
#' The Σ update matches your EEI formulas verbatim.
#'
MSTEP_coordinatewise <- function(data, gamma, params,
                                 Q_K = NULL,
                                 relative_lambda = FALSE,
                                 iterate_once = TRUE,
                                 nugget = 0,
                                 rank_deficiency = 0,
                                 tol_inner = 1e-6,
                                 max_inner = 20,
                                 verbose = FALSE) {
  # shapes
  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)

  # --- mixing proportions (same as MSTEP)
  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # --- extract current shared-diagonal Sigma from params (EEI)
  if (is.null(params$sigma) || length(params$sigma) != K)
    stop("params$sigma must be a list of length K (EEI: all identical diag matrices).")
  Sigma0 <- params$sigma[[1]]
  if (!is.matrix(Sigma0) || nrow(Sigma0) != d || ncol(Sigma0) != d)
    stop("Each params$sigma[[k]] must be a d×d matrix.")
  if (!isTRUE(all(Sigma0 == diag(diag(Sigma0), d))))
    stop("EEI required: params$sigma must be diagonal and shared across components.")

  sigma_diag <- pmax(diag(Sigma0), .Machine$double.eps)  # length d

  # Precompute sums used by both μ- and Σ-updates (match MSTEP EEI branch)
  GX  <- t(data) %*% gamma      # d × K, element (j,k) = sum_i x_ij * gamma_ik
  GX2 <- t(data^2) %*% gamma    # d × K, element (j,k) = sum_i x_ij^2 * gamma_ik

  # Helper: one μ-solve for coordinate j
  solve_mu_j <- function(j, sigma2_j, GX_row_j) {
    if (is.null(Q_K)) {
      # no penalty: μ_{jk} = (sum_i gamma_ik x_ij) / Nk[k]
      return(GX_row_j / Nk)
    } else {
      # penalized: [diag(Nk/sigma2_j) + Q_use] μ_j· = GX[j,]/sigma2_j
      Q_use <- if (relative_lambda) Q_K / sigma2_j else Q_K
      A_diag <- Nk / sigma2_j
      A <- diag(A_diag, nrow = K, ncol = K) + Q_use
      b <- GX_row_j / sigma2_j

      # Cholesky solve: A is SPD
      L <- chol(A)                            # upper-triangular
      y <- forwardsolve(t(L), b)              # solve L^T y = b
      x <- backsolve(L, y)                    # solve L x = y
      return(as.numeric(x))
    }
  }

  # Inner loop if needed (consistent with MSTEP behavior)
  mu_mat <- matrix(0, nrow = d, ncol = K)  # rows j=1..d, cols k=1..K
  inner_max <- if (!is.null(Q_K) && !iterate_once) max_inner else 1

  for (inner in seq_len(inner_max)) {
    mu_old    <- as.vector(t(mu_mat))
    sigma_old <- sigma_diag

    # ---- μ update: d independent K×K solves
    for (j in seq_len(d)) {
      mu_mat[j, ] <- solve_mu_j(j, sigma_diag[j], GX[j, ])
    }

    # ---- Σ update (EEI), verbatim match to your MSTEP EEI branch
    # V_j = sum_k [ GX2[j,k] - 2 * mu_{j,k} * GX[j,k] + mu_{j,k}^2 * Nk[k] ]
    term1 <- rowSums(GX2)                              # length d
    term2 <- 2 * rowSums(mu_mat * GX)                  # length d
    mu2Nk <- rowSums( sweep(mu_mat^2, 2, Nk, "*") )    # length d
    V <- term1 - term2 + mu2Nk                         # length d

    if (relative_lambda && !is.null(Q_K)) {
      # add prior contribution U_vec_d[j] = μ_{j·}^T Q_K μ_{j·}
      U_vec_d <- vapply(seq_len(d), function(j) {
        mj <- mu_mat[j, ]
        as.numeric(crossprod(mj, Q_K %*% mj))
      }, numeric(1))
      sigma_diag <- (V + U_vec_d) / (n + K - rank_deficiency) + nugget
    } else {
      sigma_diag <- V / n + nugget
    }
    sigma_diag <- pmax(sigma_diag, .Machine$double.eps)

    # ---- convergence check (same style as MSTEP)
    mu_diff    <- max(abs(as.vector(t(mu_mat)) - mu_old))
    sigma_diff <- max(abs(sigma_diag - sigma_old))
    if (verbose && inner_max > 1L) {
      cat(sprintf("  CW MSTEP inner %d: Δμ=%.3e, ΔΣ=%.3e\n", inner, mu_diff, sigma_diff))
    }
    if (inner_max > 1L && max(mu_diff, sigma_diff) < tol_inner) break
  }

  # Package outputs in the same structure as MSTEP
  mu_list <- lapply(seq_len(K), function(k) mu_mat[, k])
  Sigma_shared <- diag(sigma_diag, d)
  sigma_list  <- replicate(K, Sigma_shared, simplify = FALSE)

  list(pi = pi_new, mu = mu_list, sigma = sigma_list)
}





#' Single-coordinate M-step (1D), return mean/residual AND mixing proportions
#' (EEI + separable prior along K; consistent with MSTEP_coordinatewise)
#'
#' @description
#' Solve the M-step for ONE coordinate j = `coord` under:
#'   (i) EEI covariance (shared diagonal across components),
#'   (ii) separable prior along K only: Q_full = Q_K ⊗ I_d (pass Q_K here).
#' It updates μ_{j·} using the same algebra as MSTEP_coordinatewise,
#' optionally runs a small inner loop to update sigma_j^2 consistently,
#' and returns the K-vector mean for this coordinate, the n-vector residuals
#' r = x_{·j} - Gamma %*% μ_{j·}, the scalar sigma2_j, and the mixing proportions pi.
#'
#' @param data  (n × d) matrix of observations.
#' @param gamma (n × K) responsibility matrix (rows sum to 1).
#' @param params list with:
#'   - `mu`    : list length K, each a length-d mean vector (used only for shape checks).
#'   - `sigma` : list length K of identical d×d diagonal matrices (EEI); used to weight μ-update.
#' @param coord integer in {1,...,d}: which coordinate to solve.
#' @param Q_K   (K × K) precision along components for a single coordinate (NULL = no penalty).
#' @param relative_lambda logical; if TRUE, scale Q_K by 1/sigma_j^2 for this coordinate.
#' @param iterate_once    logical; if TRUE, do one μ update (weighted by current sigma_j^2)
#'                        and compute residuals; if FALSE and Q_K is not NULL,
#'                        alternate μ_j and sigma_j^2 for this coordinate until convergence.
#' @param nugget          nonnegative scalar added to sigma_j^2.
#' @param rank_deficiency integer used in (n + K - rank_deficiency) when `relative_lambda=TRUE` and `Q_K` non-NULL.
#' @param tol_inner       inner-loop tolerance for this coordinate (if iterate_once=FALSE and Q_K non-NULL).
#' @param max_inner       inner-loop maximum iterations for this coordinate.
#' @param verbose         print inner-loop deltas if TRUE.
#'
#' @return list with
#'   - `mu_j`   : length-K numeric vector (means for this coordinate),
#'   - `resid`  : length-n numeric vector, residuals r = x_j - Gamma %*% mu_j,
#'   - `sigma2_j`: scalar variance estimate for this coordinate (EEI),
#'   - `pi`     : length-K numeric vector of mixing proportions (Nk / n).
#'
MSTEP_coordinatewise_once <- function(data, gamma, params, coord,
                                      Q_K = NULL,
                                      relative_lambda = FALSE,
                                      iterate_once = TRUE,
                                      nugget = 0,
                                      rank_deficiency = 0,
                                      tol_inner = 1e-6,
                                      max_inner = 20,
                                      verbose = FALSE) {

  n <- nrow(data); d <- ncol(data); K <- ncol(gamma)
  if (length(coord) != 1L || coord < 1L || coord > d) stop("`coord` must be an integer in 1..d.")
  if (is.null(params$sigma) || length(params$sigma) != K) stop("params$sigma must be a list of length K (EEI shared diag).")

  # Mixing proportions (global; independent of coord)
  Nk <- colSums(gamma); Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # Check EEI: shared diagonal
  Sigma0 <- params$sigma[[1]]
  if (!is.matrix(Sigma0) || nrow(Sigma0) != d || ncol(Sigma0) != d)
    stop("Each params$sigma[[k]] must be a d×d matrix.")
  if (!isTRUE(all(Sigma0 == diag(diag(Sigma0), d))))
    stop("EEI required: sigma must be diagonal and shared across components.")
  for (k in seq_len(K)) {
    if (!isTRUE(all(params$sigma[[k]] == Sigma0))) stop("EEI requires all params$sigma[[k]] identical.")
  }

  # Current variance for this coordinate
  sigma_diag <- pmax(diag(Sigma0), .Machine$double.eps)
  sigma2_j <- sigma_diag[coord]

  # Sufficient stats for this coordinate
  xj  <- data[, coord]
  gx  <- as.numeric(t(gamma) %*% xj)       # length K: sum_i gamma_ik x_ij
  gx2 <- as.numeric(t(gamma) %*% (xj^2))   # length K: sum_i gamma_ik x_ij^2

  # Inner loop if penalized and iterate_once=FALSE; otherwise single shot
  inner_max <- if (!is.null(Q_K) && !iterate_once) max_inner else 1L
  mu_j <- gx / Nk  # unpenalized initializer

  for (inner in seq_len(inner_max)) {
    mu_old    <- mu_j
    sigma_old <- sigma2_j

    # μ_j update: [diag(Nk/sigma2_j) + Q_use] mu_j = gx / sigma2_j
    if (is.null(Q_K)) {
      mu_j <- gx / Nk
    } else {
      Q_use <- if (relative_lambda) Q_K / sigma2_j else Q_K
      A <- diag(Nk / sigma2_j, nrow = K, ncol = K) + Q_use
      b <- gx / sigma2_j
      # Cholesky solve (SPD)
      L <- chol(A)
      y <- forwardsolve(t(L), b)
      mu_j <- as.numeric(backsolve(L, y))
    }

    # sigma_j^2 update (EEI, 1-D)
    V_j <- sum(gx2 - 2 * mu_j * gx + (mu_j^2) * Nk)
    if (!is.null(Q_K) && relative_lambda) {
      U_j <- as.numeric(crossprod(mu_j, Q_K %*% mu_j))
      sigma2_j <- (V_j + U_j) / (n + K - rank_deficiency) + nugget
    } else {
      sigma2_j <- V_j / n + nugget
    }
    sigma2_j <- max(sigma2_j, .Machine$double.eps)

    # Convergence
    mu_diff    <- max(abs(mu_j - mu_old))
    sigma_diff <- abs(sigma2_j - sigma_old)
    if (verbose && inner_max > 1L) {
      cat(sprintf("  1D MSTEP coord %d, inner %d: Δμ=%.3e  Δσ²=%.3e\n", coord, inner, mu_diff, sigma_diff))
    }
    if (inner_max > 1L && max(mu_diff, sigma_diff) < tol_inner) break
  }

  resid <- xj - as.numeric(gamma %*% mu_j)

  list(mu_j = mu_j, resid = resid, sigma2_j = sigma2_j, pi = pi_new)
}


#' Univariate penalized ELBO for a single coordinate
#'
#' @description
#' Compute the penalized ELBO contribution for one coordinate j = `coord`,
#' mirroring the logic of `compute_penalized_ELBO`:
#'   ELBO_j = sum_{i,k} Gamma_{ik} [ log N(x_{ij}|mu_{jk}, sigma_{jk}^2) + (log pi_k if included) ]
#'            + ( - sum_{i,k} Gamma_{ik} log Gamma_{ik}  if included )
#'            + ( -1/2 * mu_j' Q_j mu_j + 1/2 * logdet(Q_j)  if prior provided )
#'
#' Notes:
#' - If you plan to SUM over j to recover the global ELBO, set
#'   `include_log_pi = FALSE` and `include_entropy = FALSE`, then add
#'   the global terms once:
#'     sum_j ELBO_j  +  sum_{i,k} Gamma_{ik} log pi_k  +  ( - sum_{i,k} Gamma_{ik} log Gamma_{ik} ).
#' - For a general (non-separable) prior Q_full, this routine uses the
#'   coordinate-specific KxK sub-block (k-major stacking). Summing across j
#'   only matches the global prior if Q_full = Q_K ⊗ I_d (separable).
#'
#' @param X      (n × d) data matrix.
#' @param coord  integer in 1..d: which coordinate to score.
#' @param Gamma  (n × K) responsibility matrix.
#' @param params list with:
#'   - mu    : list length K, each a length-d vector,
#'   - sigma : list length K, each a d×d covariance (we use [coord,coord]),
#'   - pi    : length-K mixing proportions.
#' @param Q_full_prior (optional) (dK × dK) precision matrix for the prior (already scaled if using relative-lambda).
#'                     If `NULL`, no prior term is added.
#'                     For separable prior, you may pass NULL and instead provide `Q_K`.
#' @param Q_K          (optional) (K × K) precision along components for one coordinate
#'                     (use when Q_full = Q_K ⊗ I_d; ignored if `Q_full_prior` is provided).
#' @param include_log_pi   logical; include the log pi term (default TRUE, matching global ELBO).
#' @param include_entropy  logical; include the entropy term (default TRUE, matching global ELBO).
#' @param eigen_tol        numeric; eigenvalue threshold passed to `generalized_logdet`.
#' @param rank_deficiency  integer; rank deficiency for `generalized_logdet` if the prior is improper (e.g., RW2).
#'
#' @return numeric scalar: penalized ELBO for coordinate `coord`.
#'
compute_penalized_ELBO_univariate <- function(
    X, coord, Gamma, params,
    Q_full_prior = NULL, Q_K = NULL,
    include_log_pi = TRUE,
    include_entropy = TRUE,
    eigen_tol = 0,
    rank_deficiency = 0
) {
  n <- nrow(X); d <- ncol(X); K <- ncol(Gamma)
  if (coord < 1L || coord > d) stop("coord must be in 1..d")

  # 1) Collect per-component means for this coordinate and variances
  mu_j  <- vapply(seq_len(K), function(k) params$mu[[k]][coord], numeric(1))
  sig2j <- vapply(seq_len(K), function(k) params$sigma[[k]][coord, coord], numeric(1))
  sd_j  <- sqrt(pmax(sig2j, .Machine$double.eps))

  # 2) Q-function (data) term: sum_{i,k} gamma_{ik} log N(x_ij | mu_jk, sigma_jk^2)
  xj <- X[, coord]
  log_dens_mat <- sapply(seq_len(K), function(k) {
    dnorm(xj, mean = mu_j[k], sd = sd_j[k], log = TRUE)
  }) # n × K

  log_gamma <- log_dens_mat
  if (include_log_pi) {
    log_gamma <- sweep(log_gamma, 2, log(params$pi), "+")
  }
  Qval <- sum(Gamma * log_gamma)

  # 3) Entropy term (optional)
  elbo <- Qval
  if (include_entropy) {
    entropy <- - sum(Gamma * log(Gamma + 1e-300))
    elbo <- elbo + entropy
  }

  # 4) Prior term for this coordinate (optional)
  if (!is.null(Q_full_prior) || !is.null(Q_K)) {
    Qj <- if (!is.null(Q_full_prior)) {
      # extract K×K sub-block for this coord under k-major stacking
      # indices: (k-1)*d + coord, for k=1..K
      idx <- (seq_len(K) - 1L) * d + coord
      Q_full_prior[idx, idx, drop = FALSE]
    } else {
      Q_K
    }
    # penalty and normalization constant (+ 0.5 * log|Qj|)
    penalty <- 0.5 * as.numeric(crossprod(mu_j, Qj %*% mu_j))
    logdetQ <- generalized_logdet(Q = Qj, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    elbo <- elbo - penalty + 0.5 * as.numeric(logdetQ)
  }

  return(elbo)
}










