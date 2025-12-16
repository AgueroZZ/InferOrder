# =============================================================================
# Full penalized ELBO for one ordering/subspace
# =============================================================================
compute_ordering_elbo <- function(Xsub, Gamma, params, Q_K,
                                  relative_lambda = TRUE,
                                  eigen_tol = 0, rw_order = 2) {
  dsub <- ncol(Xsub)
  if (dsub == 0L) return(NA_real_)
  Q_full <- kronecker(Q_K, diag(dsub))
  if (isTRUE(relative_lambda)) {
    sigma_diag <- diag(params$sigma[[1]])
    scale_vec  <- rep(1 / sqrt(pmax(sigma_diag, .Machine$double.eps)), times = ncol(Gamma))
    Q_full     <- Q_full * (scale_vec %o% scale_vec)
  }
  rank_def_full <- rw_order * dsub
  compute_penalized_ELBO(
    X = Xsub, Gamma = Gamma, params = params,
    Q_prior = Q_full, eigen_tol = eigen_tol, rank_deficiency = rank_def_full
  )
}



## ============================================================
## Helper 1: fit a single smooth-EM ordering on a subset of features
## ============================================================
#' Fit a single smooth-EM model on a subset of coordinates (one ordering)
#'
#' @description
#' Runs a smooth EM algorithm on a sub-matrix \code{Xsub} (n × d_sub) with:
#' - K-component Gaussian mixture (EEI covariance: shared diagonal),
#' - separable smoothness prior along components: Q_full = Q_K ⊗ I_{d_sub},
#' - coordinate-wise M-step via \code{MSTEP_coordinatewise}.
#'
#' This is used as a building block inside the backward-splitting algorithm:
#'   - First, to fit a single ordering on all features,
#'   - Then, repeatedly to refit the two orderings on their current feature sets.
#'
#' @param Xsub (n × d_sub) data matrix for a subset of coordinates.
#' @param K integer, number of mixture components.
#' @param Q_K (K × K) precision matrix along mixture components for a single coordinate.
#' @param relative_lambda logical; if TRUE, scale the prior by 1/σ_j² for each coordinate,
#'   matching the behavior of \code{MSTEP_coordinatewise}.
#' @param rank_deficiency integer; rank deficiency of Q_K (e.g., RW order) used when
#'   computing penalized ELBO.
#' @param max_iter integer; maximum number of outer EM iterations.
#' @param init_method character; one of \code{"random"} or \code{"pca"} specifying how to
#'   initialize the parameters.
#' @param nugget nonnegative scalar added to diagonal variances in the EEI covariance.
#' @param tol_inner numeric; tolerance for inner coordinate-wise M-step updates.
#' @param max_inner integer; maximum number of inner coordinate-wise updates.
#' @param eigen_tol numeric; eigenvalue threshold for generalized log-determinant in
#'   the penalized ELBO.
#' @param verbose logical; if TRUE, print ELBO per EM iteration.
#'
#' @return A list with:
#'   \item{Gamma}{(n × K) responsibility matrix at convergence.}
#'   \item{params}{List of GMM parameters \code{(pi, mu, sigma)} at convergence.}
#'   \item{elbo}{Final penalized ELBO value (including prior and normalizing constants).}
#'
local_smooth_em_one_ordering <- function(
    Xsub,
    K,
    Q_K,
    relative_lambda = TRUE,
    rank_deficiency = 0,
    max_iter = 20,
    init_method = c("random","pca"),
    nugget = 0,
    tol_inner = 1e-6,
    max_inner = 20,
    eigen_tol = 0,
    verbose = FALSE
) {
  init_method <- match.arg(init_method)
  n <- nrow(Xsub); d_sub <- ncol(Xsub)

  ## --- initialization of params ---
  if (init_method == "random") {
    params <- make_default_init(Xsub, K, ordering = TRUE)
  } else if (init_method == "pca") {
    pc1 <- prcomp(Xsub, center = TRUE, scale. = TRUE)$x[, 1]
    params <- make_init(Xsub, ordering_vec = pc1, K = K)
  }

  Gamma <- matrix(1 / K, nrow = n, ncol = K)

  ## Use penalized ELBO (including prior + constants) for convergence check
  elbo_prev <- -Inf
  elbo_curr <- NA_real_

  for (iter in seq_len(max_iter)) {
    ## E-step
    params_cached <- init_cov_cache_fast(params)
    Gamma <- ESTEP(Xsub, params_cached)

    ## M-step: coordinate-wise with separable prior Q_K
    params <- MSTEP_coordinatewise(
      data            = Xsub,
      gamma           = Gamma,
      params          = params,
      Q_K             = Q_K,
      relative_lambda = relative_lambda,
      iterate_once    = FALSE,
      nugget          = nugget,
      rank_deficiency = rank_deficiency,
      tol_inner       = tol_inner,
      max_inner       = max_inner,
      verbose         = verbose && (max_iter > 1)
    )

    ## Penalized ELBO (includes Q_K, log det(Q_K), etc.)
    elbo_curr <- compute_ordering_elbo(
      Xsub          = Xsub,
      Gamma         = Gamma,
      params        = params,
      Q_K           = Q_K,
      relative_lambda = relative_lambda,
      eigen_tol     = eigen_tol,
      rw_order      = rank_deficiency
    )

    if (verbose) {
      cat(sprintf("  [one-ordering EM] iter %d, ELBO = %.3f\n", iter, elbo_curr))
    }

    if (abs(elbo_curr - elbo_prev) < 1e-4) break
    elbo_prev <- elbo_curr
  }

  list(
    Gamma = Gamma,
    params = params,
    elbo = elbo_curr
  )
}


## ------------------------------------------------------------
## Helper 2: score ONE coordinate under ONE ordering
##   (1D M-step + 1D penalized ELBO)
## ------------------------------------------------------------
score_coord_under_ordering <- function(
    xj,             # length n
    Gamma,          # n x K responsibilities
    Q_K,            # K x K prior along components
    relative_lambda = TRUE,
    rank_deficiency = 0,
    eigen_tol = 0,
    nugget = 0,
    tol_inner = 1e-6,
    max_inner = 20
) {
  n <- length(xj)
  K <- ncol(Gamma)

  X_1d <- matrix(xj, ncol = 1)

  # 1D initialization on this coordinate
  init1d <- make_default_init(X_1d, K, ordering = TRUE)

  # Single-coordinate M-step to get optimal mu_j, sigma2_j under (Gamma, Q_K)
  one <- MSTEP_coordinatewise_once(
    data            = X_1d,
    gamma           = Gamma,
    params          = init1d,
    coord           = 1,
    Q_K             = Q_K,
    relative_lambda = relative_lambda,
    iterate_once    = FALSE,
    nugget          = nugget,
    rank_deficiency = rank_deficiency,
    tol_inner       = tol_inner,
    max_inner       = max_inner,
    verbose         = FALSE
  )

  mu_j    <- one$mu_j          # length K
  sigma2j <- one$sigma2_j      # scalar
  pi_k    <- one$pi            # length K

  # relative_lambda: scale Q_K by 1/sigma2_j
  Qj <- if (relative_lambda) Q_K / sigma2j else Q_K

  # Build params for 1D ELBO
  params_j <- list(
    mu    = lapply(seq_len(K), function(k) mu_j[k]),
    sigma = lapply(seq_len(K), function(k) matrix(sigma2j, nrow = 1, ncol = 1)),
    pi    = pi_k
  )

  elbo_j <- compute_penalized_ELBO_univariate(
    X              = X_1d,
    coord          = 1,
    Gamma          = Gamma,
    params         = params_j,
    Q_full_prior   = NULL,
    Q_K            = Qj,
    include_log_pi  = F,  # per-coordinate comparison only
    include_entropy = F,
    eigen_tol       = eigen_tol,
    rank_deficiency = rank_deficiency
  )

  return(elbo_j)
}





## ============================================================
## Main: backward split of one ordering into two orderings
## ============================================================

#' Backward split of a single ordering into two smooth orderings
#'
#' @description
#' This function implements a **backward-style** two-ordering smooth-EM strategy.
#' The idea:
#' \enumerate{
#'   \item Fit a single smooth-EM model (one ordering) on all coordinates.
#'   \item Score each coordinate by its penalized ELBO under that ordering.
#'         The worst (smallest ELBO) coordinate is used to seed a second ordering.
#'   \item Fit two separate smooth-EM models on \eqn{J_1} and \eqn{J_2} (feature sets
#'         for ordering 1 and 2).
#'   \item Perform "backfitting" iterations: at each step,
#'         \itemize{
#'           \item Refit each ordering on its current feature subset.
#'           \item For each coordinate, compute a per-coordinate penalized ELBO
#'                 under both orderings.
#'           \item Move the coordinate to the other ordering if the ELBO gain exceeds
#'                 \code{move_tol}, while ensuring neither ordering becomes empty.
#'         }
#'   \item Stop when no coordinates change ordering or when the maximum number of
#'         backfitting iterations is reached.
#' }
#'
#' The per-coordinate scoring under an ordering is done via:
#'   \itemize{
#'     \item a 1D penalized M-step (Gaussian, EEI, separable prior along K),
#'     \item followed by a 1D penalized ELBO evaluation (excluding global log-pi
#'           and entropy terms, so scores are comparable across orderings).
#'   }
#'
#' @param X (n × d) data matrix.
#' @param K integer; number of mixture components for each ordering.
#' @param Q1_K (K × K) precision along components for ordering 1.
#' @param Q2_K (K × K) precision along components for ordering 2.
#'        Defaults to \code{Q1_K} if not provided.
#' @param relative_lambda logical; if TRUE, the prior for each coordinate is scaled
#'        by 1/σ_j² (same behavior as in \code{MSTEP_coordinatewise_once}).
#' @param rank_def1 integer; rank deficiency of \code{Q1_K} (e.g., RW order).
#' @param rank_def2 integer; rank deficiency of \code{Q2_K}.
#' @param max_iter_single integer; EM iterations when fitting the single full ordering
#'        (step 0).
#' @param max_iter_per_order integer; EM iterations for each ordering during the split
#'        and backfitting steps.
#' @param max_backfit_iter integer; maximum number of backfitting iterations (outer loop).
#' @param move_tol numeric; minimum ELBO gain required to move a coordinate from one
#'        ordering to the other. Use 0 for any improvement.
#' @param init_method character; initialization method passed to
#'        \code{local_smooth_em_one_ordering}, one of \code{"random"} or \code{"pca"}.
#' @param nugget nonnegative scalar added to variances in the EEI covariance.
#' @param tol_inner numeric; tolerance for inner M-step loops (used in 1D updates).
#' @param max_inner integer; maximum inner loop iterations for coordinate-wise M-steps.
#' @param eigen_tol numeric; eigenvalue threshold for generalized log-determinant.
#' @param verbose logical; if TRUE, prints progress and details of coordinate moves.
#'
#' @return A list with:
#'   \item{coord_assign}{Integer vector of length d, with values in \{1, 2\} indicating
#'                       the final ordering assignment for each coordinate.}
#'   \item{J1}{Indices of coordinates assigned to ordering 1.}
#'   \item{J2}{Indices of coordinates assigned to ordering 2.}
#'   \item{Gamma1}{(n × K) responsibility matrix for ordering 1 at the final fit.}
#'   \item{Gamma2}{(n × K) responsibility matrix for ordering 2 at the final fit.}
#'   \item{params1}{GMM parameters for ordering 1 at the final fit.}
#'   \item{params2}{GMM parameters for ordering 2 at the final fit.}
#'   \item{elbo1}{Final penalized ELBO for ordering 1 (on \code{X[, J1]}).}
#'   \item{elbo2}{Final penalized ELBO for ordering 2 (on \code{X[, J2]}).}
#'   \item{elbo_total}{Sum \code{elbo1 + elbo2}.}
#'   \item{history}{List over backfitting iterations; each element contains:
#'                  \code{step}, \code{n_changes}, \code{moved_idx},
#'                  \code{from}, and \code{to}.}
#'
split_ordering_backward <- function(
    X,
    K,
    Q1_K,
    Q2_K = Q1_K,
    relative_lambda = TRUE,
    rank_def1 = 0,
    rank_def2 = 0,
    max_iter_single = 20,
    max_iter_per_order = 20,
    max_backfit_iter = 10,
    move_tol = 0,
    init_method = c("random","pca"),
    nugget = 0,
    tol_inner = 1e-6,
    max_inner = 20,
    eigen_tol = 0,
    verbose = TRUE
) {
  init_method <- match.arg(init_method)
  n <- nrow(X); d <- ncol(X)

  if (verbose) cat("== Backward split: step 0, fit ONE ordering with all features ==\n")

  # 0) Fit a single smooth-EM model on all features -> seed Gamma1
  fit1_full <- local_smooth_em_one_ordering(
    Xsub            = X,
    K               = K,
    Q_K             = Q1_K,
    relative_lambda = relative_lambda,
    rank_deficiency = rank_def1,
    max_iter        = max_iter_single,
    init_method     = init_method,
    nugget          = nugget,
    tol_inner       = tol_inner,
    max_inner       = max_inner,
    eigen_tol       = eigen_tol,
    verbose         = verbose
  )
  Gamma1 <- fit1_full$Gamma

  # 1) Score each coordinate j under this single ordering
  if (verbose) cat("== Step 1: score all coordinates under single ordering ==\n")
  elbo1_vec <- vapply(seq_len(d), function(j) {
    score_coord_under_ordering(
      xj              = X[, j],
      Gamma           = Gamma1,
      Q_K             = Q1_K,
      relative_lambda = relative_lambda,
      rank_deficiency = rank_def1,
      eigen_tol       = eigen_tol,
      nugget          = nugget,
      tol_inner       = tol_inner,
      max_inner       = max_inner
    )
  }, numeric(1))

  # Worst coordinate under the single ordering becomes seed for ordering 2
  j_worst <- which.min(elbo1_vec)
  if (verbose) {
    cat(sprintf("  Worst coordinate under single ordering: j = %d (ELBO = %.3f)\n",
                j_worst, elbo1_vec[j_worst]))
  }

  J2 <- j_worst
  J1 <- setdiff(seq_len(d), j_worst)

  # 2) Initial: fit two orderings on J1 and J2
  if (verbose) cat("== Step 2: fit two orderings on J1 / J2 ==\n")

  fit_ord1 <- local_smooth_em_one_ordering(
    Xsub            = X[, J1, drop = FALSE],
    K               = K,
    Q_K             = Q1_K,
    relative_lambda = relative_lambda,
    rank_deficiency = rank_def1,
    max_iter        = max_iter_per_order,
    init_method     = init_method,
    nugget          = nugget,
    tol_inner       = tol_inner,
    max_inner       = max_inner,
    eigen_tol       = eigen_tol,
    verbose         = verbose
  )
  Gamma1  <- fit_ord1$Gamma
  params1 <- fit_ord1$params

  fit_ord2 <- local_smooth_em_one_ordering(
    Xsub            = X[, J2, drop = FALSE],
    K               = K,
    Q_K             = Q2_K,
    relative_lambda = relative_lambda,
    rank_deficiency = rank_def2,
    max_iter        = max_iter_per_order,
    init_method     = init_method,
    nugget          = nugget,
    tol_inner       = tol_inner,
    max_inner       = max_inner,
    eigen_tol       = eigen_tol,
    verbose         = verbose
  )
  Gamma2  <- fit_ord2$Gamma
  params2 <- fit_ord2$params

  coord_assign <- integer(d)
  coord_assign[J1] <- 1L
  coord_assign[J2] <- 2L

  # 3) Backfitting with backward-style coordinate moves
  history <- list()
  if (verbose) cat("== Step 3: backfitting with backward-style coordinate moves ==\n")

  for (bf in seq_len(max_backfit_iter)) {
    if (verbose) cat(sprintf("[Backfitting %02d] Refit both orderings.\n", bf))

    # 3a) Refit both orderings on current J1 / J2 (fresh EM from current subsets)
    if (length(J1) > 0) {
      fit_ord1 <- local_smooth_em_one_ordering(
        Xsub            = X[, J1, drop = FALSE],
        K               = K,
        Q_K             = Q1_K,
        relative_lambda = relative_lambda,
        rank_deficiency = rank_def1,
        max_iter        = max_iter_per_order,
        init_method     = init_method,
        nugget          = nugget,
        tol_inner       = tol_inner,
        max_inner       = max_inner,
        eigen_tol       = eigen_tol,
        verbose         = verbose && FALSE
      )
      Gamma1  <- fit_ord1$Gamma
      params1 <- fit_ord1$params
    }

    if (length(J2) > 0) {
      fit_ord2 <- local_smooth_em_one_ordering(
        Xsub            = X[, J2, drop = FALSE],
        K               = K,
        Q_K             = Q2_K,
        relative_lambda = relative_lambda,
        rank_deficiency = rank_def2,
        max_iter        = max_iter_per_order,
        init_method     = init_method,
        nugget          = nugget,
        tol_inner       = tol_inner,
        max_inner       = max_inner,
        eigen_tol       = eigen_tol,
        verbose         = verbose && FALSE
      )
      Gamma2  <- fit_ord2$Gamma
      params2 <- fit_ord2$params
    }

    # 3b) For each coordinate, compute per-coordinate ELBO under both orderings
    new_assign <- coord_assign

    for (j in seq_len(d)) {
      xj <- X[, j]

      # ELBO under ordering 1 (if non-empty)
      elbo1_j <- -Inf
      if (length(J1) > 0) {
        elbo1_j <- score_coord_under_ordering(
          xj              = xj,
          Gamma           = Gamma1,
          Q_K             = Q1_K,
          relative_lambda = relative_lambda,
          rank_deficiency = rank_def1,
          eigen_tol       = eigen_tol,
          nugget          = nugget,
          tol_inner       = tol_inner,
          max_inner       = max_inner
        )
      }

      # ELBO under ordering 2 (if non-empty)
      elbo2_j <- -Inf
      if (length(J2) > 0) {
        elbo2_j <- score_coord_under_ordering(
          xj              = xj,
          Gamma           = Gamma2,
          Q_K             = Q2_K,
          relative_lambda = relative_lambda,
          rank_deficiency = rank_def2,
          eigen_tol       = eigen_tol,
          nugget          = nugget,
          tol_inner       = tol_inner,
          max_inner       = max_inner
        )
      }

      cur <- coord_assign[j]

      # Decide whether to move j, enforcing non-empty constraints
      if (cur == 1L) {
        # avoid emptying J1
        if (length(J1) > 1 && (elbo2_j > elbo1_j + move_tol)) {
          new_assign[j] <- 2L
        }
      } else if (cur == 2L) {
        # avoid emptying J2
        if (length(J2) > 1 && (elbo1_j > elbo2_j + move_tol)) {
          new_assign[j] <- 1L
        }
      } else {
        # fallback (should not really happen)
        new_assign[j] <- if (elbo1_j >= elbo2_j) 1L else 2L
      }
    }

    # Record which coordinates changed and from where to where
    moved_idx <- which(new_assign != coord_assign)
    from_vec  <- coord_assign[moved_idx]
    to_vec    <- new_assign[moved_idx]
    n_changes <- length(moved_idx)

    history[[length(history) + 1L]] <- list(
      step      = bf,
      n_changes = n_changes,
      moved_idx = moved_idx,
      from      = from_vec,
      to        = to_vec
    )

    if (verbose) {
      cat(sprintf("[Backfitting %02d] reassigned %d coords\n", bf, n_changes))
      if (n_changes > 0) {
        msg <- paste(
          sprintf("j=%d: %d->%d", moved_idx, from_vec, to_vec),
          collapse = ", "
        )
        cat("   Moves: ", msg, "\n")
      }
    }

    coord_assign <- new_assign
    J1 <- which(coord_assign == 1L)
    J2 <- which(coord_assign == 2L)

    if (n_changes == 0) break
  }

  # Final ELBOs for each ordering on their respective feature sets
  elbo1_final <- if (length(J1) > 0) {
    compute_ordering_elbo(
      Xsub          = X[, J1, drop = FALSE],
      Gamma         = Gamma1,
      params        = params1,
      Q_K           = Q1_K,
      relative_lambda = relative_lambda,
      eigen_tol     = eigen_tol,
      rw_order      = rank_def1
    )
  } else {
    NA_real_
  }

  elbo2_final <- if (length(J2) > 0) {
    compute_ordering_elbo(
      Xsub          = X[, J2, drop = FALSE],
      Gamma         = Gamma2,
      params        = params2,
      Q_K           = Q2_K,
      relative_lambda = relative_lambda,
      eigen_tol     = eigen_tol,
      rw_order      = rank_def2
    )
  } else {
    NA_real_
  }

  list(
    coord_assign = coord_assign,
    J1 = J1,
    J2 = J2,
    Gamma1 = Gamma1,
    Gamma2 = Gamma2,
    params1 = params1,
    params2 = params2,
    elbo1 = elbo1_final,
    elbo2 = elbo2_final,
    elbo_total = elbo1_final + elbo2_final,
    history = history
  )
}


