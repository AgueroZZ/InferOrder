library(Matrix)
library(mclust)

#' Match locations on an old grid to indices on a new grid
#'
#' @description
#' Given a vector of locations \code{loc_old} (typically a subset of \code{loc_new}),
#' return the indices of their nearest matches on \code{loc_new}. This is mainly used
#' for hierarchical/nested grids where coarse-grid locations appear exactly on a
#' finer grid.
#'
#' The function checks that each \code{loc_old} can be matched within tolerance
#' and that no duplicate indices occur.
#'
#' @param loc_old Numeric vector of locations to be matched (e.g. \code{u_obs}).
#' @param loc_new Numeric vector of candidate locations (e.g. \code{u_final}).
#' @param tol Nonnegative numeric tolerance. A location \code{x} in \code{loc_old}
#'   must satisfy \code{min(abs(loc_new - x)) <= tol} to be considered a valid match.
#'
#' @return Integer vector of the same length as \code{loc_old}, giving indices into \code{loc_new}.
#'
#' @examples
#' u_final <- seq(0, 1, length.out = 65)
#' u_obs   <- c(0, 0.5, 1)
#' match_locations_to_grid(u_obs, u_final)
#'
match_locations_to_grid <- function(loc_old, loc_new, tol = 1e-8) {
  idx_obs <- vapply(loc_old, function(x) {
    j <- which.min(abs(loc_new - x))
    if (abs(loc_new[j] - x) > tol) {
      stop(sprintf("Cannot match location %.17g to new grid within tol=%.3e", x, tol))
    }
    j
  }, integer(1))
  if (any(duplicated(idx_obs))) stop("Duplicated indices in matching.")
  idx_obs
}


#' Construct a hierarchy of nested grids inside a final grid
#'
#' @description
#' Build a set of *hierarchical* (nested) grid levels indexed into a final grid.
#' The final grid has \code{K_final = 2^m_max + 1} equally-spaced points on \eqn{[0,1]}.
#' Each level \code{m} has \code{K_m = 2^m + 1} points, and its indices are a subset of
#' the final grid indices.
#'
#' This is useful for progressive-resolution initialization: you can fit at a coarse
#' level and krige/interpolate to initialize the next, finer level while guaranteeing
#' exact overlap of design points.
#'
#' @param m_max Nonnegative integer; the finest level exponent.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{m_max}: the input \code{m_max}.
#'   \item \code{K_final}: \code{2^m_max + 1}.
#'   \item \code{u_final}: numeric vector of length \code{K_final} on \eqn{[0,1]}.
#'   \item \code{idx_levels}: list of length \code{m_max+1}; indices into \code{u_final}
#'         for each level \code{m=0,...,m_max}.
#'   \item \code{K_levels}: integer vector of \code{K_m} for each level.
#' }
#'
#' @examples
#' lev <- make_hierarchical_levels(m_max = 4)
#' lev$K_levels
#' lev$idx_levels[[1]]  # level m=0 -> 2 points (ends)
#'
make_hierarchical_levels <- function(m_max = 6) {
  K_final <- 2^m_max + 1
  u_final <- seq(0, 1, length.out = K_final)

  idx_levels <- lapply(0:m_max, function(m) {
    step <- 2^(m_max - m)
    seq(1, K_final, by = step)
  })
  K_levels <- vapply(idx_levels, length, integer(1))

  list(
    m_max = m_max,
    K_final = K_final,
    u_final = u_final,
    idx_levels = idx_levels,  # list of indices into final grid
    K_levels = K_levels       # K_m = 2^m + 1
  )
}


#' Scale random-walk penalty strength to account for grid spacing
#'
#' @description
#' In a random-walk prior/penalty of order \code{q} (e.g., RW2), if the precision builder
#' \code{make_random_walk_precision()} implicitly assumes unit spacing in index space,
#' then changing the physical grid spacing \code{h} changes the implied roughness penalty.
#'
#' This helper rescales the penalty strength \code{lambda} so that the *continuous*
#' roughness penalty stays comparable across nested grids on \eqn{[0,1]}.
#'
#' For nested equally-spaced grids with \code{K} points on \eqn{[0,1]}, spacing is
#' \code{h = 1/(K-1)}. With order \code{q}, a common scaling is
#' \deqn{\lambda(h) \propto h^{-(2q-1)}.}
#' Therefore, relative to a final grid \code{K_final}, the level-\code{K_level} lambda is
#' \deqn{\lambda_level = \lambda_final * ((K_level-1)/(K_final-1))^{(2q-1)}.}
#'
#' @param lambda_final Numeric scalar; the lambda used on the final grid.
#' @param K_final Integer; number of points on the final grid.
#' @param K_level Integer; number of points on the current level grid.
#' @param q Integer; RW order (e.g., 2 for RW2).
#'
#' @return Numeric scalar \code{lambda_level}.
#'
#' @examples
#' lambda_scale_for_spacing(lambda_final = 500, K_final = 65, K_level = 17, q = 2)
#'
lambda_scale_for_spacing <- function(lambda_final, K_final, K_level, q = 2) {
  if (K_level < 2 || K_final < 2) stop("K_level and K_final must be >= 2")
  ratio <- (K_level - 1) / (K_final - 1)  # h_final/h_level
  lambda_final * (ratio^(2*q - 1))
}


#' Krige/interpolate from observed nodes using a Gaussian precision matrix
#'
#' @description
#' Compute the conditional mean under a Gaussian model specified via a precision
#' matrix \code{Q_new} on a grid of size \code{K_new}.
#'
#' Let the latent vector on the full grid be partitioned into observed indices \eqn{O}
#' and unobserved indices \eqn{U}. Under a (possibly improper) Gaussian prior with
#' precision \eqn{Q}, the conditional mean satisfies
#' \deqn{E[f_U | f_O] = - Q_{UU}^{-1} Q_{UO} f_O.}
#'
#' This routine supports vector-valued observations by allowing \code{f_obs} to be a
#' matrix with rows corresponding to \eqn{O} and columns corresponding to coordinates.
#' The solve is then done once with a matrix right-hand side.
#'
#' @param f_obs Numeric vector of length \code{|O|} OR numeric matrix \code{|O| x p}.
#'   Rows correspond to observed indices \code{idx_obs}. For multivariate means, \code{p = d}.
#' @param idx_obs Integer vector of observed indices (subset of \code{1:K_new}).
#' @param Q_new \code{K_new x K_new} precision matrix (dense or sparse \code{Matrix}).
#' @param nugget Nonnegative scalar. If \code{>0}, adds \code{nugget * I} to \eqn{Q_{UU}}
#'   to stabilize the solve (useful if \eqn{Q} is rank-deficient).
#' @param keep_obs Logical. If \code{TRUE}, return a full vector/matrix of length \code{K_new}
#'   with observed entries filled by \code{f_obs} and unobserved entries filled by the
#'   conditional mean. If \code{FALSE}, only return predictions on \eqn{U}.
#'
#' @return If \code{keep_obs=TRUE}:
#' \itemize{
#'   \item a numeric vector (if \code{f_obs} was a vector), length \code{K_new}; or
#'   \item a numeric matrix \code{K_new x p} (if \code{f_obs} was a matrix).
#' }
#' If \code{keep_obs=FALSE}, returns a list with:
#' \itemize{
#'   \item \code{pred_idx}: indices \eqn{U}
#'   \item \code{f_pred}: predicted values on \eqn{U} (vector or matrix)
#' }
#'
#' @examples
#' # (Toy) RW2 precision on K=65 grid, observe two endpoints and krige the rest.
#' # See your Rmd for a full worked example.
#'
kriging_from_precision <- function(f_obs, idx_obs, Q_new,
                                   nugget = 0,
                                   keep_obs = TRUE) {
  K_new <- nrow(Q_new)
  idx_obs  <- as.integer(idx_obs)
  idx_pred <- setdiff(seq_len(K_new), idx_obs)

  if (length(idx_pred) == 0L) {
    return(if (keep_obs) f_obs else list(pred_idx = integer(0), f_pred = NULL))
  }

  if (is.null(dim(f_obs))) f_obs_mat <- matrix(f_obs, ncol = 1) else f_obs_mat <- as.matrix(f_obs)
  if (nrow(f_obs_mat) != length(idx_obs)) stop("nrow(f_obs) must equal length(idx_obs).")

  Q_UU <- Q_new[idx_pred, idx_pred, drop = FALSE]
  Q_UO <- Q_new[idx_pred, idx_obs,  drop = FALSE]

  if (nugget > 0) Q_UU <- Q_UU + nugget * Diagonal(n = nrow(Q_UU))

  rhs <- Q_UO %*% f_obs_mat
  f_pred_mat <- - Matrix::solve(Q_UU, rhs)

  if (!keep_obs) {
    out_pred <- if (ncol(f_pred_mat) == 1) as.numeric(f_pred_mat) else f_pred_mat
    return(list(pred_idx = idx_pred, f_pred = out_pred))
  }

  if (ncol(f_obs_mat) == 1) {
    f_full <- numeric(K_new)
    f_full[idx_obs]  <- f_obs_mat[, 1]
    f_full[idx_pred] <- as.numeric(f_pred_mat)
    return(f_full)
  } else {
    p <- ncol(f_obs_mat)
    f_full <- matrix(NA_real_, nrow = K_new, ncol = p)
    f_full[idx_obs, ]  <- f_obs_mat
    f_full[idx_pred, ] <- f_pred_mat
    return(f_full)
  }
}


#' Krige SmoothEM mean functions from an observed grid to the final grid
#'
#' @description
#' Convenience wrapper for SmoothEM: take a list of component means defined on an observed
#' grid (design points) and krige/interpolate them to the full final grid using a 1D
#' precision matrix \code{Q_final_1d}.
#'
#' This is used in progressive initialization, where you fit SmoothEM at a coarse set of
#' positions \code{u_obs}, then krige to the final grid \code{u_final} to obtain
#' \code{Mu_full} (a \code{K_final x d} matrix) and \code{mu_full_list} (SmoothEM format).
#'
#' @param mu_list_obs List of length \code{K_obs}; each element is a numeric vector of length \code{d}.
#' @param u_obs Numeric vector of length \code{K_obs}; design-point locations.
#' @param u_final Numeric vector of length \code{K_final}; final grid locations.
#' @param Q_final_1d \code{K_final x K_final} precision matrix acting along the grid.
#' @param nugget Nonnegative scalar passed to \code{kriging_from_precision}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Mu_full}: numeric matrix \code{K_final x d}.
#'   \item \code{mu_full_list}: list length \code{K_final}, each a length-\code{d} vector.
#'   \item \code{idx_obs}: integer indices of \code{u_obs} inside \code{u_final}.
#' }
#'
krige_mu_list_to_full_grid <- function(mu_list_obs, u_obs, u_final, Q_final_1d,
                                       nugget = 0) {
  idx_obs <- match_locations_to_grid(u_obs, u_final)
  Mu_obs  <- do.call(rbind, mu_list_obs)  # (K_obs x d)
  Mu_full <- kriging_from_precision(
    f_obs   = Mu_obs,
    idx_obs = idx_obs,
    Q_new   = Q_final_1d,
    nugget  = nugget,
    keep_obs = TRUE
  )
  mu_full_list <- lapply(seq_len(nrow(Mu_full)), function(k) as.numeric(Mu_full[k, ]))
  list(Mu_full = Mu_full, mu_full_list = mu_full_list, idx_obs = idx_obs)
}



#' Progressive-resolution initialization for SmoothEM via kriging on a final grid
#'
#' @description
#' Fit SmoothEM on a sequence of *nested* (hierarchical) grids, starting from a very coarse
#' grid (stage 1 uses \code{K=2} endpoints) and progressively refining the resolution
#' (\code{K_m = 2^m + 1} for \code{m = 1,...,m_max}). At each stage:
#' \enumerate{
#'   \item Fit a penalized GMM (SmoothEM) on the current grid locations \code{u_fit}.
#'   \item Krige/interpolate the fitted component means back to the *final* grid
#'         (size \code{K_final = 2^m_max + 1}) using a 1D random-walk precision matrix.
#'   \item Use the kriged means (restricted to the next stage’s design points) as a warm start
#'         for the next stage.
#' }
#'
#' Stage 1 (the first history entry) is special: it initializes the two endpoint means using
#' \code{mclust::Mclust(data, G = 2)} and then kriges those two means to the final grid.
#'
#' The routine is designed for your “progressive initialization” idea: it avoids PCA-based
#' initialization and instead uses a coarse-to-fine continuation strategy where the mean
#' functions evolve smoothly across grid refinements.
#'
#' @details
#' \strong{Nested grids.}
#' The function uses \code{make_hierarchical_levels(m_max)} to construct a final grid
#' \code{u_final} of length \code{K_final = 2^m_max + 1} on \eqn{[0,1]} and a list of index sets
#' \code{idx_levels} such that each coarser level is an exact subset of the final grid.
#' This guarantees that design points align exactly across stages (no floating mismatch).
#'
#' \strong{Two precisions.}
#' \itemize{
#'   \item \code{Q_final_1d}: a \code{K_final x K_final} 1D precision used only for kriging means
#'         to the final grid.
#'   \item \code{Q_next}: a \code{(K_next*d) x (K_next*d)} separable precision used in SmoothEM
#'         at the current stage with \code{d = ncol(data)} coordinates.
#' }
#'
#' \strong{Spacing-aware lambda scaling.}
#' When the grid spacing changes across stages, the discrete RW(\code{q}) penalty needs rescaling
#' to keep the implied continuous roughness penalty comparable across grids. This function uses
#' \code{lambda_scale_for_spacing(lambda_final, K_final, K_next, q)} to obtain \code{lambda_next}.
#'
#' \strong{Rank deficiency.}
#' For an RW(\code{q}) prior, each coordinate typically contributes \code{q} null-space dimensions.
#' With \code{d} coordinates and a separable prior, the overall rank deficiency is often taken as
#' \code{rank_deficiency = q * d}. This matches how you compute log-determinants / constants in
#' your penalized objectives.
#'
#' @param data Numeric matrix \code{n x d}. Each row is an observation (e.g., a cell),
#'   each column is a coordinate (e.g., loading dimension).
#' @param m_max Nonnegative integer. Finest grid uses \code{K_final = 2^m_max + 1}.
#' @param lambda_final Numeric scalar. Penalty strength on the final grid (used to build \code{Q_final_1d}
#'   and as the baseline for scaling \code{lambda_next} at each stage).
#' @param q Integer. Random-walk order (e.g., \code{q=2} for RW2).
#' @param ridge Nonnegative scalar ridge added inside \code{make_random_walk_precision()} if supported.
#' @param nugget_kriging Nonnegative scalar nugget added to \eqn{Q_{UU}} in kriging to stabilize the solve
#'   (useful if RW precisions are rank-deficient). Passed to \code{krige_mu_list_to_full_grid()}.
#' @param tol Convergence tolerance for \code{EM_algorithm()}.
#' @param max_iter Maximum iterations for \code{EM_algorithm()} at each stage.
#' @param relative_lambda Logical. Passed to \code{EM_algorithm()}. If TRUE, scales the penalty relative to
#'   the current variance estimate (your “relative lambda” option).
#' @param modelName Character. Passed to \code{EM_algorithm()} (e.g., \code{"EEI"}).
#' @param coords_show Integer vector of coordinates (columns of \code{data}) to visualize when \code{plot_each_stage=TRUE}.
#' @param plot_each_stage Logical. If TRUE, plot kriged mean curves (on the final grid) at every stage.
#' @param verbose Logical. If TRUE, print stage summaries and pass verbose to \code{EM_algorithm()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{grid}: output of \code{make_hierarchical_levels()} (final grid + nested indices).
#'   \item \code{Q_final_1d}: final-grid 1D precision used for kriging.
#'   \item \code{fits}: named list of SmoothEM fits at each stage, keyed by \code{"K_<K_next>"}.
#'   \item \code{mu_full_history}: list of length \code{m_max+1}; each element is a \code{K_final x d} matrix
#'         of kriged means on the final grid for that stage (stage 1 corresponds to mclust initialization).
#'   \item \code{mu_full_list_final}: SmoothEM-style list (length \code{K_final}) of length-\code{d} mean vectors
#'         corresponding to the last kriged final-grid means.
#'   \item \code{meta_history}: list of per-stage metadata (stage index, \code{u_obs}, \code{K_obs}, \code{lambda}, etc.).
#' }
#'
#' @examples
#' \dontrun{
#' res <- progressive_smoothEM(Loadings, m_max = 6, lambda_final = 500, q = 2,
#'                            plot_each_stage = TRUE, verbose = TRUE)
#' plot_mu_history(res$mu_full_history, res$grid$u_final, history_i = 3, res = res,
#'                 coords = c(1,2,3), add_points = TRUE)
#' compare_positions_by_history(res, h1 = 2, h2 = 3, type = "mean")
#' }
#'
progressive_smoothEM <- function(
    data,
    m_max = 6,
    lambda_final = 500,
    q = 2,
    ridge = 0,
    nugget_kriging = 0,
    tol = 1e-3,
    max_iter = 1000,
    relative_lambda = TRUE,
    modelName = "EEI",
    coords_show = c(1,2,3),
    plot_each_stage = TRUE,
    verbose = TRUE
) {
  d <- ncol(data)

  meta_history <- list()

  # 1) hierarchical grids
  grid <- make_hierarchical_levels(m_max = m_max)
  u_final <- grid$u_final
  K_final <- grid$K_final

  # 2) final-grid per-coordinate precision for kriging (1D)
  Q_final_1d <- make_random_walk_precision(K = K_final, d = 1,
                                           lambda = lambda_final, q = q, ridge = ridge)

  # 3) stage 0 init: K=2 using Mclust
  u_obs0 <- c(0, 1)
  fit0 <- mclust::Mclust(data, G = 2)
  mu_list_obs <- list(fit0$parameters$mean[,1], fit0$parameters$mean[,2])

  # krige K=2 means to full grid => mu_full_list
  kr0 <- krige_mu_list_to_full_grid(mu_list_obs, u_obs0, u_final, Q_final_1d, nugget = nugget_kriging)
  mu_full_list <- kr0$mu_full_list

  fits <- list()
  mu_full_history <- list()
  mu_full_history[[1]] <- kr0$Mu_full

  meta_history[[1]] <- list(
    stage = 1,
    u_obs = u_obs0,
    mu_obs_list = mu_list_obs,
    K_obs = length(u_obs0),
    lambda = NA_real_
  )

  if (plot_each_stage) {
    matplot(u_final, kr0$Mu_full[, coords_show, drop = FALSE], type = "l", lty = 1,
            xlab = "u", ylab = expression(mu[j](u)),
            main = "Stage 1 (K=2): kriged means on final grid")
    for (jj in coords_show) points(u_obs0, c(mu_list_obs[[1]][jj], mu_list_obs[[2]][jj]), pch = 16)
  }

  # 4) for m=1..m_max, fit at K_m = 2^m + 1 (nested by construction)
  for (m in 1:m_max) {
    idx_fit <- grid$idx_levels[[m + 1]]   # +1 because idx_levels includes m=0
    u_fit   <- u_final[idx_fit]
    K_next  <- length(u_fit)

    # spacing-aware lambda for this level
    lambda_next <- lambda_scale_for_spacing(lambda_final, K_final, K_next, q = q)

    Q_next <- make_random_walk_precision(K = K_next, d = d,
                                         lambda = lambda_next, q = q, ridge = ridge)

    # init params: default init, but overwrite mu with kriged mu_full at these nested locations
    init_params <- make_default_init(data, K = K_next)
    init_params$mu <- mu_full_list[idx_fit]  # nested => directly subset by indices

    # rank deficiency for RW separable prior
    rank_def <- q * d

    if (verbose) {
      cat(sprintf("\n[Stage %d] K=%d, lambda_next=%.3g (lambda_final=%.3g)\n",
                  m, K_next, lambda_next, lambda_final))
    }

    fit_m <- EM_algorithm(
      data = data,
      Q_prior = Q_next,
      init_params = init_params,
      rank_deficiency = rank_def,
      max_iter = max_iter,
      modelName = modelName,
      tol = tol,
      relative_lambda = relative_lambda,
      verbose = verbose
    )
    fits[[paste0("K_", K_next)]] <- fit_m

    # krige fitted means back to final grid for the next warm start
    mu_list_obs <- fit_m$params$mu
    kr_m <- krige_mu_list_to_full_grid(mu_list_obs, u_fit, u_final, Q_final_1d, nugget = nugget_kriging)
    mu_full_list <- kr_m$mu_full_list
    mu_full_history[[m + 1]] <- kr_m$Mu_full

    meta_history[[m + 1]] <- list(
      stage = (m + 1),
      u_obs = u_fit,
      mu_obs_list = mu_list_obs,
      K_obs = K_next,
      lambda = lambda_next
    )

    if (plot_each_stage) {
      matplot(u_final, kr_m$Mu_full[, coords_show, drop = FALSE], type = "l", lty = 1,
              xlab = "u", ylab = expression(mu[j](u)),
              main = sprintf("Stage %d (K=%d): kriged means on final grid", m, K_next))
      for (jj in coords_show) {
        points(u_fit, sapply(mu_list_obs, function(mu_k) mu_k[jj]), pch = 16)
      }
    }
  }

  list(
    grid = grid,
    Q_final_1d = Q_final_1d,
    fits = fits,
    mu_full_history = mu_full_history,  # each is (K_final x d)
    mu_full_list_final = mu_full_list,
    meta_history = meta_history
  )
}



#' Plot kriged mean curves from mu_full_history
#'
#' @description
#' Plot selected coordinates (columns) of a kriged mean matrix stored in
#' \code{mu_full_history[[history_i]]}. Optionally overlay the corresponding observed
#' design-point means.
#'
#' If \code{add_points=TRUE} and \code{u_obs}/\code{mu_obs_list} are not provided,
#' the function will attempt to read them from \code{res$meta_history[[history_i]]}
#' (expected fields: \code{u_obs}, \code{mu_obs_list}).
#'
#' @param mu_full_history List of matrices; each element is \code{K_final x d}.
#' @param u_final Numeric vector of length \code{K_final}.
#' @param history_i Integer; which history element to plot.
#' @param coords Integer vector of coordinates (columns) to show.
#' @param main Plot title; if \code{NULL}, a default is constructed (uses \code{res} if available).
#' @param xlab,ylab Axis labels.
#' @param lty Line type for curves.
#' @param add_points Logical; overlay observed means at design points.
#' @param u_obs Numeric vector of observed design locations (optional).
#' @param mu_obs_list List of observed means at design points (optional).
#' @param res Optional progressive-fit result object; used to auto-fill points and build title.
#' @param pch,cex Point style for overlays.
#' @param legend_loc,legend_prefix,bty Legend settings.
#'
#' @return Invisibly returns a list containing the plotted matrix and any points used.
#'
plot_mu_history <- function(mu_full_history,
                            u_final,
                            history_i = 1,
                            coords = c(1, 2, 3),
                            main = NULL,
                            xlab = "u (grid position)",
                            ylab = expression(mu[j](u)),
                            lty = 1,
                            add_points = FALSE,
                            u_obs = NULL,
                            mu_obs_list = NULL,
                            res = NULL,                 # <--- new
                            pch = 16,
                            cex = 0.8,
                            legend_loc = "topright",
                            legend_prefix = "coord ",
                            bty = "n") {

  if (history_i < 1 || history_i > length(mu_full_history)) {
    stop("history_i out of range: must be in 1..length(mu_full_history).")
  }
  Mu <- mu_full_history[[history_i]]
  if (!is.matrix(Mu)) stop("mu_full_history[[history_i]] must be a matrix (K_final x d).")

  K_final <- nrow(Mu)
  d <- ncol(Mu)
  if (length(u_final) != K_final) stop("length(u_final) must equal nrow(Mu).")

  coords <- as.integer(coords)
  if (any(coords < 1 | coords > d)) stop("coords out of range.")

  # if add_points=TRUE and not supplied, try to pull from res$meta_history
  if (add_points && (is.null(u_obs) || is.null(mu_obs_list))) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        history_i <= length(res$meta_history) &&
        !is.null(res$meta_history[[history_i]]$u_obs) &&
        !is.null(res$meta_history[[history_i]]$mu_obs_list)) {

      u_obs <- res$meta_history[[history_i]]$u_obs
      mu_obs_list <- res$meta_history[[history_i]]$mu_obs_list
    } else {
      warning("add_points=TRUE but u_obs/mu_obs_list not provided and res$meta_history missing; skipping points.")
      add_points <- FALSE
    }
  }

  # default title: show stage/K/lambda if res provided
  if (is.null(main)) {
    if (!is.null(res) && !is.null(res$meta_history) && history_i <= length(res$meta_history)) {
      mh <- res$meta_history[[history_i]]
      if (!is.null(mh$K_obs)) {
        lam_txt <- if (!is.null(mh$lambda) && is.finite(mh$lambda)) sprintf(", lambda=%.3g", mh$lambda) else ""
        main <- sprintf("history %d (K=%d%s)", history_i, mh$K_obs, lam_txt)
      } else {
        main <- sprintf("mu_full_history[[%d]] (selected coordinates)", history_i)
      }
    } else {
      main <- sprintf("mu_full_history[[%d]] (selected coordinates)", history_i)
    }
  }

  matplot(u_final,
          Mu[, coords, drop = FALSE],
          type = "l",
          lty = lty,
          xlab = xlab,
          ylab = ylab,
          main = main)

  if (add_points) {
    # mu_obs_list: list length K_obs, each length-d
    for (j in coords) {
      yj <- vapply(mu_obs_list, function(mu_k) mu_k[j], numeric(1))
      points(u_obs, yj, pch = pch, cex = cex)
    }
  }

  legend(legend_loc,
         legend = paste0(legend_prefix, coords),
         lty = rep(lty, length(coords)),
         bty = bty)

  invisible(list(Mu = Mu, coords = coords, history_i = history_i,
                 u_obs = u_obs, mu_obs_list = mu_obs_list))
}


#' Plot change in a single coordinate between two histories
#'
#' @description
#' Visualize how one coordinate \code{coord} of the kriged mean curve changes between
#' \code{history_i} and \code{history_j}. Supports:
#' \itemize{
#'   \item \code{"overlay"}: plot both curves on the same axes.
#'   \item \code{"diff"}: plot the difference curve \eqn{\Delta(u) = \mu_j(u) - \mu_i(u)}.
#'   \item \code{"both"}: produce both plots (in two plotting pages).
#' }
#'
#' Optionally highlight design points corresponding to each history (from \code{res$meta_history[[h]]$u_obs})
#' or a user-supplied set \code{highlight_u}.
#'
#' The overlay plot can auto-set \code{ylim} based on both curves to make them comparable.
#'
#' @param mu_full_history List of matrices; each element is \code{K_final x d}.
#' @param u_final Numeric vector of length \code{K_final}.
#' @param coord Integer; which coordinate (column) to plot.
#' @param history_i,history_j Integers; which two histories to compare.
#' @param type One of \code{"both"}, \code{"overlay"}, \code{"diff"}.
#' @param highlight Which design points to highlight: \code{"i"}, \code{"j"}, \code{"both"}, \code{"none"}.
#' @param highlight_u Optional numeric vector of locations to highlight (overrides \code{res}).
#' @param res Optional progressive-fit result object (for extracting \code{u_obs} per history).
#' @param main Plot title; if \code{NULL} uses \code{res$meta_history} when available.
#' @param xlab,ylab_overlay,ylab_diff Axis labels.
#' @param lty_i,lty_j,lty_diff Line types.
#' @param pch_i,pch_j,pch_diff Point styles for highlighted points.
#' @param cex_pts Point size for highlights.
#' @param add_zero_line Logical; add horizontal zero line in diff plot.
#' @param legend_loc,bty Legend settings.
#' @param ylim_overlay Optional y-limits for overlay plot. If \code{NULL}, computed automatically.
#' @param pad_ylim Nonnegative numeric; padding fraction used when auto-computing \code{ylim_overlay}.
#' @param include_highlight_in_ylim Logical; include highlighted points in y-range computation.
#' @param label_points Logical; label highlighted indices on the plot.
#'
#' @return Invisibly returns a list with extracted curves and highlighted indices.
#'
plot_coordinate_change <- function(mu_full_history,
                                   u_final,
                                   coord,
                                   history_i,
                                   history_j,
                                   type = c("both", "overlay", "diff"),
                                   highlight = c("both", "i", "j", "none"),
                                   highlight_u = NULL,
                                   res = NULL,
                                   main = NULL,
                                   xlab = "u (grid position)",
                                   ylab_overlay = expression(mu[j](u)),
                                   ylab_diff = expression(Delta*mu[j](u)),
                                   lty_i = 1,
                                   lty_j = 2,
                                   lty_diff = 1,
                                   pch_i = 16,
                                   pch_j = 17,
                                   pch_diff = 16,
                                   cex_pts = 0.9,
                                   add_zero_line = TRUE,
                                   legend_loc = "topright",
                                   bty = "n",
                                   ylim_overlay = NULL,
                                   pad_ylim = 0.04,
                                   include_highlight_in_ylim = TRUE,
                                   label_points = FALSE)
{

  type <- match.arg(type)
  highlight <- match.arg(highlight)

  H <- length(mu_full_history)
  if (history_i < 1 || history_i > H) stop("history_i out of range.")
  if (history_j < 1 || history_j > H) stop("history_j out of range.")
  if (history_i == history_j) stop("history_i and history_j must be different.")

  Mu_i <- mu_full_history[[history_i]]
  Mu_j <- mu_full_history[[history_j]]
  if (!is.matrix(Mu_i) || !is.matrix(Mu_j)) stop("mu_full_history entries must be matrices (K_final x d).")
  if (nrow(Mu_i) != nrow(Mu_j)) stop("Different K_final across histories.")
  if (length(u_final) != nrow(Mu_i)) stop("length(u_final) must equal nrow(Mu_i).")

  d <- ncol(Mu_i)
  if (coord < 1 || coord > d) stop("coord out of range.")

  yi <- Mu_i[, coord]
  yj <- Mu_j[, coord]
  dy <- yj - yi

  # --- decide highlight u's
  get_u_from_res <- function(h) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        h <= length(res$meta_history) &&
        !is.null(res$meta_history[[h]]$u_obs)) {
      return(res$meta_history[[h]]$u_obs)
    }
    NULL
  }

  u_hi_i <- NULL
  u_hi_j <- NULL

  if (!is.null(highlight_u)) {
    # user-specified: use same set for both
    u_hi_i <- highlight_u
    u_hi_j <- highlight_u
  } else if (highlight %in% c("i", "both")) {
    u_hi_i <- get_u_from_res(history_i)
  }
  if (is.null(highlight_u) && highlight %in% c("j", "both")) {
    u_hi_j <- get_u_from_res(history_j)
  }

  idx_hi_i <- if (!is.null(u_hi_i)) match_locations_to_grid(u_hi_i, u_final) else integer(0)
  idx_hi_j <- if (!is.null(u_hi_j)) match_locations_to_grid(u_hi_j, u_final) else integer(0)

  # --- default title
  if (is.null(main)) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        history_i <= length(res$meta_history) &&
        history_j <= length(res$meta_history)) {
      Ki <- res$meta_history[[history_i]]$K_obs
      Kj <- res$meta_history[[history_j]]$K_obs
      main <- sprintf("coord %d: history %d (K=%s) → history %d (K=%s)",
                      coord, history_i, ifelse(is.null(Ki), "?", Ki),
                      history_j, ifelse(is.null(Kj), "?", Kj))
    } else {
      main <- sprintf("coord %d: history %d → %d", coord, history_i, history_j)
    }
  }

  # --- overlay plot
  # --- overlay plot
  if (type %in% c("overlay", "both")) {

    # ----- compute y-limits for comparability
    if (is.null(ylim_overlay)) {
      y_all <- c(yi, yj)

      if (include_highlight_in_ylim) {
        if (length(idx_hi_i) > 0) y_all <- c(y_all, yi[idx_hi_i])
        if (length(idx_hi_j) > 0) y_all <- c(y_all, yj[idx_hi_j])
      }

      rng <- range(y_all, finite = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2])) stop("Non-finite values in curves; cannot set ylim.")
      if (rng[1] == rng[2]) {
        # flat curve: add a small buffer
        eps <- ifelse(rng[1] == 0, 1, abs(rng[1])) * 0.05
        rng <- rng + c(-eps, eps)
      } else {
        pad <- diff(rng) * pad_ylim
        rng <- rng + c(-pad, pad)
      }
      ylim_use <- rng
    } else {
      ylim_use <- ylim_overlay
    }

    plot(u_final, yi, type = "l", lty = lty_i,
         xlab = xlab, ylab = ylab_overlay,
         ylim = ylim_use,
         main = if (type == "overlay") main else paste0(main, " (overlay)"))
    lines(u_final, yj, lty = lty_j)

    # highlight points
    if (highlight %in% c("i", "both") && length(idx_hi_i) > 0) {
      points(u_final[idx_hi_i], yi[idx_hi_i], pch = pch_i, cex = cex_pts)
      if (label_points) text(u_final[idx_hi_i], yi[idx_hi_i], labels = idx_hi_i, pos = 3, cex = 0.7)
    }
    if (highlight %in% c("j", "both") && length(idx_hi_j) > 0) {
      points(u_final[idx_hi_j], yj[idx_hi_j], pch = pch_j, cex = cex_pts)
      if (label_points) text(u_final[idx_hi_j], yj[idx_hi_j], labels = idx_hi_j, pos = 3, cex = 0.7)
    }

    legend(legend_loc,
           legend = c(paste0("history ", history_i),
                      paste0("history ", history_j),
                      if (highlight %in% c("i","both") && length(idx_hi_i)>0) "design pts (i)" else NULL,
                      if (highlight %in% c("j","both") && length(idx_hi_j)>0) "design pts (j)" else NULL),
           lty = c(lty_i, lty_j, rep(NA, 2)),
           pch = c(NA, NA,
                   if (highlight %in% c("i","both") && length(idx_hi_i)>0) pch_i else NULL,
                   if (highlight %in% c("j","both") && length(idx_hi_j)>0) pch_j else NULL),
           bty = bty)
  }

  # --- difference plot
  if (type %in% c("diff", "both")) {
    if (type == "both") {
      # knitr-friendly approach: open a new plot page
      plot.new()
    }
    plot(u_final, dy, type = "l", lty = lty_diff,
         xlab = xlab, ylab = ylab_diff,
         main = if (type == "diff") main else paste0(main, " (difference)"))
    if (add_zero_line) abline(h = 0, lty = 2)

    # highlight same design points on dy
    idx_hi <- sort(unique(c(idx_hi_i, idx_hi_j)))
    if (length(idx_hi) > 0) {
      points(u_final[idx_hi], dy[idx_hi], pch = pch_diff, cex = cex_pts)
      if (label_points) text(u_final[idx_hi], dy[idx_hi], labels = idx_hi, pos = 3, cex = 0.7)
    }
  }

  invisible(list(coord = coord, history_i = history_i, history_j = history_j,
                 mu_i = yi, mu_j = yj, delta = dy,
                 idx_hi_i = idx_hi_i, idx_hi_j = idx_hi_j))
}


#' Extract a (u, gamma) block for a given progressive history index
#'
#' @description
#' Helper for progressive-resolution SmoothEM diagnostics.
#' For a given history index \code{h}, this function returns:
#' \itemize{
#'   \item \code{K}: number of design points for that history,
#'   \item \code{u}: the design-point locations,
#'   \item \code{gamma}: the responsibility matrix for that history.
#' }
#'
#' History 1 is special: it may be initialized from \code{mclust} (membership \code{z}),
#' so the function checks \code{res$mclust$z} (then several fallbacks).
#' For later histories, it expects \code{res$fits[[paste0("K_", K)]]$gamma} unless a
#' \code{res$fits_history[[h]]$gamma} override is provided.
#'
#' @param res Result object containing at least \code{meta_history}.
#' @param h Integer history index.
#' @param normalize_gamma Logical; if TRUE, row-normalize \code{gamma} to sum to 1.
#'
#' @return A list with \code{h}, \code{K}, \code{u}, \code{gamma}.
#'
get_history_block <- function(res,
                              h,
                              normalize_gamma = TRUE) {
  if (h < 1 || h > length(res$meta_history)) stop("h out of range.")

  mh <- res$meta_history[[h]]
  if (is.null(mh$K_obs) || is.null(mh$u_obs)) {
    stop("res$meta_history[[h]] must contain K_obs and u_obs.")
  }

  K <- mh$K_obs
  u <- mh$u_obs

  # ---- get gamma for this history
  gamma <- NULL

  if (h == 1) {
    # stage 0: mclust or stage0 fit
    if (!is.null(res$mclust) && !is.null(res$mclust$z)) {
      gamma <- res$mclust$z
    } else if (!is.null(res$fits$stage0) && !is.null(res$fits$stage0$gamma)) {
      gamma <- res$fits$stage0$gamma
    } else if (!is.null(res$fit_stage0) && !is.null(res$fit_stage0$z)) {
      gamma <- res$fit_stage0$z
    } else {
      stop("history 1: cannot find membership gamma. Expect res$mclust$z or res$fits$stage0$gamma.")
    }
  } else {
    key <- paste0("K_", K)
    if (!is.null(res$fits[[key]]) && !is.null(res$fits[[key]]$gamma)) {
      gamma <- res$fits[[key]]$gamma
    } else if (!is.null(res$fits_history) && !is.null(res$fits_history[[h]]) && !is.null(res$fits_history[[h]]$gamma)) {
      gamma <- res$fits_history[[h]]$gamma
    } else {
      stop(sprintf("history %d: cannot find gamma. Expected res$fits[['%s']]$gamma (or res$fits_history[[h]]$gamma).", h, key))
    }
  }

  # sanity
  if (ncol(gamma) != length(u)) {
    stop(sprintf("gamma ncol (%d) != length(u_obs) (%d) at history %d. (order mismatch?)",
                 ncol(gamma), length(u), h))
  }

  if (normalize_gamma) gamma <- gamma / rowSums(gamma)

  list(h = h, K = K, u = u, gamma = gamma)
}



#' Compare posterior position summaries across two histories
#'
#' @description
#' Compare how inferred positions change between two progressive histories.
#' For each observation \eqn{i} with responsibilities \eqn{\gamma_i} on grid \eqn{u}:
#' \itemize{
#'   \item \code{type="mean"} uses \eqn{E[u|x_i] = \sum_k \gamma_{ik} u_k}.
#'   \item \code{type="max"} uses \eqn{\arg\max_k \gamma_{ik}} mapped to \eqn{u_k}.
#' }
#'
#' Produces a scatter plot of history-2 positions vs history-1 positions and overlays:
#' \itemize{
#'   \item an identity line (optional) and
#'   \item a least-squares regression line.
#' }
#'
#' @param res Result object containing \code{meta_history} and fit memberships.
#' @param h1,h2 Integer history indices.
#' @param type One of \code{"mean"} or \code{"max"}.
#' @param add_identity Logical; add y=x reference line.
#' @param main Optional plot title.
#'
#' @return Invisibly returns a list with \code{pos1}, \code{pos2}, and their correlation.
#'
compare_positions_by_history <- function(res, h1, h2,
                                         type = c("mean", "max"),
                                         add_identity = TRUE,
                                         main = NULL) {
  type <- match.arg(type)

  b1 <- get_history_block(res, h1, normalize_gamma = TRUE)
  b2 <- get_history_block(res, h2, normalize_gamma = TRUE)

  if (type == "mean") {
    pos1 <- as.numeric(b1$gamma %*% b1$u)
    pos2 <- as.numeric(b2$gamma %*% b2$u)
  } else {
    pos1 <- b1$u[max.col(b1$gamma, ties.method = "first")]
    pos2 <- b2$u[max.col(b2$gamma, ties.method = "first")]
  }

  if (is.null(main)) {
    main <- sprintf("%s position: history %d (K=%d) vs history %d (K=%d)",
                    type, h1, b1$K, h2, b2$K)
  }

  plot(pos2 ~ pos1,
       xlab = sprintf("%s position (h=%d, K=%d)", type, h1, b1$K),
       ylab = sprintf("%s position (h=%d, K=%d)", type, h2, b2$K),
       main = main)

  if (add_identity) abline(0, 1, lty = 3)
  abline(lm(pos2 ~ pos1), lty = 2)

  invisible(list(pos1 = pos1, pos2 = pos2, cor = cor(pos1, pos2)))
}



