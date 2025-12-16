# Append a new coordinate (from 1D M-step result) into an ordering's params (EEI)
append_coord_to_params <- function(params, one_res) {
  # params: list(mu=list length K, each length d_old),
  #         sigma=list length K of diag(d_old),
  #         pi=length-K
  # one_res: returned by MSTEP_coordinatewise_once (fields: mu_j (K), sigma2_j, pi)
  K <- length(one_res$mu_j)
  if (is.null(params) || length(params$mu) == 0) {
    mu_list <- lapply(seq_len(K), function(k) one_res$mu_j[k])
    Sigma_sh <- diag(one_res$sigma2_j, 1)
    sigma_list <- replicate(K, Sigma_sh, simplify = FALSE)
    pi_vec <- one_res$pi
  } else {
    d_old <- length(params$mu[[1]])
    # append μ_jk at the tail to keep local coord order consistent with J
    mu_list <- lapply(seq_len(K), function(k) c(params$mu[[k]], one_res$mu_j[k]))
    # append σ_j^2 onto the shared diagonal
    sigma_diag <- c(diag(params$sigma[[1]]), one_res$sigma2_j)
    Sigma_sh <- diag(sigma_diag, d_old + 1)
    sigma_list <- replicate(K, Sigma_sh, simplify = FALSE)
    pi_vec <- one_res$pi
  }
  list(mu = mu_list, sigma = sigma_list, pi = pi_vec)
}

# =============================================================================
# Verbosity normalizer (0/1/2)
# =============================================================================
.normalize_verbose <- function(v) {
  if (is.logical(v)) return(if (isTRUE(v)) 1L else 0L)
  v <- as.integer(v); v[is.na(v)] <- 0L
  pmin(2L, pmax(0L, v))
}

# =============================================================================
# =============================================================================
.greedy_em_refine_once <- function(
    Xsub, Gamma, params, Q_K,
    relative_lambda = TRUE,
    rw_order = 2,                # rank deficiency = RW order
    nugget = 0,
    tol_inner = 1e-6, max_inner = 20,
    greedy_em_tol = 1e-4, greedy_em_max_iter = 20,
    verbose = 0
) {
  v <- .normalize_verbose(verbose)
  for (it in seq_len(greedy_em_max_iter)) {
    mu_old  <- unlist(params$mu)
    sig_old <- diag(params$sigma[[1]])
    G_old   <- Gamma

    params <- MSTEP_coordinatewise(
      data = Xsub, gamma = Gamma, params = params,
      Q_K = Q_K, relative_lambda = relative_lambda,
      iterate_once = FALSE, nugget = nugget,
      rank_deficiency = rw_order,
      tol_inner = tol_inner, max_inner = max_inner,
      verbose = (v >= 2)
    )
    Gamma <- ESTEP(Xsub, init_cov_cache_fast(params))
    Nk <- pmax(colSums(Gamma), 1e-8); params$pi <- Nk / sum(Nk)

    mu_new  <- unlist(params$mu)
    sig_new <- diag(params$sigma[[1]])
    d_mu  <- sqrt(sum((mu_new - mu_old)^2)) / (sqrt(sum(mu_old^2)) + 1e-12)
    d_sig <- sqrt(sum((sig_new - sig_old)^2)) / (sqrt(sum(sig_old^2)) + 1e-12)
    d_G   <- sqrt(sum((Gamma - G_old)^2))   / (sqrt(sum(G_old^2)) + 1e-12)

    if (v >= 2) cat(sprintf("    [greedy-EM] Δμ=%.2e ΔΣ=%.2e ΔΓ=%.2e\n", d_mu, d_sig, d_G))
    if (max(d_mu, d_sig, d_G) < greedy_em_tol) break
  }
  list(Gamma = Gamma, params = params)
}

# =============================================================================
# Init -> weighting params for a d-dim EEI object, only coord j_local is set
# =============================================================================
.make_weight_params_from_init_1d <- function(init1d, d, j_local) {
  K <- length(init1d$pi)
  mu1d     <- vapply(seq_len(K), function(k) init1d$mu[[k]][1], numeric(1))
  sigma2_1 <- as.numeric(init1d$sigma[[1]][1,1])
  sigma2_1 <- pmax(sigma2_1, .Machine$double.eps)

  mu_list <- lapply(seq_len(K), function(k) { v <- numeric(d); v[j_local] <- mu1d[k]; v })
  sig_diag <- rep(1.0, d); sig_diag[j_local] <- sigma2_1
  Sigma    <- diag(sig_diag, d)
  sigma_list <- replicate(K, Sigma, simplify = FALSE)
  list(mu = mu_list, sigma = sigma_list, pi = init1d$pi)
}

# =============================================================================
# Per-coordinate score using init-based weighting (NO Gamma leakage)
# =============================================================================
score_one_coord <- function(X, j, Gamma, Q_K,
                            get_init_fun, init_method, K,
                            relative_lambda = TRUE,
                            nugget = 0, rank_deficiency = 2, eigen_tol = 0,
                            include_log_pi = TRUE, include_entropy = TRUE) {
  d <- ncol(X)

  # 1) weights from init on x[, j]
  init1d    <- get_init_fun(X[, j, drop = FALSE], K)
  params_wt <- .make_weight_params_from_init_1d(init1d, d = d, j_local = j)

  # 2) 1-D penalized M-step under (Gamma, Q_K)
  one <- MSTEP_coordinatewise_once(
    data = X, gamma = Gamma, params = params_wt, coord = j,
    Q_K = Q_K, relative_lambda = relative_lambda,
    iterate_once = TRUE, nugget = nugget, rank_deficiency = rank_deficiency
  )

  # 3) univariate penalized ELBO for coord j (exclude global terms by default)
  params_j <- list(
    mu    = lapply(seq_len(K), function(k){ v <- numeric(d); v[j] <- one$mu_j[k]; v }),
    sigma = replicate(K, { S <- diag(1, d); S[j, j] <- one$sigma2_j; S }, simplify = FALSE),
    pi    = one$pi
  )
  Qj <- if (isTRUE(relative_lambda)) Q_K / one$sigma2_j else Q_K

  elbo_j <- compute_penalized_ELBO_univariate(
    X, coord = j, Gamma = Gamma, params = params_j,
    Q_full_prior = NULL, Q_K = Qj,
    include_log_pi = include_log_pi, include_entropy = include_entropy,
    eigen_tol = eigen_tol, rank_deficiency = rank_deficiency
  )
  list(elbo = elbo_j, one = one)
}

# =============================================================================
# Assemble params in a subspace by 1D updates (weights from init per coord)
# =============================================================================
assemble_params_from_1d <- function(X, J, Gamma, Q_K,
                                    get_init_fun, init_method, K,
                                    relative_lambda = TRUE,
                                    nugget = 0, rank_deficiency = 2,
                                    iterate_once = TRUE,
                                    parallel = FALSE, ncores = 2) {
  Xg <- X[, J, drop = FALSE]
  dG <- ncol(Xg); if (dG == 0) return(NULL)

  solve1 <- function(jg) {
    init1d    <- get_init_fun(Xg[, jg, drop = FALSE], K)
    params_wt <- .make_weight_params_from_init_1d(init1d, d = dG, j_local = jg)
    MSTEP_coordinatewise_once(
      data = Xg, gamma = Gamma, params = params_wt, coord = jg,
      Q_K = Q_K, relative_lambda = relative_lambda,
      iterate_once = iterate_once, nugget = nugget, rank_deficiency = rank_deficiency
    )
  }
  res_list <- if (parallel) parallel::mclapply(seq_len(dG), solve1, mc.cores = ncores)
  else          lapply(seq_len(dG), solve1)

  mu_mat     <- do.call(rbind, lapply(res_list, `[[`, "mu_j"))
  sigma_diag <- vapply(res_list, `[[`, numeric(1), "sigma2_j")
  mu_list    <- lapply(seq_len(K), function(k) mu_mat[, k])
  Sigma_sh   <- diag(sigma_diag, dG)
  sigma_list <- replicate(K, Sigma_sh, simplify = FALSE)
  Nk <- pmax(colSums(Gamma), 1e-8)
  list(mu = mu_list, sigma = sigma_list, pi = Nk / sum(Nk))
}

# =============================================================================
# Reassign all coordinates by ELBO (init-based scoring)
# =============================================================================
reassign_all_coords_by_elbo <- function(
    X, Gamma1, Gamma2, Q1, Q2,
    get_init_fun, init_method, K,
    relative_lambda = TRUE,
    nugget = 0, rank_deficiency = c(2,2), eigen_tol = 0,
    parallel = FALSE, ncores = 2
) {
  d <- ncol(X)
  score_one <- function(j) {
    s1 <- score_one_coord(X, j, Gamma1, Q1, get_init_fun, init_method, K,
                          relative_lambda, nugget, rank_deficiency[1], eigen_tol)$elbo
    s2 <- score_one_coord(X, j, Gamma2, Q2, get_init_fun, init_method, K,
                          relative_lambda, nugget, rank_deficiency[2], eigen_tol)$elbo
    c(elbo1 = s1, elbo2 = s2)
  }
  out <- if (parallel) parallel::mclapply(seq_len(d), score_one, mc.cores = ncores)
  else          lapply(seq_len(d), score_one)
  M <- do.call(rbind, out)
  new_assign <- ifelse(M[,1] >= M[,2], 1L, 2L)
  list(new_assign = new_assign, elbo_mat = M)
}

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

# =============================================================================
# Main: two_ordering_smoothEM  (now with init-based weighting everywhere)
# =============================================================================
two_ordering_smoothEM <- function(
    X,
    K = 5,
    greedy_start_index = NULL,
    greedy_start_index2 = NULL,
    make_Q_prior,
    lambda_vec = c(10, 10),
    rw_order_vec = c(2, 2),            # rank deficiency = RW order
    relative_lambda = TRUE,
    backfitting = FALSE,
    init_method = c("random","pca","tsne","warm_start"),
    max_iter_backfitting = 20,
    max_iter_EM_inner = 100,
    max_iter_EM_final = 100,
    # greedy-time EM refinement
    greedy_em_refine = TRUE,
    greedy_em_max_iter = 100,
    greedy_em_tol = 1e-4,
    # misc
    nugget = 0,
    tol_inner = 1e-6,
    max_inner = 1,
    eigen_tol = 0,
    parallel = FALSE,
    ncores = 2,
    verbose = 1              # 0: silent, 1: outer steps, 2: + internal EM
) {

  init_method <- match.arg(init_method)
  v <- .normalize_verbose(verbose)
  n <- nrow(X); d <- ncol(X)

  # Helper: init generator (as you wrote)
  get_init <- function(Xj, K) {
    if (init_method == "random") {
      make_default_init(Xj, K)
    } else if (init_method == "pca") {
      pc <- prcomp(Xj, center = TRUE, scale. = TRUE)$x[, 1]
      make_init(Xj, ordering_vec = pc, K = K)
    } else if (init_method == "tsne") {
      tsne <- Rtsne::Rtsne(
        Xj, dims = 1, perplexity = 50,
        verbose = FALSE, check_duplicates = FALSE,
        theta = 0.5
      )$Y[, 1]
      make_init(Xj, ordering_vec = tsne, K = K)
    } else {
      make_default_init(Xj, K)
    }
  }

  # Priors (K×K)
  Q1 <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1])
  Q2 <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2])

  # ----- Greedy seeds -----
  j1 <- if (is.null(greedy_start_index)) which.max(apply(X, 2, var)) else greedy_start_index
  init1_seed <- get_init(X[, j1, drop = FALSE], K)
  Gamma1 <- ESTEP(X[, j1, drop = FALSE], init_cov_cache_fast(init1_seed))
  ref <- .greedy_em_refine_once(
    Xsub = X[, j1, drop = FALSE], Gamma = Gamma1, params = init1_seed, Q_K = Q1,
    relative_lambda = relative_lambda, rw_order = rw_order_vec[1],
    nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
    greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
    verbose = v
  ); Gamma1 <- ref$Gamma; params1 <- ref$params

  if(!is.null(greedy_start_index2)){
    j2 <- greedy_start_index2
  }else{
    # pick j2 as min univariate ELBO under (Gamma1, Q1) with init-based scoring
    elbos_under_1 <- vapply(setdiff(seq_len(d), j1), function(j) {
      score_one_coord(
        X, j, Gamma1, Q1,
        get_init_fun = get_init, init_method = init_method, K = K,
        relative_lambda = relative_lambda,
        nugget = nugget, rank_deficiency = rw_order_vec[1], eigen_tol = eigen_tol
      )$elbo
    }, numeric(1))
    j2 <- setdiff(seq_len(d), j1)[ which.min(elbos_under_1) ]
  }

  init2_seed <- get_init(X[, j2, drop = FALSE], K)
  Gamma2 <- ESTEP(X[, j2, drop = FALSE], init_cov_cache_fast(init2_seed))
  ref <- .greedy_em_refine_once(
    Xsub = X[, j2, drop = FALSE], Gamma = Gamma2, params = init2_seed, Q_K = Q2,
    relative_lambda = relative_lambda, rw_order = rw_order_vec[2],
    nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
    greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
    verbose = v
  ); Gamma2 <- ref$Gamma; params2 <- ref$params


  if (v >= 1) cat(sprintf("[Greedy] seeds: j1=%d, j2=%d\n", j1, j2))

  # Build params for seeds via 1-D penalized updates (init-based weights)
  J1 <- c(j1);  J2 <- c(j2)
  Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)
  Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)

  # Optional seed refine
  if (greedy_em_refine) {
    if (length(J1) > 0) {
      ref <- .greedy_em_refine_once(
        Xsub = X[, J1, drop = FALSE], Gamma = Gamma1, params = params1, Q_K = Q1,
        relative_lambda = relative_lambda, rw_order = rw_order_vec[1],
        nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
        greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
        verbose = v
      ); Gamma1 <- ref$Gamma; params1 <- ref$params
    }
    if (length(J2) > 0) {
      ref <- .greedy_em_refine_once(
        Xsub = X[, J2, drop = FALSE], Gamma = Gamma2, params = params2, Q_K = Q2,
        relative_lambda = relative_lambda, rw_order = rw_order_vec[2],
        nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
        greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
        verbose = v
      ); Gamma2 <- ref$Gamma; params2 <- ref$params
    }
  }

  # ----- Greedy assignment for remaining coords -----
  remaining <- setdiff(seq_len(d), c(j1, j2))
  for (t in seq_along(remaining)) {
    j <- remaining[t]

    s1 <- score_one_coord(
      X, j, Gamma1, Q1,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, rank_deficiency = rw_order_vec[1], eigen_tol = eigen_tol
    )
    s2 <- score_one_coord(
      X, j, Gamma2, Q2,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, rank_deficiency = rw_order_vec[2], eigen_tol = eigen_tol
    )

    if (s1$elbo >= s2$elbo) {
      params1 <- append_coord_to_params(params1, s1$one)
      J1 <- c(J1, j)
      X1 <- X[, J1, drop = FALSE]
      Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1)); Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)

      if (greedy_em_refine) {
        ref <- .greedy_em_refine_once(
          Xsub = X1, Gamma = Gamma1, params = params1, Q_K = Q1,
          relative_lambda = relative_lambda, rw_order = rw_order_vec[1],
          nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
          greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
          verbose = v
        ); Gamma1 <- ref$Gamma; params1 <- ref$params
      }

      if (v >= 1) cat(sprintf("[Greedy] coord %d -> ord1 (Δ=%.3f)\n", j, s1$elbo - s2$elbo))
    } else {
      params2 <- append_coord_to_params(params2, s2$one)
      J2 <- c(J2, j)
      X2 <- X[, J2, drop = FALSE]
      Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2)); Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)

      if (greedy_em_refine) {
        ref <- .greedy_em_refine_once(
          Xsub = X2, Gamma = Gamma2, params = params2, Q_K = Q2,
          relative_lambda = relative_lambda, rw_order = rw_order_vec[2],
          nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
          greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
          verbose = v
        ); Gamma2 <- ref$Gamma; params2 <- ref$params
      }

      if (v >= 1) cat(sprintf("[Greedy] coord %d -> ord2 (Δ=%.3f)\n", j, s2$elbo - s1$elbo))
    }
  }
  coord_assign <- integer(d); coord_assign[J1] <- 1L; coord_assign[J2] <- 2L
  if (v >= 1) cat(sprintf("[Greedy] done: |J1|=%d, |J2|=%d\n", length(J1), length(J2)))

  # ----- Greedy-only return -----
  if (!backfitting) {
    elbo1 <- if (length(J1) > 0) compute_ordering_elbo(
      Xsub = X[, J1, drop = FALSE], Gamma = Gamma1, params = params1, Q_K = Q1,
      relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[1]
    ) else NA_real_
    elbo2 <- if (length(J2) > 0) compute_ordering_elbo(
      Xsub = X[, J2, drop = FALSE], Gamma = Gamma2, params = params2, Q_K = Q2,
      relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[2]
    ) else NA_real_

    return(list(
      backfitting   = FALSE,
      coord_assign  = coord_assign,
      J1 = J1, J2 = J2,
      elbo1 = elbo1, elbo2 = elbo2, elbo_total = elbo1 + elbo2,
      Gamma1 = Gamma1, Gamma2 = Gamma2,
      params1 = params1, params2 = params2
    ))
  }

  # =========================
  # Backfitting with reassignment by ELBO
  # =========================
  history <- list()
  X1 <- X[, J1, drop = FALSE]; X2 <- X[, J2, drop = FALSE]

  for (bf in seq_len(max_iter_backfitting)) {
    if (length(J1) > 0) {
      for (r in seq_len(max_iter_EM_inner)) {
        params1 <- MSTEP_coordinatewise(
          data = X1, gamma = Gamma1, params = params1,
          Q_K = Q1, relative_lambda = relative_lambda,
          iterate_once = FALSE, nugget = nugget,
          rank_deficiency = rw_order_vec[1],
          tol_inner = tol_inner, max_inner = max_inner,
          verbose = (v >= 2)
        )
        Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1))
      }
    }
    if (length(J2) > 0) {
      for (r in seq_len(max_iter_EM_inner)) {
        params2 <- MSTEP_coordinatewise(
          data = X2, gamma = Gamma2, params = params2,
          Q_K = Q2, relative_lambda = relative_lambda,
          iterate_once = FALSE, nugget = nugget,
          rank_deficiency = rw_order_vec[2],
          tol_inner = tol_inner, max_inner = max_inner,
          verbose = (v >= 2)
        )
        Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2))
      }
    }

    # Reassign by ELBO (init-based)
    reas <- reassign_all_coords_by_elbo(
      X = X,
      Gamma1 = if (length(J1) > 0) Gamma1 else matrix(1/K, n, K),
      Gamma2 = if (length(J2) > 0) Gamma2 else matrix(1/K, n, K),
      Q1 = Q1, Q2 = Q2,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, eigen_tol = eigen_tol,
      rank_deficiency = rw_order_vec,
      parallel = parallel, ncores = ncores
    )

    new_assign <- reas$new_assign
    n_changes <- sum(new_assign != coord_assign)
    history[[length(history)+1L]] <- list(step = bf, changes = n_changes)
    if (v >= 1) cat(sprintf("[Backfitting %02d] reassigned %d coords\n", bf, n_changes))
    if (n_changes == 0) break

    # Refresh groups
    coord_assign <- new_assign
    J1 <- which(coord_assign == 1L); J2 <- which(coord_assign == 2L)
    X1 <- X[, J1, drop = FALSE];     X2 <- X[, J2, drop = FALSE]
    if (length(J1) > 0) {
      params1 <- get_init(X1, K)
      Gamma1  <- ESTEP(X1, init_cov_cache_fast(params1))
      Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)
    } else { Gamma1 <- NULL; params1 <- NULL }
    if (length(J2) > 0) {
      params2 <- get_init(X2, K)
      Gamma2  <- ESTEP(X2, init_cov_cache_fast(params2))
      Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)
    } else { Gamma2 <- NULL; params2 <- NULL }
  }

  # Final polish
  if (length(J1) > 0) {
    for (r in seq_len(max_iter_EM_final)) {
      params1 <- MSTEP_coordinatewise(
        data = X1, gamma = Gamma1, params = params1,
        Q_K = Q1, relative_lambda = relative_lambda,
        iterate_once = FALSE, nugget = nugget,
        rank_deficiency = rw_order_vec[1],
        tol_inner = tol_inner, max_inner = max_inner,
        verbose = (v >= 2)
      )
      Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1))
    }
  }
  if (length(J2) > 0) {
    for (r in seq_len(max_iter_EM_final)) {
      params2 <- MSTEP_coordinatewise(
        data = X2, gamma = Gamma2, params = params2,
        Q_K = Q2, relative_lambda = relative_lambda,
        iterate_once = FALSE, nugget = nugget,
        rank_deficiency = rw_order_vec[2],
        tol_inner = tol_inner, max_inner = max_inner,
        verbose = (v >= 2)
      )
      Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2))
    }
  }

  elbo1 <- if (length(J1) > 0) compute_ordering_elbo(
    Xsub = X[, J1, drop = FALSE], Gamma = Gamma1, params = params1, Q_K = Q1,
    relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[1]
  ) else NA_real_
  elbo2 <- if (length(J2) > 0) compute_ordering_elbo(
    Xsub = X[, J2, drop = FALSE], Gamma = Gamma2, params = params2, Q_K = Q2,
    relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[2]
  ) else NA_real_

  list(
    backfitting  = TRUE,
    coord_assign = coord_assign,
    J1 = J1, J2 = J2,
    elbo1 = elbo1, elbo2 = elbo2, elbo_total = elbo1 + elbo2,
    Gamma1 = Gamma1, Gamma2 = Gamma2,
    params1 = if (exists("params1")) params1 else NULL,
    params2 = if (exists("params2")) params2 else NULL,
    history = history
  )
}



# =============================================================================
# Main: two_ordering_smoothEM_v2 (max-gap greedy)
# =============================================================================

#' Two-ordering smooth EM (v2): max-gap greedy with global rescoring
#'
#' @description
#' Partition variables into two orderings and fit a K-component Gaussian mixture
#' with EEI covariance and a separable smoothing prior along components.
#' Compared to the standard greedy, this version is **more conservative**:
#' at each greedy iteration it computes univariate penalized ELBOs (excluding
#' log-mixture and entropy) for **all remaining** coordinates under both
#' orderings, picks the coordinate with the **largest absolute ELBO gap**, assigns
#' it to the higher-ELBO ordering, updates that ordering's responsibilities and
#' parameters (optionally short EM refine), and then **rescoring is repeated**
#' over the remaining coordinates. This improves stability on real data at the
#' cost of higher compute.
#'
#' @param X  (n × d) data matrix.
#' @param K  number of mixture components.
#' @param greedy_start_index,greedy_start_index2 optional seeds for the two orderings;
#'        if `greedy_start_index2` is NULL, it is chosen as the min-ELBO coord under (Γ₁, Q₁).
#' @param make_Q_prior function(K, d, lambda, q) returning the (K×K) prior block for d=1
#'        (assumed separable: Q_full = Q_K ⊗ I_d).
#' @param lambda_vec length-2 prior strengths for the two orderings.
#' @param rw_order_vec length-2 RW orders; **rank deficiency = RW order** per ordering.
#' @param relative_lambda logical; if TRUE, 1-D updates scale Q_K by σ_j^{-2}, and full ELBO
#'        uses Q_K ⊗ diag(σ_j^{-2}).
#' @param backfitting logical; if TRUE, runs EM + coordinate reassignment after greedy.
#' @param init_method character; one of c("random","pca","tsne","warm_start"=not implemented)
#'        for per-coordinate initializers used in scoring (prevents Γ leakage).
#' @param max_iter_backfitting,max_iter_EM_inner,max_iter_EM_final integer loop counts for EM.
#' @param greedy_em_refine logical; after each greedy addition, run short EM refinement.
#' @param greedy_em_max_iter,greedy_em_tol controls for the short EM refine.
#' @param nugget,tol_inner,max_inner,eigen_tol numeric controls passed to M-steps and logdet.
#' @param parallel,ncores kept for API symmetry; the v2 greedy loop itself rescans sequentially.
#' @param verbose 0=silent, 1=greedy/backfitting steps, 2=+internal EM progress.
#'
#' @return A list with:
#'   - `coord_assign` (length-d, values in {1,2}),
#'   - `J1`, `J2` (index sets),
#'   - `Gamma1`, `Gamma2` (n×K),
#'   - `params1`, `params2` (mu/sigma/pi),
#'   - `elbo1`, `elbo2`, `elbo_total`,
#'   - `backfitting` (logical),
#'   - `history` (reassignment counts per outer backfitting loop).
#'
two_ordering_smoothEM_v2 <- function(
    X,
    K = 5,
    greedy_start_index = NULL,
    greedy_start_index2 = NULL,
    make_Q_prior,
    lambda_vec = c(10, 10),
    rw_order_vec = c(2, 2),
    relative_lambda = TRUE,
    backfitting = FALSE,
    init_method = c("random","pca","tsne","warm_start"),
    max_iter_backfitting = 20,
    max_iter_EM_inner = 100,
    max_iter_EM_final = 100,
    # greedy-time EM refinement
    greedy_em_refine = TRUE,
    greedy_em_max_iter = 100,
    greedy_em_tol = 1e-4,
    # misc
    nugget = 0,
    tol_inner = 1e-6,
    max_inner = 1,
    eigen_tol = 0,
    parallel = FALSE,
    ncores = 2,
    verbose = 1
) {
  init_method <- match.arg(init_method)
  v <- .normalize_verbose(verbose)
  n <- nrow(X); d <- ncol(X)

  # Initializer for scoring (init-based weighting everywhere)
  get_init <- function(Xj, K) {
    if (init_method == "random") {
      make_default_init(Xj, K)
    } else if (init_method == "pca") {
      pc <- prcomp(Xj, center = TRUE, scale. = TRUE)$x[, 1]
      make_init(Xj, ordering_vec = pc, K = K)
    } else if (init_method == "tsne") {
      tsne <- Rtsne::Rtsne(
        Xj, dims = 1, perplexity = 50,
        verbose = FALSE, check_duplicates = FALSE,
        theta = 0.5
      )$Y[, 1]
      make_init(Xj, ordering_vec = tsne, K = K)
    } else {
      stop("Warm start initialization not implemented yet.")
    }
  }

  # Priors (K×K)
  Q1 <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[1], q = rw_order_vec[1])
  Q2 <- make_Q_prior(K = K, d = 1, lambda = lambda_vec[2], q = rw_order_vec[2])

  # ----- Seeds -----
  j1 <- if (is.null(greedy_start_index)) which.max(apply(X, 2, var)) else greedy_start_index
  init1_seed <- get_init(X[, j1, drop = FALSE], K)
  Gamma1 <- ESTEP(X[, j1, drop = FALSE], init_cov_cache_fast(init1_seed))
  ref <- .greedy_em_refine_once(
    Xsub = X[, j1, drop = FALSE], Gamma = Gamma1, params = init1_seed, Q_K = Q1,
    relative_lambda = relative_lambda, rw_order = rw_order_vec[1],
    nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
    greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
    verbose = v
  ); Gamma1 <- ref$Gamma; params1 <- ref$params
  Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)

  if (!is.null(greedy_start_index2) && greedy_start_index2 == j1) {
    warning("greedy_start_index2 equals greedy_start_index; selecting j2 automatically.")
    greedy_start_index2 <- NULL
  }
  if (!is.null(greedy_start_index2)) {
    j2 <- greedy_start_index2
  } else {
    elbos_under_1 <- vapply(setdiff(seq_len(d), j1), function(j) {
      score_one_coord(
        X, j, Gamma1, Q1,
        get_init_fun = get_init, init_method = init_method, K = K,
        relative_lambda = relative_lambda,
        nugget = nugget, rank_deficiency = rw_order_vec[1], eigen_tol = eigen_tol
      )$elbo
    }, numeric(1))
    j2 <- setdiff(seq_len(d), j1)[ which.min(elbos_under_1) ]
  }

  init2_seed <- get_init(X[, j2, drop = FALSE], K)
  Gamma2 <- ESTEP(X[, j2, drop = FALSE], init_cov_cache_fast(init2_seed))
  ref <- .greedy_em_refine_once(
    Xsub = X[, j2, drop = FALSE], Gamma = Gamma2, params = init2_seed, Q_K = Q2,
    relative_lambda = relative_lambda, rw_order = rw_order_vec[2],
    nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
    greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
    verbose = v
  ); Gamma2 <- ref$Gamma; params2 <- ref$params
  Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)

  if (v >= 1) cat(sprintf("[Greedy-v2] seeds: j1=%d, j2=%d\n", j1, j2))

  J1 <- c(j1); J2 <- c(j2)
  remaining <- setdiff(seq_len(d), c(j1, j2))

  # ----- v2 Greedy loop: global rescoring with ELBO vectors -----
  while (length(remaining) > 0) {
    idxs <- remaining

    elbo_vec1 <- vapply(idxs, function(j) {
      score_one_coord(
        X, j, Gamma1, Q1,
        get_init_fun = get_init, init_method = init_method, K = K,
        relative_lambda = relative_lambda,
        nugget = nugget, rank_deficiency = rw_order_vec[1], eigen_tol = eigen_tol
      )$elbo
    }, numeric(1))

    elbo_vec2 <- vapply(idxs, function(j) {
      score_one_coord(
        X, j, Gamma2, Q2,
        get_init_fun = get_init, init_method = init_method, K = K,
        relative_lambda = relative_lambda,
        nugget = nugget, rank_deficiency = rw_order_vec[2], eigen_tol = eigen_tol
      )$elbo
    }, numeric(1))

    gap_vec  <- elbo_vec1 - elbo_vec2
    gap_vec[!is.finite(gap_vec)] <- -Inf
    pick_pos <- which.max(abs(gap_vec))
    j_star   <- idxs[pick_pos]
    gap_star <- gap_vec[pick_pos]

    # Recompute for j_star to obtain the 1-D "one" object used to plug in
    s1_star <- score_one_coord(
      X, j_star, Gamma1, Q1,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, rank_deficiency = rw_order_vec[1], eigen_tol = eigen_tol
    )
    s2_star <- score_one_coord(
      X, j_star, Gamma2, Q2,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, rank_deficiency = rw_order_vec[2], eigen_tol = eigen_tol
    )

    if (gap_star >= 0) {
      # assign to ordering 1
      params1 <- append_coord_to_params(params1, s1_star$one)
      J1 <- c(J1, j_star)
      X1 <- X[, J1, drop = FALSE]
      Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1))
      Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)

      if (greedy_em_refine) {
        ref <- .greedy_em_refine_once(
          Xsub = X1, Gamma = Gamma1, params = params1, Q_K = Q1,
          relative_lambda = relative_lambda, rw_order = rw_order_vec[1],
          nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
          greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
          verbose = v
        ); Gamma1 <- ref$Gamma; params1 <- ref$params
      }
      if (v >= 1) cat(sprintf("[Greedy-v2] pick %d -> ord1 (|Δ|=%.3f), remaining=%d\n",
                              j_star, abs(gap_star), length(remaining)-1))
    } else {
      # assign to ordering 2
      params2 <- append_coord_to_params(params2, s2_star$one)
      J2 <- c(J2, j_star)
      X2 <- X[, J2, drop = FALSE]
      Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2))
      Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)

      if (greedy_em_refine) {
        ref <- .greedy_em_refine_once(
          Xsub = X2, Gamma = Gamma2, params = params2, Q_K = Q2,
          relative_lambda = relative_lambda, rw_order = rw_order_vec[2],
          nugget = nugget, tol_inner = tol_inner, max_inner = max_inner,
          greedy_em_tol = greedy_em_tol, greedy_em_max_iter = greedy_em_max_iter,
          verbose = v
        ); Gamma2 <- ref$Gamma; params2 <- ref$params
      }
      if (v >= 1) cat(sprintf("[Greedy-v2] pick %d -> ord2 (|Δ|=%.3f), remaining=%d\n",
                              j_star, abs(gap_star), length(remaining)-1))
    }

    remaining <- remaining[-pick_pos]  # global rescoring next round
  }

  coord_assign <- integer(d); coord_assign[J1] <- 1L; coord_assign[J2] <- 2L
  if (v >= 1) cat(sprintf("[Greedy-v2] done: |J1|=%d, |J2|=%d\n", length(J1), length(J2)))

  # ----- Greedy-only return -----
  if (!backfitting) {
    elbo1 <- if (length(J1) > 0) compute_ordering_elbo(
      Xsub = X[, J1, drop = FALSE], Gamma = Gamma1, params = params1, Q_K = Q1,
      relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[1]
    ) else NA_real_
    elbo2 <- if (length(J2) > 0) compute_ordering_elbo(
      Xsub = X[, J2, drop = FALSE], Gamma = Gamma2, params = params2, Q_K = Q2,
      relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[2]
    ) else NA_real_

    return(list(
      backfitting   = FALSE,
      coord_assign  = coord_assign,
      J1 = J1, J2 = J2,
      elbo1 = elbo1, elbo2 = elbo2, elbo_total = elbo1 + elbo2,
      Gamma1 = Gamma1, Gamma2 = Gamma2,
      params1 = params1, params2 = params2
    ))
  }

  # =========================
  # Backfitting (same as v1)
  # =========================
  history <- list()
  X1 <- X[, J1, drop = FALSE]; X2 <- X[, J2, drop = FALSE]

  for (bf in seq_len(max_iter_backfitting)) {
    if (length(J1) > 0) {
      for (r in seq_len(max_iter_EM_inner)) {
        params1 <- MSTEP_coordinatewise(
          data = X1, gamma = Gamma1, params = params1,
          Q_K = Q1, relative_lambda = relative_lambda,
          iterate_once = FALSE, nugget = nugget,
          rank_deficiency = rw_order_vec[1],
          tol_inner = tol_inner, max_inner = max_inner,
          verbose = (v >= 2)
        )
        Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1))
      }
    }
    if (length(J2) > 0) {
      for (r in seq_len(max_iter_EM_inner)) {
        params2 <- MSTEP_coordinatewise(
          data = X2, gamma = Gamma2, params = params2,
          Q_K = Q2, relative_lambda = relative_lambda,
          iterate_once = FALSE, nugget = nugget,
          rank_deficiency = rw_order_vec[2],
          tol_inner = tol_inner, max_inner = max_inner,
          verbose = (v >= 2)
        )
        Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2))
      }
    }

    reas <- reassign_all_coords_by_elbo(
      X = X,
      Gamma1 = if (length(J1) > 0) Gamma1 else matrix(1/K, n, K),
      Gamma2 = if (length(J2) > 0) Gamma2 else matrix(1/K, n, K),
      Q1 = Q1, Q2 = Q2,
      get_init_fun = get_init, init_method = init_method, K = K,
      relative_lambda = relative_lambda,
      nugget = nugget, eigen_tol = eigen_tol,
      rank_deficiency = rw_order_vec,
      parallel = FALSE, ncores = ncores
    )

    new_assign <- reas$new_assign
    n_changes <- sum(new_assign != coord_assign)
    history[[length(history)+1L]] <- list(step = bf, changes = n_changes)
    if (v >= 1) cat(sprintf("[Backfitting %02d] reassigned %d coords\n", bf, n_changes))
    if (n_changes == 0) break

    coord_assign <- new_assign
    J1 <- which(coord_assign == 1L); J2 <- which(coord_assign == 2L)
    X1 <- X[, J1, drop = FALSE];     X2 <- X[, J2, drop = FALSE]
    if (length(J1) > 0) {
      params1 <- get_init(X1, K)
      Gamma1  <- ESTEP(X1, init_cov_cache_fast(params1))
      Nk1 <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk1 / sum(Nk1)
    } else { Gamma1 <- NULL; params1 <- NULL }
    if (length(J2) > 0) {
      params2 <- get_init(X2, K)
      Gamma2  <- ESTEP(X2, init_cov_cache_fast(params2))
      Nk2 <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk2 / sum(Nk2)
    } else { Gamma2 <- NULL; params2 <- NULL }
  }

  # Final polish
  if (length(J1) > 0) {
    for (r in seq_len(max_iter_EM_final)) {
      params1 <- MSTEP_coordinatewise(
        data = X1, gamma = Gamma1, params = params1,
        Q_K = Q1, relative_lambda = relative_lambda,
        iterate_once = FALSE, nugget = nugget,
        rank_deficiency = rw_order_vec[1],
        tol_inner = tol_inner, max_inner = max_inner,
        verbose = (v >= 2)
      )
      Gamma1 <- ESTEP(X1, init_cov_cache_fast(params1))
    }
  }
  if (length(J2) > 0) {
    for (r in seq_len(max_iter_EM_final)) {
      params2 <- MSTEP_coordinatewise(
        data = X2, gamma = Gamma2, params = params2,
        Q_K = Q2, relative_lambda = relative_lambda,
        iterate_once = FALSE, nugget = nugget,
        rank_deficiency = rw_order_vec[2],
        tol_inner = tol_inner, max_inner = max_inner,
        verbose = (v >= 2)
      )
      Gamma2 <- ESTEP(X2, init_cov_cache_fast(params2))
    }
  }

  elbo1 <- if (length(J1) > 0) compute_ordering_elbo(
    Xsub = X[, J1, drop = FALSE], Gamma = Gamma1, params = params1, Q_K = Q1,
    relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[1]
  ) else NA_real_
  elbo2 <- if (length(J2) > 0) compute_ordering_elbo(
    Xsub = X[, J2, drop = FALSE], Gamma = Gamma2, params = params2, Q_K = Q2,
    relative_lambda = relative_lambda, eigen_tol = eigen_tol, rw_order = rw_order_vec[2]
  ) else NA_real_

  list(
    backfitting  = TRUE,
    coord_assign = coord_assign,
    J1 = J1, J2 = J2,
    elbo1 = elbo1, elbo2 = elbo2, elbo_total = elbo1 + elbo2,
    Gamma1 = Gamma1, Gamma2 = Gamma2,
    params1 = if (exists("params1")) params1 else NULL,
    params2 = if (exists("params2")) params2 else NULL,
    history = history
  )
}

