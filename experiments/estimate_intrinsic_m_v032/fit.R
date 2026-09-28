# Resumable candidate fitting. Source common.R before this file.
fit_score <- function(fit, M) {
  trace <- if (M == 1L) fit$fit$elbo_trace else fit$fit$objective_history
  as.numeric(tail(trace, 1L))
}
compact_fit <- function(fit, M, prior) {
  raw <- fit$fit
  stopifnot(identical(raw$control$convergence, "normalized"))
  stopifnot(identical(if (M == 1L) raw$control$tol else raw$control$tol_outer, design$tolerance))
  if (M == 1L) {
    W <- matrix(1, ncol(fit$data), 1L, dimnames = list(colnames(fit$data), "A"))
    omega <- 1
    positions <- matrix(fit$locations$mean$pseudotime, ncol = 1L,
      dimnames = list(rownames(fit$data), "A"))
    trace <- raw$elbo_trace
    temperature <- rep(1, length(trace))
  } else {
    W <- fit$partition$pi_weights
    omega <- as.numeric(MPCurver::fitted_prior(fit, type = "partition")$omega)
    positions <- vapply(fit$locations, function(x) x$mean$pseudotime, numeric(nrow(fit$data)))
    trace <- raw$objective_history
    temperature <- raw$temperature_history
    stopifnot(raw$control$partition_prior == prior, raw$control$partition_init == "similarity",
      length(omega) == M, tail(temperature, 1L) == 1)
    if (prior == "fixed") stopifnot(max(abs(omega - 1 / M)) < 1e-10)
    if (prior == "adaptive") stopifnot(max(abs(omega - colMeans(W))) < 1e-10)
  }
  stopifnot(fit$K == design$K, all(is.finite(positions)), all(is.finite(W)),
    all(W >= 0), max(abs(rowSums(W) - 1)) < 1e-10,
    all(is.finite(fit$params$sigma2)), all(fit$params$sigma2 > 0), all(is.finite(trace)))
  if (M == 1L) getFromNamespace(".mpcurve_dimension_single_score_check", "MPCurver")(fit, tail(trace, 1))
  t1 <- trace[temperature == 1]
  stopifnot(length(t1) > 1L, min(diff(t1)) > -1e-5)
  active <- if (prior == "adaptive") which(omega > design$effective_weight_tol) else seq_len(M)
  position_sd <- apply(positions[, active, drop = FALSE], 2, sd)
  correlations <- suppressWarnings(cor(positions[, active, drop = FALSE], method = "spearman"))
  off_diagonal <- abs(correlations[upper.tri(correlations)])
  list(M = M, effective_M = if (prior == "adaptive") length(active) else M,
    omega = omega, active = active, assignments = max.col(W, ties.method = "first"), W = W,
    positions = positions, score = tail(trace, 1), converged = isTRUE(raw$converged),
    iterations = raw$iter, objective = trace, temperature = temperature,
    convergence = raw$control$convergence, tolerance = design$tolerance,
    last_normalized_increment = abs(tail(diff(trace), 1)) / (design$n * design$D),
    min_T1_increment = min(diff(t1)),
    last_relative_increment = tail(diff(t1), 1) / (abs(tail(t1, 2)[1]) + 1e-12),
    threshold_counts = vapply(c(1e-12, 1e-9, 1e-6, 1e-3), function(tol) sum(omega > tol), integer(1)),
    constant_orderings = sum(position_sd < 1e-10),
    near_duplicate_orderings = any(off_diagonal > 0.995, na.rm = TRUE),
    init_info = raw$init_info, sigma2 = fit$params$sigma2)
}
run_candidate <- function(dataset, method, M) {
  key <- sprintf("%s_%s_M%02d", dataset$row$id, method, M)
  output_path <- file.path(study_dir, "candidates", paste0(key, ".rds"))
  checkpoint_path <- file.path(study_dir, "checkpoints", paste0(key, ".rds"))
  if (file.exists(output_path)) {
    saved <- readRDS(output_path)
    stopifnot(saved$design_hash == design_hash, saved$input_hash == dataset$input_hash)
    return(saved)
  }
  prior <- if (method == "adaptive") "adaptive" else "fixed"
  state <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else
    list(fit = NULL, elapsed_seconds = 0, warnings = character(), continuations = 0L,
         design_hash = design_hash, input_hash = dataset$input_hash)
  stopifnot(state$design_hash == design_hash, state$input_hash == dataset$input_hash)
  status <- "success"
  error <- ""
  fit_started <- proc.time()[["elapsed"]]
  block_running <- FALSE
  cat(sprintf("%s %s start; resumed=%s\n", Sys.time(), key, !is.null(state$fit)))
  flush.console()
  value <- tryCatch(withCallingHandlers({
    repeat {
      fit_started <- proc.time()[["elapsed"]]
      block_running <- TRUE
      if (is.null(state$fit)) {
        set.seed(dataset$row$seed + 1000000L)
        args <- c(list(X = dataset$X, intrinsic_dim = as.integer(M),
          partition_init = "similarity", partition_prior = prior,
          effective_weight_tol = design$effective_weight_tol),
          if (M > 1L) list(hard_assign_final = FALSE) else list(), fit_options())
        state$fit <- do.call(MPCurver::fit_mpcurve, args)
      } else {
        remaining <- design$maximum_sweeps - state$fit$fit$iter - 1L
        if (remaining <= 0L || isTRUE(state$fit$fit$converged)) {
          block_running <- FALSE
          break
        }
        state$fit <- MPCurver::do_mpcurve(state$fit,
          iter = min(design$continuation_budget, remaining),
          tol = design$tolerance, tol_outer = design$tolerance,
          convergence = design$convergence, verbose = FALSE)
        state$continuations <- state$continuations + 1L
      }
      state$elapsed_seconds <- state$elapsed_seconds + proc.time()[["elapsed"]] - fit_started
      block_running <- FALSE
      stopifnot(identical(state$fit$data, dataset$X))
      atomic_save(state, checkpoint_path)
      cat(sprintf("%s %s sweeps=%d score=%.6f converged=%s elapsed=%.2f\n", Sys.time(),
        key, state$fit$fit$iter, fit_score(state$fit, M), state$fit$fit$converged, state$elapsed_seconds))
      flush.console()
      if (isTRUE(state$fit$fit$converged) || state$fit$fit$iter >= design$maximum_sweeps - 1L) break
    }
    compact_fit(state$fit, M, prior)
  }, warning = function(w) {
    state$warnings <<- c(state$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    if (block_running) state$elapsed_seconds <<- state$elapsed_seconds + proc.time()[["elapsed"]] - fit_started
    NULL
  })
  if (is.null(value)) status <- "error" else if (!value$converged) status <- "nonconverged"
  saved <- list(key = key, dataset_id = dataset$row$id, method = method, M = M,
    status = status, error = error, warnings = unique(state$warnings), fit = value,
    elapsed_seconds = state$elapsed_seconds, continuations = state$continuations,
    design_hash = design_hash, input_hash = dataset$input_hash)
  atomic_save(saved, output_path)
  keep_full <- dataset$row$phase == "pilot" || dataset$row$replicate == 1L || status != "success"
  if (keep_full && !is.null(state$fit)) atomic_save(state,
    file.path(study_dir, "full_fits", paste0(key, ".rds")))
  if (status == "success" && file.exists(checkpoint_path)) unlink(checkpoint_path)
  cat(sprintf("%s %s finished: %s %s\n", Sys.time(), key, status, error))
  flush.console()
  saved
}
