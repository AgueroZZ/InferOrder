# Resumable fixed-dimension fitting with the original study's settings.
# Source common.R before this file.

execution_provenance <- function() {
  scripts <- c("common.R", "fit.R", "evaluate.R", "run_dataset.R", "reporting_metrics.R")
  script_paths <- file.path(study_dir, scripts)
  archive <- file.path(study_dir, "source", "MPCurver_0.3.2.tar.gz")
  source_commit <- file.path(study_dir, "source", "source_commit.txt")
  stopifnot(all(file.exists(script_paths)), file.exists(archive), file.exists(source_commit))
  hashes <- setNames(vapply(script_paths, digest::digest, character(1),
    algo = "sha256", file = TRUE), scripts)
  record <- list(code_hashes = hashes,
    package_archive_hash = digest::digest(file = archive, algo = "sha256"),
    package_source_commit = paste(readLines(source_commit, warn = FALSE), collapse = "\n"),
    package_version = as.character(utils::packageVersion("MPCurver")),
    R_version = as.character(getRversion()))
  stopifnot(record$package_archive_hash == design$package_archive_sha256,
    record$package_source_commit == design$package_source_commit,
    record$package_version == design$package_version)
  record$fingerprint <- digest::digest(record, algo = "sha256")
  record
}

fit_score <- function(fit, M) {
  trace <- if (M == 1L) fit$fit$elbo_trace else fit$fit$objective_history
  as.numeric(tail(trace, 1L))
}

# These tolerances match the package's numerical nondecrease safeguards.
# A single-ordering fit requires a nonnegative last increment to converge.
objective_nondecrease_guard <- function(previous, M, stopping = FALSE) {
  if (M == 1L) {
    rep(if (stopping) 0 else 1e-8, length(previous))
  } else {
    1e-8 * (abs(previous) + 1)
  }
}

fit_trace_diagnostics <- function(fit, M) {
  if (is.null(fit)) return(NULL)
  raw <- fit$fit
  trace <- if (M == 1L) raw$elbo_trace else raw$objective_history
  temperature <- if (M == 1L) rep(1, length(trace)) else raw$temperature_history
  T1 <- trace[temperature == 1]
  increments <- diff(T1)
  previous <- head(T1, -1L)
  valid_trace <- length(trace) == length(temperature) && length(T1) > 1L &&
    all(is.finite(trace)) && all(is.finite(temperature))
  nondecrease <- valid_trace && all(increments >= -objective_nondecrease_guard(previous, M))
  final_increment <- if (length(increments)) as.numeric(tail(increments, 1L)) else NA_real_
  final_previous <- if (length(previous)) as.numeric(tail(previous, 1L)) else NA_real_
  normalized_increment <- abs(final_increment) / (design$n * design$D)
  stopping_passed <- valid_trace && tail(temperature, 1L) == 1 &&
    final_increment >= -objective_nondecrease_guard(final_previous, M, stopping = TRUE) &&
    normalized_increment < design$tolerance
  list(objective = as.numeric(trace), temperature = as.numeric(temperature),
    iterations = as.integer(raw$iter), converged = isTRUE(raw$converged),
    convergence = raw$control$convergence,
    tolerance = if (M == 1L) raw$control$tol else raw$control$tol_outer,
    final_temperature = if (length(temperature)) as.numeric(tail(temperature, 1L)) else NA_real_,
    score = if (length(trace)) as.numeric(tail(trace, 1L)) else NA_real_,
    valid_trace = valid_trace, T1_nondecrease_passed = nondecrease,
    stopping_passed = stopping_passed,
    last_normalized_increment = normalized_increment,
    min_T1_increment = if (length(increments)) min(increments) else NA_real_,
    last_relative_increment = final_increment / (abs(final_previous) + 1e-12))
}

compact_fit <- function(fit, M, prior) {
  raw <- fit$fit
  diagnostics <- fit_trace_diagnostics(fit, M)
  stopifnot(identical(raw$control$convergence, "normalized"),
    identical(diagnostics$tolerance, design$tolerance))
  if (M == 1L) {
    W <- matrix(1, ncol(fit$data), 1L, dimnames = list(colnames(fit$data), "A"))
    omega <- 1
    positions <- matrix(fit$locations$mean$pseudotime, ncol = 1L,
      dimnames = list(rownames(fit$data), "A"))
    getFromNamespace(".mpcurve_dimension_single_score_check", "MPCurver")(
      fit, diagnostics$score)
  } else {
    W <- fit$partition$pi_weights
    omega <- as.numeric(MPCurver::fitted_prior(fit, type = "partition")$omega)
    positions <- vapply(fit$locations, function(x) x$mean$pseudotime,
      numeric(nrow(fit$data)))
    stopifnot(raw$control$partition_prior == prior,
      raw$control$partition_init == "similarity", length(omega) == M)
    if (prior == "fixed") stopifnot(max(abs(omega - 1 / M)) < 1e-10)
    if (prior == "adaptive") stopifnot(max(abs(omega - colMeans(W))) < 1e-10)
  }
  stopifnot(fit$K == design$K, all(is.finite(positions)), all(is.finite(W)),
    all(W >= 0), max(abs(rowSums(W) - 1)) < 1e-10,
    all(is.finite(fit$params$sigma2)), all(fit$params$sigma2 > 0),
    diagnostics$valid_trace, diagnostics$final_temperature == 1,
    diagnostics$T1_nondecrease_passed)
  if (diagnostics$converged) stopifnot(diagnostics$stopping_passed)

  # Occupancy uses posterior assignment weights for both fitting strategies.
  occupancy <- colMeans(W)
  active <- which(occupancy > design$effective_weight_tol)
  assignments <- max.col(W, ties.method = "first")
  stopifnot(length(active) > 0L, all(assignments %in% active))
  position_sd <- apply(positions[, active, drop = FALSE], 2, stats::sd)
  correlations <- suppressWarnings(stats::cor(positions[, active, drop = FALSE],
    method = "spearman"))
  off_diagonal <- abs(correlations[upper.tri(correlations)])
  list(M = as.integer(M), selected_model_M = as.integer(M),
    effective_M = as.integer(length(active)), omega = omega, occupancy = occupancy,
    active = active, assignments = assignments, W = W, positions = positions,
    score = diagnostics$score, converged = diagnostics$converged,
    iterations = diagnostics$iterations, objective = diagnostics$objective,
    temperature = diagnostics$temperature, convergence = diagnostics$convergence,
    tolerance = diagnostics$tolerance,
    last_normalized_increment = diagnostics$last_normalized_increment,
    min_T1_increment = diagnostics$min_T1_increment,
    last_relative_increment = diagnostics$last_relative_increment,
    threshold_counts = vapply(c(1e-12, 1e-9, 1e-6, 1e-3),
      function(tol) as.integer(sum(occupancy > tol)), integer(1)),
    constant_orderings = as.integer(sum(position_sd < 1e-10)),
    near_duplicate_orderings = any(off_diagonal > 0.995, na.rm = TRUE),
    init_info = raw$init_info, sigma2 = fit$params$sigma2)
}

empty_block_history <- function() {
  data.frame(block = integer(), kind = character(), requested_budget = integer(),
    start_iteration = integer(), end_iteration = integer(), score = numeric(),
    converged = logical(), elapsed_seconds = numeric())
}

run_candidate <- function(dataset, method, M) {
  stopifnot(method %in% c("adaptive", "forward"), M >= 1L, M <= design$max_M,
    identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
  provenance <- execution_provenance()
  key <- sprintf("%s_%s_M%02d", dataset$row$id, method, M)
  output_path <- file.path(study_dir, "candidates", paste0(key, ".rds"))
  checkpoint_path <- file.path(study_dir, "checkpoints", paste0(key, ".rds"))
  if (file.exists(output_path)) {
    saved <- readRDS(output_path)
    stopifnot(saved$design_hash == design_hash, saved$input_hash == dataset$input_hash,
      identical(saved$provenance, provenance))
    return(saved)
  }
  prior <- if (method == "adaptive") "adaptive" else "fixed"
  state <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else
    list(fit = NULL, elapsed_seconds = 0, warnings = character(), continuations = 0L,
      block_history = empty_block_history(), rng_state = NULL,
      design_hash = design_hash, input_hash = dataset$input_hash, provenance = provenance)
  stopifnot(state$design_hash == design_hash, state$input_hash == dataset$input_hash,
    identical(state$provenance, provenance))
  if (!is.null(state$rng_state)) assign(".Random.seed", state$rng_state, envir = .GlobalEnv)
  status <- "success"
  error <- ""
  fit_started <- proc.time()[["elapsed"]]
  block_running <- FALSE
  cat(sprintf("%s %s start; resumed=%s\n", Sys.time(), key, !is.null(state$fit)))
  flush.console()
  value <- tryCatch(withCallingHandlers({
    repeat {
      fit_started <- proc.time()[["elapsed"]]
      start_iteration <- if (is.null(state$fit)) 0L else state$fit$fit$iter
      block_running <- TRUE
      if (is.null(state$fit)) {
        set.seed(dataset$row$seed + 1000000L)
        block_kind <- "initial"
        requested_budget <- design$initial_budget
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
        block_kind <- "continuation"
        requested_budget <- min(design$continuation_budget, remaining)
        state$fit <- MPCurver::do_mpcurve(state$fit, iter = requested_budget,
          tol = design$tolerance, tol_outer = design$tolerance,
          convergence = design$convergence, verbose = FALSE)
        state$continuations <- state$continuations + 1L
      }
      elapsed <- proc.time()[["elapsed"]] - fit_started
      state$elapsed_seconds <- state$elapsed_seconds + elapsed
      block_running <- FALSE
      state$block_history <- rbind(state$block_history,
        data.frame(block = nrow(state$block_history) + 1L, kind = block_kind,
          requested_budget = as.integer(requested_budget),
          start_iteration = as.integer(start_iteration),
          end_iteration = as.integer(state$fit$fit$iter), score = fit_score(state$fit, M),
          converged = isTRUE(state$fit$fit$converged), elapsed_seconds = elapsed))
      state$rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      stopifnot(identical(state$fit$data, dataset$X))
      atomic_save(state, checkpoint_path)
      cat(sprintf("%s %s sweeps=%d score=%.6f converged=%s elapsed=%.2f\n", Sys.time(),
        key, state$fit$fit$iter, fit_score(state$fit, M), state$fit$fit$converged,
        state$elapsed_seconds))
      flush.console()
      if (isTRUE(state$fit$fit$converged) ||
          state$fit$fit$iter >= design$maximum_sweeps - 1L) break
    }
    compact_fit(state$fit, M, prior)
  }, warning = function(w) {
    state$warnings <<- c(state$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    if (block_running) {
      state$elapsed_seconds <<- state$elapsed_seconds + proc.time()[["elapsed"]] - fit_started
    }
    NULL
  })
  if (is.null(value)) status <- "error" else if (!value$converged) status <- "nonconverged"
  saved <- list(key = key, dataset_id = dataset$row$id, method = method, M = as.integer(M),
    status = status, error = error, warnings = unique(state$warnings), fit = value,
    diagnostics = fit_trace_diagnostics(state$fit, M),
    elapsed_seconds = state$elapsed_seconds, continuations = state$continuations,
    block_history = state$block_history, design_hash = design_hash,
    input_hash = dataset$input_hash, provenance = provenance)
  atomic_save(saved, output_path)
  keep_full <- dataset$row$replicate == 1L || status != "success"
  if (keep_full && !is.null(state$fit)) {
    atomic_save(state, file.path(study_dir, "full_fits", paste0(key, ".rds")))
  }
  if (status == "success" && file.exists(checkpoint_path)) unlink(checkpoint_path)
  cat(sprintf("%s %s finished: %s %s\n", Sys.time(), key, status, error))
  flush.console()
  saved
}
