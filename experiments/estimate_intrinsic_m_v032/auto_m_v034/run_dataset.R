#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/auto_m_v034/common.R")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript run_dataset.R <one-based manifest index or dataset id>")
}
if (grepl("^[0-9]+$", arguments[1L])) {
  index <- as.integer(arguments[1L])
  stopifnot(!is.na(index), index >= 1L, index <= nrow(manifest))
  row <- manifest[index, , drop = FALSE]
} else {
  row <- manifest[manifest$id == arguments[1L], , drop = FALSE]
}
stopifnot(nrow(row) == 1L)
dataset <- load_fixed_dataset(row)
provenance <- execution_provenance()
output_path <- file.path(extension_dir, "results", paste0(row$id, "_auto_adaptive.rds"))
checkpoint_path <- file.path(extension_dir, "checkpoints", paste0(row$id, ".rds"))

if (file.exists(output_path)) {
  saved <- readRDS(output_path)
  stopifnot(saved$design_hash == design_hash,
    saved$input_hash == dataset$input_hash,
    identical(saved$provenance, provenance))
  cat("RESULT ALREADY COMPLETE\n")
  print(saved$row, row.names = FALSE)
  quit(status = 0L)
}

empty_block_history <- function() {
  data.frame(block = integer(), kind = character(), requested_budget = integer(),
    start_iteration = integer(), end_iteration = integer(), score = numeric(),
    converged = logical(), elapsed_seconds = numeric())
}

state <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else list(
  fit = NULL,
  elapsed_seconds = 0,
  warnings = character(),
  continuations = 0L,
  block_history = empty_block_history(),
  rng_state = NULL,
  design_hash = design_hash,
  input_hash = dataset$input_hash,
  provenance = provenance
)
stopifnot(state$design_hash == design_hash,
  state$input_hash == dataset$input_hash,
  identical(state$provenance, provenance))
if (!is.null(state$rng_state)) {
  assign(".Random.seed", state$rng_state, envir = .GlobalEnv)
}

status <- "success"
error <- ""
block_started <- proc.time()[["elapsed"]]
block_running <- FALSE
cat(sprintf("%s %s start; resumed=%s\n", Sys.time(), row$id, !is.null(state$fit)))
flush.console()

compact <- tryCatch(withCallingHandlers({
  repeat {
    start_iteration <- if (is.null(state$fit)) 0L else state$fit$fit$iter
    block_started <- proc.time()[["elapsed"]]
    block_running <- TRUE
    if (is.null(state$fit)) {
      set.seed(row$seed + design$seed_offset)
      block_kind <- "initial"
      requested_budget <- design$initial_budget
      state$fit <- MPCurver::fit_mpcurve(
        X = dataset$X,
        intrinsic_dim = "auto",
        max_intrinsic_dim = design$max_intrinsic_dim,
        similarity_metric = design$similarity_metric,
        spline_r2_df = design$spline_r2_df,
        similarity_min_cluster_size = design$similarity_min_cluster_size,
        partition_init = "similarity",
        partition_prior = design$partition_prior,
        effective_weight_tol = design$effective_weight_tol,
        algorithm = "cavi",
        method = "isomap",
        K = design$K,
        rw_q = design$rw_q,
        ridge = design$ridge,
        lambda = 1,
        fix_lambda = FALSE,
        S = NULL,
        position_prior = design$position_prior,
        cluster_linkage = design$cluster_linkage,
        discretization = design$discretization,
        num_cores = 1L,
        iter = design$initial_budget,
        tol = design$tolerance,
        convergence = design$convergence,
        T_start = design$T_start,
        T_end = design$T_end,
        n_outer = design$n_outer,
        inner_iter = design$inner_iter,
        max_converge_iter = design$initial_budget,
        tol_outer = design$tolerance,
        verbose = FALSE
      )
    } else {
      remaining <- design$maximum_sweeps - state$fit$fit$iter - 1L
      if (remaining <= 0L || isTRUE(state$fit$fit$converged)) {
        block_running <- FALSE
        break
      }
      block_kind <- "continuation"
      requested_budget <- min(design$continuation_budget, remaining)
      state$fit <- MPCurver::do_mpcurve(
        state$fit,
        iter = requested_budget,
        tol = design$tolerance,
        tol_outer = design$tolerance,
        convergence = design$convergence,
        verbose = FALSE
      )
      state$continuations <- state$continuations + 1L
    }

    elapsed <- proc.time()[["elapsed"]] - block_started
    state$elapsed_seconds <- state$elapsed_seconds + elapsed
    block_running <- FALSE
    diagnostics <- fit_trace_diagnostics(state$fit)
    state$block_history <- rbind(state$block_history, data.frame(
      block = nrow(state$block_history) + 1L,
      kind = block_kind,
      requested_budget = as.integer(requested_budget),
      start_iteration = as.integer(start_iteration),
      end_iteration = as.integer(state$fit$fit$iter),
      score = diagnostics$score,
      converged = diagnostics$converged,
      elapsed_seconds = elapsed
    ))
    state$rng_state <- if (exists(".Random.seed", envir = .GlobalEnv,
      inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
    stopifnot(identical(state$fit$data, dataset$X))
    atomic_save(state, checkpoint_path)
    cat(sprintf("%s %s sweeps=%d selected_M=%d score=%.6f converged=%s elapsed=%.2f\n",
      Sys.time(), row$id, state$fit$fit$iter, state$fit$intrinsic_dim,
      diagnostics$score, diagnostics$converged, state$elapsed_seconds))
    flush.console()
    if (diagnostics$converged ||
        state$fit$fit$iter >= design$maximum_sweeps - 1L) break
  }
  compact_auto_fit(state$fit, dataset)
}, warning = function(w) {
  state$warnings <<- c(state$warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) {
  error <<- conditionMessage(e)
  if (block_running) {
    state$elapsed_seconds <<- state$elapsed_seconds +
      proc.time()[["elapsed"]] - block_started
  }
  NULL
})

if (is.null(compact)) {
  status <- "error"
} else if (!compact$converged) {
  status <- "nonconverged"
}
evaluation <- if (!is.null(compact)) evaluate_compact_fit(compact, dataset) else NULL
result_row <- data.frame(
  id = row$id,
  phase = row$phase,
  true_M = row$true_M,
  snr = row$snr,
  replicate = row$replicate,
  method = design$method,
  status = status,
  selected_initial_M = if (is.null(compact)) NA_integer_ else compact$selected_initial_M,
  selected_model_M = if (is.null(compact)) NA_integer_ else compact$M,
  estimated_M = if (is.null(compact)) NA_integer_ else compact$effective_M,
  effective_M = if (is.null(compact)) NA_integer_ else compact$effective_M,
  initial_partition_ARI = if (is.null(compact)) NA_real_ else compact$initial_partition_ARI,
  ARI = if (is.null(evaluation)) NA_real_ else evaluation$ARI,
  ordering_recovery = if (is.null(evaluation)) NA_real_ else evaluation$ordering_recovery,
  matched_ordering_recovery = if (is.null(evaluation)) NA_real_ else
    evaluation$matched_ordering_recovery,
  true_ordering_coverage = if (is.null(evaluation)) NA_real_ else
    evaluation$true_ordering_coverage,
  minimum_initial_cluster_size = if (is.null(compact)) NA_integer_ else
    min(compact$dimension_initialization$cluster_sizes),
  elapsed_seconds = state$elapsed_seconds,
  candidate_count = 1L,
  warning_count = length(unique(state$warnings)),
  continuations = state$continuations,
  stop_reason = if (status == "success") "automatic_similarity_dimension" else status,
  error = error
)

result <- list(
  row = result_row,
  design_hash = design_hash,
  input_hash = dataset$input_hash,
  provenance = provenance,
  compact = compact,
  evaluation = evaluation,
  warnings = unique(state$warnings),
  block_history = state$block_history
)
atomic_save(result, output_path)
if (row$replicate == 1L || status != "success") {
  atomic_save(state, file.path(extension_dir, "full_fits", paste0(row$id, ".rds")))
}
if (status == "success" && file.exists(checkpoint_path)) unlink(checkpoint_path)

cat("DATASET COMPLETE\n")
print(result_row, row.names = FALSE)
if (status != "success") quit(status = 1L)
