#!/usr/bin/env Rscript

# Experiment-only refit using a fixed-df directional spline-R-squared similarity.
# The installed MPCurver namespace is overridden only within this R process;
# package source, installed files, and the study's published results are unchanged.

source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "fit.R"))
source(file.path(study_dir, "evaluate.R"))

arguments <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(arguments) >= 1L) arguments[[1L]] else "main_M5_S1_r001"
spline_df <- if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 5L
stopifnot(length(spline_df) == 1L, is.finite(spline_df), spline_df >= 1L)

manifest <- active_manifest()
row <- manifest[manifest$id == dataset_id, , drop = FALSE]
stopifnot(nrow(row) == 1L, row$phase == "main")
dataset <- load_dataset(row)

metric_name <- paste0("spline_r2_df", spline_df)
output_dir <- file.path(study_dir, "exploratory_similarity", dataset_id,
  paste0(metric_name, "_fits"))
candidate_dir <- file.path(output_dir, "candidates")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

experiment_design <- list(
  baseline_design_hash = design_hash,
  dataset_id = dataset_id,
  input_hash = dataset$input_hash,
  similarity_metric = "directional natural-cubic-spline training R-squared",
  spline_df = spline_df,
  predictor = "normalized feature rank",
  symmetrization = "maximum of two directional R-squared values",
  cluster_linkage = "single",
  fitting_settings = design)
experiment_hash <- digest::digest(experiment_design, algo = "sha256")
script_path <- file.path(study_dir, "run_spline_r2_initialization.R")
provenance <- list(
  experiment_hash = experiment_hash,
  script_sha256 = digest::digest(file = script_path, algo = "sha256"),
  input_hash = dataset$input_hash,
  package_version = as.character(utils::packageVersion("MPCurver")),
  package_source_commit = design$package_source_commit,
  R_version = as.character(getRversion()))

compute_spline_r2_similarity <- function(X, min_feature_sd = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  feature_names <- colnames(X)
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(D))
  feature_sd <- apply(X, 2L, stats::sd)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd

  rank_position <- (seq_len(n) - 0.5) / n
  spline_basis <- cbind("(Intercept)" = 1,
    splines::ns(rank_position, df = spline_df, intercept = FALSE))
  basis_qr <- qr(spline_basis)
  Q <- qr.Q(basis_qr)
  stopifnot(basis_qr$rank == ncol(spline_basis))

  centered_X <- sweep(X, 2L, colMeans(X), FUN = "-")
  total_sum_squares <- colSums(centered_X^2)
  directional_r2 <- matrix(0, D, D, dimnames = list(feature_names, feature_names))
  valid_response <- !low_variance & total_sum_squares > sqrt(.Machine$double.eps)
  for (predictor in seq_len(D)) {
    if (low_variance[predictor]) next
    ordered_X <- X[order(X[, predictor], method = "radix"), , drop = FALSE]
    projected_coordinates <- crossprod(Q, ordered_X)
    residual_sum_squares <- pmax(
      colSums(ordered_X^2) - colSums(projected_coordinates^2), 0)
    directional_r2[predictor, valid_response] <-
      1 - residual_sum_squares[valid_response] / total_sum_squares[valid_response]
  }
  directional_r2 <- pmin(pmax(directional_r2, 0), 1)
  similarity <- pmax(directional_r2, t(directional_r2))
  if (any(low_variance)) {
    similarity[low_variance, ] <- 0
    similarity[, low_variance] <- 0
  }
  diag(similarity) <- 1
  list(
    S = similarity,
    distance = 1 - similarity,
    metric = metric_name,
    feature_info = data.frame(feature = feature_names, sd = as.numeric(feature_sd),
      low_variance = low_variance, stringsAsFactors = FALSE),
    directional_r2 = directional_r2,
    spline_df = spline_df)
}

similarity_cache <- NULL
similarity_override <- function(X,
                                metric = c("spearman", "pearson"),
                                use = "pairwise.complete.obs",
                                abs_value = TRUE,
                                min_feature_sd = 1e-8) {
  metric <- match.arg(metric)
  if (is.null(similarity_cache)) {
    similarity_cache <<- compute_spline_r2_similarity(X, min_feature_sd)
    stopifnot(identical(digest::digest(X, algo = "sha256"), dataset$input_hash))
  } else {
    stopifnot(identical(dim(similarity_cache$S), c(ncol(X), ncol(X))))
  }
  similarity_cache[c("S", "distance", "metric", "feature_info")]
}

original_similarity_cor <- getFromNamespace(
  ".compute_same_ordering_similarity_cor", "MPCurver")
assignInNamespace(".compute_same_ordering_similarity_cor", similarity_override,
  ns = "MPCurver")
on.exit(assignInNamespace(".compute_same_ordering_similarity_cor",
  original_similarity_cor, ns = "MPCurver"), add = TRUE)

annotate_similarity <- function(fit) {
  if (is.null(fit) || !inherits(fit$fit, "soft_partition_cavi")) return(fit)
  fit$fit$control$similarity_metric <- metric_name
  if (!is.null(fit$fit$similarity_init)) {
    fit$fit$similarity_init$similarity_metric <- metric_name
    fit$fit$similarity_init$directional_r2 <- similarity_cache$directional_r2
    fit$fit$similarity_init$spline_df <- spline_df
  }
  fit
}

run_spline_candidate <- function(dataset, method, M) {
  stopifnot(method %in% c("adaptive", "forward"), M >= 2L, M <= design$max_M,
    identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
  key <- sprintf("%s_%s_M%02d", dataset$row$id, method, M)
  output_path <- file.path(candidate_dir, paste0(key, ".rds"))
  checkpoint_path <- file.path(checkpoint_dir, paste0(key, ".rds"))
  if (file.exists(output_path)) {
    saved <- readRDS(output_path)
    reusable <- identical(saved$experiment_hash, experiment_hash) &&
      identical(saved$input_hash, dataset$input_hash) &&
      identical(saved$provenance, provenance) &&
      saved$status %in% c("success", "nonconverged")
    if (reusable) return(saved)
  }

  prior <- if (method == "adaptive") "adaptive" else "fixed"
  state <- if (file.exists(checkpoint_path)) readRDS(checkpoint_path) else
    list(fit = NULL, elapsed_seconds = 0, warnings = character(), continuations = 0L,
      block_history = empty_block_history(), rng_state = NULL,
      experiment_hash = experiment_hash, input_hash = dataset$input_hash,
      provenance = provenance)
  stopifnot(state$experiment_hash == experiment_hash,
    state$input_hash == dataset$input_hash, identical(state$provenance, provenance))
  if (!is.null(state$rng_state)) {
    assign(".Random.seed", state$rng_state, envir = .GlobalEnv)
  }

  status <- "success"
  error <- ""
  fit_started <- proc.time()[["elapsed"]]
  block_running <- FALSE
  cat(sprintf("%s %s %s start; resumed=%s\n",
    Sys.time(), metric_name, key, !is.null(state$fit)))
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
          effective_weight_tol = design$effective_weight_tol,
          hard_assign_final = FALSE), fit_options())
        state$fit <- annotate_similarity(do.call(MPCurver::fit_mpcurve, args))
      } else {
        remaining <- design$maximum_sweeps - state$fit$fit$iter - 1L
        if (remaining <= 0L || isTRUE(state$fit$fit$converged)) {
          block_running <- FALSE
          break
        }
        block_kind <- "continuation"
        requested_budget <- min(design$continuation_budget, remaining)
        state$fit <- annotate_similarity(MPCurver::do_mpcurve(state$fit,
          iter = requested_budget, tol = design$tolerance,
          tol_outer = design$tolerance, convergence = design$convergence,
          verbose = FALSE))
        state$continuations <- state$continuations + 1L
      }
      elapsed <- proc.time()[["elapsed"]] - fit_started
      state$elapsed_seconds <- state$elapsed_seconds + elapsed
      block_running <- FALSE
      state$block_history <- rbind(state$block_history,
        data.frame(block = nrow(state$block_history) + 1L, kind = block_kind,
          requested_budget = as.integer(requested_budget),
          start_iteration = as.integer(start_iteration),
          end_iteration = as.integer(state$fit$fit$iter),
          score = fit_score(state$fit, M),
          converged = isTRUE(state$fit$fit$converged), elapsed_seconds = elapsed))
      state$rng_state <- if (exists(".Random.seed", envir = .GlobalEnv,
          inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
      stopifnot(identical(state$fit$data, dataset$X))
      atomic_save(state, checkpoint_path)
      cat(sprintf("%s %s sweeps=%d score=%.6f converged=%s elapsed=%.2f\n",
        Sys.time(), key, state$fit$fit$iter, fit_score(state$fit, M),
        state$fit$fit$converged, state$elapsed_seconds))
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
      state$elapsed_seconds <<- state$elapsed_seconds +
        proc.time()[["elapsed"]] - fit_started
    }
    NULL
  })

  if (is.null(value)) {
    status <- "error"
  } else if (!value$converged) {
    status <- "nonconverged"
  }
  saved <- list(key = key, dataset_id = dataset$row$id, method = method,
    M = as.integer(M), status = status, error = error,
    warnings = unique(state$warnings), fit = value,
    diagnostics = fit_trace_diagnostics(state$fit, M),
    elapsed_seconds = state$elapsed_seconds, continuations = state$continuations,
    block_history = state$block_history, experiment_hash = experiment_hash,
    input_hash = dataset$input_hash, provenance = provenance)
  atomic_save(saved, output_path)
  if (status == "success" && file.exists(checkpoint_path)) unlink(checkpoint_path)
  cat(sprintf("%s %s finished: %s %s\n", Sys.time(), key, status, error))
  flush.console()
  saved
}

# Adaptive EB retains the original fixed maximum dimension M = 8.
adaptive_candidate <- run_spline_candidate(dataset, "adaptive", design$max_M)
adaptive_selected <- if (adaptive_candidate$status == "success") adaptive_candidate else NULL

# M = 1 does not use feature-similarity initialization, so reuse the exact
# published candidate and refit only M >= 2 under the spline-R-squared metric.
original_M1_path <- file.path(study_dir, "candidates",
  paste0(dataset_id, "_forward_M01.rds"))
stopifnot(file.exists(original_M1_path))
forward_candidates <- list(readRDS(original_M1_path))
stopifnot(forward_candidates[[1L]]$status == "success",
  forward_candidates[[1L]]$M == 1L,
  forward_candidates[[1L]]$input_hash == dataset$input_hash)
forward_selected <- forward_candidates[[1L]]
forward_history <- data.frame(current_M = integer(), candidate_M = integer(),
  current_score = numeric(), candidate_score = numeric(), delta = numeric(),
  accepted = logical())
forward_stop_reason <- "maximum_dimension"

for (M in 2:design$max_M) {
  candidate <- run_spline_candidate(dataset, "forward", M)
  forward_candidates[[M]] <- candidate
  if (candidate$status != "success") {
    forward_selected <- NULL
    forward_stop_reason <- paste0("candidate_", candidate$status)
    break
  }
  delta <- candidate$fit$score - forward_selected$fit$score
  accepted <- is.finite(delta) && delta > 0
  forward_history <- rbind(forward_history,
    data.frame(current_M = forward_selected$M, candidate_M = M,
      current_score = forward_selected$fit$score,
      candidate_score = candidate$fit$score, delta = delta, accepted = accepted))
  if (!accepted) {
    forward_stop_reason <- "first_nonimprovement"
    break
  }
  forward_selected <- candidate
}

summarize_selected <- function(method, selected, candidates, stop_reason) {
  successful <- !is.null(selected) && selected$status == "success"
  evaluation <- if (successful) evaluate_fit(selected$fit, dataset) else NULL
  data.frame(
    id = dataset_id,
    method = method,
    status = if (successful) "success" else "unresolved",
    true_M = dataset$row$true_M,
    selected_model_M = if (successful) selected$M else NA_integer_,
    effective_M = if (successful) selected$fit$effective_M else NA_integer_,
    exact_M = if (successful) selected$fit$effective_M == dataset$row$true_M else FALSE,
    ARI = if (successful) evaluation$ARI else NA_real_,
    ordering_recovery = if (successful) evaluation$ordering_recovery else NA_real_,
    matched_ordering_recovery = if (successful) evaluation$matched_ordering_recovery else NA_real_,
    true_ordering_coverage = if (successful) evaluation$true_ordering_coverage else NA_real_,
    candidate_count = length(candidates),
    total_elapsed_seconds = sum(vapply(candidates, function(x) x$elapsed_seconds,
      numeric(1))),
    warning_count = sum(vapply(candidates, function(x) length(x$warnings), integer(1))),
    stop_reason = stop_reason,
    stringsAsFactors = FALSE)
}

adaptive_summary <- summarize_selected("adaptive", adaptive_selected,
  list(adaptive_candidate),
  if (adaptive_candidate$status == "success") "fixed_maximum_dimension" else
    paste0("candidate_", adaptive_candidate$status))
forward_summary <- summarize_selected("forward", forward_selected,
  forward_candidates, forward_stop_reason)
summary_table <- rbind(adaptive_summary, forward_summary)

stopifnot(all(summary_table$status == "success"),
  all(vapply(c(list(adaptive_candidate), forward_candidates[-1L]),
    function(x) x$status == "success", logical(1))),
  all(vapply(c(list(adaptive_candidate), forward_candidates[-1L]),
    function(x) isTRUE(x$fit$converged), logical(1))))

utils::write.csv(summary_table, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(forward_history, file.path(output_dir, "forward_history.csv"),
  row.names = FALSE)
jsonlite::write_json(list(
  experiment_design = experiment_design,
  experiment_hash = experiment_hash,
  provenance = provenance,
  validation = list(
    package_files_modified = FALSE,
    published_results_modified = FALSE,
    M1_candidate_reused_because_similarity_independent = TRUE,
    all_new_candidates_converged = TRUE,
    vectorized_similarity_cached_within_process = TRUE)),
  file.path(output_dir, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
saveRDS(list(summary = summary_table, forward_history = forward_history,
  adaptive = adaptive_selected, forward = forward_selected,
  experiment_design = experiment_design, provenance = provenance),
  file.path(output_dir, "selected_results.rds"), compress = FALSE)

cat("\nSPLINE INITIALIZATION EXPERIMENT COMPLETE\n")
print(summary_table, row.names = FALSE)
print(forward_history, row.names = FALSE)
