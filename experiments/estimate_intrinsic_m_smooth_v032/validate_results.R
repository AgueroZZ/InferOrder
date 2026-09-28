#!/usr/bin/env Rscript
# Verify every saved outcome, including scientific failures and stopping steps.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "fit.R"))
source(file.path(study_dir, "evaluate.R"))
source(file.path(study_dir, "reporting_metrics.R"))
manifest <- active_manifest()
manifest <- manifest[manifest$phase == "main", , drop = FALSE]
stopifnot(nrow(manifest) == 90L, all(table(manifest$true_M, manifest$snr) == 10L))
provenance <- execution_provenance()
outcome_rows <- list()
candidate_rows <- list()

validate_candidate <- function(candidate, dataset) {
  stopifnot(candidate$dataset_id == dataset$row$id,
    candidate$method %in% c("adaptive", "forward"),
    candidate$M >= 1L, candidate$M <= design$max_M,
    candidate$design_hash == design_hash, candidate$input_hash == dataset$input_hash,
    identical(candidate$provenance, provenance),
    candidate$status %in% c("success", "error", "nonconverged"),
    is.character(candidate$warnings), is.character(candidate$error),
    length(candidate$error) == 1L, is.finite(candidate$elapsed_seconds),
    candidate$elapsed_seconds >= 0)
  blocks <- candidate$block_history
  stopifnot(nrow(blocks) == candidate$continuations + as.integer(nrow(blocks) > 0L))
  if (nrow(blocks)) {
    stopifnot(identical(blocks$block, seq_len(nrow(blocks))),
      blocks$kind[1L] == "initial", blocks$requested_budget[1L] == design$initial_budget,
      blocks$start_iteration[1L] == 0L, all(blocks$elapsed_seconds >= 0),
      all(blocks$end_iteration <= design$maximum_sweeps - 1L),
      all(blocks$end_iteration > blocks$start_iteration))
    if (nrow(blocks) > 1L) {
      stopifnot(all(blocks$kind[-1L] == "continuation"),
        all(blocks$start_iteration[-1L] == head(blocks$end_iteration, -1L)),
        all(blocks$requested_budget[-1L] == pmin(design$continuation_budget,
          design$maximum_sweeps - blocks$start_iteration[-1L] - 1L)),
        !any(head(blocks$converged, -1L)))
    }
    stopifnot(candidate$elapsed_seconds + 1e-8 >= sum(blocks$elapsed_seconds))
  }
  diagnostics <- candidate$diagnostics
  if (!is.null(diagnostics)) {
    stopifnot(diagnostics$iterations <= design$maximum_sweeps - 1L,
      diagnostics$convergence == design$convergence, diagnostics$tolerance == design$tolerance)
    trace <- diagnostics$objective
    temperature <- diagnostics$temperature
    T1 <- trace[temperature == 1]
    valid_trace <- length(trace) == length(temperature) && length(T1) > 1L &&
      all(is.finite(trace)) && all(is.finite(temperature))
    nondecrease <- valid_trace &&
      all(diff(T1) >= -objective_nondecrease_guard(head(T1, -1L), candidate$M))
    final_increment <- if (length(T1) > 1L) as.numeric(tail(diff(T1), 1L)) else NA_real_
    final_previous <- if (length(T1) > 1L) as.numeric(tail(head(T1, -1L), 1L)) else NA_real_
    normalized_increment <- abs(final_increment) / (design$n * design$D)
    stopping_passed <- valid_trace && tail(temperature, 1L) == 1 &&
      final_increment >= -objective_nondecrease_guard(final_previous, candidate$M, stopping = TRUE) &&
      normalized_increment < design$tolerance
    stopifnot(identical(diagnostics$valid_trace, valid_trace),
      identical(diagnostics$T1_nondecrease_passed, nondecrease),
      identical(diagnostics$stopping_passed, stopping_passed),
      identical(diagnostics$last_normalized_increment, normalized_increment))
    if (nrow(blocks)) {
      stopifnot(tail(blocks$end_iteration, 1L) == diagnostics$iterations,
        tail(blocks$converged, 1L) == diagnostics$converged)
    }
    if (candidate$status != "error") {
      stopifnot(valid_trace, nondecrease, tail(temperature, 1L) == 1,
        identical(diagnostics$objective, candidate$fit$objective),
        identical(diagnostics$temperature, candidate$fit$temperature),
        identical(diagnostics$iterations, candidate$fit$iterations))
    }
  }
  if (candidate$status == "error") {
    stopifnot(is.null(candidate$fit), nzchar(candidate$error))
  } else {
    fit <- candidate$fit
    stopifnot(!is.null(fit), fit$M == candidate$M, fit$selected_model_M == candidate$M,
      fit$score == tail(fit$objective, 1L), all(is.finite(fit$W)), all(fit$W >= 0),
      max(abs(rowSums(fit$W) - 1)) < 1e-10,
      all(is.finite(fit$positions)), all(is.finite(fit$sigma2)), all(fit$sigma2 > 0),
      identical(fit$active, which(colMeans(fit$W) > design$effective_weight_tol)),
      fit$effective_M == length(fit$active),
      identical(fit$occupancy, colMeans(fit$W)),
      all(fit$assignments %in% fit$active))
    if (candidate$method == "adaptive") {
      stopifnot(candidate$M == design$max_M,
        max(abs(fit$omega - colMeans(fit$W))) < 1e-10)
    } else {
      stopifnot(max(abs(fit$omega - 1 / candidate$M)) < 1e-10)
    }
    if (candidate$status == "success") {
      stopifnot(fit$converged, diagnostics$converged, diagnostics$stopping_passed,
        diagnostics$last_normalized_increment < design$tolerance,
        !nzchar(candidate$error))
    } else {
      stopifnot(!fit$converged, !diagnostics$converged,
        diagnostics$iterations >= design$maximum_sweeps - 1L, !nzchar(candidate$error))
    }
  }
  data.frame(key = candidate$key, dataset_id = candidate$dataset_id,
    method = candidate$method, M = candidate$M, status = candidate$status,
    iterations = if (is.null(diagnostics)) NA_integer_ else diagnostics$iterations,
    converged = if (is.null(diagnostics)) FALSE else diagnostics$converged,
    T1_nondecrease_passed = if (is.null(diagnostics)) NA else diagnostics$T1_nondecrease_passed,
    stopping_passed = if (is.null(diagnostics)) NA else diagnostics$stopping_passed,
    continuations = candidate$continuations, warning_count = length(candidate$warnings),
    elapsed_seconds = candidate$elapsed_seconds, error = candidate$error)
}

for (i in seq_len(nrow(manifest))) {
  stopifnot(file.exists(file.path(study_dir, "data", paste0(manifest$id[i], ".rds"))))
  dataset <- load_dataset(manifest[i, , drop = FALSE])
  stopifnot(identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
  for (method in c("adaptive", "forward")) {
    result_path <- file.path(study_dir, "results", paste0(manifest$id[i], "_", method, ".rds"))
    stopifnot(file.exists(result_path))
    result <- readRDS(result_path)
    row <- result$row
    stopifnot(result$design_hash == design_hash, result$input_hash == dataset$input_hash,
      identical(result$provenance, provenance), row$id == dataset$row$id,
      row$method == method, row$status %in% c("success", "error", "nonconverged"),
      row$true_M == dataset$row$true_M, row$snr == dataset$row$snr,
      row$replicate == dataset$row$replicate,
      length(result$candidate_keys) == row$candidate_count,
      !anyDuplicated(result$candidate_keys))
    candidates <- lapply(result$candidate_keys, function(key) {
      path <- file.path(study_dir, "candidates", paste0(key, ".rds"))
      stopifnot(file.exists(path))
      candidate <- readRDS(path)
      stopifnot(candidate$key == key, candidate$method == method)
      candidate_rows[[length(candidate_rows) + 1L]] <<- validate_candidate(candidate, dataset)
      candidate
    })
    stopifnot(abs(row$elapsed_seconds -
      sum(vapply(candidates, function(x) x$elapsed_seconds, numeric(1)))) < 1e-8,
      row$warning_count == sum(vapply(candidates, function(x) length(x$warnings), integer(1))))
    terminal <- candidates[[length(candidates)]]
    selected_index <- NA_integer_
    if (method == "adaptive") {
      stopifnot(length(candidates) == 1L, terminal$M == design$max_M,
        is.null(result$history))
      if (terminal$status == "success") selected_index <- 1L
      stopifnot(row$stop_reason == if (terminal$status == "success")
        "fixed_maximum_dimension" else paste0("candidate_", terminal$status))
    } else {
      candidate_M <- vapply(candidates, function(x) x$M, integer(1))
      stopifnot(identical(candidate_M, seq_along(candidates)))
      expected_history <- data.frame(current_M = integer(), candidate_M = integer(),
        current_score = numeric(), candidate_score = numeric(), delta = numeric(), accepted = logical())
      current_index <- 1L
      stopped <- FALSE
      for (j in seq_along(candidates)) {
        candidate <- candidates[[j]]
        if (candidate$status != "success") {
          stopifnot(j == length(candidates))
          current_index <- NA_integer_
          stopped <- TRUE
          break
        }
        if (j == 1L) next
        previous <- candidates[[current_index]]
        delta <- candidate$fit$score - previous$fit$score
        accepted <- is.finite(delta) && delta > 0
        expected_history <- rbind(expected_history,
          data.frame(current_M = previous$M, candidate_M = candidate$M,
            current_score = previous$fit$score, candidate_score = candidate$fit$score,
            delta = delta, accepted = accepted))
        if (!accepted) {
          stopifnot(j == length(candidates))
          stopped <- TRUE
          break
        }
        current_index <- j
      }
      stopifnot(isTRUE(all.equal(result$history, expected_history, tolerance = 0)),
        stopped || length(candidates) == design$max_M)
      selected_index <- current_index
      expected_stop_reason <- if (terminal$status != "success") paste0("candidate_", terminal$status) else
        if (stopped) "first_nonimprovement" else "maximum_dimension"
      stopifnot(row$stop_reason == expected_stop_reason)
    }
    if (is.na(selected_index)) {
      stopifnot(row$status == terminal$status, row$status != "success",
        is.null(result$selected), is.null(result$evaluation),
        is.na(result$selected_candidate_key), is.na(row$estimated_M),
        is.na(row$effective_M), is.na(row$selected_model_M), is.na(row$ARI),
        is.na(row$ordering_recovery), row$error == terminal$error)
    } else {
      selected <- candidates[[selected_index]]
      stopifnot(row$status == "success", selected$status == "success",
        identical(result$selected, selected$fit), result$selected_candidate_key == selected$key,
        row$estimated_M == reported_effective_M(result), row$effective_M == row$estimated_M,
        row$selected_model_M == selected$M)
      evaluation <- evaluate_fit(selected$fit, dataset)
      stopifnot(isTRUE(all.equal(result$evaluation, evaluation, tolerance = 1e-12)),
        isTRUE(all.equal(row$ARI, evaluation$ARI, tolerance = 1e-12)),
        isTRUE(all.equal(row$ordering_recovery, evaluation$ordering_recovery, tolerance = 1e-12)),
        isTRUE(all.equal(row$matched_ordering_recovery, evaluation$matched_ordering_recovery, tolerance = 1e-12)),
        isTRUE(all.equal(row$true_ordering_coverage, evaluation$true_ordering_coverage, tolerance = 1e-12)),
        row$map_groups == length(unique(selected$fit$assignments)), !nzchar(row$error))
    }
    outcome_rows[[length(outcome_rows) + 1L]] <- reported_result_row(result)
  }
}
outcomes <- do.call(rbind, outcome_rows)
candidates <- do.call(rbind, candidate_rows)
stopifnot(nrow(outcomes) == 180L, !anyDuplicated(paste(outcomes$id, outcomes$method)),
  !anyDuplicated(candidates$key))
write.csv(candidates, file.path(study_dir, "candidate_validation.csv"), row.names = FALSE)
write.csv(outcomes, file.path(study_dir, "outcome_validation.csv"), row.names = FALSE)
jsonlite::write_json(list(validated = TRUE, datasets = nrow(manifest),
  method_outcomes = nrow(outcomes), candidates = nrow(candidates),
  successful_outcomes = sum(outcomes$status == "success"),
  unresolved_outcomes = sum(outcomes$status != "success"),
  error_outcomes = sum(outcomes$status == "error"),
  nonconverged_outcomes = sum(outcomes$status == "nonconverged"),
  candidate_errors = sum(candidates$status == "error"),
  candidate_nonconverged = sum(candidates$status == "nonconverged"),
  objective_guard_violations = sum(candidates$T1_nondecrease_passed == FALSE, na.rm = TRUE),
  candidate_warnings = sum(candidates$warning_count),
  all_saved_outcomes_accounted_for = TRUE, converged_scores_at_temperature_one = TRUE,
  posterior_occupancy_and_ordering_matching = TRUE, forward_decisions_consistent = TRUE,
  design_hash = design_hash, fitting_fingerprint = provenance$fingerprint),
  file.path(study_dir, "main_validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Validated all 180 method outcomes; unresolved outcomes:",
  sum(outcomes$status != "success"), "\n")
