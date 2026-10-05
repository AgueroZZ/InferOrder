source("experiments/isomap_elbo_screen_v040/common.R")
dir.create(file.path(experiment_dir, "full_fits"), showWarnings = FALSE)
warnings <- character()
result <- withCallingHandlers({
  screening <- timed(score_candidates(design$dense_neighbors))
  choices <- screening$value
  # Choose using one-sweep ELBO before calculating any converged endpoints.
  selected_k <- choices$selected_k
  cat(sprintf("One-sweep selection: k=%d; screening %.3f seconds\n",
    selected_k, screening$seconds))
  selected_run <- timed(continue_candidate(choices$states[[choices$selected_index]]))
  final_fits <- choices$states
  continuation_seconds <- rep(NA_real_, length(final_fits))
  final_fits[[choices$selected_index]] <- selected_run$value
  continuation_seconds[choices$selected_index] <- selected_run$seconds
  for (index in seq_along(final_fits)) {
    if (is.null(final_fits[[index]]) || index == choices$selected_index) next
    run <- timed(continue_candidate(final_fits[[index]]))
    final_fits[[index]] <- run$value
    continuation_seconds[index] <- run$seconds
  }
  candidates <- lapply(seq_along(final_fits), function(index) {
    k <- choices$neighbors[index]
    if (is.null(final_fits[[index]])) return(list(k = k, status = "failed",
      failure = choices$failures[index]))
    raw <- do.call(isomap_ordering, c(list(X = observations), initial_control(k)$method_args))
    stopifnot(raw$n_components == 1L, all(is.finite(raw$t)),
      length(raw$keep_idx) == nrow(observations))
    early <- choices$states[[index]]
    final <- final_fits[[index]]
    cat(sprintf("k=%2d: first ELBO=%10.3f, final ELBO=%10.3f, rho=%.6f, sweeps=%d\n",
      k, choices$scores[index], last_elbo(final), recovery(position(final)), final$fit$iter))
    list(k = k, status = "converged", raw_position = as.numeric(raw$t),
      graph_components = raw$n_components, early = compact_fit(early),
      final = compact_fit(final), continuation_seconds = continuation_seconds[index])
  })
  # Uninterrupted public fits check that selection genuinely resumes a state.
  uninterrupted_checks <- lapply(unique(c(15L, selected_k)), function(k) {
    uninterrupted <- fit_candidate(k, design$total_budget)
    resumed <- final_fits[[match(k, choices$neighbors)]]
    checks <- c(position = max(abs(position(uninterrupted) - position(resumed))),
      responsibilities = max(abs(uninterrupted$gamma - resumed$gamma)),
      trajectory = max(abs(uninterrupted$params$mu - resumed$params$mu)),
      elbo_trace = max(abs(uninterrupted$elbo_trace - resumed$elbo_trace)))
    stopifnot(uninterrupted$converged, uninterrupted$fit$iter == resumed$fit$iter,
      all(checks < 1e-8))
    list(k = k, max_absolute_differences = checks)
  })
  saveRDS(list(early = choices$states, final = final_fits),
    file.path(experiment_dir, "full_fits", "candidates.rds"))
  list(selected_k = selected_k, candidates = candidates,
    screening_seconds = screening$seconds,
    selected_continuation_seconds = selected_run$seconds,
    uninterrupted_checks = uninterrupted_checks)
}, warning = function(warning) {
  warnings <<- c(warnings, conditionMessage(warning))
  invokeRestart("muffleWarning")
})
result$warnings <- warnings
result$provenance <- provenance()
saveRDS(result, file.path(experiment_dir, "results.rds"))
stopifnot(length(warnings) == 0L)
