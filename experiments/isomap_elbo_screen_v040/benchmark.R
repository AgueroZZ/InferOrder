source("experiments/isomap_elbo_screen_v040/common.R")

screen_pipeline <- function(neighbors) {
  scoring <- timed(score_candidates(neighbors))
  choices <- scoring$value
  continuation <- timed(continue_candidate(choices$states[[choices$selected_index]]))
  list(fit = continuation$value, selected_k = choices$selected_k,
    screening_seconds = scoring$seconds, continuation_seconds = continuation$seconds,
    failures = choices$failures)
}

full_search_pipeline <- function() {
  # Endpoints remain concealed until every complete candidate has been fitted.
  fits <- lapply(design$dense_neighbors, function(k) {
    tryCatch(continue_candidate(fit_candidate(k, 2000L)), error = function(error) NULL)
  })
  scores <- vapply(fits, function(fit) if (is.null(fit)) -Inf else last_elbo(fit), numeric(1))
  best <- order(-scores, design$dense_neighbors)[1L]
  stopifnot(is.finite(scores[best]))
  list(fit = fits[[best]], selected_k = design$dense_neighbors[best],
    screening_seconds = NA_real_, continuation_seconds = NA_real_,
    failures = vapply(fits, is.null, logical(1)))
}

pipelines <- list(
  fixed_k15 = function() list(fit = continue_candidate(fit_candidate(15L, 2000L)),
    selected_k = 15L, screening_seconds = 0, continuation_seconds = NA_real_, failures = ""),
  dense_screen = function() screen_pipeline(design$dense_neighbors),
  sparse_screen = function() screen_pipeline(design$sparse_neighbors),
  dense_full_search = full_search_pipeline)
reference <- readRDS(file.path(experiment_dir, "results.rds"))
warnings <- character()
rows <- withCallingHandlers({
  for (method in names(pipelines)) {
    cat(sprintf("Warming %s\n", method))
    invisible(pipelines[[method]]())
  }
  set.seed(design$timing_seed)
  # Precompute the schedule because initializer seeds reset the global RNG.
  schedule <- lapply(seq_len(design$timing_repeats), function(repeat_index)
    sample(names(pipelines)))
  records <- list()
  for (repeat_index in seq_along(schedule)) {
    for (order_index in seq_along(schedule[[repeat_index]])) {
      method <- schedule[[repeat_index]][order_index]
      gc()
      run <- timed(pipelines[[method]]())
      endpoint <- run$value
      stopifnot(endpoint$fit$converged)
      expected <- reference$candidates[[match(endpoint$selected_k, design$dense_neighbors)]]$final
      stopifnot(max(abs(position(endpoint$fit) - expected$position)) < 1e-8,
        abs(last_elbo(endpoint$fit) - tail(expected$elbo_trace, 1L)) < 1e-8,
        !any(nzchar(as.character(endpoint$failures)) & as.character(endpoint$failures) != "FALSE"))
      records[[length(records) + 1L]] <- data.frame(repeat_index = repeat_index,
        order_index = order_index, method = method, selected_k = endpoint$selected_k,
        total_seconds = run$seconds, screening_seconds = endpoint$screening_seconds,
        continuation_seconds = endpoint$continuation_seconds,
        final_elbo = last_elbo(endpoint$fit), final_rho = recovery(position(endpoint$fit)),
        iterations = endpoint$fit$fit$iter, converged = endpoint$fit$converged)
      cat(sprintf("Repeat %d %s: %.3f seconds, k=%d\n",
        repeat_index, method, run$seconds, endpoint$selected_k))
    }
  }
  do.call(rbind, records)
}, warning = function(warning) {
  warnings <<- c(warnings, conditionMessage(warning))
  invokeRestart("muffleWarning")
})
write.csv(rows, file.path(experiment_dir, "timings.csv"), row.names = FALSE)
saveRDS(list(timings = rows, warnings = warnings, provenance = provenance()),
  file.path(experiment_dir, "benchmark.rds"))
stopifnot(length(warnings) == 0L)
