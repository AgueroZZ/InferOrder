# Test the connected, first-violation interval using the unchanged fitting API.
source("experiments/isomap_elbo_screen_v040/common.R")
bounds <- readRDS(file.path(experiment_dir, "graph_bounds.rds"))
interval <- bounds$summary[bounds$summary$interpretation ==
  "first_violation_from_connectivity", , drop = FALSE]
stopifnot(nrow(interval) == 1L, !interval$empty_interval,
  interval$k_min == 2L, interval$k_max == 4L)
neighbors <- seq.int(interval$k_min, interval$k_max)
reference <- readRDS(file.path(experiment_dir, "results.rds"))
stopifnot(identical(digest::digest(observations, algo = "sha256"),
  reference$provenance$group_input_sha256),
  identical(package_source_hashes, reference$provenance$package_source_hashes))

warnings <- character()
result <- withCallingHandlers({
  choices <- score_candidates(neighbors)
  stopifnot(all(choices$failures == ""))
  # The choice is fixed using early ELBO before any converged fit is available.
  cat(sprintf("Narrow interval [%d,%d]: one-sweep selection k=%d\n",
    min(neighbors), max(neighbors), choices$selected_k))
  final_fits <- choices$states
  run_order <- c(choices$selected_index,
    setdiff(seq_along(neighbors), choices$selected_index))
  candidates <- vector("list", length(neighbors))
  for (index in run_order) {
    k <- neighbors[index]
    final <- continue_candidate(choices$states[[index]])
    raw <- do.call(isomap_ordering,
      c(list(X = observations), initial_control(k)$method_args))
    stopifnot(raw$n_components == 1L,
      length(raw$keep_idx) == nrow(observations), all(is.finite(raw$t)),
      final$K == design$num_bins, final$converged,
      all(is.finite(position(final))),
      identical(final$elbo_trace[1:2], choices$states[[index]]$elbo_trace),
      min(diff(final$elbo_trace)) >= -1e-8,
      abs(tail(diff(final$elbo_trace), 1L)) / length(observations) < design$tolerance,
      identical(final$fit$init_info, choices$states[[index]]$fit$init_info),
      identical(final$fit$control$method, "isomap"))
    final_fits[[index]] <- final
    candidates[[index]] <- list(k = k, status = "converged",
      raw_position = as.numeric(raw$t), graph_components = raw$n_components,
      early = compact_fit(choices$states[[index]]), final = compact_fit(final))
    cat(sprintf("k=%d: raw rho=%.6f, one-step ELBO=%.6f, final ELBO=%.6f, final rho=%.6f, sweeps=%d\n",
      k, recovery(raw$t), choices$scores[index], last_elbo(final),
      recovery(position(final)), final$fit$iter))
  }
  uninterrupted <- fit_candidate(choices$selected_k, design$total_budget)
  resumed <- final_fits[[choices$selected_index]]
  continuation_check <- c(
    positions = max(abs(position(uninterrupted) - position(resumed))),
    responsibilities = max(abs(uninterrupted$gamma - resumed$gamma)),
    trajectories = max(abs(uninterrupted$params$mu - resumed$params$mu)),
    elbo_trace = max(abs(uninterrupted$elbo_trace - resumed$elbo_trace)))
  stopifnot(uninterrupted$converged,
    uninterrupted$fit$iter == resumed$fit$iter, all(continuation_check < 1e-8))
  saveRDS(list(early = choices$states, final = final_fits),
    file.path(experiment_dir, "full_fits", "narrow_candidates.rds"))
  list(candidate_neighbors = neighbors, selected_k = choices$selected_k,
    range_convention = "first violation of degree rule after connectivity",
    candidates = candidates, continuation_check = continuation_check)
}, warning = function(warning) {
  warnings <<- c(warnings, conditionMessage(warning))
  invokeRestart("muffleWarning")
})

comparison <- c(result$candidates,
  reference$candidates[match(c(10L, 15L), design$dense_neighbors)])
summary <- do.call(rbind, lapply(comparison, function(candidate) data.frame(
  k = candidate$k,
  source = if (candidate$k %in% neighbors) "narrow-range follow-up" else "original 5:30 experiment",
  selected_in_narrow_range = candidate$k == result$selected_k,
  raw_rho = recovery(candidate$raw_position),
  one_step_elbo = candidate$early$elbo_trace[2L],
  one_step_rho = recovery(candidate$early$position),
  final_elbo = tail(candidate$final$elbo_trace, 1L),
  final_rho = recovery(candidate$final$position), iterations = candidate$final$iter,
  converged = candidate$final$converged)))
result$summary <- summary
result$warnings <- warnings
result$provenance <- provenance()
# Describe this follow-up's actual candidate set separately from the parent design.
result$provenance$design$dense_neighbors <- neighbors
result$provenance$design$sparse_neighbors <- NULL
result$provenance$design$timing_seed <- NULL
result$provenance$design$timing_repeats <- NULL
result$provenance$range_source_sha256 <- digest::digest(
  file = file.path(experiment_dir, "graph_bounds.rds"), algo = "sha256")
result$provenance$followup_script_hashes <- setNames(vapply(
  file.path(experiment_dir, c("fit_narrow_range.R", "plot_narrow_range.R")),
  function(path) digest::digest(file = path, algo = "sha256"), character(1)),
  c("fit_narrow_range.R", "plot_narrow_range.R"))
write.csv(summary, file.path(experiment_dir, "narrow_range_summary.csv"), row.names = FALSE)
saveRDS(result, file.path(experiment_dir, "narrow_range.rds"))
stopifnot(length(warnings) == 0L)
print(summary, row.names = FALSE)
print(result$continuation_check)
