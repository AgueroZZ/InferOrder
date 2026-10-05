source("experiments/isomap_elbo_screen_v040/common.R")
result <- readRDS(file.path(experiment_dir, "results.rds"))
benchmark <- readRDS(file.path(experiment_dir, "benchmark.rds"))
current <- provenance()
stopifnot(identical(current$group_input_sha256, result$provenance$group_input_sha256),
  identical(current$truth_sha256, result$provenance$truth_sha256),
  identical(current$package_source_hashes, result$provenance$package_source_hashes),
  identical(current$script_hashes, result$provenance$script_hashes),
  identical(current$group_input_sha256, benchmark$provenance$group_input_sha256),
  identical(current$script_hashes, benchmark$provenance$script_hashes),
  identical(design$dense_neighbors, vapply(result$candidates, `[[`, integer(1), "k")),
  length(result$warnings) == 0L, length(benchmark$warnings) == 0L)
scores <- rep(-Inf, length(result$candidates))
for (index in seq_along(result$candidates)) {
  candidate <- result$candidates[[index]]
  if (candidate$status != "converged") next
  stopifnot(candidate$early$iter == 1L, length(candidate$early$elbo_trace) == 2L,
    candidate$graph_components == 1L, candidate$final$converged,
    all(is.finite(candidate$final$position)),
    length(candidate$final$elbo_trace) == candidate$final$iter + 1L,
    identical(candidate$early$elbo_trace, candidate$final$elbo_trace[1:2]),
    min(diff(candidate$final$elbo_trace)) >= -1e-8,
    abs(tail(diff(candidate$final$elbo_trace), 1L)) / length(observations) < design$tolerance,
    identical(candidate$early$initialization$method_used, "isomap"),
    identical(candidate$early$initialization, candidate$final$initialization))
  scores[index] <- candidate$early$elbo_trace[2L]
}
stopifnot(result$selected_k == design$dense_neighbors[order(-scores, design$dense_neighbors)[1L]],
  nrow(benchmark$timings) == 4L * design$timing_repeats,
  all(benchmark$timings$converged), all(benchmark$timings$total_seconds > 0))

# Recompute representative public one-sweep scores independently of saved states.
score_checks <- lapply(unique(c(5L, 10L, 15L, result$selected_k, 30L)), function(k) {
  fit <- fit_candidate(k, 1L, tolerance = 0)
  expected <- result$candidates[[match(k, design$dense_neighbors)]]$early
  differences <- c(elbo = max(abs(fit$elbo_trace - expected$elbo_trace)),
    positions = max(abs(position(fit) - expected$position)))
  stopifnot(all(differences < 1e-10))
  list(k = k, max_absolute_differences = differences)
})

# Verify the installed source matches the recorded adjacent package source.
namespace <- asNamespace("MPCurver")
comparison_environment <- new.env(parent = namespace)
for (path in list.files(file.path(package_source, "R"), pattern = "\\.R$", full.names = TRUE))
  sys.source(path, envir = comparison_environment)
function_names <- ls(comparison_environment, all.names = TRUE)
function_names <- function_names[vapply(function_names, function(name)
  is.function(get(name, envir = comparison_environment, inherits = FALSE)), logical(1))]
stopifnot(all(vapply(function_names, function(name) {
  installed <- get(name, envir = namespace, inherits = FALSE)
  source <- get(name, envir = comparison_environment, inherits = FALSE)
  identical(formals(installed), formals(source)) && identical(body(installed), body(source))
}, logical(1))))
verification <- list(verified_at_utc = format(Sys.time(), tz = "UTC"),
  candidate_count = sum(vapply(result$candidates, function(x) x$status == "converged", logical(1))),
  timing_count = nrow(benchmark$timings), score_checks = score_checks,
  installed_function_count = length(function_names),
  uninterrupted_checks = result$uninterrupted_checks, provenance = current)
saveRDS(verification, file.path(experiment_dir, "verification.rds"))
writeLines(capture.output(str(verification[c("verified_at_utc", "candidate_count", "timing_count",
  "score_checks", "installed_function_count", "uninterrupted_checks")])),
  file.path(experiment_dir, "verification.txt"))
cat(sprintf("Verified %d candidate fits, %d timed pipelines, and %d installed function bodies.\n",
  verification$candidate_count, verification$timing_count, verification$installed_function_count))
