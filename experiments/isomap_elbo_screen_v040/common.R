# Public-API helpers for the fixed-case internal experiment.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(MPCurver))
stopifnot(as.character(packageVersion("MPCurver")) == "0.4.0")
experiment_dir <- "experiments/isomap_elbo_screen_v040"
input_path <- "experiments/estimate_intrinsic_m_smooth_v032/data/main_M5_S4_r001.rds"
dataset <- readRDS(input_path)
stopifnot(identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
feature_indices <- which(dataset$true_assign == "B")
observations <- dataset$X[, feature_indices, drop = FALSE]
truth <- dataset$latent_positions[, match("B", dataset$ordering_labels)]
stopifnot(identical(dim(observations), c(300L, 12L)), all(is.finite(observations)))

design <- list(dataset = dataset$row$id, num_bins = 50L,
  dense_neighbors = 5:30, sparse_neighbors = c(5L, 10L, 15L, 20L, 30L),
  screening_sweeps = 1L, tolerance = 1e-6, total_budget = 10000L,
  fit_seed = as.integer(dataset$row$seed + 1000000L),
  timing_seed = 20261002L, timing_repeats = 5L,
  control = mpcurve_control(rw_order = 2L, lambda_init = 1, ridge = 0,
    convergence = "normalized"),
  data_orientation = "samples by features", feature_standardization = FALSE)

initial_control <- function(k) {
  mpcurve_init_control(discretization = "quantile", on_failure = "error",
    method_args = list(num_neighbors = as.integer(k), seed = design$fit_seed,
      control = list(num_landmarks = nrow(observations), keep = "all")))
}

fit_candidate <- function(k, sweeps, tolerance = design$tolerance) {
  set.seed(design$fit_seed)
  fit_mpcurve(observations, num_bins = design$num_bins, intrinsic_dim = 1L,
    initial_method = "isomap", position_prior = "adaptive", max_iter = sweeps,
    tol = tolerance, init_control = initial_control(k), control = design$control)
}

continue_candidate <- function(fit) {
  # Resume the actual fitted posterior/parameters, without regenerating Isomap.
  while (!isTRUE(fit$converged) && fit$fit$iter < design$total_budget) {
    fit <- do_mpcurve(fit, max_iter = min(2000L, design$total_budget - fit$fit$iter),
      tol = design$tolerance)
  }
  if (!isTRUE(fit$converged)) stop("Candidate exhausted the convergence budget.")
  fit
}

score_candidates <- function(neighbors) {
  states <- vector("list", length(neighbors))
  scores <- rep(-Inf, length(neighbors))
  failures <- rep("", length(neighbors))
  for (index in seq_along(neighbors)) {
    states[index] <- list(tryCatch(
      fit_candidate(neighbors[index], design$screening_sweeps, tolerance = 0),
      error = function(error) { failures[index] <<- conditionMessage(error); NULL }))
    if (!is.null(states[[index]])) {
      stopifnot(states[[index]]$fit$iter == 1L,
        length(states[[index]]$elbo_trace) == 2L,
        states[[index]]$K == design$num_bins,
        states[[index]]$fit$control$method == "isomap")
      scores[index] <- tail(states[[index]]$elbo_trace, 1L)
    }
  }
  if (!any(is.finite(scores))) stop("No eligible Isomap candidate.")
  # Sorting makes exact-tie resolution explicit and independent of loop order.
  best <- order(-scores, neighbors)[1L]
  list(neighbors = neighbors, states = states, scores = scores,
    failures = failures, selected_index = best, selected_k = neighbors[best])
}

position <- function(fit) as.numeric(fitted_positions(fit)[, 1L])
recovery <- function(values) abs(cor(truth, values, method = "spearman"))
last_elbo <- function(fit) tail(fit$elbo_trace, 1L)
compact_fit <- function(fit) list(position = position(fit),
  elbo_trace = fit$elbo_trace, iter = fit$fit$iter, converged = fit$converged,
  sigma2 = fit$params$sigma2, lambda = fit$fit$lambda_vec,
  position_prior = fit$params$pi, initialization = fit$fit$init_info)

timed <- function(expression) {
  start <- proc.time()[["elapsed"]]
  value <- force(expression)
  list(value = value, seconds = proc.time()[["elapsed"]] - start)
}

package_source <- normalizePath("../MPCurver")
source_paths <- c(file.path(package_source, "DESCRIPTION"),
  file.path(package_source, "NAMESPACE"), list.files(file.path(package_source, "R"),
    pattern = "\\.R$", full.names = TRUE))
package_source_hashes <- setNames(vapply(source_paths, function(path)
  digest::digest(file = path, algo = "sha256"), character(1)),
  substring(source_paths, nchar(package_source) + 2L))
script_paths <- file.path(experiment_dir, c("common.R", "run.R", "benchmark.R",
  "report.R", "verify.R", "run_r.sh", "DESIGN.md"))
provenance <- function() list(design = design,
  input_path = input_path, input_file_sha256 = digest::digest(file = input_path, algo = "sha256"),
  full_input_sha256 = dataset$input_hash,
  group_input_sha256 = digest::digest(observations, algo = "sha256"),
  truth_sha256 = digest::digest(truth, algo = "sha256"),
  feature_indices = feature_indices, feature_names = colnames(observations),
  generating_seed = dataset$row$seed, shape_seed = dataset$shape_seed,
  package_version = as.character(packageVersion("MPCurver")),
  package_source_commit = system2("git", c("-C", package_source, "rev-parse", "HEAD"), stdout = TRUE),
  inferorder_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  package_source_hashes = package_source_hashes,
  script_hashes = setNames(vapply(script_paths, function(path)
    digest::digest(file = path, algo = "sha256"), character(1)), basename(script_paths)),
  session = sessionInfo())
