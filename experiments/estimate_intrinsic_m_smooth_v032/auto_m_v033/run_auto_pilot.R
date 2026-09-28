#!/usr/bin/env Rscript

# Pilot the package-level automatic ordering-count initialization on fixed
# datasets from the one-monotone-anchor study.
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

study_dir <- normalizePath("experiments/estimate_intrinsic_m_smooth_v032")
extension_dir <- file.path(study_dir, "auto_m_v033")
pilot_dir <- file.path(extension_dir, "pilot")
dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)

candidate_library <- Sys.getenv("AUTO_MPCURVE_LIBRARY",
  file.path(study_dir, "library_auto_candidate"))
dependency_library <- file.path(study_dir, "library")
.libPaths(c(candidate_library, dependency_library, .libPaths()))

suppressPackageStartupMessages({
  library(MPCurver)
  library(clue)
  library(digest)
  library(mclust)
})

manifest <- read.csv(file.path(study_dir, "active_manifest.csv"))
default_ids <- sprintf("main_M%d_S%d_r001", rep(3:5, each = 3L), rep(c(1, 4, 16), 3L))
arguments <- commandArgs(trailingOnly = TRUE)
dataset_ids <- if (length(arguments)) arguments else default_ids
stopifnot(!anyDuplicated(dataset_ids), all(dataset_ids %in% manifest$id))

maximum_sweeps <- 10000L
block_sweeps <- 1500L
effective_weight_tol <- 1e-12

objective_trace <- function(fit) {
  if (fit$intrinsic_dim == 1L) fit$fit$elbo_trace else fit$fit$objective_history
}

fit_one <- function(dataset) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  set.seed(dataset$row$seed + 1000000L)
  fit <- withCallingHandlers(MPCurver::fit_mpcurve(
    X = dataset$X,
    intrinsic_dim = "auto",
    max_intrinsic_dim = 8L,
    similarity_metric = "spline_r2",
    spline_r2_df = 5L,
    similarity_min_cluster_size = 2L,
    partition_init = "similarity",
    partition_prior = "adaptive",
    effective_weight_tol = effective_weight_tol,
    algorithm = "cavi",
    method = "isomap",
    K = 50L,
    rw_q = 2L,
    ridge = 0,
    lambda = 1,
    fix_lambda = FALSE,
    S = NULL,
    position_prior = "adaptive",
    cluster_linkage = "single",
    discretization = "quantile",
    num_cores = 1L,
    iter = block_sweeps,
    tol = 1e-6,
    convergence = "normalized",
    T_start = 5,
    T_end = 1,
    n_outer = 25L,
    inner_iter = 1L,
    max_converge_iter = block_sweeps,
    tol_outer = 1e-6,
    verbose = FALSE
  ), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  while (!isTRUE(fit$fit$converged) && fit$fit$iter < maximum_sweeps - 1L) {
    remaining <- maximum_sweeps - fit$fit$iter - 1L
    requested <- min(block_sweeps, remaining)
    fit <- withCallingHandlers(MPCurver::do_mpcurve(
      fit,
      iter = requested,
      tol = 1e-6,
      tol_outer = 1e-6,
      convergence = "normalized",
      verbose = FALSE
    ), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - started
  initialization <- fit$dimension_initialization
  selected_M <- as.integer(initialization$selected_M)
  stopifnot(isTRUE(initialization$automatic), selected_M == fit$intrinsic_dim,
    initialization$similarity_metric == "spline_r2",
    initialization$min_cluster_size == 2L,
    initialization$cluster_linkage == "single")
  if (selected_M >= 2L) {
    stopifnot(min(initialization$cluster_sizes) >= 2L,
      max(abs(initialization$initial_partition_probabilities -
        initialization$cluster_sizes / ncol(dataset$X))) < 1e-12)
  }

  if (selected_M == 1L) {
    weights <- matrix(1, nrow = ncol(dataset$X), ncol = 1L)
    positions <- matrix(fit$locations$mean$pseudotime, ncol = 1L)
  } else {
    weights <- fit$partition$pi_weights
    positions <- vapply(fit$locations, function(x) x$mean$pseudotime,
      numeric(nrow(dataset$X)))
  }
  stopifnot(all(is.finite(weights)), all(weights >= 0),
    max(abs(rowSums(weights) - 1)) < 1e-10,
    all(is.finite(positions)))
  occupancy <- colMeans(weights)
  active <- which(occupancy > effective_weight_tol)
  assignments <- max.col(weights, ties.method = "first")

  truth <- dataset$latent_positions
  active_positions <- positions[, active, drop = FALSE]
  correlation <- abs(suppressWarnings(stats::cor(truth, active_positions,
    method = "spearman")))
  correlation[!is.finite(correlation)] <- 0
  assignment_size <- max(ncol(truth), ncol(active_positions))
  matching_score <- matrix(0, assignment_size, assignment_size)
  matching_score[seq_len(ncol(truth)), seq_len(ncol(active_positions))] <- correlation
  matching <- as.integer(clue::solve_LSAP(matching_score, maximum = TRUE))[
    seq_len(ncol(truth))]
  matched <- matching <= ncol(active_positions)
  matched_correlations <- numeric(ncol(truth))
  matched_correlations[matched] <- correlation[cbind(which(matched), matching[matched])]

  trace <- objective_trace(fit)
  row <- data.frame(
    id = dataset$row$id,
    true_M = dataset$row$true_M,
    snr = dataset$row$snr,
    replicate = dataset$row$replicate,
    selected_initial_M = selected_M,
    effective_M = length(active),
    exact_initial_M = selected_M == dataset$row$true_M,
    exact_effective_M = length(active) == dataset$row$true_M,
    ARI = mclust::adjustedRandIndex(dataset$true_assign, assignments),
    ordering_recovery = mean(matched_correlations),
    true_ordering_coverage = mean(matched),
    converged = isTRUE(fit$fit$converged),
    iterations = as.integer(fit$fit$iter),
    objective = tail(trace, 1L),
    elapsed_seconds = elapsed_seconds,
    warning_count = length(unique(warnings)),
    minimum_initial_cluster_size = min(initialization$cluster_sizes),
    stringsAsFactors = FALSE
  )
  compact <- list(
    row = row,
    package_version = as.character(utils::packageVersion("MPCurver")),
    input_hash = digest::digest(dataset$X, algo = "sha256"),
    initialization = initialization,
    occupancy = occupancy,
    assignments = assignments,
    positions = positions,
    warnings = unique(warnings)
  )
  saveRDS(compact, file.path(pilot_dir, paste0(dataset$row$id, ".rds")),
    compress = FALSE)
  row
}

rows <- vector("list", length(dataset_ids))
for (i in seq_along(dataset_ids)) {
  id <- dataset_ids[i]
  manifest_row <- manifest[manifest$id == id, , drop = FALSE]
  dataset <- readRDS(file.path(study_dir, "data", paste0(id, ".rds")))
  stopifnot(nrow(manifest_row) == 1L, dataset$row$id == id,
    dataset$row$seed == manifest_row$seed,
    identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
  cat(sprintf("%s: starting %s\n", Sys.time(), id))
  rows[[i]] <- fit_one(dataset)
  print(rows[[i]], row.names = FALSE)
  flush.console()
}

results <- do.call(rbind, rows)
write.csv(results, file.path(pilot_dir, "runs.csv"), row.names = FALSE)
summary <- data.frame(
  package_version = as.character(utils::packageVersion("MPCurver")),
  datasets = nrow(results),
  initial_M_exact = sum(results$exact_initial_M),
  effective_M_exact = sum(results$exact_effective_M),
  converged = sum(results$converged),
  warnings = sum(results$warning_count),
  mean_ARI = mean(results$ARI),
  mean_ordering_recovery = mean(results$ordering_recovery),
  total_seconds = sum(results$elapsed_seconds)
)
write.csv(summary, file.path(pilot_dir, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
