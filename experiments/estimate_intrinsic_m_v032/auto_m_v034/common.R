# Shared controls and validation helpers for the MPCurver 0.3.4 automatic-M
# extension of the all-monotone study.
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

if (!dir.exists("experiments/estimate_intrinsic_m_v032") ||
    !dir.exists("analysis")) {
  stop("Run this study from the InferOrder repository root.")
}

study_dir <- normalizePath("experiments/estimate_intrinsic_m_v032")
extension_dir <- file.path(study_dir, "auto_m_v034")
dependency_library <- file.path(study_dir, "library")
shared_dependency_library <- normalizePath(
  "experiments/estimate_intrinsic_m_smooth_v032/library"
)
extension_library <- file.path(extension_dir, "library")
.libPaths(c(extension_library, dependency_library, shared_dependency_library,
  .libPaths()))

stopifnot(as.character(utils::packageVersion("MPCurver")) == "0.3.4")
suppressPackageStartupMessages({
  library(clue)
  library(digest)
  library(mclust)
})

design <- list(
  method = "auto_adaptive",
  n = 300L,
  D = 60L,
  true_M = 3:5,
  snr = c(1, 4, 16),
  repetitions = 10L,
  max_intrinsic_dim = 8L,
  similarity_metric = "spline_r2",
  spline_r2_df = 5L,
  cluster_linkage = "single",
  similarity_min_cluster_size = 2L,
  K = 50L,
  rw_q = 2L,
  ridge = 0,
  tolerance = 1e-6,
  convergence = "normalized",
  initial_budget = 1500L,
  continuation_budget = 1500L,
  maximum_sweeps = 10000L,
  effective_weight_tol = 1e-12,
  position_prior = "adaptive",
  partition_prior = "adaptive",
  discretization = "quantile",
  T_start = 5,
  T_end = 1,
  n_outer = 25L,
  inner_iter = 1L,
  seed_offset = 1000000L,
  input_study = "estimate_intrinsic_m_v032",
  package_version = "0.3.4",
  package_source_commit = "15f2b0bbe5dfa61cd46da5160b2bc251e75a0475",
  package_archive_sha256 =
    "58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6"
)
design_hash <- digest::digest(design, algo = "sha256")

for (directory in c("results", "checkpoints", "full_fits", "logs", "slurm", "summary")) {
  dir.create(file.path(extension_dir, directory), recursive = TRUE,
    showWarnings = FALSE)
}

manifest <- read.csv(file.path(study_dir, "active_manifest.csv"))
manifest <- manifest[manifest$phase == "main", , drop = FALSE]
stopifnot(nrow(manifest) == 90L, !anyDuplicated(manifest$id),
  all(table(manifest$true_M, manifest$snr) == design$repetitions))

atomic_save <- function(object, path, compress = FALSE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  stopifnot(file.rename(temporary, path))
}

load_fixed_dataset <- function(row) {
  path <- file.path(study_dir, "data", paste0(row$id, ".rds"))
  dataset <- readRDS(path)
  stopifnot(dataset$row$id == row$id, dataset$row$seed == row$seed,
    nrow(dataset$X) == design$n, ncol(dataset$X) == design$D,
    identical(dataset$input_hash, digest::digest(dataset$X, algo = "sha256")))
  dataset
}

execution_provenance <- function() {
  script_names <- c("common.R", "run_dataset.R")
  script_paths <- file.path(extension_dir, script_names)
  archive <- "experiments/mpcurve_v034_site/source/MPCurver_0.3.4.tar.gz"
  source_commit <- design$package_source_commit
  hashes <- setNames(vapply(script_paths, digest::digest, character(1),
    algo = "sha256", file = TRUE), script_names)
  value <- list(
    code_hashes = hashes,
    package_archive_hash = digest::digest(archive, algo = "sha256", file = TRUE),
    package_source_commit = source_commit,
    package_version = as.character(utils::packageVersion("MPCurver")),
    R_version = as.character(getRversion()),
    design_hash = design_hash
  )
  stopifnot(value$package_archive_hash == design$package_archive_sha256,
    value$package_source_commit == design$package_source_commit,
    value$package_version == design$package_version)
  value$fingerprint <- digest::digest(value, algo = "sha256")
  value
}

objective_nondecrease_guard <- function(previous, M, stopping = FALSE) {
  if (M == 1L) {
    rep(if (stopping) 0 else 1e-8, length(previous))
  } else {
    1e-8 * (abs(previous) + 1)
  }
}

fit_trace_diagnostics <- function(fit) {
  M <- fit$intrinsic_dim
  raw <- fit$fit
  trace <- if (M == 1L) raw$elbo_trace else raw$objective_history
  temperature <- if (M == 1L) rep(1, length(trace)) else raw$temperature_history
  temperature_one <- trace[temperature == 1]
  increments <- diff(temperature_one)
  previous <- head(temperature_one, -1L)
  valid_trace <- length(trace) == length(temperature) && length(temperature_one) > 1L &&
    all(is.finite(trace)) && all(is.finite(temperature))
  nondecrease <- valid_trace &&
    all(increments >= -objective_nondecrease_guard(previous, M))
  final_increment <- if (length(increments)) tail(increments, 1L) else NA_real_
  final_previous <- if (length(previous)) tail(previous, 1L) else NA_real_
  normalized_increment <- abs(final_increment) / (design$n * design$D)
  stopping_passed <- valid_trace && tail(temperature, 1L) == 1 &&
    final_increment >= -objective_nondecrease_guard(final_previous, M, stopping = TRUE) &&
    normalized_increment < design$tolerance
  list(
    objective = as.numeric(trace),
    temperature = as.numeric(temperature),
    iterations = as.integer(raw$iter),
    converged = isTRUE(raw$converged),
    score = tail(trace, 1L),
    valid_trace = valid_trace,
    T1_nondecrease_passed = nondecrease,
    stopping_passed = stopping_passed,
    last_normalized_increment = normalized_increment,
    min_T1_increment = if (length(increments)) min(increments) else NA_real_
  )
}

compact_auto_fit <- function(fit, dataset) {
  M <- as.integer(fit$intrinsic_dim)
  initialization <- fit$dimension_initialization
  diagnostics <- fit_trace_diagnostics(fit)
  stopifnot(isTRUE(initialization$automatic), initialization$selected_M == M,
    initialization$similarity_metric == design$similarity_metric,
    initialization$spline_r2_df == design$spline_r2_df,
    initialization$cluster_linkage == design$cluster_linkage,
    initialization$min_cluster_size == design$similarity_min_cluster_size,
    diagnostics$valid_trace, diagnostics$T1_nondecrease_passed,
    fit$fit$control$convergence == design$convergence)
  if (diagnostics$converged) stopifnot(diagnostics$stopping_passed)

  if (M == 1L) {
    weights <- matrix(1, nrow = design$D, ncol = 1L,
      dimnames = list(colnames(dataset$X), "A"))
    positions <- matrix(fit$locations$mean$pseudotime, ncol = 1L,
      dimnames = list(rownames(dataset$X), "A"))
    omega <- 1
  } else {
    weights <- fit$partition$pi_weights
    positions <- vapply(fit$locations, function(x) x$mean$pseudotime,
      numeric(design$n))
    omega <- as.numeric(MPCurver::fitted_prior(fit, type = "partition")$omega)
    stopifnot(max(abs(omega - colMeans(weights))) < 1e-10)
  }
  stopifnot(all(is.finite(weights)), all(weights >= 0),
    max(abs(rowSums(weights) - 1)) < 1e-10,
    all(is.finite(positions)), all(is.finite(fit$params$sigma2)),
    all(fit$params$sigma2 > 0))

  if (M >= 2L) {
    expected_initial <- initialization$cluster_sizes / design$D
    stopifnot(min(initialization$cluster_sizes) >=
        design$similarity_min_cluster_size,
      max(abs(unname(initialization$initial_partition_probabilities) -
        expected_initial)) < 1e-12)
  }

  occupancy <- colMeans(weights)
  active <- which(occupancy > design$effective_weight_tol)
  assignments <- max.col(weights, ties.method = "first")
  stopifnot(length(active) > 0L, all(assignments %in% active))
  active_positions <- positions[, active, drop = FALSE]
  position_sd <- apply(active_positions, 2, stats::sd)
  correlations <- suppressWarnings(stats::cor(active_positions, method = "spearman"))
  off_diagonal <- abs(correlations[upper.tri(correlations)])

  list(
    M = M,
    selected_initial_M = M,
    effective_M = as.integer(length(active)),
    omega = omega,
    occupancy = occupancy,
    active = active,
    assignments = assignments,
    weights = weights,
    positions = positions,
    score = diagnostics$score,
    converged = diagnostics$converged,
    iterations = diagnostics$iterations,
    objective = diagnostics$objective,
    temperature = diagnostics$temperature,
    last_normalized_increment = diagnostics$last_normalized_increment,
    min_T1_increment = diagnostics$min_T1_increment,
    threshold_counts = vapply(c(1e-12, 1e-9, 1e-6, 1e-3),
      function(tol) as.integer(sum(occupancy > tol)), integer(1)),
    constant_orderings = as.integer(sum(position_sd < 1e-10)),
    near_duplicate_orderings = any(off_diagonal > 0.995, na.rm = TRUE),
    dimension_initialization = initialization,
    initial_partition_ARI = mclust::adjustedRandIndex(
      dataset$true_assign, initialization$feature_cluster),
    sigma2 = fit$params$sigma2
  )
}

evaluate_compact_fit <- function(compact, dataset) {
  positions <- compact$positions[, compact$active, drop = FALSE]
  correlation <- abs(suppressWarnings(stats::cor(dataset$latent_positions,
    positions, method = "spearman")))
  correlation[!is.finite(correlation)] <- 0
  size <- max(ncol(dataset$latent_positions), ncol(positions))
  score <- matrix(0, size, size)
  score[seq_len(ncol(dataset$latent_positions)), seq_len(ncol(positions))] <-
    correlation
  matching <- as.integer(clue::solve_LSAP(score, maximum = TRUE))[
    seq_len(ncol(dataset$latent_positions))]
  matched <- matching <= ncol(positions)
  rho <- numeric(ncol(dataset$latent_positions))
  rho[matched] <- correlation[cbind(which(matched), matching[matched])]
  list(
    ARI = mclust::adjustedRandIndex(dataset$true_assign, compact$assignments),
    ordering_recovery = mean(rho),
    matched_ordering_recovery = if (any(matched)) mean(rho[matched]) else NA_real_,
    true_ordering_coverage = mean(matched)
  )
}
