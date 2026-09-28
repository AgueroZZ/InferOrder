# Shared settings, paired smooth-trajectory generation, and checkpoint utilities.
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

if (!dir.exists("experiments/estimate_intrinsic_m_v032") || !dir.exists("analysis")) {
  stop("Run this study from the InferOrder repository root.")
}
study_dir <- normalizePath("experiments/estimate_intrinsic_m_smooth_v032", mustWork = FALSE)
baseline_dir <- normalizePath("experiments/estimate_intrinsic_m_v032")
.libPaths(c(file.path(study_dir, "library"), .libPaths()))
stopifnot(as.character(utils::packageVersion("MPCurver")) == "0.3.2")

baseline_archive <- file.path(baseline_dir, "source", "MPCurver_0.3.2.tar.gz")
baseline_archive_sha256 <- digest::digest(file = baseline_archive, algo = "sha256")
baseline_source_commit <- trimws(readLines(file.path(baseline_dir, "source", "source_commit.txt")))
baseline_design <- jsonlite::read_json(file.path(baseline_dir, "design.json"), simplifyVector = TRUE)
planned_repetitions <- 10L
design <- list(n = 300L, D = 60L, true_M = 3:5, snr = c(1, 4, 16),
  repetitions = planned_repetitions, max_M = 8L, K = 50L, tolerance = 1e-6,
  convergence = "normalized", initial_budget = 1500L, continuation_budget = 1500L,
  maximum_sweeps = 10000L, effective_weight_tol = 1e-12,
  method = "isomap", similarity_metric = "spearman", cluster_linkage = "single",
  rw_q = 2L, ridge = 0, T_start = 5, T_end = 1, n_outer = 25L, inner_iter = 1L,
  position_prior = "adaptive", generator = "paired_one_monotone_fourier_v1",
  family = "random_sine_cosine_mixture", frequencies = 2:4,
  angular_frequency = "pi * k", spectral_attenuation = "1 / k^2",
  coefficient_distribution = "independent standard normal",
  monotone_features_per_ordering = 1L,
  anchor_selection = "first current feature column within each true ordering",
  derivative_grid_size = 1001L, derivative_relative_tolerance = 1e-6,
  maximum_curve_draws = 10000L, shape_seed_offset = 2000000L,
  rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
  pairing = "same latent positions, assignments, feature order, monotone anchor, and residual noise",
  baseline_study = "estimate_intrinsic_m_v032",
  baseline_design_hash = baseline_design$design_hash,
  package_version = "0.3.2", package_source_commit = baseline_source_commit,
  package_archive_sha256 = baseline_archive_sha256)
design_hash <- digest::digest(design, algo = "sha256")

for (directory in c("data", "results", "candidates", "checkpoints", "logs", "full_fits", "figures")) {
  dir.create(file.path(study_dir, directory), showWarnings = FALSE, recursive = TRUE)
}

atomic_save <- function(object, path, compress = FALSE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  stopifnot(file.rename(temporary, path))
}

make_manifest <- function() {
  registry <- read.csv(file.path(baseline_dir, "active_manifest.csv"))
  main <- registry[registry$phase == "main" & registry$replicate <= planned_repetitions, , drop = FALSE]
  rownames(main) <- NULL
  stopifnot(nrow(main) == 90L, !anyDuplicated(main$id), !anyDuplicated(main$seed),
    all(table(main$true_M, main$snr) == planned_repetitions))
  main
}

active_manifest <- function() {
  path <- file.path(study_dir, "manifest.csv")
  if (file.exists(path)) read.csv(path) else make_manifest()
}

# Two signs of the analytic derivative certify a nonmonotone smooth curve.
draw_smooth_curve <- function(t) {
  frequencies <- design$frequencies
  attenuation <- 1 / frequencies^2
  derivative_grid <- seq(0, 1, length.out = design$derivative_grid_size)
  sample_angles <- outer(t, pi * frequencies)
  grid_angles <- outer(derivative_grid, pi * frequencies)
  for (draw in seq_len(design$maximum_curve_draws)) {
    sine_coefficients <- stats::rnorm(length(frequencies))
    cosine_coefficients <- stats::rnorm(length(frequencies))
    derivative <- as.numeric(
      cos(grid_angles) %*% (sine_coefficients * attenuation * pi * frequencies) -
      sin(grid_angles) %*% (cosine_coefficients * attenuation * pi * frequencies))
    derivative_scale <- max(abs(derivative))
    threshold <- design$derivative_relative_tolerance * derivative_scale
    certified <- min(derivative) < -threshold && max(derivative) > threshold
    if (!certified) next
    raw <- as.numeric(sin(sample_angles) %*% (sine_coefficients * attenuation) +
      cos(sample_angles) %*% (cosine_coefficients * attenuation))
    raw_center <- mean(raw)
    raw_scale <- stats::sd(raw)
    if (!is.finite(raw_scale) || raw_scale <= 0) next
    return(list(signal = (raw - raw_center) / raw_scale,
      sine_coefficients = sine_coefficients, cosine_coefficients = cosine_coefficients,
      raw_center = raw_center, raw_scale = raw_scale, draws = draw,
      derivative_min = min(derivative), derivative_max = max(derivative),
      derivative_scale = derivative_scale, derivative_threshold = threshold))
  }
  stop("Could not generate a certified nonmonotone trajectory within the declared draw limit.")
}

baseline_input_hashes <- function(baseline, path) {
  list(file_sha256 = digest::digest(file = path, algo = "sha256"),
    input_hash = digest::digest(baseline$X, algo = "sha256"),
    signal_hash = digest::digest(baseline$signal, algo = "sha256"),
    latent_hash = digest::digest(baseline$latent_positions, algo = "sha256"),
    assignment_hash = digest::digest(baseline$true_assign, algo = "sha256"),
    residual_noise_hash = digest::digest(baseline$X - baseline$signal, algo = "sha256"))
}

generate_paired_dataset <- function(row) {
  baseline_path <- file.path(baseline_dir, "data", paste0(row$id, ".rds"))
  baseline <- readRDS(baseline_path)
  stopifnot(baseline$row$seed == row$seed, baseline$row$id == row$id,
    baseline$design_hash == design$baseline_design_hash,
    identical(baseline$input_hash, digest::digest(baseline$X, algo = "sha256")))
  signal <- baseline$signal
  noise <- baseline$X - baseline$signal
  anchor_indices <- vapply(baseline$ordering_labels, function(label) {
    which(baseline$true_assign == label)[1L]
  }, integer(1))
  shape_seed <- as.integer(row$seed + design$shape_seed_offset)
  set.seed(shape_seed, kind = design$rng_kind[1L],
    normal.kind = design$rng_kind[2L], sample.kind = design$rng_kind[3L])
  coefficients <- data.frame(feature_index = seq_len(design$D),
    feature = colnames(signal), ordering = baseline$true_assign,
    is_anchor = seq_len(design$D) %in% anchor_indices,
    sine_k2 = NA_real_, sine_k3 = NA_real_, sine_k4 = NA_real_,
    cosine_k2 = NA_real_, cosine_k3 = NA_real_, cosine_k4 = NA_real_,
    raw_center = NA_real_, raw_scale = NA_real_)
  diagnostics <- coefficients[, c("feature_index", "feature", "ordering", "is_anchor")]
  diagnostics$draws <- 0L
  diagnostics$derivative_min <- NA_real_
  diagnostics$derivative_max <- NA_real_
  diagnostics$derivative_scale <- NA_real_
  diagnostics$derivative_threshold <- NA_real_
  diagnostics$nonmonotone_certified <- FALSE
  for (j in which(!coefficients$is_anchor)) {
    ordering <- match(baseline$true_assign[j], baseline$ordering_labels)
    curve <- draw_smooth_curve(baseline$latent_positions[, ordering])
    signal[, j] <- curve$signal
    coefficients[j, c("sine_k2", "sine_k3", "sine_k4")] <- curve$sine_coefficients
    coefficients[j, c("cosine_k2", "cosine_k3", "cosine_k4")] <- curve$cosine_coefficients
    coefficients$raw_center[j] <- curve$raw_center
    coefficients$raw_scale[j] <- curve$raw_scale
    diagnostics$draws[j] <- curve$draws
    diagnostics$derivative_min[j] <- curve$derivative_min
    diagnostics$derivative_max[j] <- curve$derivative_max
    diagnostics$derivative_scale[j] <- curve$derivative_scale
    diagnostics$derivative_threshold[j] <- curve$derivative_threshold
    diagnostics$nonmonotone_certified[j] <- TRUE
  }
  X <- signal + noise
  list(row = row, design_hash = design_hash, X = X, signal = signal, noise = noise,
    true_assign = baseline$true_assign, latent_positions = baseline$latent_positions,
    ordering_labels = baseline$ordering_labels, noise_variance = baseline$noise_variance,
    realized_noise_variance = apply(noise, 2, stats::var),
    input_hash = digest::digest(X, algo = "sha256"), shape_seed = shape_seed,
    anchor_indices = anchor_indices, anchor_features = colnames(signal)[anchor_indices],
    coefficients = coefficients, derivative_diagnostics = diagnostics,
    baseline_hashes = baseline_input_hashes(baseline, baseline_path),
    baseline_archive_sha256 = baseline_archive_sha256,
    package_source_commit = baseline_source_commit)
}

load_dataset <- function(row) {
  path <- file.path(study_dir, "data", paste0(row$id, ".rds"))
  if (file.exists(path)) {
    value <- readRDS(path)
    stopifnot(identical(value$design_hash, design_hash), value$row$seed == row$seed,
      identical(value$input_hash, digest::digest(value$X, algo = "sha256")))
    return(value)
  }
  value <- generate_paired_dataset(row)
  atomic_save(value, path)
  value
}

# Keep fitting options identical to the monotone study.
fit_options <- function() list(algorithm = "cavi", method = design$method, K = design$K,
  rw_q = design$rw_q, ridge = design$ridge, lambda = 1, fix_lambda = FALSE, S = NULL,
  position_prior = design$position_prior, similarity_metric = design$similarity_metric,
  cluster_linkage = design$cluster_linkage, discretization = "quantile", num_cores = 1L,
  iter = design$initial_budget, tol = design$tolerance, convergence = design$convergence,
  T_start = design$T_start, T_end = design$T_end, n_outer = design$n_outer,
  inner_iter = design$inner_iter, max_converge_iter = design$initial_budget,
  tol_outer = design$tolerance, verbose = FALSE)
