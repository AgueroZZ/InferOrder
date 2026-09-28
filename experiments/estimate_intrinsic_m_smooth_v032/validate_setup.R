#!/usr/bin/env Rscript
# Validate data invariants, exact pairing, and independent deterministic regeneration.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
manifest <- active_manifest()
stopifnot(identical(manifest, make_manifest()), nrow(manifest) == 90L,
  all(manifest$phase == "main"), all(table(manifest$true_M, manifest$snr) == 10L),
  !anyDuplicated(manifest$id), !anyDuplicated(manifest$seed))
archive_reference <- jsonlite::read_json(file.path(baseline_dir, "execution_checksums.json"))
stopifnot(identical(baseline_archive_sha256,
  archive_reference[["source/MPCurver_0.3.2.tar.gz"]]))

derivative_grid <- seq(0, 1, length.out = design$derivative_grid_size)
grid_angles <- outer(derivative_grid, pi * design$frequencies)
sampled_nonmonotone <- 0L
total_redraws <- 0L
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  data <- load_dataset(row)
  baseline_path <- file.path(baseline_dir, "data", paste0(row$id, ".rds"))
  baseline <- readRDS(baseline_path)
  stopifnot(identical(dim(data$X), c(design$n, design$D)), all(is.finite(data$X)),
    max(abs(colMeans(data$signal))) < 1e-10,
    max(abs(apply(data$signal, 2, stats::var) - 1)) < 1e-10,
    all(table(data$true_assign) == design$D / row$true_M),
    data$noise_variance == 1 / row$snr,
    identical(data$true_assign, baseline$true_assign),
    identical(data$latent_positions, baseline$latent_positions),
    identical(data$ordering_labels, baseline$ordering_labels),
    identical(dimnames(data$X), dimnames(baseline$X)),
    identical(data$noise, baseline$X - baseline$signal),
    isTRUE(all.equal(data$X - data$signal, data$noise, tolerance = 1e-14)),
    identical(data$baseline_hashes, baseline_input_hashes(baseline, baseline_path)),
    identical(data$baseline_archive_sha256, baseline_archive_sha256),
    identical(data$package_source_commit, baseline_source_commit),
    identical(data$input_hash, digest::digest(data$X, algo = "sha256")),
    data$shape_seed == as.integer(row$seed + design$shape_seed_offset),
    length(data$anchor_indices) == row$true_M,
    !anyDuplicated(data$anchor_indices),
    all(table(data$true_assign[data$anchor_indices]) == 1L),
    identical(data$signal[, data$anchor_indices, drop = FALSE],
      baseline$signal[, data$anchor_indices, drop = FALSE]))
  expected_anchors <- vapply(data$ordering_labels, function(label) {
    which(data$true_assign == label)[1L]
  }, integer(1))
  stopifnot(identical(data$anchor_indices, expected_anchors),
    identical(data$anchor_features, colnames(data$signal)[expected_anchors]))
  for (m in seq_along(data$ordering_labels)) {
    j <- data$anchor_indices[m]
    differences <- diff(data$signal[order(data$latent_positions[, m]), j])
    stopifnot(all(differences > 0) || all(differences < 0))
  }
  nonanchors <- which(!data$coefficients$is_anchor)
  stopifnot(length(nonanchors) == design$D - row$true_M,
    all(data$derivative_diagnostics$nonmonotone_certified[nonanchors]),
    all(data$derivative_diagnostics$draws[nonanchors] >= 1L))
  for (j in nonanchors) {
    sine <- unlist(data$coefficients[j, c("sine_k2", "sine_k3", "sine_k4")], use.names = FALSE)
    cosine <- unlist(data$coefficients[j, c("cosine_k2", "cosine_k3", "cosine_k4")], use.names = FALSE)
    derivative <- as.numeric(cos(grid_angles) %*% (sine * pi / design$frequencies) -
      sin(grid_angles) %*% (cosine * pi / design$frequencies))
    threshold <- design$derivative_relative_tolerance * max(abs(derivative))
    stopifnot(min(derivative) < -threshold, max(derivative) > threshold,
      isTRUE(all.equal(min(derivative), data$derivative_diagnostics$derivative_min[j])),
      isTRUE(all.equal(max(derivative), data$derivative_diagnostics$derivative_max[j])))
    m <- match(data$true_assign[j], data$ordering_labels)
    angles <- outer(data$latent_positions[, m], pi * design$frequencies)
    reconstructed <- as.numeric(sin(angles) %*% (sine / design$frequencies^2) +
      cos(angles) %*% (cosine / design$frequencies^2))
    reconstructed <- (reconstructed - data$coefficients$raw_center[j]) /
      data$coefficients$raw_scale[j]
    stopifnot(isTRUE(all.equal(data$signal[, j], reconstructed, check.attributes = FALSE,
      tolerance = 1e-13)))
    differences <- diff(data$signal[order(data$latent_positions[, m]), j])
    if (any(differences > 0) && any(differences < 0)) sampled_nonmonotone <- sampled_nonmonotone + 1L
  }
  # This call ignores saved new data and regenerates all coefficients from the shape seed.
  regenerated <- generate_paired_dataset(row)
  stopifnot(identical(data, regenerated))
  total_redraws <- total_redraws + sum(data$derivative_diagnostics$draws[nonanchors] - 1L)
}
stopifnot(sampled_nonmonotone == sum(design$D - manifest$true_M))

# Pair against the current posterior-occupancy baseline, including its reporting correction.
baseline_status <- jsonlite::read_json(file.path(baseline_dir, "main_summary", "status.json"))
baseline_runs <- read.csv(file.path(baseline_dir, "main_summary", "runs.csv"))
baseline_summary <- read.csv(file.path(baseline_dir, "main_summary", "summary.csv"))
stopifnot(isTRUE(baseline_status$complete), baseline_status$reported_dimension == "posterior_occupancy",
  baseline_status$completed_method_runs == 180L, nrow(baseline_runs) == 180L,
  all(baseline_runs$status == "success"),
  all(table(baseline_runs$id) == 2L), setequal(baseline_runs$id, manifest$id),
  all(baseline_runs$estimated_M == baseline_runs$effective_M),
  sum(baseline_runs$method == "adaptive" & baseline_runs$effective_M == baseline_runs$true_M) == 87L,
  sum(baseline_runs$method == "forward" & baseline_runs$effective_M == baseline_runs$true_M) == 90L,
  sum(baseline_summary$exact[baseline_summary$method == "adaptive"]) == 87L,
  sum(baseline_summary$exact[baseline_summary$method == "forward"]) == 90L)
jsonlite::write_json(list(data_generation = TRUE, datasets_checked = nrow(manifest),
  exact_latent_assignment_noise_pairing = TRUE, unchanged_monotone_anchors = TRUE,
  analytic_nonmonotonicity_certificates = TRUE, sampled_nonmonotone_features = sampled_nonmonotone,
  deterministic_regeneration = TRUE, source_archive_hash_verified = TRUE,
  total_curve_redraws = total_redraws,
  baseline_effective_m_exact = list(adaptive = 87L, forward = 90L),
  design_hash = design_hash), file.path(study_dir, "setup_validation.json"),
  auto_unbox = TRUE, pretty = TRUE)
cat("All 90 paired datasets passed generation, provenance, and regeneration checks.\n")
