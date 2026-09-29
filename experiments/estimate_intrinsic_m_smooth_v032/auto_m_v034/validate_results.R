#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R")

validation_rows <- vector("list", nrow(manifest))
provenance_reference <- execution_provenance()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  path <- file.path(extension_dir, "results",
    paste0(row$id, "_auto_adaptive.rds"))
  stopifnot(file.exists(path))
  result <- readRDS(path)
  dataset <- load_fixed_dataset(row)
  compact <- result$compact
  stopifnot(
    result$design_hash == design_hash,
    result$input_hash == dataset$input_hash,
    identical(result$provenance, provenance_reference),
    result$row$status == "success",
    compact$converged,
    result$row$warning_count == 0L,
    result$row$selected_initial_M == compact$selected_initial_M,
    result$row$effective_M == compact$effective_M,
    result$row$minimum_initial_cluster_size >= design$similarity_min_cluster_size,
    compact$dimension_initialization$selected_M == compact$M,
    min(compact$dimension_initialization$cluster_sizes) >=
      design$similarity_min_cluster_size,
    max(abs(rowSums(compact$weights) - 1)) < 1e-10,
    all(is.finite(compact$positions)),
    all(is.finite(compact$sigma2)), all(compact$sigma2 > 0)
  )
  validation_rows[[i]] <- data.frame(
    id = row$id,
    true_M = row$true_M,
    snr = row$snr,
    replicate = row$replicate,
    input_hash = result$input_hash,
    selected_initial_M = compact$selected_initial_M,
    effective_M = compact$effective_M,
    exact_initial_M = compact$selected_initial_M == row$true_M,
    exact_effective_M = compact$effective_M == row$true_M,
    initial_partition_ARI = compact$initial_partition_ARI,
    final_ARI = result$evaluation$ARI,
    ordering_recovery = result$evaluation$ordering_recovery,
    minimum_initial_cluster_size =
      min(compact$dimension_initialization$cluster_sizes),
    converged = compact$converged,
    iterations = compact$iterations,
    warning_count = result$row$warning_count,
    elapsed_seconds = result$row$elapsed_seconds,
    M_at_1e_12 = compact$threshold_counts[1],
    M_at_1e_9 = compact$threshold_counts[2],
    M_at_1e_6 = compact$threshold_counts[3],
    M_at_1e_3 = compact$threshold_counts[4]
  )
}

validation <- do.call(rbind, validation_rows)
stopifnot(nrow(validation) == 90L, !anyDuplicated(validation$id),
  all(validation$converged), sum(validation$warning_count) == 0L,
  all(validation$minimum_initial_cluster_size >=
    design$similarity_min_cluster_size))
write.csv(validation, file.path(extension_dir, "candidate_validation.csv"),
  row.names = FALSE)

replicate_one <- manifest[manifest$replicate == 1L, , drop = FALSE]
for (id in replicate_one$id) {
  state <- readRDS(file.path(extension_dir, "full_fits", paste0(id, ".rds")))
  init_info <- state$fit$fit$init_info
  stopifnot(length(init_info) == state$fit$intrinsic_dim,
    all(vapply(init_info, function(x) x$method_requested == "isomap", logical(1))),
    all(vapply(init_info, function(x) x$method_used == "isomap", logical(1))),
    all(vapply(init_info, function(x) is.na(x$pca_component_used), logical(1))))
}

record <- list(
  complete = TRUE,
  design_hash = design_hash,
  package_version = design$package_version,
  package_source_commit = design$package_source_commit,
  package_archive_sha256 = design$package_archive_sha256,
  provenance_fingerprint = provenance_reference$fingerprint,
  datasets = nrow(validation),
  converged = sum(validation$converged),
  warnings = sum(validation$warning_count),
  independent_isomap_initialization_checks = nrow(replicate_one),
  exact_initial_M = sum(validation$exact_initial_M),
  exact_effective_M = sum(validation$exact_effective_M),
  under_initial_M = sum(validation$selected_initial_M < validation$true_M),
  over_initial_M = sum(validation$selected_initial_M > validation$true_M),
  under_effective_M = sum(validation$effective_M < validation$true_M),
  over_effective_M = sum(validation$effective_M > validation$true_M),
  mean_initial_partition_ARI = mean(validation$initial_partition_ARI),
  mean_final_ARI = mean(validation$final_ARI),
  mean_ordering_recovery = mean(validation$ordering_recovery),
  total_seconds = sum(validation$elapsed_seconds),
  threshold_stability = list(
    M_1e_12_equals_1e_9 = sum(validation$M_at_1e_12 == validation$M_at_1e_9),
    M_1e_12_equals_1e_6 = sum(validation$M_at_1e_12 == validation$M_at_1e_6),
    M_1e_12_equals_1e_3 = sum(validation$M_at_1e_12 == validation$M_at_1e_3)
  )
)
jsonlite::write_json(record, file.path(extension_dir, "validation.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(record)
