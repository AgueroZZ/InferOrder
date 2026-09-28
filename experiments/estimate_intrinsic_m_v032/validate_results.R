#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
manifest <- active_manifest()
manifest <- manifest[manifest$phase == "main", ]
candidate_count <- 0L
failures <- 0L
for (i in seq_len(nrow(manifest))) {
  data <- load_dataset(manifest[i, ])
  for (method in c("adaptive", "forward")) {
    result <- readRDS(file.path(study_dir, "results",
      paste0(manifest$id[i], "_", method, ".rds")))
    stopifnot(result$design_hash == design_hash, result$input_hash == data$input_hash,
      result$row$status %in% c("success", "error", "nonconverged"))
    failures <- failures + as.integer(result$row$status != "success")
    for (key in result$candidate_keys) {
      candidate <- readRDS(file.path(study_dir, "candidates", paste0(key, ".rds")))
      stopifnot(candidate$design_hash == design_hash, candidate$input_hash == data$input_hash)
      candidate_count <- candidate_count + 1L
      if (candidate$status == "success") {
        fit <- candidate$fit
        stopifnot(fit$converged, fit$convergence == "normalized", fit$tolerance == 1e-6,
          tail(fit$temperature, 1) == 1, fit$last_normalized_increment < 1e-6)
      }
    }
    if (method == "forward" && nrow(result$history)) {
      stopifnot(all(result$history$accepted == (result$history$delta > 0)))
    }
  }
}
jsonlite::write_json(list(validated = TRUE, datasets = nrow(manifest),
  method_outcomes = 2L * nrow(manifest), candidates = candidate_count,
  unresolved_outcomes = failures, design_hash = design_hash),
  file.path(study_dir, "main_validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Validated all 180 method outcomes; unresolved outcomes:", failures, "\n")
