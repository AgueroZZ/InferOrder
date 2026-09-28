#!/usr/bin/env Rscript
# Prepare the fixed 90-dataset paired comparison from the InferOrder root.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
manifest <- make_manifest()
write.csv(manifest, file.path(study_dir, "manifest.csv"), row.names = FALSE)
write.csv(manifest, file.path(study_dir, "active_manifest.csv"), row.names = FALSE)
dispatch <- manifest
dispatch$task_index <- seq_len(nrow(manifest)) - 1L
stopifnot(identical(dispatch$task_index, 0:89),
  all(table(dispatch$true_M, dispatch$snr) == planned_repetitions))
write.csv(dispatch, file.path(study_dir, "dispatch_manifest.csv"), row.names = FALSE)
jsonlite::write_json(c(design, list(planned_repetitions = planned_repetitions,
  planned_main_datasets = nrow(manifest), design_hash = design_hash)),
  file.path(study_dir, "design.json"), pretty = TRUE, auto_unbox = TRUE)
atomic_save(design, file.path(study_dir, "design.rds"))
hashes <- lapply(seq_len(nrow(manifest)), function(i) {
  data <- load_dataset(manifest[i, , drop = FALSE])
  data.frame(id = manifest$id[i], input_hash = data$input_hash,
    baseline_input_hash = data$baseline_hashes$input_hash,
    baseline_signal_hash = data$baseline_hashes$signal_hash,
    baseline_file_sha256 = data$baseline_hashes$file_sha256,
    latent_hash = data$baseline_hashes$latent_hash,
    noise_hash = data$baseline_hashes$residual_noise_hash,
    shape_seed = data$shape_seed)
})
write.csv(do.call(rbind, hashes), file.path(study_dir, "input_hashes.csv"), row.names = FALSE)
jsonlite::write_json(list(package_version = design$package_version,
  package_archive_sha256 = baseline_archive_sha256,
  source_commit = baseline_source_commit,
  baseline_design_hash = design$baseline_design_hash,
  baseline_manifest_sha256 = digest::digest(
    file = file.path(baseline_dir, "active_manifest.csv"), algo = "sha256"),
  new_design_hash = design_hash),
  file.path(study_dir, "source_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(c(paste("Package path:", find.package("MPCurver")),
  paste("Source commit:", baseline_source_commit),
  paste("Source archive SHA-256:", baseline_archive_sha256),
  capture.output(sessionInfo())), file.path(study_dir, "session_info.txt"))
cat("Prepared 90 paired datasets with one retained monotone feature per ordering.\n")
