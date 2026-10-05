source("experiments/isomap_kmin_m1_p4_v040/common.R")
metadata <- provenance()
metadata_path <- file.path(study, "provenance.rds")
if (file.exists(metadata_path)) {
  existing <- readRDS(metadata_path)
  stopifnot(identical(existing$design_hash, metadata$design_hash),
    identical(existing$package_source_hashes, metadata$package_source_hashes),
    identical(existing$archive_sha256, metadata$archive_sha256),
    identical(existing$script_hashes, metadata$script_hashes))
} else save_atomic(metadata, metadata_path)
manifest <- lapply(seq_len(design$replications), function(replication) {
  input <- simulate_replicate(replication)
  path <- input_path(replication)
  if (file.exists(path)) stopifnot(identical(input, readRDS(path))) else save_atomic(input, path)
  stopifnot(all(apply(input$dense_signal, 2L, is_nonmonotone)),
    abs(mean(input$dense_signal^2) - 1) < 1e-12)
  data.frame(replication = replication, generating_seed = input$seed, parent_seed = input$parent_seed,
    parent_input_sha256 = input$parent_input_sha256,
    fit_seed = design$fit_seed_base + replication, shape_attempts = input$shape_attempts,
    samples = nrow(input$X), features = ncol(input$X), noise_sd = design$noise_sd,
    average_variance_snr = 1 / design$noise_sd^2, input_sha256 = input$input_sha256)
})
write.csv(do.call(rbind, manifest), file.path(study, "manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(study, "R_session.txt"))
cat(sprintf("Prepared %d independent nonmonotone M=1, P=4 paired inputs.\n", design$replications))
