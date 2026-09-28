#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
registry <- read.csv(file.path(study_dir, "manifest.csv"))
expected <- make_manifest()
stopifnot(identical(names(registry), names(expected)), nrow(registry) == nrow(expected),
  all(vapply(names(registry), function(name)
    isTRUE(all.equal(registry[[name]], expected[[name]])), logical(1))))
active <- active_manifest()
main <- active[active$phase == "main", ]
stopifnot(nrow(main) == 90L, !anyDuplicated(main$seed),
  all(table(main$true_M, main$snr) == 10L))
write.csv(active, file.path(study_dir, "active_manifest.csv"), row.names = FALSE)
dispatch <- main
dispatch$worker <- (seq_len(nrow(main)) - 1L) %% 2L + 1L
stopifnot(all(table(dispatch$worker) == 45L),
  all(table(dispatch$worker, dispatch$true_M, dispatch$snr) == 5L))
write.csv(dispatch, file.path(study_dir, "dispatch_manifest.csv"), row.names = FALSE)
jsonlite::write_json(c(design, list(planned_repetitions = planned_repetitions,
  planned_main_datasets = 90L, design_hash = design_hash)),
  file.path(study_dir, "design.json"), pretty = TRUE, auto_unbox = TRUE)
atomic_save(design, file.path(study_dir, "design.rds"))
hashes <- lapply(seq_len(nrow(active)), function(i) {
  data <- load_dataset(active[i, ])
  if (active$phase[i] == "pilot") {
    reference <- read.csv(file.path(study_dir, "source/pilot_reference_hashes.csv"))
    reference <- reference[reference$id == active$id[i], ]
    stopifnot(nrow(reference) == 1L, data$input_hash == reference$input_hash,
      digest::digest(data$latent_positions, algo = "sha256") == reference$latent_hash)
  }
  data.frame(id = active$id[i], input_hash = data$input_hash)
})
write.csv(do.call(rbind, hashes), file.path(study_dir, "input_hashes.csv"), row.names = FALSE)
writeLines(c(paste("Package path:", find.package("MPCurver")),
  paste("Source commit:", readLines(file.path(study_dir, "source/source_commit.txt"))),
  capture.output(sessionInfo())), file.path(study_dir, "session_info.txt"))
jsonlite::write_json(list(approved_for_main = TRUE,
  authorization = "User approved two CPU threads and starting the 90-dataset study.",
  convergence = "normalized", tolerance = 1e-6, workers = 2L,
  cpu_threads_per_worker = 1L, design_hash = design_hash,
  pilot_input_reproduction = "All nine pinned pilot matrices and latent positions match exactly.",
  package_regression = "MPCurver experiments/v032_convergence_regression: 47 fits passed.",
  prior_pilots_in_main_totals = FALSE),
  file.path(study_dir, "execution_authorization.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Prepared 90 main datasets, nine input-reproduction checks, and two disjoint workers.\n")
