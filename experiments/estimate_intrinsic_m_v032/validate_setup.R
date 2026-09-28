#!/usr/bin/env Rscript
# Scientific checks for truth generation and label/reversal-invariant scoring.
source("experiments/estimate_intrinsic_m_v032/common.R")
source(file.path(study_dir, "evaluate.R"))
manifest <- read.csv(file.path(study_dir, "manifest.csv"))
stopifnot(sum(manifest$phase == "main") == 450L,
          sum(manifest$phase == "pilot") == 9L, !anyDuplicated(manifest$seed))
for (i in which(manifest$phase == "pilot")) {
  data <- load_dataset(manifest[i, ])
  stopifnot(all(dim(data$X) == c(300L, 60L)),
    max(abs(colMeans(data$signal))) < 1e-10,
    max(abs(apply(data$signal, 2, var) - 1)) < 1e-10,
    data$noise_variance == 1 / manifest$snr[i],
    ncol(data$latent_positions) == manifest$true_M[i],
    all(data$latent_positions >= 0 & data$latent_positions <= 1),
    identical(data$input_hash, digest::digest(data$X, algo = "sha256")))
}
truth <- cbind(seq_len(6), c(3, 1, 4, 6, 2, 5), c(4, 6, 1, 3, 5, 2)) / 7
example <- list(latent_positions = truth, true_assign = rep(LETTERS[1:3], each = 2),
                ordering_labels = LETTERS[1:3])
perfect <- list(active = 1:3, positions = cbind(truth[, 3], 1 - truth[, 1], truth[, 2]),
                assignments = c(2, 2, 3, 3, 1, 1))
result <- evaluate_fit(perfect, example)
stopifnot(abs(result$ARI - 1) < 1e-12, abs(result$ordering_recovery - 1) < 1e-12,
          all(result$alignment$estimated_slot == c(2, 3, 1)))
missing <- perfect
missing$active <- 1:2
missing$assignments <- c(2, 2, 2, 2, 1, 1)
result <- evaluate_fit(missing, example)
stopifnot(abs(result$ordering_recovery - 2 / 3) < 1e-12,
          abs(result$true_ordering_coverage - 2 / 3) < 1e-12,
          result$ARI < 1, sum(!result$alignment$matched) == 1L)
constant <- perfect
constant$positions[,] <- 0.5
stopifnot(evaluate_fit(constant, example)$ordering_recovery == 0)
jsonlite::write_json(list(data_generation = TRUE, paired_input_hashes = TRUE,
  label_and_reversal_invariance = TRUE, unmatched_ordering_penalty = TRUE,
  constant_ordering_score = TRUE, design_hash = design_hash),
  file.path(study_dir, "setup_validation.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Data-generation and scoring checks passed.\n")

active <- active_manifest()
stopifnot(sum(active$phase == "main") == 90L,
          all(table(active$true_M[active$phase == "main"], active$snr[active$phase == "main"]) == 10L))
