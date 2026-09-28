#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
source(file.path(study_dir, "fit.R"))
source(file.path(study_dir, "evaluate.R"))
id <- commandArgs(trailingOnly = TRUE)[1]
manifest <- active_manifest()
row <- manifest[manifest$id == id, , drop = FALSE]
stopifnot(nrow(row) == 1L)
dataset <- load_dataset(row)
adaptive <- run_candidate(dataset, "adaptive", design$max_M)
eb_result <- method_result(dataset, "adaptive",
  if (adaptive$status == "success") adaptive else NULL, list(adaptive))
history <- data.frame(current_M = integer(), candidate_M = integer(),
  current_score = numeric(), candidate_score = numeric(), delta = numeric(), accepted = logical())
candidates <- list()
selected <- NULL
for (M in seq_len(design$max_M)) {
  candidate <- run_candidate(dataset, "forward", M)
  candidates[[M]] <- candidate
  if (candidate$status != "success") {
    selected <- NULL
    break
  }
  if (M == 1L) {
    selected <- candidate
    next
  }
  delta <- candidate$fit$score - selected$fit$score
  accepted <- is.finite(delta) && delta > 0
  history <- rbind(history, data.frame(current_M = selected$M, candidate_M = M,
    current_score = selected$fit$score, candidate_score = candidate$fit$score,
    delta = delta, accepted = accepted))
  if (!accepted) break
  selected <- candidate
}
forward_result <- method_result(dataset, "forward", selected, candidates, history)
cat("DATASET COMPLETE\n")
print(rbind(eb_result$row, forward_result$row))
