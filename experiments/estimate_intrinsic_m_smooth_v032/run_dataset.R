#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "fit.R"))
source(file.path(study_dir, "evaluate.R"))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript run_dataset.R <one-based manifest index or dataset id>")
}
manifest <- active_manifest()
manifest <- manifest[manifest$phase == "main", , drop = FALSE]
stopifnot(nrow(manifest) == 90L, !anyDuplicated(manifest$id))
if (grepl("^[0-9]+$", arguments[1L])) {
  index <- as.integer(arguments[1L])
  stopifnot(!is.na(index), index >= 1L, index <= nrow(manifest))
  row <- manifest[index, , drop = FALSE]
} else {
  row <- manifest[manifest$id == arguments[1L], , drop = FALSE]
}
stopifnot(nrow(row) == 1L)
dataset <- load_dataset(row)

adaptive <- run_candidate(dataset, "adaptive", design$max_M)
eb_result <- method_result(dataset, "adaptive",
  if (adaptive$status == "success") adaptive else NULL, list(adaptive),
  stop_reason = if (adaptive$status == "success") "fixed_maximum_dimension" else
    paste0("candidate_", adaptive$status))

history <- data.frame(current_M = integer(), candidate_M = integer(),
  current_score = numeric(), candidate_score = numeric(), delta = numeric(), accepted = logical())
candidates <- list()
selected <- NULL
stop_reason <- "maximum_dimension"
for (M in seq_len(design$max_M)) {
  candidate <- run_candidate(dataset, "forward", M)
  candidates[[M]] <- candidate
  if (candidate$status != "success") {
    selected <- NULL
    stop_reason <- paste0("candidate_", candidate$status)
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
  if (!accepted) {
    stop_reason <- "first_nonimprovement"
    break
  }
  selected <- candidate
}
forward_result <- method_result(dataset, "forward", selected, candidates, history, stop_reason)
cat("DATASET COMPLETE\n")
print(rbind(eb_result$row, forward_result$row))
