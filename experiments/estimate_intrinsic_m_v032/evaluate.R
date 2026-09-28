# Truth-based evaluation with one-to-one, reversal-invariant ordering matching.
evaluate_fit <- function(compact, dataset) {
  active <- compact$active
  positions <- compact$positions[, active, drop = FALSE]
  truth <- dataset$latent_positions
  correlation <- abs(suppressWarnings(cor(truth, positions, method = "spearman")))
  correlation[!is.finite(correlation)] <- 0
  size <- max(ncol(truth), ncol(positions))
  score <- matrix(0, size, size)
  score[seq_len(ncol(truth)), seq_len(ncol(positions))] <- correlation
  assignment <- as.integer(clue::solve_LSAP(score, maximum = TRUE))[seq_len(ncol(truth))]
  matched <- assignment <= ncol(positions)
  rho <- numeric(ncol(truth))
  rho[matched] <- correlation[cbind(which(matched), assignment[matched])]
  alignment <- data.frame(true_ordering = dataset$ordering_labels,
    estimated_slot = ifelse(matched, active[pmin(assignment, length(active))], NA_integer_),
    abs_spearman = rho, matched = matched)
  list(ARI = mclust::adjustedRandIndex(dataset$true_assign, compact$assignments),
    ordering_recovery = mean(rho), matched_ordering_recovery = mean(rho[matched]),
    true_ordering_coverage = mean(matched), alignment = alignment)
}
method_result <- function(dataset, method, selected, candidates, history = NULL) {
  successful <- !is.null(selected) && selected$status == "success"
  evaluation <- if (successful) evaluate_fit(selected$fit, dataset) else NULL
  failed <- candidates[[length(candidates)]]
  row <- data.frame(id = dataset$row$id, phase = dataset$row$phase,
    true_M = dataset$row$true_M, snr = dataset$row$snr,
    replicate = dataset$row$replicate, method = method,
    status = if (successful) "success" else failed$status,
    estimated_M = if (successful) selected$fit$effective_M else NA_integer_,
    ARI = if (successful) evaluation$ARI else NA_real_,
    ordering_recovery = if (successful) evaluation$ordering_recovery else NA_real_,
    matched_ordering_recovery = if (successful) evaluation$matched_ordering_recovery else NA_real_,
    true_ordering_coverage = if (successful) evaluation$true_ordering_coverage else NA_real_,
    map_groups = if (successful) length(unique(selected$fit$assignments)) else NA_integer_,
    near_duplicate_orderings = if (successful) selected$fit$near_duplicate_orderings else NA,
    constant_orderings = if (successful) selected$fit$constant_orderings else NA_integer_,
    elapsed_seconds = sum(vapply(candidates, function(x) x$elapsed_seconds, numeric(1))),
    candidate_count = length(candidates),
    warning_count = sum(vapply(candidates, function(x) length(x$warnings), integer(1))),
    error = if (successful) "" else failed$error)
  value <- list(row = row, design_hash = design_hash, input_hash = dataset$input_hash,
    selected = if (successful) selected$fit else NULL, evaluation = evaluation,
    candidate_keys = vapply(candidates, function(x) x$key, character(1)), history = history)
  atomic_save(value, file.path(study_dir, "results", paste0(dataset$row$id, "_", method, ".rds")))
  value
}
