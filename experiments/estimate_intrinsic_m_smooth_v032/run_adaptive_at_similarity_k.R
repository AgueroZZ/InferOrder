#!/usr/bin/env Rscript

# Fit adaptive EB at the data-only cluster count selected from a saved
# spline-R-squared similarity diagnostic. This experiment changes only the
# fitted dimension; all initialization and fitting controls remain unchanged.

source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "fit.R"))
source(file.path(study_dir, "evaluate.R"))

arguments <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(arguments) >= 1L) arguments[[1L]] else "main_M5_S1_r001"
spline_df <- if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 5L
assignment_initialization <- if (length(arguments) >= 3L) arguments[[3L]] else "default"
assignment_initialization <- match.arg(assignment_initialization,
  c("default", "cluster_concentrated", "cluster_proportion_prior"))
assignment_concentration <- if (length(arguments) >= 4L) {
  as.numeric(arguments[[4L]])
} else {
  0.9
}
cluster_k_override <- if (length(arguments) >= 5L) {
  as.integer(arguments[[5L]])
} else {
  NA_integer_
}
stopifnot(length(spline_df) == 1L, is.finite(spline_df), spline_df >= 1L)
stopifnot(length(assignment_concentration) == 1L,
  is.finite(assignment_concentration), assignment_concentration > 0,
  assignment_concentration < 1)

manifest <- active_manifest()
row <- manifest[manifest$id == dataset_id, , drop = FALSE]
stopifnot(nrow(row) == 1L, row$phase == "main")
dataset <- load_dataset(row)

metric_name <- paste0("spline_r2_df", spline_df)
diagnostic_dir <- file.path(study_dir, "exploratory_similarity", dataset_id)
diagnostic_path <- file.path(diagnostic_dir,
  paste0(metric_name, "_direct_M_diagnostics.csv"))
stopifnot(file.exists(diagnostic_path))
diagnostics <- utils::read.csv(diagnostic_path, stringsAsFactors = FALSE)
silhouette_k <- diagnostics$k[which.max(diagnostics$mean_silhouette)]
cluster_k <- if (is.na(cluster_k_override)) silhouette_k else cluster_k_override
stopifnot(length(cluster_k) == 1L, cluster_k >= 2L, cluster_k <= design$max_M)
selection_rule <- if (is.na(cluster_k_override)) {
  "spline R-squared mean silhouette"
} else {
  "explicit spline R-squared dendrogram cut"
}

output_suffix <- if (assignment_initialization == "default" && cluster_k == silhouette_k) {
  paste0(metric_name, "_adaptive_at_silhouette_k")
} else if (assignment_initialization == "cluster_concentrated" &&
    cluster_k == silhouette_k) {
  paste0(metric_name, "_adaptive_at_silhouette_k_cluster_weights")
} else {
  paste0(metric_name, "_adaptive_at_k", cluster_k, "_",
    assignment_initialization)
}
output_dir <- file.path(diagnostic_dir, output_suffix)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

experiment_design <- list(
  baseline_design_hash = design_hash,
  dataset_id = dataset_id,
  input_hash = dataset$input_hash,
  similarity_metric = "directional natural-cubic-spline training R-squared",
  spline_df = spline_df,
  cluster_count_rule = selection_rule,
  selected_cluster_k = as.integer(cluster_k),
  fitted_dimension = as.integer(cluster_k),
  partition_prior = "adaptive",
  assignment_initialization = assignment_initialization,
  assignment_concentration = if (assignment_initialization ==
      "cluster_concentrated") assignment_concentration else NULL,
  fitting_settings = design)
experiment_hash <- digest::digest(experiment_design, algo = "sha256")
script_path <- file.path(study_dir, "run_adaptive_at_similarity_k.R")
provenance <- list(
  experiment_hash = experiment_hash,
  script_sha256 = digest::digest(file = script_path, algo = "sha256"),
  diagnostic_sha256 = digest::digest(file = diagnostic_path, algo = "sha256"),
  input_hash = dataset$input_hash,
  package_version = as.character(utils::packageVersion("MPCurver")),
  package_source_commit = design$package_source_commit,
  R_version = as.character(getRversion()))

compute_spline_r2_similarity <- function(X, min_feature_sd = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  feature_names <- colnames(X)
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(D))
  feature_sd <- apply(X, 2L, stats::sd)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd
  rank_position <- (seq_len(n) - 0.5) / n
  spline_basis <- cbind("(Intercept)" = 1,
    splines::ns(rank_position, df = spline_df, intercept = FALSE))
  basis_qr <- qr(spline_basis)
  Q <- qr.Q(basis_qr)
  stopifnot(basis_qr$rank == ncol(spline_basis))
  centered_X <- sweep(X, 2L, colMeans(X), FUN = "-")
  total_sum_squares <- colSums(centered_X^2)
  directional_r2 <- matrix(0, D, D,
    dimnames = list(feature_names, feature_names))
  valid_response <- !low_variance & total_sum_squares > sqrt(.Machine$double.eps)
  for (predictor in seq_len(D)) {
    if (low_variance[predictor]) next
    ordered_X <- X[order(X[, predictor], method = "radix"), , drop = FALSE]
    projected_coordinates <- crossprod(Q, ordered_X)
    residual_sum_squares <- pmax(
      colSums(ordered_X^2) - colSums(projected_coordinates^2), 0)
    directional_r2[predictor, valid_response] <-
      1 - residual_sum_squares[valid_response] / total_sum_squares[valid_response]
  }
  directional_r2 <- pmin(pmax(directional_r2, 0), 1)
  similarity <- pmax(directional_r2, t(directional_r2))
  if (any(low_variance)) {
    similarity[low_variance, ] <- 0
    similarity[, low_variance] <- 0
  }
  diag(similarity) <- 1
  list(
    S = similarity,
    distance = 1 - similarity,
    metric = metric_name,
    feature_info = data.frame(feature = feature_names, sd = as.numeric(feature_sd),
      low_variance = low_variance, stringsAsFactors = FALSE))
}

similarity_cache <- NULL
similarity_override <- function(X,
                                metric = c("spearman", "pearson"),
                                use = "pairwise.complete.obs",
                                abs_value = TRUE,
                                min_feature_sd = 1e-8) {
  metric <- match.arg(metric)
  if (is.null(similarity_cache)) {
    similarity_cache <<- compute_spline_r2_similarity(X, min_feature_sd)
    stopifnot(identical(digest::digest(X, algo = "sha256"), dataset$input_hash))
  }
  similarity_cache
}

original_similarity_cor <- getFromNamespace(
  ".compute_same_ordering_similarity_cor", "MPCurver")
assignInNamespace(".compute_same_ordering_similarity_cor", similarity_override,
  ns = "MPCurver")
on.exit(assignInNamespace(".compute_same_ordering_similarity_cor",
  original_similarity_cor, ns = "MPCurver"), add = TRUE)

# Experiment-only option: replace the structural model's uniform initial q(Z)
# with probabilities concentrated on the same canonical feature clusters used
# for the trajectory initialization. No subsequent update is changed.
initial_prior <- rep(1 / cluster_k, cluster_k)
cluster_sizes <- rep(ncol(dataset$X) / cluster_k, cluster_k)
if (assignment_initialization != "default") {
  saved_similarity_path <- file.path(diagnostic_dir,
    paste0(metric_name, "_similarity.csv"))
  saved_similarity <- as.matrix(utils::read.csv(saved_similarity_path,
    row.names = 1L, check.names = FALSE))
  storage.mode(saved_similarity) <- "double"
  cluster_tree <- stats::hclust(stats::as.dist(1 - saved_similarity),
    method = design$cluster_linkage)
  raw_cluster <- stats::cutree(cluster_tree, k = cluster_k)
  cluster_info <- getFromNamespace(
    ".cavi_canonicalize_feature_clusters", "MPCurver")(raw_cluster)
  feature_cluster <- cluster_info$feature_cluster
  cluster_sizes <- as.numeric(cluster_info$cluster_sizes)
  initial_prior <- cluster_sizes / sum(cluster_sizes)
  if (assignment_initialization == "cluster_concentrated") {
    off_cluster_probability <- (1 - assignment_concentration) / (cluster_k - 1L)
    initial_assignment_weights <- matrix(off_cluster_probability,
      nrow = ncol(dataset$X), ncol = cluster_k,
      dimnames = list(colnames(dataset$X), LETTERS[seq_len(cluster_k)]))
    initial_assignment_weights[cbind(seq_len(ncol(dataset$X)), feature_cluster)] <-
      assignment_concentration
  } else {
    initial_assignment_weights <- matrix(rep(initial_prior,
      each = ncol(dataset$X)), nrow = ncol(dataset$X), ncol = cluster_k,
      dimnames = list(colnames(dataset$X), LETTERS[seq_len(cluster_k)]))
  }

  original_initial_state <- getFromNamespace(
    ".structural_partition_initial_state", "MPCurver")
  replacement_count <- 0L
  replace_uniform_weights <- function(node) {
    if (!is.call(node)) return(node)
    if (length(node) == 3L && identical(node[[1L]], as.name("<-")) &&
        identical(node[[2L]], as.name("weights")) && is.call(node[[3L]]) &&
        identical(node[[3L]][[1L]], as.name("matrix")) &&
        grepl("1/M", paste(deparse(node[[3L]]), collapse = ""), fixed = TRUE)) {
      replacement_count <<- replacement_count + 1L
      node[[3L]] <- as.name(".cluster_assignment_weights")
      return(node)
    }
    for (element in seq_along(node)) {
      if (is.call(node[[element]])) {
        node[[element]] <- replace_uniform_weights(node[[element]])
      }
    }
    node
  }
  patched_initial_state <- original_initial_state
  body(patched_initial_state) <- replace_uniform_weights(body(original_initial_state))
  stopifnot(replacement_count == 1L)
  patched_environment <- new.env(parent = environment(original_initial_state))
  patched_environment$.cluster_assignment_weights <- initial_assignment_weights
  environment(patched_initial_state) <- patched_environment
  assignInNamespace(".structural_partition_initial_state", patched_initial_state,
    ns = "MPCurver")
  on.exit(assignInNamespace(".structural_partition_initial_state",
    original_initial_state, ns = "MPCurver"), add = TRUE)
}

utils::write.csv(data.frame(component = seq_len(cluster_k),
  cluster_size = cluster_sizes, initial_prior = initial_prior),
  file.path(output_dir, "initial_prior.csv"), row.names = FALSE)

set.seed(dataset$row$seed + 1000000L)
warnings <- character()
elapsed <- system.time({
  fitted <- withCallingHandlers({
    args <- c(list(X = dataset$X, intrinsic_dim = as.integer(cluster_k),
      partition_init = "similarity", partition_prior = "adaptive",
      effective_weight_tol = design$effective_weight_tol,
      hard_assign_final = FALSE), fit_options())
    value <- do.call(MPCurver::fit_mpcurve, args)
    repeat {
      remaining <- design$maximum_sweeps - value$fit$iter - 1L
      if (remaining <= 0L || isTRUE(value$fit$converged)) break
      value <- MPCurver::do_mpcurve(value,
        iter = min(design$continuation_budget, remaining),
        tol = design$tolerance, tol_outer = design$tolerance,
        convergence = design$convergence, verbose = FALSE)
    }
    value
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
})

compact <- compact_fit(fitted, cluster_k, "adaptive")
evaluation <- evaluate_fit(compact, dataset)
stopifnot(compact$converged)

summary <- data.frame(
  id = dataset_id,
  selection_rule = selection_rule,
  assignment_initialization = assignment_initialization,
  assignment_concentration = if (assignment_initialization ==
      "cluster_concentrated") assignment_concentration else NA_real_,
  selected_cluster_k = as.integer(cluster_k),
  fitted_model_M = compact$selected_model_M,
  effective_M = compact$effective_M,
  true_M = dataset$row$true_M,
  exact_M = compact$effective_M == dataset$row$true_M,
  ARI = evaluation$ARI,
  ordering_recovery = evaluation$ordering_recovery,
  matched_ordering_recovery = evaluation$matched_ordering_recovery,
  true_ordering_coverage = evaluation$true_ordering_coverage,
  objective = compact$score,
  iterations = compact$iterations,
  elapsed_seconds = unname(elapsed[["elapsed"]]),
  warning_count = length(unique(warnings)),
  stringsAsFactors = FALSE)

utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(component = seq_along(compact$occupancy),
  occupancy = compact$occupancy, prior_weight = compact$omega),
  file.path(output_dir, "occupancy.csv"), row.names = FALSE)
jsonlite::write_json(list(experiment_design = experiment_design,
  provenance = provenance, warnings = unique(warnings)),
  file.path(output_dir, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
saveRDS(list(summary = summary, fit = compact, evaluation = evaluation,
  experiment_design = experiment_design, provenance = provenance),
  file.path(output_dir, "result.rds"), compress = FALSE)

print(summary, row.names = FALSE)
print(data.frame(component = seq_along(compact$occupancy),
  occupancy = compact$occupancy, prior_weight = compact$omega), row.names = FALSE)
