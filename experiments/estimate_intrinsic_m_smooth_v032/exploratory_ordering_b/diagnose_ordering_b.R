#!/usr/bin/env Rscript

# Diagnose the imperfect ordering-B recovery in main_M3_S4_r001.
# This is a read-only analysis of saved study fits plus deterministic
# reconstruction of their Isomap initializations.

source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(igraph)
})

dataset_id <- "main_M3_S4_r001"
output_dir <- file.path(study_dir, "exploratory_ordering_b")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset <- readRDS(file.path(study_dir, "data", paste0(dataset_id, ".rds")))
adaptive <- readRDS(file.path(study_dir, "results", paste0(dataset_id, "_adaptive.rds")))
forward <- readRDS(file.path(study_dir, "results", paste0(dataset_id, "_forward.rds")))
auto <- readRDS(file.path(study_dir, "auto_m_v034", "results",
  paste0(dataset_id, "_auto_adaptive.rds")))

stopifnot(
  identical(dataset$ordering_labels, c("A", "B", "C")),
  adaptive$row$effective_M == 3L,
  forward$row$effective_M == 3L,
  auto$row$effective_M == 3L,
  adaptive$evaluation$ARI == 1,
  forward$evaluation$ARI == 1,
  auto$evaluation$ARI == 1
)

grid <- seq(0, 1, length.out = design$K)
expected_position <- function(fit) as.numeric(fit$gamma %*% grid)
absolute_spearman <- function(x, y) abs(stats::cor(x, y, method = "spearman"))

reconstruct_initialization <- function(M) {
  set.seed(dataset$row$seed + 1000000L)
  MPCurver:::init_m_trajectories_cavi(
    X = dataset$X,
    S = NULL,
    M = M,
    methods = rep(design$method, M),
    K = design$K,
    rw_q = design$rw_q,
    ridge = design$ridge,
    lambda_sd_prior_rate = NULL,
    smooth_fit_lambda_mode = "optimize",
    smooth_fit_lambda_value = 1,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    sigma_min = 1e-10,
    sigma_max = 1e10,
    discretization = "quantile",
    partition_init = "similarity",
    similarity_metric = "spearman",
    cluster_linkage = "single",
    similarity_min_feature_sd = 1e-8,
    num_iter = 2L,
    verbose = FALSE
  )
}

initialization_table <- function(initialization, M) {
  do.call(rbind, lapply(seq_len(M), function(m) {
    feature_idx <- initialization$init_info[[m]]$feature_idx
    truth_counts <- sort(table(dataset$true_assign[feature_idx]), decreasing = TRUE)
    majority_truth <- names(truth_counts)[1L]
    truth_index <- match(majority_truth, dataset$ordering_labels)
    position <- expected_position(initialization$fits[[m]])
    data.frame(
      fitted_M = M,
      slot = m,
      slot_label = names(initialization$init_info)[m],
      majority_truth = majority_truth,
      cluster_size = length(feature_idx),
      majority_count = unname(truth_counts[1L]),
      contains_anchor = any(feature_idx %in%
        dataset$anchor_indices[dataset$ordering_labels == majority_truth]),
      initial_abs_spearman = absolute_spearman(
        dataset$latent_positions[, truth_index], position),
      feature_indices = paste(feature_idx, collapse = ";"),
      feature_names = paste(colnames(dataset$X)[feature_idx], collapse = ";")
    )
  }))
}

initial_M8 <- reconstruct_initialization(8L)
initial_M3 <- reconstruct_initialization(3L)
initialization_clusters <- rbind(
  initialization_table(initial_M8, 8L),
  initialization_table(initial_M3, 3L)
)
write.csv(initialization_clusters,
  file.path(output_dir, "initialization_clusters.csv"), row.names = FALSE)

adaptive_B_slot <- adaptive$evaluation$alignment$estimated_slot[
  adaptive$evaluation$alignment$true_ordering == "B"]
forward_B_slot <- forward$evaluation$alignment$estimated_slot[
  forward$evaluation$alignment$true_ordering == "B"]
auto_B_slot <- as.integer(names(which.max(table(
  auto$compact$assignments[dataset$true_assign == "B"]))))
true_B_position <- dataset$latent_positions[, match("B", dataset$ordering_labels)]

core_slot <- adaptive_B_slot
core_features <- initial_M8$init_info[[core_slot]]$feature_idx
all_B_features <- which(dataset$true_assign == "B")
excluded_B_features <- setdiff(all_B_features, core_features)
stopifnot(length(core_features) == 16L, length(excluded_B_features) == 4L)

fit_subset_ordering <- function(feature_idx) {
  set.seed(dataset$row$seed + 1000000L)
  subset_fit <- MPCurver:::.cavi_similarity_subset_fit(
    X_sub = dataset$X[, feature_idx, drop = FALSE],
    S = NULL,
    K = design$K,
    method = design$method,
    pca_component = NA_integer_,
    rw_q = design$rw_q,
    ridge = design$ridge,
    lambda_sd_prior_rate = NULL,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    sigma_min = 1e-10,
    sigma_max = 1e10,
    max_iter = 2L,
    tol = 1e-6,
    discretization = "quantile",
    cluster_label = "ordering-B diagnostic",
    verbose = FALSE
  )
  expected_position(subset_fit$fit)
}

addition_sets <- unlist(lapply(0:length(excluded_B_features), function(k) {
  combn(excluded_B_features, k, simplify = FALSE)
}), recursive = FALSE)
subset_sensitivity <- do.call(rbind, lapply(addition_sets, function(added) {
  position <- fit_subset_ordering(c(core_features, added))
  data.frame(
    added_indices = if (length(added)) paste(added, collapse = ";") else "none",
    added_features = if (length(added)) {
      paste(colnames(dataset$X)[added], collapse = ";")
    } else {
      "none"
    },
    contains_anchor = any(added %in% dataset$anchor_indices),
    feature_count = length(core_features) + length(added),
    initial_abs_spearman = absolute_spearman(true_B_position, position)
  )
}))
subset_sensitivity <- subset_sensitivity[
  order(subset_sensitivity$feature_count, -subset_sensitivity$initial_abs_spearman), ]
write.csv(subset_sensitivity,
  file.path(output_dir, "ordering_b_subset_sensitivity.csv"), row.names = FALSE)

geodesic_fidelity <- function(feature_idx, subset_label) {
  X_subset <- dataset$X[, feature_idx, drop = FALSE]
  neighbors <- RANN::nn2(X_subset, X_subset, k = 16L)
  edges <- data.frame(
    from = rep(seq_len(nrow(X_subset)), each = 15L),
    to = as.vector(t(neighbors$nn.idx[, -1L, drop = FALSE])),
    weight = as.vector(t(neighbors$nn.dists[, -1L, drop = FALSE]))
  )
  graph <- igraph::graph_from_data_frame(
    edges, directed = FALSE, vertices = seq_len(nrow(X_subset)))
  geodesic_distance <- igraph::distances(graph, weights = igraph::E(graph)$weight)
  true_distance <- as.matrix(stats::dist(true_B_position))
  upper <- upper.tri(true_distance)
  data.frame(
    subset = subset_label,
    feature_count = length(feature_idx),
    spearman = stats::cor(
      geodesic_distance[upper], true_distance[upper], method = "spearman"),
    pearson = stats::cor(geodesic_distance[upper], true_distance[upper]),
    connected_components = igraph::components(graph)$no
  )
}

single_addition <- subset_sensitivity[
  subset_sensitivity$feature_count == length(core_features) + 1L, ]
best_added_name <- single_addition$added_features[
  which.max(single_addition$initial_abs_spearman)]
best_added_index <- match(best_added_name, colnames(dataset$X))
stopifnot(length(best_added_index) == 1L, !is.na(best_added_index))

geodesic_summary <- rbind(
  geodesic_fidelity(core_features, "M=8 core (16)"),
  geodesic_fidelity(c(core_features, best_added_index),
    paste0("Core + ", best_added_name, " (17)")),
  geodesic_fidelity(all_B_features, "All B features (20)")
)
write.csv(geodesic_summary,
  file.path(output_dir, "ordering_b_geodesic_fidelity.csv"), row.names = FALSE)

feature_diagnostics <- data.frame(
  feature_index = all_B_features,
  feature = colnames(dataset$X)[all_B_features],
  is_anchor = all_B_features %in% dataset$anchor_indices,
  M8_initial_slot = vapply(all_B_features, function(j) {
    which(vapply(initial_M8$init_info, function(info) j %in% info$feature_idx, logical(1)))
  }, integer(1)),
  noiseless_abs_spearman = vapply(all_B_features, function(j) {
    absolute_spearman(true_B_position, dataset$signal[, j])
  }, numeric(1)),
  observed_abs_spearman = vapply(all_B_features, function(j) {
    absolute_spearman(true_B_position, dataset$X[, j])
  }, numeric(1)),
  realized_noise_variance = dataset$realized_noise_variance[all_B_features]
)
write.csv(feature_diagnostics,
  file.path(output_dir, "ordering_b_feature_diagnostics.csv"), row.names = FALSE)

noise_summary <- aggregate(dataset$realized_noise_variance,
  list(ordering = dataset$true_assign),
  function(x) c(mean = mean(x), minimum = min(x), maximum = max(x)))
noise_summary <- data.frame(
  ordering = noise_summary$ordering,
  mean = noise_summary$x[, "mean"],
  minimum = noise_summary$x[, "minimum"],
  maximum = noise_summary$x[, "maximum"]
)
write.csv(noise_summary,
  file.path(output_dir, "noise_variance_by_ordering.csv"), row.names = FALSE)

initial_B_position <- expected_position(initial_M8$fits[[core_slot]])
final_adaptive_B_position <- adaptive$selected$positions[, adaptive_B_slot]
final_forward_B_position <- forward$selected$positions[, forward_B_slot]
final_auto_B_position <- auto$compact$positions[, auto_B_slot]

orient_position <- function(position) {
  if (stats::cor(true_B_position, position, method = "spearman") < 0) 1 - position else position
}

position_frame <- rbind(
  data.frame(
    true_position = true_B_position,
    inferred_position = orient_position(initial_B_position),
    fit = "M=8 initialization\n(16-feature core)"
  ),
  data.frame(
    true_position = true_B_position,
    inferred_position = orient_position(final_adaptive_B_position),
    fit = "M=8 adaptive EB\n(final)"
  ),
  data.frame(
    true_position = true_B_position,
    inferred_position = orient_position(final_auto_B_position),
    fit = "Auto-selected M=3\n(final)"
  )
)
position_frame$fit <- factor(position_frame$fit,
  levels = c("M=8 initialization\n(16-feature core)",
    "M=8 adaptive EB\n(final)", "Auto-selected M=3\n(final)"))
position_labels <- do.call(rbind, lapply(levels(position_frame$fit), function(label) {
  frame <- position_frame[position_frame$fit == label, ]
  data.frame(
    fit = factor(label, levels = levels(position_frame$fit)),
    label = sprintf("|rho| = %.3f",
      absolute_spearman(frame$true_position, frame$inferred_position)),
    x = 0.04,
    y = 0.96
  )
}))

position_plot <- ggplot(position_frame,
  aes(true_position, inferred_position)) +
  geom_point(alpha = 0.38, size = 0.75, color = "#287D8E") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35,
    linetype = "dashed", color = "grey50") +
  geom_text(data = position_labels,
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3.1) +
  facet_wrap(~fit, nrow = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "The M=8 fit preserves its imperfect ordering-B initialization",
    x = "True ordering-B position",
    y = "Inferred position"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

single_addition$label <- paste0("+ ", single_addition$added_features)
single_addition$label <- factor(single_addition$label,
  levels = single_addition$label[order(single_addition$initial_abs_spearman)])
core_rho <- subset_sensitivity$initial_abs_spearman[
  subset_sensitivity$feature_count == 16L]
all_rho <- subset_sensitivity$initial_abs_spearman[
  subset_sensitivity$feature_count == 20L]

sensitivity_plot <- ggplot(single_addition,
  aes(label, initial_abs_spearman,
    fill = added_features == best_added_name)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = core_rho, linetype = "dashed", color = "#D55E00") +
  geom_hline(yintercept = all_rho, linetype = "dotted", color = "#0072B2") +
  annotate("text", x = 0.55, y = core_rho + 0.018,
    label = sprintf("16-feature core: %.3f", core_rho), hjust = 0, size = 3) +
  annotate("text", x = 0.55, y = all_rho - 0.005,
    label = sprintf("All 20 features: %.3f", all_rho), hjust = 1, vjust = 1, size = 3) +
  scale_fill_manual(values = c("TRUE" = "#009E73", "FALSE" = "grey65"),
    guide = "none") +
  coord_flip(ylim = c(0.65, 1.01)) +
  labs(
    title = "One excluded trajectory restores the Isomap ordering",
    subtitle = "Each bar adds one of the four B features excluded from the 16-feature core",
    x = NULL,
    y = "Initial ordering |Spearman correlation|"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

geodesic_summary$subset <- factor(geodesic_summary$subset,
  levels = geodesic_summary$subset)
geodesic_plot <- ggplot(geodesic_summary, aes(subset, spearman)) +
  geom_col(fill = c("#D55E00", "#009E73", "#0072B2"), width = 0.68) +
  geom_text(aes(label = sprintf("%.3f", spearman)), vjust = -0.35, size = 3.2) +
  coord_cartesian(ylim = c(0.70, 1.01)) +
  labs(
    title = "The missing trajectory repairs Isomap's geodesic geometry",
    subtitle = "Correlation between 15-NN graph distance and true latent distance",
    x = NULL,
    y = "Spearman correlation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

figure <- position_plot / (sensitivity_plot | geodesic_plot) +
  plot_annotation(
    title = "Why ordering B is imperfect in main_M3_S4_r001",
    subtitle = paste0(
      "The old M=8 absolute-Spearman cut splits B into 3 + 16 + 1 features.\n",
      "The 16-feature component survives adaptive EB but starts from a folded Isomap ordering."
    ),
    caption = paste0(
      "MPCurver 0.3.2; true M=3, variance SNR=4, replicate 1. ",
      "The highlighted feature is selected only for this truth-aware diagnostic; ",
      "it is not a proposed fitting rule."
    ),
    theme = theme(
      plot.title = element_text(size = 17, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 11, lineheight = 1.05,
        margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, margin = margin(t = 7)),
      plot.margin = margin(12, 12, 12, 12)
    )
  )

ggsave(file.path(output_dir, "ordering_b_diagnostic.png"), figure,
  width = 12, height = 9.2, dpi = 180, bg = "white")
ggsave(file.path(output_dir, "ordering_b_diagnostic.pdf"), figure,
  width = 12, height = 9.2, bg = "white")

summary <- data.frame(
  dataset_id = dataset_id,
  true_M = dataset$row$true_M,
  snr = dataset$row$snr,
  replicate = dataset$row$replicate,
  old_adaptive_ARI = adaptive$evaluation$ARI,
  old_adaptive_B_abs_spearman = adaptive$evaluation$alignment$abs_spearman[
    adaptive$evaluation$alignment$true_ordering == "B"],
  old_adaptive_B_initial_abs_spearman = absolute_spearman(
    true_B_position, initial_B_position),
  forward_B_abs_spearman = absolute_spearman(
    true_B_position, final_forward_B_position),
  auto_B_abs_spearman = absolute_spearman(true_B_position, final_auto_B_position),
  M8_B_cluster_sizes = paste(sort(vapply(initial_M8$init_info, function(info) {
    sum(dataset$true_assign[info$feature_idx] == "B")
  }, integer(1))[vapply(initial_M8$init_info, function(info) {
    any(dataset$true_assign[info$feature_idx] == "B")
  }, logical(1))]), collapse = "+"),
  core_feature_count = length(core_features),
  best_single_added_feature = best_added_name,
  core_initial_abs_spearman = core_rho,
  core_plus_best_initial_abs_spearman = max(single_addition$initial_abs_spearman),
  all_B_initial_abs_spearman = all_rho,
  core_geodesic_spearman = geodesic_summary$spearman[1L],
  core_plus_best_geodesic_spearman = geodesic_summary$spearman[2L],
  all_B_geodesic_spearman = geodesic_summary$spearman[3L],
  mean_noise_variance_A = noise_summary$mean[noise_summary$ordering == "A"],
  mean_noise_variance_B = noise_summary$mean[noise_summary$ordering == "B"],
  mean_noise_variance_C = noise_summary$mean[noise_summary$ordering == "C"]
)
write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)

cat("Ordering-B diagnostic complete.\n")
print(summary)
