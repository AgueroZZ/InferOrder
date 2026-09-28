#!/usr/bin/env Rscript

# Explore whether the feature-similarity initialization separates true orderings.
# This script reads saved inputs and result summaries without modifying fitted results.

if (!file.exists("_workflowr.yml")) {
  stop("Run this script from the InferOrder repository root.")
}

study_dir <- file.path("experiments", "estimate_intrinsic_m_smooth_v032")
arguments <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(arguments)) arguments[[1L]] else "main_M5_S1_r001"
output_dir <- file.path(study_dir, "exploratory_similarity", dataset_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

runs <- read.csv(file.path(study_dir, "main_summary", "runs.csv"),
  stringsAsFactors = FALSE)
selected_runs <- runs[runs$id == dataset_id & runs$method %in% c("adaptive", "forward"), ]
stopifnot(nrow(selected_runs) == 2L, all(selected_runs$status == "success"),
  length(unique(selected_runs$true_M)) == 1L,
  all(selected_runs$effective_M != selected_runs$true_M))

dataset_path <- file.path(study_dir, "data", paste0(dataset_id, ".rds"))
dataset <- readRDS(dataset_path)
X <- as.matrix(dataset$X)
true_assignment <- as.character(dataset$true_assign)
true_levels <- if (!is.null(dataset$ordering_labels)) {
  as.character(dataset$ordering_labels)
} else {
  unique(true_assignment)
}
true_assignment <- factor(true_assignment, levels = true_levels)
stopifnot(ncol(X) == length(true_assignment), !anyNA(true_assignment))

# Reproduce MPCurver's absolute-Spearman initialization similarity.
feature_sd <- apply(X, 2L, stats::sd)
low_variance <- !is.finite(feature_sd) | feature_sd < 1e-8
similarity <- suppressWarnings(stats::cor(X,
  method = "spearman", use = "pairwise.complete.obs"))
similarity[!is.finite(similarity)] <- 0
similarity <- abs(similarity)
similarity <- pmin(pmax(similarity, 0), 1)
similarity <- 0.5 * (similarity + t(similarity))
if (any(low_variance)) {
  similarity[low_variance, ] <- 0
  similarity[, low_variance] <- 0
}
diag(similarity) <- 1
dissimilarity <- 1 - similarity

cluster_tree <- stats::hclust(stats::as.dist(dissimilarity), method = "single")
cluster_order <- cluster_tree$order

adjusted_rand_index <- function(first, second) {
  contingency <- table(first, second)
  choose_two <- function(values) values * (values - 1) / 2
  agreement <- sum(choose_two(contingency))
  first_pairs <- sum(choose_two(rowSums(contingency)))
  second_pairs <- sum(choose_two(colSums(contingency)))
  all_pairs <- choose_two(sum(contingency))
  expected <- first_pairs * second_pairs / all_pairs
  denominator <- 0.5 * (first_pairs + second_pairs) - expected
  if (denominator == 0) return(as.numeric(agreement == expected))
  (agreement - expected) / denominator
}

mean_silhouette <- function(cluster, distance_matrix) {
  cluster <- as.integer(factor(cluster))
  widths <- numeric(length(cluster))
  for (feature_index in seq_along(cluster)) {
    same_cluster <- which(cluster == cluster[feature_index])
    same_cluster <- setdiff(same_cluster, feature_index)
    if (!length(same_cluster)) {
      widths[feature_index] <- 0
      next
    }
    within_distance <- mean(distance_matrix[feature_index, same_cluster])
    other_clusters <- setdiff(unique(cluster), cluster[feature_index])
    between_distance <- min(vapply(other_clusters, function(other_cluster) {
      mean(distance_matrix[feature_index, cluster == other_cluster])
    }, numeric(1)))
    scale <- max(within_distance, between_distance)
    widths[feature_index] <- if (scale > 0) {
      (between_distance - within_distance) / scale
    } else {
      0
    }
  }
  mean(widths)
}

candidate_k <- 2:8
diagnostics <- do.call(rbind, lapply(candidate_k, function(k) {
  cluster <- stats::cutree(cluster_tree, k = k)
  cluster_sizes <- tabulate(cluster, nbins = k)
  merge_index <- ncol(X) - k
  height_gap <- cluster_tree$height[merge_index + 1L] - cluster_tree$height[merge_index]
  data.frame(k = k,
    mean_silhouette = mean_silhouette(cluster, dissimilarity),
    height_gap = height_gap,
    ARI_to_truth = adjusted_rand_index(true_assignment, cluster),
    minimum_cluster_size = min(cluster_sizes),
    maximum_cluster_size = max(cluster_sizes),
    singleton_clusters = sum(cluster_sizes == 1L))
}))
silhouette_k <- diagnostics$k[which.max(diagnostics$mean_silhouette)]
height_gap_k <- diagnostics$k[which.max(diagnostics$height_gap)]
true_M <- unique(selected_runs$true_M)

fit_partition_comparison <- do.call(rbind, lapply(c("adaptive", "forward"), function(method) {
  result <- readRDS(file.path(study_dir, "results", paste0(dataset_id, "_", method, ".rds")))
  final_assignment <- result$selected$assignments
  effective_M <- selected_runs$effective_M[selected_runs$method == method]
  tree_assignment <- if (effective_M == 1L) {
    rep(1L, ncol(X))
  } else {
    stats::cutree(cluster_tree, k = effective_M)
  }
  ARI_to_truth <- adjusted_rand_index(true_assignment, final_assignment)
  stopifnot(isTRUE(all.equal(ARI_to_truth,
    selected_runs$ARI[selected_runs$method == method], tolerance = 1e-12)))
  data.frame(method = method, effective_M = effective_M,
    ARI_to_truth = ARI_to_truth,
    ARI_to_single_linkage_cut_at_effective_M =
      adjusted_rand_index(final_assignment, tree_assignment))
}))

feature_metadata <- data.frame(
  feature_index = seq_len(ncol(X)),
  feature = colnames(X),
  true_ordering = as.character(true_assignment),
  monotone_anchor = seq_len(ncol(X)) %in% dataset$anchor_indices,
  dendrogram_position = match(seq_len(ncol(X)), cluster_order),
  cluster_at_silhouette_k = stats::cutree(cluster_tree, k = silhouette_k),
  cluster_at_true_M = stats::cutree(cluster_tree, k = true_M),
  stringsAsFactors = FALSE)

utils::write.csv(diagnostics,
  file.path(output_dir, "direct_M_diagnostics.csv"), row.names = FALSE)
utils::write.csv(feature_metadata,
  file.path(output_dir, "feature_order_and_clusters.csv"), row.names = FALSE)
utils::write.csv(fit_partition_comparison,
  file.path(output_dir, "fit_partition_comparison.csv"), row.names = FALSE)
utils::write.csv(similarity,
  file.path(output_dir, "absolute_spearman_similarity.csv"), row.names = TRUE)

ordering_colors <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
  "#56B4E9", "#F0E442", "#000000")[seq_along(true_levels)], true_levels)
side_colors <- unname(ordering_colors[as.character(true_assignment)])
heatmap_colors <- grDevices::colorRampPalette(c(
  "#081D58", "#225EA8", "#41B6C4", "#C7E9B4", "#FFFFD9"))(256L)

adaptive_M <- selected_runs$effective_M[selected_runs$method == "adaptive"]
forward_M <- selected_runs$effective_M[selected_runs$method == "forward"]
heatmap_path <- file.path(output_dir, "clustered_absolute_spearman_similarity.png")
grDevices::png(heatmap_path, width = 2600, height = 2300, res = 220)
graphics::par(oma = c(1, 1, 8, 1))
stats::heatmap(similarity,
  Rowv = stats::as.dendrogram(cluster_tree), Colv = "Rowv",
  scale = "none", symm = TRUE, revC = TRUE,
  col = heatmap_colors, margins = c(6, 6),
  labRow = seq_len(ncol(X)), labCol = seq_len(ncol(X)),
  RowSideColors = side_colors, ColSideColors = side_colors,
  xlab = "Feature index", ylab = "Feature index",
  main = "")
graphics::mtext(paste0("Clustered feature similarity: ", dataset_id),
  outer = TRUE, side = 3, line = 5.0, cex = 1.2, font = 2)
graphics::mtext("Absolute Spearman correlation; single linkage",
  outer = TRUE, side = 3, line = 3.3, cex = 1.0)
graphics::mtext(sprintf(
  "True M = %d; Adaptive EB effective M = %d; Uniform + forward effective M = %d",
  true_M, adaptive_M, forward_M), outer = TRUE, side = 3, line = 1.7, cex = 0.9)
graphics::legend("topright", inset = 0.01, title = "True ordering",
  legend = true_levels, fill = unname(ordering_colors[true_levels]),
  border = NA, bty = "n", cex = 0.75)
similarity_breaks <- seq(0, 1, by = 0.25)
similarity_color_index <- 1L + round(similarity_breaks * (length(heatmap_colors) - 1L))
graphics::legend("bottomright", inset = 0.01, title = "Similarity",
  legend = format(similarity_breaks, trim = TRUE),
  fill = heatmap_colors[similarity_color_index], border = NA, bty = "n", cex = 0.65)
grDevices::dev.off()

diagnostic_path <- file.path(output_dir, "direct_M_diagnostics.png")
grDevices::png(diagnostic_path, width = 2200, height = 1000, res = 200)
old_par <- graphics::par(mfrow = c(1, 2), mar = c(9.0, 4.8, 3, 4.8))
graphics::plot(diagnostics$k, diagnostics$mean_silhouette, type = "b", pch = 16,
  lwd = 2, col = "#0072B2", ylim = range(c(diagnostics$mean_silhouette,
    diagnostics$ARI_to_truth)), xlab = "Number of clusters k", ylab = "Score",
  main = "Single-linkage cuts")
graphics::lines(diagnostics$k, diagnostics$ARI_to_truth, type = "b", pch = 17,
  lwd = 2, col = "#D55E00")
graphics::abline(v = true_M, lty = 3, col = "#333333")
graphics::legend("bottomleft", inset = c(0, -0.62), xpd = NA,
  legend = c("Mean silhouette (data-only)", "ARI to truth (evaluation)", "True M"),
  col = c("#0072B2", "#D55E00", "#333333"), lty = c(1, 1, 3),
  pch = c(16, 17, NA), bty = "n")
graphics::plot(diagnostics$k, diagnostics$height_gap, type = "b", pch = 16,
  lwd = 2, col = "#009E73", xlab = "Number of clusters k",
  ylab = "Next merge-height gap", main = "Dendrogram gap heuristic")
graphics::abline(v = true_M, lty = 3, col = "#333333")
graphics::mtext(sprintf("Data-only estimates: silhouette k = %d; height-gap k = %d",
  silhouette_k, height_gap_k), side = 3, line = 0.25, cex = 0.85)
graphics::par(old_par)
grDevices::dev.off()

summary_lines <- c(
  paste0("Dataset: ", dataset_id),
  sprintf("True M: %d", true_M),
  sprintf("Adaptive EB effective M: %d", adaptive_M),
  sprintf("Uniform + forward effective M: %d", forward_M),
  sprintf("Single-linkage silhouette estimate: %d", silhouette_k),
  sprintf("Single-linkage merge-height-gap estimate: %d", height_gap_k),
  sprintf("ARI of the k = true M tree cut: %.6f",
    diagnostics$ARI_to_truth[diagnostics$k == true_M]),
  sprintf("ARI between the forward partition and the single-linkage k = %d cut: %.6f",
    forward_M,
    fit_partition_comparison$ARI_to_single_linkage_cut_at_effective_M[
      fit_partition_comparison$method == "forward"]),
  "Similarity: absolute Spearman correlation of noisy observed features",
  "Clustering: single linkage on distance 1 - similarity",
  paste0("Input: ", dataset_path),
  "The two direct-M estimates use no truth labels; ARI is reported only for evaluation.")
writeLines(summary_lines, file.path(output_dir, "README.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Saved exploratory outputs under", output_dir, "\n")
