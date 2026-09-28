#!/usr/bin/env Rscript

# Explore a fast directional natural-spline R-squared feature similarity.
# This script reads saved inputs and result summaries without modifying fitted results.

if (!file.exists("_workflowr.yml")) {
  stop("Run this script from the InferOrder repository root.")
}

study_dir <- file.path("experiments", "estimate_intrinsic_m_smooth_v032")
arguments <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(arguments) >= 1L) arguments[[1L]] else "main_M5_S1_r001"
spline_df <- if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 5L
stopifnot(length(spline_df) == 1L, is.finite(spline_df), spline_df >= 1L)

output_dir <- file.path(study_dir, "exploratory_similarity", dataset_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_prefix <- paste0("spline_r2_df", spline_df)

runs <- read.csv(file.path(study_dir, "main_summary", "runs.csv"),
  stringsAsFactors = FALSE)
selected_runs <- runs[runs$id == dataset_id & runs$method %in% c("adaptive", "forward"), ]
stopifnot(nrow(selected_runs) == 2L, all(selected_runs$status == "success"),
  length(unique(selected_runs$true_M)) == 1L,
  all(selected_runs$effective_M != selected_runs$true_M))

dataset_path <- file.path(study_dir, "data", paste0(dataset_id, ".rds"))
dataset <- readRDS(dataset_path)
X <- as.matrix(dataset$X)
n <- nrow(X)
D <- ncol(X)
stopifnot(n > spline_df + 2L, D >= 2L, all(is.finite(X)))

true_assignment <- as.character(dataset$true_assign)
true_levels <- if (!is.null(dataset$ordering_labels)) {
  as.character(dataset$ordering_labels)
} else {
  unique(true_assignment)
}
true_assignment <- factor(true_assignment, levels = true_levels)
stopifnot(length(true_assignment) == D, !anyNA(true_assignment))

# Every predictor uses the same normalized-rank design after its samples are
# sorted, so the natural-spline basis and its QR decomposition are computed once.
rank_position <- (seq_len(n) - 0.5) / n
spline_basis <- cbind("(Intercept)" = 1,
  splines::ns(rank_position, df = spline_df, intercept = FALSE))
basis_qr <- qr(spline_basis)
orthonormal_basis <- qr.Q(basis_qr)
basis_rank <- basis_qr$rank
stopifnot(basis_rank == ncol(spline_basis))

centered_X <- sweep(X, 2L, colMeans(X), FUN = "-")
total_sum_squares <- colSums(centered_X^2)
directional_r2 <- matrix(0, D, D,
  dimnames = list(colnames(X), colnames(X)))

runtime <- system.time({
  for (predictor in seq_len(D)) {
    ordered_X <- X[order(X[, predictor], method = "radix"), , drop = FALSE]
    projected_coordinates <- crossprod(orthonormal_basis, ordered_X)
    residual_sum_squares <- colSums(ordered_X^2) - colSums(projected_coordinates^2)
    residual_sum_squares <- pmax(residual_sum_squares, 0)
    valid <- total_sum_squares > sqrt(.Machine$double.eps)
    directional_r2[predictor, valid] <-
      1 - residual_sum_squares[valid] / total_sum_squares[valid]
  }
})
directional_r2 <- pmin(pmax(directional_r2, 0), 1)

# Validate the vectorized projection against an explicit least-squares fit.
validation_pairs <- rbind(c(1L, 2L), c(1L, D), c(D, 1L))
for (pair_index in seq_len(nrow(validation_pairs))) {
  predictor <- validation_pairs[pair_index, 1L]
  response <- validation_pairs[pair_index, 2L]
  ordered_response <- X[order(X[, predictor], method = "radix"), response]
  explicit_fit <- stats::lm.fit(x = spline_basis, y = ordered_response)
  explicit_r2 <- 1 - sum(explicit_fit$residuals^2) / total_sum_squares[response]
  stopifnot(abs(directional_r2[predictor, response] - explicit_r2) < 1e-10)
}

similarity <- pmax(directional_r2, t(directional_r2))
diag(similarity) <- 1
dissimilarity <- 1 - similarity
cluster_tree <- stats::hclust(stats::as.dist(dissimilarity), method = "single")

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
    same_cluster <- setdiff(which(cluster == cluster[feature_index]), feature_index)
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
  merge_index <- D - k
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

same_group <- outer(as.character(true_assignment), as.character(true_assignment), "==")
upper <- upper.tri(similarity, diag = FALSE)
within_similarity <- similarity[upper & same_group]
between_similarity <- similarity[upper & !same_group]
anchor_similarity <- unlist(lapply(dataset$anchor_indices, function(anchor) {
  peers <- which(true_assignment == true_assignment[anchor])
  similarity[anchor, setdiff(peers, anchor)]
}), use.names = FALSE)
similarity_summary <- data.frame(
  relation = c("within true ordering", "between true orderings", "anchor to same-ordering feature"),
  pairs = c(length(within_similarity), length(between_similarity), length(anchor_similarity)),
  mean = c(mean(within_similarity), mean(between_similarity), mean(anchor_similarity)),
  median = c(stats::median(within_similarity), stats::median(between_similarity),
    stats::median(anchor_similarity)),
  q25 = c(stats::quantile(within_similarity, 0.25),
    stats::quantile(between_similarity, 0.25), stats::quantile(anchor_similarity, 0.25)),
  q75 = c(stats::quantile(within_similarity, 0.75),
    stats::quantile(between_similarity, 0.75), stats::quantile(anchor_similarity, 0.75)),
  stringsAsFactors = FALSE)

feature_metadata <- data.frame(
  feature_index = seq_len(D),
  feature = colnames(X),
  true_ordering = as.character(true_assignment),
  monotone_anchor = seq_len(D) %in% dataset$anchor_indices,
  dendrogram_position = match(seq_len(D), cluster_tree$order),
  cluster_at_silhouette_k = stats::cutree(cluster_tree, k = silhouette_k),
  cluster_at_true_M = stats::cutree(cluster_tree, k = true_M),
  stringsAsFactors = FALSE)

utils::write.csv(directional_r2,
  file.path(output_dir, paste0(output_prefix, "_directional_r2.csv")), row.names = TRUE)
utils::write.csv(similarity,
  file.path(output_dir, paste0(output_prefix, "_similarity.csv")), row.names = TRUE)
utils::write.csv(diagnostics,
  file.path(output_dir, paste0(output_prefix, "_direct_M_diagnostics.csv")), row.names = FALSE)
utils::write.csv(similarity_summary,
  file.path(output_dir, paste0(output_prefix, "_similarity_summary.csv")), row.names = FALSE)
utils::write.csv(feature_metadata,
  file.path(output_dir, paste0(output_prefix, "_feature_order_and_clusters.csv")),
  row.names = FALSE)

ordering_colors <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00",
  "#56B4E9", "#F0E442", "#000000")[seq_along(true_levels)], true_levels)
side_colors <- unname(ordering_colors[as.character(true_assignment)])
heatmap_colors <- grDevices::colorRampPalette(c(
  "#081D58", "#225EA8", "#41B6C4", "#C7E9B4", "#FFFFD9"))(256L)
adaptive_M <- selected_runs$effective_M[selected_runs$method == "adaptive"]
forward_M <- selected_runs$effective_M[selected_runs$method == "forward"]

heatmap_path <- file.path(output_dir, paste0(output_prefix, "_clustered_similarity.png"))
grDevices::png(heatmap_path, width = 2600, height = 2300, res = 220)
graphics::par(oma = c(1, 1, 8, 1))
stats::heatmap(similarity,
  Rowv = stats::as.dendrogram(cluster_tree), Colv = "Rowv",
  scale = "none", symm = TRUE, revC = TRUE,
  col = heatmap_colors, margins = c(6, 6),
  labRow = seq_len(D), labCol = seq_len(D),
  RowSideColors = side_colors, ColSideColors = side_colors,
  xlab = "Feature index", ylab = "Feature index", main = "")
graphics::mtext(paste0("Clustered spline R-squared similarity: ", dataset_id),
  outer = TRUE, side = 3, line = 5.0, cex = 1.2, font = 2)
graphics::mtext(sprintf(
  "Directional natural cubic spline on normalized ranks; fixed df = %d; symmetric maximum",
  spline_df), outer = TRUE, side = 3, line = 3.3, cex = 1.0)
graphics::mtext(sprintf(
  "True M = %d; Adaptive EB effective M = %d; Uniform + forward effective M = %d",
  true_M, adaptive_M, forward_M), outer = TRUE, side = 3, line = 1.7, cex = 0.9)
graphics::legend("topright", inset = 0.01, title = "True ordering",
  legend = true_levels, fill = unname(ordering_colors[true_levels]),
  border = NA, bty = "n", cex = 0.75)
similarity_breaks <- seq(0, 1, by = 0.25)
similarity_color_index <- 1L + round(similarity_breaks * (length(heatmap_colors) - 1L))
graphics::legend("bottomright", inset = 0.01, title = "Proportion variance explained",
  legend = format(similarity_breaks, trim = TRUE),
  fill = heatmap_colors[similarity_color_index], border = NA, bty = "n", cex = 0.6)
grDevices::dev.off()

diagnostic_path <- file.path(output_dir, paste0(output_prefix, "_direct_M_diagnostics.png"))
grDevices::png(diagnostic_path, width = 2200, height = 1000, res = 200)
old_par <- graphics::par(mfrow = c(1, 2), mar = c(9.0, 4.8, 3, 4.8))
graphics::plot(diagnostics$k, diagnostics$mean_silhouette, type = "b", pch = 16,
  lwd = 2, col = "#0072B2", ylim = range(c(diagnostics$mean_silhouette,
    diagnostics$ARI_to_truth)), xlab = "Number of clusters k", ylab = "Score",
  main = sprintf("Spline R-squared cuts (df = %d)", spline_df))
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
  sprintf("Natural cubic spline df: %d", spline_df),
  sprintf("Vectorized similarity runtime in seconds: %.6f", unname(runtime[["elapsed"]])),
  sprintf("Single-linkage silhouette estimate: %d", silhouette_k),
  sprintf("Single-linkage merge-height-gap estimate: %d", height_gap_k),
  sprintf("ARI of the k = true M tree cut: %.6f",
    diagnostics$ARI_to_truth[diagnostics$k == true_M]),
  "Directional score: training proportion variance explained by a fixed-df natural cubic spline on normalized predictor ranks",
  "Symmetrization: maximum of the two directional R-squared values",
  "Clustering: single linkage on distance 1 - similarity",
  paste0("Input: ", dataset_path),
  "Truth labels are used only for side colors and evaluation summaries.")
writeLines(summary_lines, file.path(output_dir, paste0(output_prefix, "_README.txt")))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Saved exploratory outputs under", output_dir, "\n")
