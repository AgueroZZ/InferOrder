#!/usr/bin/env Rscript

# Compare the single-linkage dendrograms induced by two feature similarities.
# This reads previously saved exploratory matrices and does not refit any model.

if (!file.exists("_workflowr.yml")) {
  stop("Run this script from the InferOrder repository root.")
}

study_dir <- file.path("experiments", "estimate_intrinsic_m_smooth_v032")
arguments <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(arguments)) arguments[[1L]] else "main_M5_S1_r001"
output_dir <- file.path(study_dir, "exploratory_similarity", dataset_id)

read_similarity <- function(filename) {
  similarity <- as.matrix(utils::read.csv(
    file.path(output_dir, filename), row.names = 1L, check.names = FALSE))
  storage.mode(similarity) <- "double"
  stopifnot(nrow(similarity) == ncol(similarity), all(is.finite(similarity)))
  similarity
}

read_diagnostics <- function(filename) {
  utils::read.csv(file.path(output_dir, filename), stringsAsFactors = FALSE)
}

spearman_similarity <- read_similarity("absolute_spearman_similarity.csv")
spline_similarity <- read_similarity("spline_r2_df5_similarity.csv")
stopifnot(identical(dim(spearman_similarity), dim(spline_similarity)))

spearman_tree <- stats::hclust(
  stats::as.dist(1 - spearman_similarity), method = "single")
spline_tree <- stats::hclust(
  stats::as.dist(1 - spline_similarity), method = "single")
spearman_diagnostics <- read_diagnostics("direct_M_diagnostics.csv")
spline_diagnostics <- read_diagnostics("spline_r2_df5_direct_M_diagnostics.csv")

favored_k <- function(diagnostics) {
  c(
    silhouette = diagnostics$k[which.max(diagnostics$mean_silhouette)],
    height_gap = diagnostics$k[which.max(diagnostics$height_gap)]
  )
}

cut_height <- function(tree, k) {
  feature_count <- length(tree$order)
  merge_index <- feature_count - k
  mean(tree$height[c(merge_index, merge_index + 1L)])
}

plot_tree <- function(tree, diagnostics, title, subtitle) {
  selected_k <- favored_k(diagnostics)
  graphics::plot(tree, labels = FALSE, hang = -1,
    xlab = "Features", ylab = "Dissimilarity at merge (1 - similarity)",
    main = title, sub = subtitle, lwd = 1.1)
  heights <- vapply(selected_k, function(k) cut_height(tree, k), numeric(1))
  graphics::abline(h = heights[["silhouette"]], col = "#0072B2", lwd = 2.5)
  graphics::abline(h = heights[["height_gap"]], col = "#D55E00", lwd = 2.5,
    lty = 2)
  graphics::legend("topleft", inset = 0.01, bty = "o", bg = "white", lwd = 2.5,
    lty = c(1, 2), col = c("#0072B2", "#D55E00"),
    legend = c(
      sprintf("Best mean silhouette: k = %d", selected_k[["silhouette"]]),
      sprintf("Largest merge-height gap: k = %d", selected_k[["height_gap"]])
    ))
}

output_path <- file.path(output_dir,
  "comparison_spearman_vs_spline_r2_df5_dendrograms.png")
grDevices::png(output_path, width = 3000, height = 1350, res = 220)
old_par <- graphics::par(mfrow = c(1, 2), mar = c(5, 5, 5, 2), oma = c(0, 0, 3, 0))
plot_tree(spearman_tree, spearman_diagnostics,
  "Absolute Spearman similarity", "Single linkage")
plot_tree(spline_tree, spline_diagnostics,
  "Spline R-squared similarity", "Single linkage; natural cubic spline df = 5")
graphics::mtext(sprintf("Feature-similarity dendrograms: %s", dataset_id),
  outer = TRUE, side = 3, line = 1, font = 2, cex = 1.25)
graphics::par(old_par)
grDevices::dev.off()

cat("Saved", output_path, "\n")
cat(sprintf("Spearman favors k = %d by silhouette and k = %d by height gap.\n",
  favored_k(spearman_diagnostics)[["silhouette"]],
  favored_k(spearman_diagnostics)[["height_gap"]]))
cat(sprintf("Spline R-squared favors k = %d by silhouette and k = %d by height gap.\n",
  favored_k(spline_diagnostics)[["silhouette"]],
  favored_k(spline_diagnostics)[["height_gap"]]))
