source("experiments/isomap_elbo_screen_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
result <- readRDS(file.path(experiment_dir, "narrow_range.rds"))
stopifnot(identical(result$provenance$group_input_sha256,
  digest::digest(observations, algo = "sha256")))
normalized_rank <- function(values)
  (rank(values, ties.method = "average") - 1) / (length(values) - 1)
aligned_rank <- function(values) {
  ranks <- normalized_rank(values)
  if (cor(truth, values, method = "spearman") < 0) 1 - ranks else ranks
}
points <- list()
for (candidate in result$candidates) {
  stages <- list("Raw Isomap" = candidate$raw_position,
    "After one CAVI sweep" = candidate$early$position,
    "At convergence" = candidate$final$position)
  label <- sprintf("k=%d%s", candidate$k,
    if (candidate$k == result$selected_k) " (one-sweep selected)" else "")
  for (stage in names(stages)) {
    points[[length(points) + 1L]] <- data.frame(k_label = label, stage = stage,
      true_rank = normalized_rank(truth), inferred_rank = aligned_rank(stages[[stage]]),
      rho_label = sprintf("|Spearman rho| = %.4f", recovery(stages[[stage]])))
  }
}
points <- do.call(rbind, points)
points$stage <- factor(points$stage,
  levels = c("Raw Isomap", "After one CAVI sweep", "At convergence"))
annotations <- unique(points[c("k_label", "stage", "rho_label")])
caption <- paste("Figure 1. All 300 sample ranks under Isomap k=2,3,4 on the original 12-feature B group.",
  "Rows show raw Isomap, one MPCurve CAVI sweep, and continuation to convergence.",
  "The column marked selected has the highest one-sweep ELBO within [2,4].",
  "Orientation is aligned using truth for display only; the gray diagonal is perfect rank recovery.",
  "All candidates use 50 position bins and identical fitting settings.")
plot <- ggplot(points, aes(true_rank, inferred_rank)) +
  geom_abline(slope = 1, intercept = 0, color = "#999999", linewidth = .4) +
  geom_point(color = "#0072B2", size = .8, alpha = .6) +
  geom_text(data = annotations, aes(x = .035, y = .965, label = rho_label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_grid(stage ~ k_label) + coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "True latent-position rank (normalized to [0,1])",
    y = "Inferred position rank (orientation aligned)",
    title = "Testing the graph-rule interval [2,4]",
    subtitle = sprintf("One-sweep ELBO selects k=%d | Same observations and MPCurver 0.4.0 controls",
      result$selected_k), caption = paste(strwrap(caption, 130), collapse = "\n")) +
  theme_minimal(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 9))
ggsave(file.path(experiment_dir, "narrow_positions.png"), plot,
  width = 12, height = 11, dpi = 160, bg = "white")
ggsave(file.path(experiment_dir, "narrow_positions.pdf"), plot, width = 12, height = 11)
write.csv(points, file.path(experiment_dir, "narrow_plotted_positions.csv"), row.names = FALSE)
