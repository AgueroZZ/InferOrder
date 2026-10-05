# Compare raw Isomap orderings by color on the same observed feature plane.
source("experiments/isomap_kmin_m1_p2_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
output <- file.path(study, "initialization_diagnostic_rep06")
samples <- read.csv(file.path(output, "sample_coordinates.csv"))
scores <- read.csv(file.path(output, "scores.csv"))
input <- readRDS(input_path(6L))
stopifnot(identical(input$input_sha256, hash_object(input$X)))
labels <- c("True positions\n(reference)", sprintf("%s (k=%d)\nRaw |Spearman| = %.3f",
  c("Auto kmin", "Fixed", "Fixed"), scores$k_used, scores$raw_absolute_spearman))
reference <- samples[samples$method == "auto_kmin", ]
points <- data.frame(sample = reference$sample, panel = labels[1L],
  feature_1 = reference$feature_1, feature_2 = reference$feature_2,
  color_position = reference$truth)
for (index in seq_along(design$methods)) {
  method <- design$methods[index]
  selected <- samples[samples$method == method, ]
  fit <- readRDS(result_path(6L, method))
  stopifnot(identical(selected$sample, rownames(input$X)),
    max(abs(selected$feature_1 - input$X[, 1L])) < 1e-12,
    max(abs(selected$feature_2 - input$X[, 2L])) < 1e-12,
    max(abs(selected$truth - input$truth)) < 1e-12,
    max(abs(selected$raw_position - fit$raw_positions)) < 1e-12)
  expected <- if (scores$initial_display_reversed[index]) 1 - fit$raw_positions else fit$raw_positions
  stopifnot(max(abs(selected$raw_display_position - expected)) < 1e-12)
  points <- rbind(points, data.frame(sample = selected$sample, panel = labels[index + 1L],
    feature_1 = selected$feature_1, feature_2 = selected$feature_2,
    color_position = selected$raw_display_position))
}
points$panel <- factor(points$panel, levels = labels)
write.csv(points, file.path(output, "plotted_feature_plane_points.csv"), row.names = FALSE)
caption <- paste("Figure 4. The same 200 observed samples with feature 1 on the horizontal axis and feature 2 on the vertical axis.",
  "Colors show true positions in the reference panel and raw Isomap positions in the other panels, before quantile binning",
  "or MPCurve fitting. Automatic k=3 positions are globally reflected for display; fixed-k positions retain their original direction.",
  "All panels share axis limits and a [0,1] color scale. No rank transformation or nonlinear alignment is applied.")
plot <- ggplot(points, aes(feature_1, feature_2, color = color_position)) +
  geom_point(size = 1.7, alpha = .85) + facet_wrap(~ panel, nrow = 1) +
  coord_equal() + scale_color_viridis_c(name = "Position along ordering", limits = c(0, 1),
    breaks = seq(0, 1, .25), option = "D", end = .95) +
  labs(x = "Observed feature 1", y = "Observed feature 2",
    title = "Same feature scatter, different initial orderings: replicate 6",
    subtitle = "X and Y are the observed features in every panel | Color shows true or initial Isomap position",
    caption = paste(strwrap(caption, 142), collapse = "\n")) + theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9), legend.position = "bottom",
    legend.key.width = grid::unit(1.5, "cm"))
for (extension in c("png", "pdf")) ggsave(
  file.path(output, paste0("feature_plane_comparison.", extension)), plot,
  width = 15, height = 5.8, dpi = 180, bg = "white")
record <- list(replication = 6L, created_at_utc = format(Sys.time(), tz = "UTC"),
  script_sha256 = hash_file(file.path(study, "plot_feature_plane_initializations.R")),
  coordinate_csv_sha256 = hash_file(file.path(output, "sample_coordinates.csv")),
  score_csv_sha256 = hash_file(file.path(output, "scores.csv")),
  parent_provenance = readRDS(file.path(output, "provenance.rds")),
  displayed_reversal = scores[, c("method", "initial_display_reversed")],
  plotted_rows = nrow(points), shared_axes = TRUE, color_limits = c(0, 1))
save_atomic(record, file.path(output, "feature_plane_provenance.rds"))
cat("Verified all 800 plotted rows against saved observations and original Isomap positions.\n")
