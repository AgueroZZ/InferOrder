# Inspect a selected automatic-kmin failure using the saved, pre-fit Isomap coordinates.
source("experiments/isomap_kmin_m1_p2_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
metadata <- readRDS(file.path(study, "provenance.rds"))
current <- provenance()
stopifnot(identical(metadata$design_hash, current$design_hash),
  identical(metadata$package_source_hashes, current$package_source_hashes),
  identical(metadata$script_hashes, current$script_hashes))

# This is an outcome-selected illustration, separate from the complete paired study.
metrics <- read.csv(file.path(study, "metrics.csv"))
selection <- do.call(rbind, lapply(seq_len(design$replications), function(replication) {
  rows <- metrics[metrics$replication == replication, ]
  automatic <- rows$spearman[rows$method == "auto_kmin"]
  fixed_mean <- mean(rows$spearman[rows$method != "auto_kmin"])
  data.frame(replication = replication, automatic_final_spearman = automatic,
    fixed_mean_final_spearman = fixed_mean, fixed_minus_auto = fixed_mean - automatic,
    eligible = is.finite(automatic) && automatic < .2)
}))
eligible <- selection[selection$eligible, ]
stopifnot(nrow(eligible) > 0L)
replication <- eligible$replication[order(-eligible$fixed_minus_auto, eligible$replication)][1L]
output <- file.path(study, sprintf("initialization_diagnostic_rep%02d", replication))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(selection, file.path(output, "case_selection.csv"), row.names = FALSE)
input <- readRDS(input_path(replication))
fits <- setNames(lapply(design$methods, function(method)
  readRDS(result_path(replication, method))), design$methods)
stopifnot(identical(input$input_sha256, hash_object(input$X)),
  identical(input, simulate_replicate(replication)))
summaries <- list()
coordinates <- list()
for (method in design$methods) {
  fit <- fits[[method]]
  stopifnot(identical(fit$input_sha256, input$input_sha256),
    identical(fit$package_source_hashes, current$package_source_hashes),
    identical(fit$metrics, position_metrics(input$truth, fit$positions)),
    identical(fit$raw_metrics, position_metrics(input$truth, fit$raw_positions)),
    !fit$initialization$fallback, fit$initialization$method_used == "isomap")
  # Independently reproduce the initializer, with no MPCurve fitting or optimization.
  raw <- do.call(isomap_ordering, c(list(X = input$X), method_arguments(method, fit$seed)))
  stopifnot(max(abs(raw$t - fit$raw_positions)) < 1e-12,
    raw$k_used == fit$k_used, raw$n_components == 1L,
    length(raw$keep_idx) == design$N)
  initial_signed <- cor(input$truth, fit$raw_positions, method = "spearman")
  final_signed <- cor(input$truth, fit$positions, method = "spearman")
  summaries[[method]] <- data.frame(method = method, k_used = fit$k_used,
    raw_signed_spearman = initial_signed, raw_absolute_spearman = abs(initial_signed),
    final_signed_spearman = final_signed, final_absolute_spearman = abs(final_signed),
    initial_display_reversed = initial_signed < 0, final_display_reversed = final_signed < 0,
    graph_components = raw$n_components, retained_samples = length(raw$keep_idx))
  coordinates[[method]] <- data.frame(sample = rownames(input$X), replication = replication,
    method = method, truth = input$truth, raw_position = fit$raw_positions,
    raw_display_position = if (initial_signed < 0) 1 - fit$raw_positions else fit$raw_positions,
    final_position = fit$positions,
    final_display_position = if (final_signed < 0) 1 - fit$positions else fit$positions,
    feature_1 = input$X[, 1L], feature_2 = input$X[, 2L])
}
scores <- do.call(rbind, summaries)
row.names(scores) <- NULL
samples <- do.call(rbind, coordinates)
row.names(samples) <- NULL
write.csv(scores, file.path(output, "scores.csv"), row.names = FALSE)
write.csv(samples, file.path(output, "sample_coordinates.csv"), row.names = FALSE)
dense <- data.frame(truth = input$grid, feature_1 = input$dense_signal[, 1L],
  feature_2 = input$dense_signal[, 2L])
write.csv(dense, file.path(output, "true_curve.csv"), row.names = FALSE)

truth_colors <- scale_color_viridis_c(name = "True sample position", limits = c(0, 1),
  breaks = seq(0, 1, .25), option = "D", end = .95)
plot_theme <- theme_minimal(base_size = 12) + theme(
  panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
  strip.text = element_text(size = 10, face = "bold"),
  legend.position = "bottom", legend.key.width = grid::unit(1.5, "cm"))
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(output, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 180, bg = "white")
  ggsave(file.path(output, paste0(stem, ".pdf")), plot, width = width, height = height)
}
caption <- function(text, width = 108) paste(strwrap(text, width), collapse = "\n")
geometry_points <- data.frame(truth = input$truth,
  feature_1 = input$X[, 1L], feature_2 = input$X[, 2L])
geometry <- ggplot(geometry_points, aes(feature_1, feature_2)) +
  geom_path(data = dense, color = "#333333", linewidth = .65) +
  geom_point(aes(color = truth), size = 2, alpha = .85) + truth_colors +
  coord_equal() + labs(x = "Observed feature 1", y = "Observed feature 2",
    title = sprintf("Observed feature geometry: replicate %d", replication),
    subtitle = "Points: noisy observations | Black line: noiseless generating curve",
    caption = caption(paste("Figure 1. Two observed features for the selected N=200, P=2 dataset.",
      "Colors identify true positions along the shared generating curve; the black line traces",
      "that curve in increasing true-position order. Independent observation-noise SD is 0.25."), 83)) + plot_theme
save_figure(geometry, "observed_feature_scatter", 8, 7)

method_labels <- setNames(sprintf("%s (k=%d)\nRaw |Spearman| = %.3f",
  c("Auto kmin", "Fixed", "Fixed"), scores$k_used, scores$raw_absolute_spearman), design$methods)
choice_levels <- c("True positions\n(reference)", unname(method_labels))
feature_points <- do.call(rbind, lapply(seq_len(design$P), function(feature) {
  reference <- data.frame(sample = rownames(input$X), feature = paste("Feature", feature),
    choice = choice_levels[1L], position = input$truth, value = input$X[, feature], truth = input$truth)
  ordered <- do.call(rbind, lapply(design$methods, function(method) {
    points <- samples[samples$method == method, ]
    data.frame(sample = points$sample, feature = paste("Feature", feature),
      choice = method_labels[[method]], position = points$raw_display_position,
      value = input$X[, feature], truth = input$truth)
  }))
  rbind(reference, ordered)
}))
feature_points$choice <- factor(feature_points$choice, levels = choice_levels)
feature_lines <- do.call(rbind, lapply(seq_len(design$P), function(feature)
  data.frame(feature = paste("Feature", feature),
    choice = factor(choice_levels[1L], levels = choice_levels),
    position = input$grid, value = input$dense_signal[, feature])))
write.csv(feature_points, file.path(output, "plotted_feature_points.csv"), row.names = FALSE)
feature_plot <- ggplot(feature_points, aes(position, value)) +
  geom_line(data = feature_lines, color = "#333333", linewidth = .65) +
  geom_point(aes(color = truth), size = 1.5, alpha = .8) +
  facet_grid(feature ~ choice, scales = "free_y") + truth_colors +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  labs(x = "True position (first column) or raw Isomap position (remaining columns)",
    y = "Observed feature value", title = sprintf("Feature-by-feature Isomap initialization: replicate %d", replication),
    subtitle = "Same samples and colors in every column | No MPCurve iterations in the Isomap columns",
    caption = caption(paste("Figure 2. Each row shows one observed feature versus true positions or saved raw Isomap positions.",
      "The true-position column overlays the noiseless generating trajectory. Raw coordinates are before quantile",
      "binning and MPCurve optimization, scaled to [0,1] by the initializer. A global reflection, if needed,",
      "makes rank correlation positive for display; no rank transformation or nonlinear alignment is applied."), 133)) + plot_theme
save_figure(feature_plot, "feature_by_position", 14, 8.2)

samples$choice <- factor(method_labels[samples$method], levels = unname(method_labels))
initial_plot <- ggplot(samples, aes(truth, raw_display_position)) +
  geom_abline(slope = 1, intercept = 0, color = "#777777", linetype = "dashed") +
  geom_point(aes(color = truth), size = 1.7, alpha = .85) +
  facet_wrap(~ choice, nrow = 1) + truth_colors +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  coord_fixed() + labs(x = "True sample position", y = "Raw Isomap position",
    title = sprintf("Original Isomap estimates: replicate %d", replication),
    subtitle = "Coordinates before MPCurve fitting | All 200 samples retained in every graph",
    caption = caption(paste("Figure 3. Saved initial Isomap coordinates against true positions for automatic kmin, k=15, and k=10.",
      "A global reflection allows ordering reversal for display. The dashed diagonal marks equal coordinates;",
      "Spearman assesses rank agreement rather than distance from that diagonal. Colors identify true positions."), 128)) + plot_theme
save_figure(initial_plot, "initial_positions_vs_truth", 13, 5.5)

record <- list(replication = replication,
  selection_rule = "Among automatic final |Spearman| < 0.2, largest mean fixed-minus-auto final score; ties by replicate index.",
  selection = selection, scores = scores, input_sha256 = input$input_sha256,
  generating_seed = input$seed, fitting_seed = fits[[1L]]$seed,
  input_file_sha256 = hash_file(input_path(replication)),
  fit_file_sha256 = setNames(vapply(design$methods, function(method)
    hash_file(result_path(replication, method)), character(1)), design$methods),
  diagnostic_script_sha256 = hash_file(file.path(study, "diagnose_initialization.R")),
  parent_provenance = metadata, package_source_hashes = current$package_source_hashes,
  created_at_utc = format(Sys.time(), tz = "UTC"))
save_atomic(record, file.path(output, "provenance.rds"))
writeLines(c("Verified saved inputs and all three initializer/final metrics.",
  "Reproduced all three raw Isomap coordinate vectors within 1e-12.",
  "All three graphs are connected and retain all 200 samples; no PCA fallback.",
  "No MPCurve fitting was rerun; parent results and plots were preserved."),
  file.path(output, "verification.txt"))
print(scores, digits = 6)
cat("Saved diagnostic figures and plotted data in ", output, ".\n", sep = "")
