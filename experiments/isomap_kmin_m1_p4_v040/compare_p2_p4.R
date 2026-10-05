# Compare the paired feature additions and inspect the previously selected replicate.
source("experiments/isomap_kmin_m1_p4_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
rows <- do.call(rbind, lapply(c(2L, 4L), function(P) {
  root <- if (P == 2L) parent_study else study
  do.call(rbind, lapply(seq_len(design$replications), function(replication)
    do.call(rbind, lapply(design$methods, function(method) {
      fit <- readRDS(file.path(root, "results", sprintf("rep%02d_%s.rds", replication, method)))
      input <- readRDS(file.path(root, "inputs", sprintf("rep%02d.rds", replication)))
      stopifnot(identical(fit$input_sha256, input$input_sha256),
        identical(fit$metrics, position_metrics(input$truth, fit$positions)),
        identical(fit$raw_metrics, position_metrics(input$truth, fit$raw_positions)))
      data.frame(replication = replication, P = P, method = method, k_used = fit$k_used,
        converged = fit$converged, status = fit$status,
        raw_spearman = fit$raw_metrics[["spearman"]], final_spearman = fit$metrics[["spearman"]])
    }))))
}))
write.csv(rows, file.path(study, "paired_dimension_metrics.csv"), row.names = FALSE)
differences <- do.call(rbind, lapply(design$methods, function(method) {
  old <- rows[rows$P == 2L & rows$method == method, ]
  added <- rows[rows$P == 4L & rows$method == method, ]
  matched <- match(old$replication, added$replication)
  do.call(rbind, lapply(c("raw", "final"), function(stage) {
    column <- paste0(stage, "_spearman")
    data.frame(replication = old$replication, method = method, stage = stage,
      p2_score = old[[column]], p4_score = added[[column]][matched],
      difference = added[[column]][matched] - old[[column]])
  }))
}))
write.csv(differences, file.path(study, "paired_dimension_differences.csv"), row.names = FALSE)
summaries <- do.call(rbind, lapply(split(differences, interaction(differences$method, differences$stage)), function(data) {
  finite <- is.finite(data$difference)
  scores <- data$difference[finite]
  data.frame(method = data$method[1L], stage = data$stage[1L],
    complete_pairs = length(scores), missing_pairs = sum(!finite),
    p2_median = median(data$p2_score, na.rm = TRUE), p4_median = median(data$p4_score, na.rm = TRUE),
    p2_mean = mean(data$p2_score, na.rm = TRUE), p4_mean = mean(data$p4_score, na.rm = TRUE),
    mean_difference = mean(scores), median_difference = median(scores),
    monte_carlo_se = sd(scores) / sqrt(length(scores)),
    improvements = sum(scores > 1e-10), ties = sum(abs(scores) <= 1e-10), declines = sum(scores < -1e-10))
}))
write.csv(summaries, file.path(study, "paired_dimension_summary.csv"), row.names = FALSE)

save_plot <- function(plot, stem, width, height) {
  for (extension in c("png", "pdf")) ggsave(file.path(study, paste0(stem, ".", extension)),
    plot, width = width, height = height, dpi = 180, bg = "white")
}
caption <- function(text, width = 120) paste(strwrap(text, width), collapse = "\n")
theme <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(size = 11, face = "bold"),
  plot.caption = element_text(hjust = 0, size = 9))
labels <- c(auto_kmin = "Auto kmin", fixed_k15 = "Fixed k=15", fixed_k10 = "Fixed k=10")
plot_data <- rows
plot_data$method <- factor(plot_data$method, levels = design$methods, labels = labels[design$methods])
plot_data$features <- factor(plot_data$P, levels = c(2L, 4L), labels = c("P=2", "P=4"))
write.csv(plot_data, file.path(study, "plotted_dimension_scores.csv"), row.names = FALSE)
boxplot <- ggplot(plot_data, aes(features, final_spearman)) +
  geom_line(aes(group = replication), color = "#888888", alpha = .45, linewidth = .4, na.rm = TRUE) +
  geom_boxplot(aes(fill = features), alpha = .4, width = .45, outlier.shape = NA, na.rm = TRUE) +
  geom_point(aes(color = features), size = 1.8, alpha = .85, na.rm = TRUE) +
  facet_wrap(~ method, nrow = 1) +
  scale_color_manual(values = c("P=2" = "#0072B2", "P=4" = "#D55E00"), guide = "none") +
  scale_fill_manual(values = c("P=2" = "#0072B2", "P=4" = "#D55E00"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "Number of observed signal features", y = "Final absolute Spearman correlation",
    title = "Position recovery after adding two trajectories",
    subtitle = "20 paired replicates | Original two features, sample positions, and noise preserved",
    caption = caption(paste("Figure 1. Converged MPCurve position recovery with P=2 or P=4, N=200 and M=1.",
      "P=4 appends two independent nonmonotone cubic B-spline features; average signal variance remains one and noise SD=0.25.",
      "Gray lines connect the same replication. Scores allow global reversal; boxes show median/quartiles and whiskers extend",
      "to the most extreme values within 1.5 IQR. All 20 endpoints per feature-count/method combination are included."))) + theme
save_plot(boxplot, "paired_p2_p4_spearman", 12, 6.5)

# The inspection case stays replicate 6, selected before any P=4 outcomes were available.
replication <- 6L
point_rows <- list()
case_scores <- rows[rows$replication == replication, ]
for (P in c(2L, 4L)) {
  root <- if (P == 2L) parent_study else study
  input <- readRDS(file.path(root, "inputs", sprintf("rep%02d.rds", replication)))
  if (P == 4L) stopifnot(identical(input$X[, 1:2, drop = FALSE], parent_input(replication)$X))
  reference <- sprintf("P=%d: True positions\n(reference)", P)
  point_rows[[length(point_rows) + 1L]] <- data.frame(sample = rownames(input$X),
    panel = reference, feature_1 = input$X[, 1L], feature_2 = input$X[, 2L], color_position = input$truth)
  for (method in design$methods) {
    fit <- readRDS(file.path(root, "results", sprintf("rep%02d_%s.rds", replication, method)))
    raw <- fit$raw_positions
    if (cor(input$truth, raw, method = "spearman") < 0) raw <- 1 - raw
    label <- sprintf("P=%d: %s (k=%d)\nRaw |Spearman| = %.3f", P, labels[[method]],
      fit$k_used, fit$raw_metrics[["spearman"]])
    point_rows[[length(point_rows) + 1L]] <- data.frame(sample = rownames(input$X), panel = label,
      feature_1 = input$X[, 1L], feature_2 = input$X[, 2L], color_position = raw)
  }
}
case_points <- do.call(rbind, point_rows)
case_points$panel <- factor(case_points$panel, levels = unique(case_points$panel))
write.csv(case_points, file.path(study, "plotted_rep06_feature_plane.csv"), row.names = FALSE)
write.csv(case_scores, file.path(study, "rep06_scores.csv"), row.names = FALSE)
position_colors <- scale_color_viridis_c(name = "Position along ordering", limits = c(0, 1),
  option = "D", end = .95, breaks = seq(0, 1, .25))
case_plot <- ggplot(case_points, aes(feature_1, feature_2, color = color_position)) +
  geom_point(size = 1.5, alpha = .85) + facet_wrap(~ panel, ncol = 4) + coord_equal() + position_colors +
  labs(x = "Observed feature 1", y = "Observed feature 2", title = "Same feature plane, more information for Isomap: replicate 6",
    subtitle = "Top: Isomap uses two features | Bottom: Isomap uses four features | X and Y remain features 1 and 2",
    caption = caption(paste("Figure 2. Identical observed feature-1/feature-2 scatter in all eight panels.",
      "Colors show true positions or raw Isomap positions computed from the indicated total number of features, before MPCurve fitting.",
      "Only a global reflection makes signed rank correlation positive for display. All panels share axes and color limits.",
      "The case was chosen in the preceding P=2 diagnostic; no new case selection uses P=4 outcomes."), 140)) +
  theme + theme(legend.position = "bottom", legend.key.width = grid::unit(1.5, "cm"))
save_plot(case_plot, "rep06_feature_plane_p2_p4", 15, 9)

input <- readRDS(input_path(replication))
pairs <- combn(seq_len(design$P), 2L, simplify = FALSE)
projections <- do.call(rbind, lapply(pairs, function(pair) data.frame(
  pair = sprintf("Feature %d (x) / Feature %d (y)", pair[1L], pair[2L]),
  x = input$X[, pair[1L]], y = input$X[, pair[2L]], truth = input$truth)))
curves <- do.call(rbind, lapply(pairs, function(pair) data.frame(
  pair = sprintf("Feature %d (x) / Feature %d (y)", pair[1L], pair[2L]),
  x = input$dense_signal[, pair[1L]], y = input$dense_signal[, pair[2L]], truth = input$grid)))
write.csv(projections, file.path(study, "plotted_rep06_projections.csv"), row.names = FALSE)
projection_plot <- ggplot(projections, aes(x, y)) +
  geom_path(data = curves, color = "#555555", linewidth = .45) +
  geom_point(aes(color = truth), size = 1.5, alpha = .85) + facet_wrap(~ pair, ncol = 3) +
  coord_equal() + scale_color_viridis_c(name = "True sample position", limits = c(0, 1), option = "D", end = .95) +
  labs(x = "First feature named in panel", y = "Second feature named in panel",
    title = "Six two-feature projections of the four-feature curve: replicate 6",
    subtitle = "Colors: true positions | Gray lines: noiseless generating trajectories",
    caption = caption(paste("Figure 3. All six pairwise projections of the same P=4 observations, with named horizontal/vertical features.",
      "The first two features and their noise are unchanged from P=2. The added independent trajectories provide extra joint information",
      "even when individual two-feature projections have crossings or close approaches. Axes share observed-feature units."))) +
  theme + theme(legend.position = "bottom", legend.key.width = grid::unit(1.5, "cm"))
save_plot(projection_plot, "rep06_six_feature_projections", 12, 9)
save_atomic(list(primary_metric = "spearman", summaries = summaries, case_scores = case_scores,
  paired_metric_rows = nrow(rows), comparison_rows = nrow(differences),
  evaluated_at_utc = format(Sys.time(), tz = "UTC"), script_sha256 = hash_file(file.path(study, "compare_p2_p4.R"))),
  file.path(study, "paired_dimension_evaluation.rds"))
print(summaries, digits = 6)
print(case_scores, digits = 6)
