#!/usr/bin/env Rscript

# Render public-report figures using only the MPCurver 0.3.4 automatic-M fits.

source("experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

summary_dir <- file.path(extension_dir, "summary")
runs <- read.csv(file.path(summary_dir, "auto_runs.csv"))
stopifnot(nrow(runs) == 90L, all(runs$status == "success"),
  all(runs$method == "auto_adaptive"), sum(runs$warning_count) == 0L)

theme_set(theme_minimal(base_size = 12))
save_plot <- function(plot, name, width, height) {
  for (extension in c("png", "pdf")) {
    ggsave(file.path(summary_dir, paste0(name, ".", extension)), plot,
      width = width, height = height, dpi = 180, bg = "white")
  }
}

wilson <- function(x, n) {
  z <- stats::qnorm(0.975)
  p <- x / n
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, center - half), upper = min(1, center + half))
}

condition_summary <- do.call(rbind, lapply(split(runs,
  interaction(runs$true_M, runs$snr, drop = TRUE)), function(x) {
  exact <- sum(x$effective_M == x$true_M)
  interval <- wilson(exact, nrow(x))
  data.frame(
    true_M = x$true_M[1L],
    snr = x$snr[1L],
    datasets = nrow(x),
    initial_exact = sum(x$selected_initial_M == x$true_M),
    final_exact = exact,
    final_under = sum(x$effective_M < x$true_M),
    final_over = sum(x$effective_M > x$true_M),
    exact_rate = exact / nrow(x),
    lower = interval[1L],
    upper = interval[2L],
    mean_ARI = mean(x$ARI),
    mean_ordering_recovery = mean(x$ordering_recovery),
    median_seconds = median(x$elapsed_seconds)
  )
}))
rownames(condition_summary) <- NULL
write.csv(condition_summary,
  file.path(summary_dir, "current_condition_summary.csv"), row.names = FALSE)

condition_summary$snr_label <- factor(condition_summary$snr, levels = design$snr)
accuracy <- ggplot(condition_summary, aes(snr_label, exact_rate, group = 1)) +
  geom_line(color = "#009E73") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12,
    color = "#009E73") +
  geom_point(size = 2.7, color = "#009E73") +
  facet_wrap(~true_M, nrow = 1,
    labeller = label_bquote(M[true] == .(true_M))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Variance signal-to-noise ratio",
    y = "Exact effective-M recovery rate",
    title = "Ordering-count recovery with MPCurver 0.3.4",
    subtitle = "Each point summarizes 10 automatic-M adaptive fits; bars are 95% Wilson intervals"
  ) +
  theme(panel.grid.minor = element_blank())
save_plot(accuracy, "current_exact_recovery", 10.5, 4.5)

stage_frame <- rbind(
  data.frame(true_M = runs$true_M, snr = runs$snr,
    stage = "Similarity cut", exact = runs$selected_initial_M == runs$true_M),
  data.frame(true_M = runs$true_M, snr = runs$snr,
    stage = "Adaptive EB final", exact = runs$effective_M == runs$true_M)
)
stage_summary <- aggregate(exact ~ true_M + snr + stage, stage_frame, mean)
stage_summary$snr_label <- factor(stage_summary$snr, levels = design$snr)
stage_summary$stage <- factor(stage_summary$stage,
  levels = c("Similarity cut", "Adaptive EB final"))
stage_plot <- ggplot(stage_summary,
  aes(snr_label, exact, color = stage, group = stage)) +
  geom_line() + geom_point(size = 2.6) +
  facet_wrap(~true_M, nrow = 1,
    labeller = label_bquote(M[true] == .(true_M))) +
  scale_color_manual(values = c("Similarity cut" = "#0072B2",
    "Adaptive EB final" = "#009E73")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Variance signal-to-noise ratio",
    y = "Exact-M recovery rate",
    color = NULL,
    title = "Similarity initialization and final adaptive estimate"
  ) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(stage_plot, "current_initial_vs_final", 10.5, 4.5)

distribution_groups <- split(runs,
  interaction(runs$true_M, runs$snr, drop = TRUE))
distribution <- do.call(rbind, lapply(distribution_groups, function(x) {
  levels <- as.character(seq_len(design$max_intrinsic_dim))
  counts <- table(factor(as.character(x$effective_M), levels = levels))
  data.frame(
    true_M = x$true_M[1L],
    snr = x$snr[1L],
    estimated = levels,
    count = as.integer(counts),
    probability = as.integer(counts) / nrow(x)
  )
}))
distribution$estimated <- factor(distribution$estimated,
  levels = as.character(seq_len(design$max_intrinsic_dim)))
distribution$snr_label <- factor(paste0("SNR = ", distribution$snr),
  levels = paste0("SNR = ", design$snr))
distribution$truth <- factor(distribution$true_M, levels = rev(design$true_M))
heatmap <- ggplot(distribution, aes(estimated, truth, fill = probability)) +
  geom_tile(color = "white") +
  geom_tile(data = distribution[
    as.character(distribution$estimated) == as.character(distribution$true_M), ],
    fill = NA, color = "#333333", linewidth = 0.7) +
  geom_text(aes(label = ifelse(count > 0, count, "")), size = 3.6) +
  facet_wrap(~snr_label, nrow = 1) +
  scale_fill_gradient(low = "#FFFFFF", high = "#56B4E9", limits = c(0, 1)) +
  labs(
    x = "Estimated effective M",
    y = "True M",
    fill = "Proportion",
    title = "Distribution of estimated ordering counts",
    subtitle = "Cell labels give counts out of 10; outlined cells recover the true M"
  ) +
  theme(panel.grid = element_blank())
save_plot(heatmap, "current_estimated_M_distribution", 11, 4.4)
write.csv(distribution,
  file.path(summary_dir, "current_estimated_M_distribution.csv"), row.names = FALSE)

metrics <- rbind(
  data.frame(condition_summary,
    metric = "Feature partition: mean ARI", value = condition_summary$mean_ARI),
  data.frame(condition_summary,
    metric = "Ordering recovery: mean |Spearman|",
    value = condition_summary$mean_ordering_recovery)
)
quality <- ggplot(metrics, aes(snr_label, value, group = 1)) +
  geom_line(color = "#009E73") +
  geom_point(size = 2.6, color = "#009E73") +
  facet_grid(metric ~ true_M) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Variance signal-to-noise ratio",
    y = "Mean recovery score",
    title = "Recovery of feature groups and sample orderings"
  ) +
  theme(panel.grid.minor = element_blank())
save_plot(quality, "current_structural_recovery", 11.5, 6.2)

timing <- ggplot(condition_summary, aes(snr_label, median_seconds, group = 1)) +
  geom_line(color = "#009E73") +
  geom_point(size = 2.6, color = "#009E73") +
  facet_wrap(~true_M, nrow = 1,
    labeller = label_bquote(M[true] == .(true_M))) +
  labs(
    x = "Variance signal-to-noise ratio",
    y = "Median fitting time (seconds)",
    title = "Runtime of automatic-M adaptive fitting"
  ) +
  theme(panel.grid.minor = element_blank())
save_plot(timing, "current_runtime", 10.5, 4.3)

representatives <- manifest[manifest$snr == 4 & manifest$replicate == 1L, ]
stopifnot(nrow(representatives) == 3L)

for (i in seq_len(nrow(representatives))) {
  row <- representatives[i, , drop = FALSE]
  dataset <- load_fixed_dataset(row)
  result <- readRDS(file.path(extension_dir, "results",
    paste0(row$id, "_auto_adaptive.rds")))
  stopifnot(result$row$status == "success", result$provenance$package_version == "0.3.4")

  positions <- result$compact$positions[, result$compact$active, drop = FALSE]
  correlation <- abs(stats::cor(dataset$latent_positions, positions,
    method = "spearman"))
  size <- max(ncol(dataset$latent_positions), ncol(positions))
  score <- matrix(0, size, size)
  score[seq_len(ncol(dataset$latent_positions)), seq_len(ncol(positions))] <-
    correlation
  matching <- as.integer(clue::solve_LSAP(score, maximum = TRUE))[
    seq_len(ncol(dataset$latent_positions))]

  panels <- list()
  for (m in seq_len(row$true_M)) {
    columns <- head(which(dataset$true_assign == dataset$ordering_labels[m]), 3L)
    frame <- do.call(rbind, lapply(columns, function(j) data.frame(
      position = dataset$latent_positions[, m],
      signal = dataset$signal[, j],
      feature = factor(which(columns == j), levels = 1:3,
        labels = c("Monotone anchor", "Smooth trajectory 1",
          "Smooth trajectory 2"))
    )))
    panels[[length(panels) + 1L]] <- ggplot(frame,
      aes(position, signal, color = feature)) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = c("#0072B2", "#D55E00", "#7B4F9D")) +
      labs(title = paste("True ordering", dataset$ordering_labels[m]),
        x = "True position", y = "Noiseless feature signal") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
  }

  for (m in seq_len(row$true_M)) {
    matched <- matching[m] <= ncol(positions)
    if (matched) {
      x <- dataset$latent_positions[, m]
      y <- positions[, matching[m]]
      signed_rho <- stats::cor(x, y, method = "spearman")
      if (signed_rho < 0) y <- 1 - y
      panel <- ggplot(data.frame(x = x, y = y), aes(x, y)) +
        geom_point(alpha = 0.4, size = 0.8, color = "#009E73") +
        labs(title = sprintf("MPCurver 0.3.4 auto-M\n|rho| = %.3f",
          abs(signed_rho)))
    } else {
      panel <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Unmatched ordering") +
        labs(title = "MPCurver 0.3.4 auto-M")
    }
    panels[[length(panels) + 1L]] <- panel +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "True position", y = "Inferred position") +
      theme_minimal(base_size = 10) +
      theme(panel.grid.minor = element_blank())
  }

  figure <- wrap_plots(panels, ncol = row$true_M) +
    plot_annotation(
      title = sprintf("Current-package example: true M = %d", row$true_M),
      subtitle = sprintf(
        "MPCurver 0.3.4 auto-M adaptive fit; SNR = 4, replicate 1; initial M = %d, effective M = %d.",
        result$row$selected_initial_M, result$row$effective_M),
      caption = paste0(
        "Top row: one monotone anchor (blue) and two nonmonotone smooth trajectories per true ordering.\n",
        "Bottom row: one-to-one matched inferred positions, with reversal allowed."
      )
    )
  save_plot(figure, sprintf("current_example_M%d", row$true_M),
    max(10.5, 2.8 * row$true_M), 6.7)
}

cat("Rendered current-package report figures.\n")
print(condition_summary, row.names = FALSE)
