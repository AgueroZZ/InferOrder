#!/usr/bin/env Rscript
# Illustrate every simulation condition from saved replicate-one inputs.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
suppressPackageStartupMessages(library(ggplot2))
manifest <- active_manifest()
examples <- manifest[manifest$phase == "main" & manifest$replicate == 1L, ]
stopifnot(nrow(examples) == 9L)
out_dir <- file.path(study_dir, "main_summary")
dir.create(out_dir, showWarnings = FALSE)
frames <- list()
selection <- list()
for (i in seq_len(nrow(examples))) {
  row <- examples[i, ]
  dataset <- readRDS(file.path(study_dir, "data", paste0(row$id, ".rds")))
  stopifnot(dataset$input_hash == digest::digest(dataset$X, algo = "sha256"))
  for (m in seq_len(row$true_M)) {
    columns <- head(which(dataset$true_assign == dataset$ordering_labels[m]), 3L)
    stopifnot(length(columns) == 3L)
    for (k in seq_along(columns)) {
      j <- columns[k]
      frames[[length(frames) + 1L]] <- data.frame(
        dataset_id = row$id, true_M = row$true_M, snr = row$snr,
        ordering = dataset$ordering_labels[m], example_feature = k,
        feature = colnames(dataset$X)[j],
        role = if (k == 1L) "Monotone anchor" else paste("Smooth trajectory", k - 1L),
        position = dataset$latent_positions[, m],
        signal = dataset$signal[, j], observed = dataset$X[, j])
      selection[[length(selection) + 1L]] <- data.frame(
        dataset_id = row$id, seed = row$seed, true_M = row$true_M,
        snr = row$snr, replicate = row$replicate,
        ordering = dataset$ordering_labels[m], example_feature = k,
        feature_index = j, feature = colnames(dataset$X)[j],
        role = if (k == 1L) "Monotone anchor" else "Smooth nonmonotone",
        input_hash = dataset$input_hash)
    }
  }
}
frame <- do.call(rbind, frames)
frame$role <- factor(frame$role,
  levels = c("Monotone anchor", "Smooth trajectory 1", "Smooth trajectory 2"))
frame$snr_label <- factor(paste0("SNR = ", frame$snr),
                         levels = paste0("SNR = ", design$snr))
frame$ordering_label <- paste("Ordering", frame$ordering)
limit <- ceiling(max(abs(c(frame$signal, frame$observed))))
palette <- c("Monotone anchor" = "#0072B2", "Smooth trajectory 1" = "#D55E00",
             "Smooth trajectory 2" = "#7B4F9D")
for (M in design$true_M) {
  selected <- frame[frame$true_M == M, ]
  plot <- ggplot(selected, aes(position, color = role)) +
    geom_point(aes(y = observed), size = 0.55, alpha = 0.20, stroke = 0) +
    geom_line(aes(y = signal, group = role), linewidth = 0.9) +
    facet_grid(ordering_label ~ snr_label) +
    scale_color_manual(values = palette, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1),
                       expand = expansion(mult = 0.02)) +
    coord_cartesian(ylim = c(-limit, limit)) +
    labs(title = sprintf("One monotone anchor and smooth trajectories: true M = %d", M),
      subtitle = sprintf("300 samples; %d features per ordering. Lines: true signals. Points: noisy observations.", design$D / M),
      x = "True latent position within each ordering", y = "Feature value",
      caption = "Replicate 1 in each condition. The first feature is the monotone anchor; the next two are smooth nonmonotone trajectories.\nLines show noiseless signals; points show all 300 observations. Axes share the same scales across all figures.") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "#F2F4F5", color = NA),
      legend.position = "bottom", plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      panel.spacing = grid::unit(0.8, "lines"))
  for (extension in c("png", "pdf")) {
    ggsave(file.path(out_dir, sprintf("design_M%d.%s", M, extension)), plot,
      width = 11, height = 2 * M + 1.6, dpi = 180, bg = "white")
  }
}
write.csv(do.call(rbind, selection), file.path(out_dir, "design_illustration_features.csv"),
          row.names = FALSE)
cat("Illustrated all nine conditions from saved replicate-one inputs.\n")
