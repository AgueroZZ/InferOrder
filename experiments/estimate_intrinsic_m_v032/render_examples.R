#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
source(file.path(study_dir, "reporting_metrics.R"))
suppressPackageStartupMessages(library(ggplot2))
manifest <- read.csv(file.path(study_dir, "manifest.csv"))
representatives <- manifest[manifest$phase == "main" & manifest$snr == 4 & manifest$replicate == 1L, ]
out_dir <- file.path(study_dir, "main_summary")
for (i in seq_len(nrow(representatives))) {
  row <- representatives[i, ]
  dataset <- readRDS(file.path(study_dir, "data", paste0(row$id, ".rds")))
  results <- lapply(c("adaptive", "forward"), function(method)
    readRDS(file.path(study_dir, "results", paste0(row$id, "_", method, ".rds"))))
  names(results) <- c("adaptive", "forward")
  panels <- list()
  for (m in seq_len(row$true_M)) {
    columns <- head(which(dataset$true_assign == dataset$ordering_labels[m]), 3L)
    frame <- do.call(rbind, lapply(columns, function(j) data.frame(
      position = dataset$latent_positions[, m], signal = dataset$signal[, j], feature = factor(j))))
    panels[[length(panels) + 1L]] <- ggplot(frame, aes(position, signal, color = feature)) +
      geom_line(linewidth = 0.7) + scale_color_manual(values = c("#4477AA", "#EE6677", "#228833")) +
      labs(title = paste("True ordering", dataset$ordering_labels[m]),
        x = "True position", y = "Noiseless feature signal") + theme_minimal(base_size = 10) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
  }
  for (method in names(results)) for (m in seq_len(row$true_M)) {
    result <- results[[method]]
    label <- if (method == "adaptive") "Adaptive EB" else "Uniform + forward"
    if (result$row$status == "success") {
      alignment <- result$evaluation$alignment[m, ]
    } else alignment <- data.frame(matched = FALSE)
    if (isTRUE(alignment$matched)) {
      x <- dataset$latent_positions[, m]
      y <- result$selected$positions[, alignment$estimated_slot]
      if (is.finite(cor(x, y, method = "spearman")) && cor(x, y, method = "spearman") < 0) y <- 1 - y
      p <- ggplot(data.frame(x = x, y = y), aes(x, y)) +
        geom_point(alpha = 0.4, size = 0.8, color = if (method == "adaptive") "#287D8E" else "#D87443") +
        labs(title = sprintf("%s: |rho| = %.2f", label, alignment$abs_spearman))
    } else {
      p <- ggplot() + annotate("text", x = 0.5, y = 0.5,
        label = if (result$row$status == "success") "Unmatched ordering" else "Unresolved fit", size = 3) +
        labs(title = label)
    }
    panels[[length(panels) + 1L]] <- p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "True position", y = "Inferred position") + theme_minimal(base_size = 10) +
      theme(panel.grid.minor = element_blank())
  }
  estimates <- vapply(results, function(x) if (x$row$status == "success")
    as.character(reported_effective_M(x)) else "unresolved", character(1))
  figure <- patchwork::wrap_plots(panels, ncol = row$true_M) +
    patchwork::plot_annotation(title = sprintf("True M = %d: a preselected example", row$true_M),
      subtitle = sprintf("SNR = 4, replicate 1. Effective M: adaptive EB %s; uniform + forward %s.",
                         estimates[1], estimates[2]),
      caption = "The same dataset is used by both methods. Estimated orderings are matched one-to-one to truth and may be reversed for evaluation.")
  for (extension in c("png", "pdf")) ggsave(file.path(out_dir, sprintf("example_M%d.%s", row$true_M, extension)),
    figure, width = max(11, 2.8 * row$true_M), height = 9, dpi = 170, bg = "white")
}
