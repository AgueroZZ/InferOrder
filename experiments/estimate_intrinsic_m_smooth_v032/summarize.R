#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "reporting_metrics.R"))
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "main"
stopifnot(phase %in% c("pilot", "main"))
manifest <- active_manifest()
manifest <- manifest[manifest$phase == phase, , drop = FALSE]
if (phase == "main") stopifnot(nrow(manifest) == 90L)
rows <- list()
input_hashes <- list()
threshold_rows <- list()
index <- 0L
for (i in seq_len(nrow(manifest))) for (method in c("adaptive", "forward")) {
  index <- index + 1L
  path <- file.path(study_dir, "results", paste0(manifest$id[i], "_", method, ".rds"))
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(saved$design_hash == design_hash)
    id <- manifest$id[i]
    if (!is.null(input_hashes[[id]])) stopifnot(identical(input_hashes[[id]], saved$input_hash))
    input_hashes[[id]] <- saved$input_hash
    if (method == "adaptive" && saved$row$status == "success") {
      threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(id = id,
        true_M = manifest$true_M[i], snr = manifest$snr[i],
        threshold = c(1e-12, 1e-9, 1e-6, 1e-3),
        effective_M = saved$selected$threshold_counts)
    }
    rows[[index]] <- reported_result_row(saved)
  } else {
    rows[[index]] <- data.frame(id = manifest$id[i], phase = phase,
      true_M = manifest$true_M[i], snr = manifest$snr[i], replicate = manifest$replicate[i],
      method = method, status = "pending", estimated_M = NA_integer_, ARI = NA_real_,
      ordering_recovery = NA_real_, matched_ordering_recovery = NA_real_,
      true_ordering_coverage = NA_real_, map_groups = NA_integer_,
      near_duplicate_orderings = NA, constant_orderings = NA_integer_,
      elapsed_seconds = NA_real_, candidate_count = NA_integer_, warning_count = NA_integer_,
      stop_reason = "pending", error = "",
      original_reported_M = NA_integer_, selected_model_M = NA_integer_, effective_M = NA_integer_)
  }
}
runs <- do.call(rbind, rows)
stopifnot(all(runs$status %in% c("pending", "success", "error", "nonconverged")))
out_dir <- file.path(study_dir, paste0(phase, "_summary"))
dir.create(out_dir, showWarnings = FALSE)
write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)
if (length(threshold_rows)) write.csv(do.call(rbind, threshold_rows),
  file.path(out_dir, "adaptive_threshold_sensitivity.csv"), row.names = FALSE)
wilson <- function(x, n) {
  z <- qnorm(0.975)
  p <- x / n
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, center - half), upper = min(1, center + half))
}
mean_available <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_available <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
groups <- split(runs, interaction(runs$true_M, runs$snr, runs$method, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) {
  success <- x$status == "success"
  exact <- sum(success & x$estimated_M == x$true_M, na.rm = TRUE)
  ci <- wilson(exact, nrow(x))
  data.frame(true_M = x$true_M[1], snr = x$snr[1], method = x$method[1],
    planned = nrow(x), completed = sum(x$status != "pending"), successful = sum(success),
    failed = sum(x$status %in% c("error", "nonconverged")),
    exact = exact, under = sum(success & x$estimated_M < x$true_M, na.rm = TRUE),
    over = sum(success & x$estimated_M > x$true_M, na.rm = TRUE),
    exact_rate = exact / nrow(x), lower = ci[1], upper = ci[2],
    mean_ARI = mean_available(x$ARI), mean_ordering_recovery = mean_available(x$ordering_recovery),
    median_seconds = median_available(x$elapsed_seconds),
    total_seconds = sum(x$elapsed_seconds, na.rm = TRUE))
}))
rownames(summary) <- NULL
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
complete <- !any(runs$status == "pending")
jsonlite::write_json(list(phase = phase, design_hash = design_hash,
  reported_dimension = "posterior_occupancy", effective_weight_tol = design$effective_weight_tol,
  complete = complete, planned_method_runs = nrow(runs),
  completed_method_runs = sum(runs$status != "pending"), successful_method_runs = sum(runs$status == "success"),
  fitting_seconds = sum(runs$elapsed_seconds, na.rm = TRUE)),
  file.path(out_dir, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary[, c("true_M", "snr", "method", "completed", "planned", "exact", "failed", "total_seconds")], row.names = FALSE)
if (!complete) {
  cat("Phase incomplete; final figures are deferred.\n")
  quit(status = 0L)
}
suppressPackageStartupMessages(library(ggplot2))
labels <- c(adaptive = "Adaptive EB", forward = "Uniform + forward")
summary$method_label <- factor(labels[summary$method], levels = unname(labels))
summary$snr_label <- factor(summary$snr, levels = design$snr)
palette <- c("Adaptive EB" = "#287D8E", "Uniform + forward" = "#D87443")
theme_set(theme_minimal(base_size = 12))
save_plot <- function(plot, name, width, height) {
  for (extension in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", extension)),
    plot, width = width, height = height, dpi = 170, bg = "white")
}
accuracy <- ggplot(summary, aes(snr_label, exact_rate, color = method_label, group = method_label)) +
  geom_line(position = position_dodge(width = 0.12)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, position = position_dodge(width = 0.12)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.12)) +
  facet_wrap(~true_M, nrow = 1, labeller = label_bquote(M[true] == .(true_M))) +
  scale_color_manual(values = palette) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Variance signal-to-noise ratio", y = "Exact effective-M recovery rate", color = NULL,
    title = "Recovery of the intrinsic number of orderings",
    subtitle = "One monotone anchor per ordering with smooth nonmonotone trajectories") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(accuracy, "exact_recovery", 11, 4.6)
distribution <- do.call(rbind, lapply(groups, function(x) {
  values <- ifelse(x$status == "success", as.character(x$estimated_M), "Unresolved")
  levels <- c(as.character(seq_len(design$max_M)), "Unresolved")
  counts <- table(factor(values, levels = levels))
  data.frame(true_M = x$true_M[1], snr = x$snr[1], method = x$method[1],
    estimated = levels, count = as.integer(counts), probability = as.integer(counts) / nrow(x))
}))
distribution$estimated <- factor(distribution$estimated,
  levels = c(as.character(seq_len(design$max_M)), "Unresolved"))
distribution$method_label <- factor(labels[distribution$method], levels = unname(labels))
distribution$snr_label <- factor(paste0("SNR = ", distribution$snr), levels = paste0("SNR = ", design$snr))
distribution$truth <- factor(distribution$true_M, levels = rev(design$true_M))
heatmap <- ggplot(distribution, aes(estimated, truth, fill = probability)) +
  geom_tile(color = "white") +
  geom_tile(data = distribution[as.character(distribution$estimated) == as.character(distribution$true_M), ],
    fill = NA, color = "#333333", linewidth = 0.65) +
  geom_text(aes(label = ifelse(count > 0, count, "")), size = 3.4) +
  facet_grid(snr_label ~ method_label) +
  scale_fill_gradient(low = "#FFFFFF", high = "#62A6C3", limits = c(0, 1)) +
  labs(x = "Estimated effective M", y = "True M", fill = "Proportion",
    title = "Distribution of estimated ordering counts", subtitle = "Cell labels give counts; outlined cells recover the true M") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(heatmap, "estimated_M_distribution", 12, 8)
write.csv(distribution, file.path(out_dir, "estimated_M_distribution.csv"), row.names = FALSE)
metrics <- rbind(transform(summary, metric = "Feature partition: mean ARI", value = mean_ARI),
  transform(summary, metric = "Ordering recovery: mean absolute Spearman", value = mean_ordering_recovery))
quality <- ggplot(metrics, aes(snr_label, value, color = method_label, group = method_label)) +
  geom_line() + geom_point(size = 2.3) + facet_grid(metric ~ true_M) +
  scale_color_manual(values = palette) +
  labs(x = "Variance signal-to-noise ratio", y = "Recovery among completed fits", color = NULL) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(quality, "structural_recovery", 12, 7)
timing <- ggplot(summary, aes(snr_label, median_seconds, color = method_label, group = method_label)) +
  geom_line() + geom_point(size = 2.5) + facet_wrap(~true_M, nrow = 1) +
  scale_color_manual(values = palette) + labs(x = "Variance signal-to-noise ratio",
    y = "Median fitting time (seconds)", color = NULL,
    title = "Runtime includes all candidates and continuations") + theme(legend.position = "bottom")
save_plot(timing, "runtime", 11, 4.5)
cat("Completed summaries and figures for", phase, "\n")
