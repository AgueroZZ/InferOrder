#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R")

summary_dir <- file.path(extension_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

new_rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  path <- file.path(extension_dir, "results",
    paste0(manifest$id[i], "_auto_adaptive.rds"))
  if (file.exists(path)) {
    result <- readRDS(path)
    stopifnot(result$design_hash == design_hash,
      result$input_hash == load_fixed_dataset(manifest[i, , drop = FALSE])$input_hash,
      identical(result$provenance, execution_provenance()))
    new_rows[[i]] <- result$row
  } else {
    new_rows[[i]] <- data.frame(
      id = manifest$id[i], phase = "main", true_M = manifest$true_M[i],
      snr = manifest$snr[i], replicate = manifest$replicate[i],
      method = design$method, status = "pending", selected_initial_M = NA_integer_,
      selected_model_M = NA_integer_, estimated_M = NA_integer_,
      effective_M = NA_integer_, initial_partition_ARI = NA_real_, ARI = NA_real_,
      ordering_recovery = NA_real_, matched_ordering_recovery = NA_real_,
      true_ordering_coverage = NA_real_, minimum_initial_cluster_size = NA_integer_,
      elapsed_seconds = NA_real_, candidate_count = 1L, warning_count = NA_integer_,
      continuations = NA_integer_, stop_reason = "pending", error = ""
    )
  }
}
new_runs <- do.call(rbind, new_rows)
write.csv(new_runs, file.path(summary_dir, "auto_runs.csv"), row.names = FALSE)

old_runs <- read.csv(file.path(study_dir, "main_summary", "runs.csv"))
old_runs <- old_runs[old_runs$method %in% c("adaptive", "forward"), , drop = FALSE]
old_runs$selected_initial_M <- NA_integer_
old_runs$initial_partition_ARI <- NA_real_

common_columns <- c("id", "true_M", "snr", "replicate", "method", "status",
  "selected_initial_M", "effective_M", "initial_partition_ARI", "ARI",
  "ordering_recovery", "elapsed_seconds", "warning_count")
combined <- rbind(old_runs[, common_columns], new_runs[, common_columns])
combined$method <- factor(combined$method,
  levels = c("adaptive", "auto_adaptive", "forward"))
write.csv(combined, file.path(summary_dir, "combined_runs.csv"), row.names = FALSE)

wilson <- function(x, n) {
  z <- stats::qnorm(0.975)
  p <- x / n
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, center - half), upper = min(1, center + half))
}
mean_available <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_available <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

condition_groups <- split(combined,
  interaction(combined$true_M, combined$snr, combined$method, drop = TRUE))
condition_summary <- do.call(rbind, lapply(condition_groups, function(x) {
  success <- x$status == "success"
  exact <- sum(success & x$effective_M == x$true_M, na.rm = TRUE)
  interval <- wilson(exact, nrow(x))
  data.frame(
    true_M = x$true_M[1], snr = x$snr[1], method = as.character(x$method[1]),
    planned = nrow(x), successful = sum(success), exact = exact,
    under = sum(success & x$effective_M < x$true_M, na.rm = TRUE),
    over = sum(success & x$effective_M > x$true_M, na.rm = TRUE),
    exact_rate = exact / nrow(x), lower = interval[1], upper = interval[2],
    mean_ARI = mean_available(x$ARI),
    mean_ordering_recovery = mean_available(x$ordering_recovery),
    median_seconds = median_available(x$elapsed_seconds),
    total_seconds = sum(x$elapsed_seconds, na.rm = TRUE),
    warnings = sum(x$warning_count, na.rm = TRUE)
  )
}))
rownames(condition_summary) <- NULL
write.csv(condition_summary, file.path(summary_dir, "condition_summary.csv"),
  row.names = FALSE)

method_groups <- split(combined, combined$method, drop = TRUE)
overall <- do.call(rbind, lapply(method_groups, function(x) {
  success <- x$status == "success"
  data.frame(
    method = as.character(x$method[1]),
    planned = nrow(x), successful = sum(success),
    exact = sum(success & x$effective_M == x$true_M, na.rm = TRUE),
    under = sum(success & x$effective_M < x$true_M, na.rm = TRUE),
    over = sum(success & x$effective_M > x$true_M, na.rm = TRUE),
    mean_ARI = mean_available(x$ARI),
    mean_ordering_recovery = mean_available(x$ordering_recovery),
    median_seconds = median_available(x$elapsed_seconds),
    total_seconds = sum(x$elapsed_seconds, na.rm = TRUE),
    warnings = sum(x$warning_count, na.rm = TRUE)
  )
}))
rownames(overall) <- NULL
write.csv(overall, file.path(summary_dir, "overall_summary.csv"), row.names = FALSE)

auto_success <- new_runs[new_runs$status == "success", , drop = FALSE]
initialization_summary <- do.call(rbind, lapply(split(auto_success,
  interaction(auto_success$true_M, auto_success$snr, drop = TRUE)), function(x) {
  data.frame(
    true_M = x$true_M[1], snr = x$snr[1], datasets = nrow(x),
    exact_initial_M = sum(x$selected_initial_M == x$true_M),
    exact_effective_M = sum(x$effective_M == x$true_M),
    initial_under = sum(x$selected_initial_M < x$true_M),
    initial_over = sum(x$selected_initial_M > x$true_M),
    final_under = sum(x$effective_M < x$true_M),
    final_over = sum(x$effective_M > x$true_M),
    mean_initial_partition_ARI = mean(x$initial_partition_ARI),
    mean_final_ARI = mean(x$ARI)
  )
}))
rownames(initialization_summary) <- NULL
write.csv(initialization_summary,
  file.path(summary_dir, "initialization_summary.csv"), row.names = FALSE)

old_adaptive <- old_runs[old_runs$method == "adaptive", c("id", "true_M", "snr",
  "replicate", "effective_M", "ARI", "ordering_recovery")]
names(old_adaptive)[5:7] <- c("old_effective_M", "old_ARI", "old_ordering_recovery")
paired <- merge(old_adaptive, auto_success[, c("id", "selected_initial_M",
  "effective_M", "initial_partition_ARI", "ARI", "ordering_recovery")],
  by = "id", all.x = TRUE, sort = FALSE)
paired$old_exact <- paired$old_effective_M == paired$true_M
paired$initial_exact <- paired$selected_initial_M == paired$true_M
paired$auto_exact <- paired$effective_M == paired$true_M
paired$gained_exact <- !paired$old_exact & paired$auto_exact
paired$lost_exact <- paired$old_exact & !paired$auto_exact
paired$delta_ARI <- paired$ARI - paired$old_ARI
paired$delta_ordering_recovery <- paired$ordering_recovery - paired$old_ordering_recovery
write.csv(paired, file.path(summary_dir, "paired_vs_original_adaptive.csv"),
  row.names = FALSE)

comparison_summary <- data.frame(
  paired_datasets = nrow(paired),
  old_adaptive_exact = sum(paired$old_exact, na.rm = TRUE),
  similarity_cut_exact = sum(paired$initial_exact, na.rm = TRUE),
  auto_adaptive_exact = sum(paired$auto_exact, na.rm = TRUE),
  gained_exact = sum(paired$gained_exact, na.rm = TRUE),
  lost_exact = sum(paired$lost_exact, na.rm = TRUE),
  mean_delta_ARI = mean(paired$delta_ARI, na.rm = TRUE),
  mean_delta_ordering_recovery = mean(paired$delta_ordering_recovery, na.rm = TRUE)
)
write.csv(comparison_summary, file.path(summary_dir, "comparison_summary.csv"),
  row.names = FALSE)

if (any(new_runs$status == "pending")) {
  cat("Automatic-M runs remain pending; wrote partial tabular summaries.\n")
  print(overall, row.names = FALSE)
  quit(status = 0L)
}

suppressPackageStartupMessages(library(ggplot2))
method_labels <- c(
  adaptive = "Original adaptive EB (M = 8)",
  auto_adaptive = "Automatic-M adaptive EB",
  forward = "Uniform + forward"
)
method_colors <- c(
  "Original adaptive EB (M = 8)" = "#6A51A3",
  "Automatic-M adaptive EB" = "#009E73",
  "Uniform + forward" = "#D55E00"
)
theme_set(theme_minimal(base_size = 12))
save_plot <- function(plot, name, width, height) {
  for (extension in c("png", "pdf")) {
    ggsave(file.path(summary_dir, paste0(name, ".", extension)), plot,
      width = width, height = height, dpi = 180, bg = "white")
  }
}

condition_summary$method_label <- factor(method_labels[condition_summary$method],
  levels = unname(method_labels))
condition_summary$snr_label <- factor(condition_summary$snr, levels = design$snr)
accuracy <- ggplot(condition_summary,
  aes(snr_label, exact_rate, color = method_label, group = method_label)) +
  geom_line(position = position_dodge(width = 0.12)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12,
    position = position_dodge(width = 0.12)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.12)) +
  facet_wrap(~true_M, nrow = 1, labeller = label_bquote(M[true] == .(true_M))) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Variance signal-to-noise ratio",
    y = "Exact effective-M recovery rate", color = NULL,
    title = "Automatic initialization improves adaptive ordering-count recovery") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(accuracy, "exact_recovery_comparison", 12, 4.8)

distribution_groups <- split(combined,
  interaction(combined$true_M, combined$snr, combined$method, drop = TRUE))
distribution <- do.call(rbind, lapply(distribution_groups, function(x) {
  values <- ifelse(x$status == "success", as.character(x$effective_M), "Unresolved")
  levels <- c(as.character(seq_len(design$max_intrinsic_dim)), "Unresolved")
  counts <- table(factor(values, levels = levels))
  data.frame(true_M = x$true_M[1], snr = x$snr[1],
    method = as.character(x$method[1]), estimated = levels,
    count = as.integer(counts), probability = as.integer(counts) / nrow(x))
}))
distribution$estimated <- factor(distribution$estimated,
  levels = c(as.character(seq_len(design$max_intrinsic_dim)), "Unresolved"))
distribution$method_label <- factor(method_labels[distribution$method],
  levels = unname(method_labels))
distribution$snr_label <- factor(paste0("SNR = ", distribution$snr),
  levels = paste0("SNR = ", design$snr))
distribution$truth <- factor(distribution$true_M, levels = rev(design$true_M))
heatmap <- ggplot(distribution, aes(estimated, truth, fill = probability)) +
  geom_tile(color = "white") +
  geom_tile(data = distribution[
    as.character(distribution$estimated) == as.character(distribution$true_M), ],
    fill = NA, color = "#333333", linewidth = 0.65) +
  geom_text(aes(label = ifelse(count > 0, count, "")), size = 3.2) +
  facet_grid(snr_label ~ method_label) +
  scale_fill_gradient(low = "#FFFFFF", high = "#56B4E9", limits = c(0, 1)) +
  labs(x = "Estimated effective M", y = "True M", fill = "Proportion",
    title = "Distribution of estimated ordering counts",
    subtitle = "Cell labels give counts out of 10; outlined cells recover the true M") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35,
    hjust = 1))
save_plot(heatmap, "estimated_M_distribution_comparison", 15, 8)
write.csv(distribution, file.path(summary_dir,
  "estimated_M_distribution_comparison.csv"), row.names = FALSE)

metrics <- rbind(
  transform(condition_summary, metric = "Feature partition: mean ARI",
    value = mean_ARI),
  transform(condition_summary,
    metric = "Ordering recovery: mean |Spearman|",
    value = mean_ordering_recovery)
)
quality <- ggplot(metrics,
  aes(snr_label, value, color = method_label, group = method_label)) +
  geom_line() + geom_point(size = 2.3) + facet_grid(metric ~ true_M) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Variance signal-to-noise ratio",
    y = "Mean recovery score", color = NULL) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(quality, "structural_recovery_comparison", 13, 7)

auto_stage <- rbind(
  data.frame(true_M = auto_success$true_M, snr = auto_success$snr,
    replicate = auto_success$replicate, stage = "Similarity cut",
    exact = auto_success$selected_initial_M == auto_success$true_M),
  data.frame(true_M = auto_success$true_M, snr = auto_success$snr,
    replicate = auto_success$replicate, stage = "Adaptive EB final",
    exact = auto_success$effective_M == auto_success$true_M)
)
stage_summary <- aggregate(exact ~ true_M + snr + stage, auto_stage, mean)
stage_summary$snr_label <- factor(stage_summary$snr, levels = design$snr)
stage_summary$stage <- factor(stage_summary$stage,
  levels = c("Similarity cut", "Adaptive EB final"))
stage_plot <- ggplot(stage_summary,
  aes(snr_label, exact, color = stage, group = stage)) +
  geom_line() + geom_point(size = 2.5) +
  facet_wrap(~true_M, nrow = 1, labeller = label_bquote(M[true] == .(true_M))) +
  scale_color_manual(values = c("Similarity cut" = "#0072B2",
    "Adaptive EB final" = "#009E73")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Variance signal-to-noise ratio", y = "Exact-M recovery rate",
    color = NULL, title = "Data-only initialization and final adaptive estimate") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(stage_plot, "auto_initial_vs_final", 11, 4.7)

timing <- ggplot(condition_summary,
  aes(snr_label, median_seconds, color = method_label, group = method_label)) +
  geom_line() + geom_point(size = 2.3) + facet_wrap(~true_M, nrow = 1) +
  scale_color_manual(values = method_colors) +
  labs(x = "Variance signal-to-noise ratio", y = "Median fitting time (seconds)",
    color = NULL, title = "Runtime per dataset") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_plot(timing, "runtime_comparison", 12, 4.7)

print(overall, row.names = FALSE)
print(comparison_summary, row.names = FALSE)
