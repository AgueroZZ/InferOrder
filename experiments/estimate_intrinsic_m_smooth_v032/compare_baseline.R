#!/usr/bin/env Rscript
# Compare matched datasets with a common occupancy-based structural evaluation.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
source(file.path(study_dir, "evaluate.R"))

baseline_dir <- "experiments/estimate_intrinsic_m_v032"
out_dir <- file.path(study_dir, "main_summary")
baseline <- read.csv(file.path(baseline_dir, "main_summary", "runs.csv"))
smooth <- read.csv(file.path(out_dir, "runs.csv"))
manifest <- active_manifest()
manifest <- manifest[manifest$phase == "main", , drop = FALSE]
stopifnot(nrow(manifest) == 90L, nrow(smooth) == 180L,
          nrow(baseline) == 180L, !any(smooth$status == "pending"),
          all(baseline$status == "success"))

result_key <- function(rows) paste(rows$id, rows$method, sep = "_")
stopifnot(!anyDuplicated(result_key(baseline)),
          !anyDuplicated(result_key(smooth)),
          setequal(result_key(baseline), result_key(smooth)))
smooth <- smooth[match(result_key(baseline), result_key(smooth)), , drop = FALSE]
for (column in c("id", "true_M", "snr", "replicate", "method")) {
  stopifnot(identical(baseline[[column]], smooth[[column]]))
}

pairing_checks <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  id <- manifest$id[i]
  original <- readRDS(file.path(baseline_dir, "data", paste0(id, ".rds")))
  modified <- readRDS(file.path(study_dir, "data", paste0(id, ".rds")))
  anchors <- vapply(original$ordering_labels, function(ordering)
    which(original$true_assign == ordering)[1], integer(1))
  noise_difference <- max(abs((original$X - original$signal) -
                               (modified$X - modified$signal)))
  stopifnot(identical(original$latent_positions, modified$latent_positions),
            identical(original$true_assign, modified$true_assign),
            identical(original$ordering_labels, modified$ordering_labels),
            identical(original$signal[, anchors, drop = FALSE],
                      modified$signal[, anchors, drop = FALSE]),
            noise_difference < 1e-12)
  pairing_checks[[i]] <- data.frame(id = id, true_M = manifest$true_M[i],
    snr = manifest$snr[i], replicate = manifest$replicate[i],
    same_latent_positions = TRUE, same_feature_groups = TRUE,
    same_monotone_anchors = TRUE, maximum_noise_difference = noise_difference,
    baseline_input_hash = original$input_hash,
    smooth_input_hash = modified$input_hash)

  for (method in c("adaptive", "forward")) {
    row_index <- which(baseline$id == id & baseline$method == method)
    saved <- readRDS(file.path(baseline_dir, "results",
                              paste0(id, "_", method, ".rds")))
    compact <- saved$selected
    compact$active <- which(colMeans(compact$W) > design$effective_weight_tol)
    stopifnot(baseline$estimated_M[row_index] == length(compact$active))
    evaluation <- evaluate_fit(compact, original)
    baseline$ARI[row_index] <- evaluation$ARI
    baseline$ordering_recovery[row_index] <- evaluation$ordering_recovery
    baseline$matched_ordering_recovery[row_index] <- evaluation$matched_ordering_recovery
    baseline$true_ordering_coverage[row_index] <- evaluation$true_ordering_coverage
  }
}
write.csv(do.call(rbind, pairing_checks),
          file.path(out_dir, "baseline_pairing_checks.csv"), row.names = FALSE)
write.csv(baseline, file.path(out_dir, "baseline_runs_aligned.csv"), row.names = FALSE)

paired <- baseline[, c("id", "true_M", "snr", "replicate", "method")]
paired$baseline_status <- baseline$status
paired$smooth_status <- smooth$status
paired$baseline_M <- baseline$estimated_M
paired$smooth_M <- smooth$estimated_M
paired$baseline_exact <- baseline$status == "success" &
  baseline$estimated_M == baseline$true_M
paired$smooth_exact <- smooth$status == "success" &
  !is.na(smooth$estimated_M) & smooth$estimated_M == smooth$true_M
paired$both_exact <- paired$baseline_exact & paired$smooth_exact
paired$lost_exact <- paired$baseline_exact & !paired$smooth_exact
paired$gained_exact <- !paired$baseline_exact & paired$smooth_exact
paired$neither_exact <- !paired$baseline_exact & !paired$smooth_exact
paired$completed_pair <- baseline$status == "success" & smooth$status == "success"
for (metric in c("ARI", "ordering_recovery")) {
  paired[[paste0("baseline_", metric)]] <- baseline[[metric]]
  paired[[paste0("smooth_", metric)]] <- smooth[[metric]]
  paired[[paste0("difference_", metric)]] <- ifelse(paired$completed_pair,
    smooth[[metric]] - baseline[[metric]], NA_real_)
}
write.csv(paired, file.path(out_dir, "baseline_comparison_pairs.csv"), row.names = FALSE)

available_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summarize_pairs <- function(rows) {
  data.frame(planned = nrow(rows), baseline_exact = sum(rows$baseline_exact),
    smooth_exact = sum(rows$smooth_exact), both_exact = sum(rows$both_exact),
    lost_exact = sum(rows$lost_exact), gained_exact = sum(rows$gained_exact),
    neither_exact = sum(rows$neither_exact),
    unresolved = sum(rows$smooth_status != "success"),
    completed_pairs = sum(rows$completed_pair),
    baseline_mean_ARI = available_mean(rows$baseline_ARI),
    smooth_mean_ARI = available_mean(rows$smooth_ARI),
    mean_difference_ARI = available_mean(rows$difference_ARI),
    baseline_mean_ordering_recovery = available_mean(rows$baseline_ordering_recovery),
    smooth_mean_ordering_recovery = available_mean(rows$smooth_ordering_recovery),
    mean_difference_ordering_recovery = available_mean(rows$difference_ordering_recovery))
}
condition_groups <- split(paired,
  interaction(paired$true_M, paired$snr, paired$method, drop = TRUE))
comparison <- do.call(rbind, lapply(condition_groups, function(rows)
  cbind(rows[1, c("true_M", "snr", "method")], summarize_pairs(rows))))
rownames(comparison) <- NULL
write.csv(comparison, file.path(out_dir, "baseline_comparison_summary.csv"), row.names = FALSE)
overall <- do.call(rbind, lapply(split(paired, paired$method), function(rows)
  cbind(method = rows$method[1], summarize_pairs(rows))))
rownames(overall) <- NULL
write.csv(overall, file.path(out_dir, "baseline_comparison_overall.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
method_labels <- c(adaptive = "Adaptive EB", forward = "Uniform + forward")
comparison$method_label <- factor(method_labels[comparison$method],
                                  levels = unname(method_labels))
comparison$snr_label <- factor(comparison$snr, levels = design$snr)
wilson <- function(x, n) {
  z <- qnorm(0.975)
  p <- x / n
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, center - half), upper = min(1, center + half))
}
recovery <- rbind(
  transform(comparison, family = "All monotone", exact = baseline_exact),
  transform(comparison, family = "One monotone anchor", exact = smooth_exact))
intervals <- t(vapply(seq_len(nrow(recovery)), function(i)
  wilson(recovery$exact[i], recovery$planned[i]), numeric(2)))
recovery$rate <- recovery$exact / recovery$planned
recovery$lower <- intervals[, "lower"]
recovery$upper <- intervals[, "upper"]
recovery$family <- factor(recovery$family,
                         levels = c("All monotone", "One monotone anchor"))
plot_recovery <- ggplot(recovery,
  aes(snr_label, rate, color = family, group = family)) +
  geom_line(position = position_dodge(width = 0.12)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12,
                position = position_dodge(width = 0.12)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.12)) +
  facet_grid(method_label ~ true_M,
             labeller = labeller(true_M = function(x) paste("True M =", x))) +
  scale_color_manual(values = c("All monotone" = "#0072B2",
                                "One monotone anchor" = "#D55E00")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Variance signal-to-noise ratio", y = "Exact effective-M recovery rate",
       color = "Trajectory design") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
structure <- rbind(
  transform(comparison, metric = "Feature partition: ARI",
            value = mean_difference_ARI),
  transform(comparison, metric = "Ordering: absolute Spearman",
            value = mean_difference_ordering_recovery))
plot_structure <- ggplot(structure,
  aes(snr_label, value, color = method_label, group = method_label)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.4) +
  geom_line() + geom_point(size = 2.5) +
  facet_grid(metric ~ true_M, scales = "free_y",
             labeller = labeller(true_M = function(x) paste("True M =", x))) +
  scale_color_manual(values = c("Adaptive EB" = "#287D8E",
                                "Uniform + forward" = "#D87443")) +
  labs(x = "Variance signal-to-noise ratio",
       y = "Mean paired change: smooth design minus all monotone", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
for (extension in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("baseline_comparison_recovery.", extension)),
         plot_recovery, width = 11.5, height = 7, dpi = 180, bg = "white")
  ggsave(file.path(out_dir, paste0("baseline_comparison_structure.", extension)),
         plot_structure, width = 11.5, height = 7, dpi = 180, bg = "white")
}
jsonlite::write_json(list(validated = TRUE, paired_datasets = nrow(manifest),
  paired_method_outcomes = nrow(paired), occupancy_threshold = design$effective_weight_tol,
  structural_evaluation = "common occupied-slot one-to-one matching",
  maximum_noise_difference = max(vapply(pairing_checks,
    function(x) x$maximum_noise_difference, numeric(1)))),
  file.path(out_dir, "baseline_comparison_validation.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(overall, row.names = FALSE)
