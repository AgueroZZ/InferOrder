source("experiments/isomap_elbo_screen_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
result <- readRDS(file.path(experiment_dir, "results.rds"))
timings <- read.csv(file.path(experiment_dir, "timings.csv"))
candidates <- do.call(rbind, lapply(result$candidates, function(candidate) {
  if (candidate$status != "converged") return(data.frame(k = candidate$k,
    status = candidate$status, initial_elbo = NA_real_, one_step_elbo = NA_real_,
    final_elbo = NA_real_, raw_rho = NA_real_, one_step_rho = NA_real_,
    final_rho = NA_real_, iterations = NA_integer_))
  data.frame(k = candidate$k, status = candidate$status,
    initial_elbo = candidate$early$elbo_trace[1L],
    one_step_elbo = candidate$early$elbo_trace[2L],
    final_elbo = tail(candidate$final$elbo_trace, 1L),
    raw_rho = recovery(candidate$raw_position),
    one_step_rho = recovery(candidate$early$position),
    final_rho = recovery(candidate$final$position), iterations = candidate$final$iter)
}))
write.csv(candidates, file.path(experiment_dir, "candidates.csv"), row.names = FALSE)
eligible <- candidates[candidates$status == "converged", ]
full_choice <- eligible$k[order(-eligible$final_elbo, eligible$k)[1L]]
sparse <- eligible[eligible$k %in% design$sparse_neighbors, ]
sparse_choice <- sparse$k[order(-sparse$one_step_elbo, sparse$k)[1L]]

budget_selection <- do.call(rbind, lapply(c(0L, 1L, 2L, 3L, 5L, 10L), function(sweeps) {
  scores <- vapply(result$candidates, function(candidate) {
    if (candidate$status != "converged" || length(candidate$final$elbo_trace) < sweeps + 1L)
      return(-Inf)
    candidate$final$elbo_trace[sweeps + 1L]
  }, numeric(1))
  best <- order(-scores, design$dense_neighbors)[1L]
  chosen <- candidates[best, ]
  data.frame(sweeps = sweeps, selected_k = chosen$k, selection_elbo = scores[best],
    final_elbo = chosen$final_elbo, final_rho = chosen$final_rho,
    final_elbo_loss = max(eligible$final_elbo) - chosen$final_elbo)
}))
write.csv(budget_selection, file.path(experiment_dir, "selection_by_budget.csv"), row.names = FALSE)

timing_summary <- do.call(rbind, lapply(split(timings, timings$method), function(rows) {
  default <- timings[timings$method == "fixed_k15", ]
  paired <- match(rows$repeat_index, default$repeat_index)
  data.frame(method = rows$method[1L], selected_k = rows$selected_k[1L],
    median_seconds = median(rows$total_seconds), min_seconds = min(rows$total_seconds),
    max_seconds = max(rows$total_seconds),
    median_screening_seconds = median(rows$screening_seconds),
    median_continuation_seconds = median(rows$continuation_seconds),
    median_paired_ratio_to_default = median(rows$total_seconds / default$total_seconds[paired]),
    median_paired_extra_seconds = median(rows$total_seconds - default$total_seconds[paired]),
    final_rho = rows$final_rho[1L], final_elbo = rows$final_elbo[1L])
}))
write.csv(timing_summary, file.path(experiment_dir, "timing_summary.csv"), row.names = FALSE)

score_plot_data <- rbind(
  data.frame(k = eligible$k, value = eligible$one_step_elbo, panel = "ELBO after one CAVI sweep"),
  data.frame(k = eligible$k, value = eligible$final_elbo, panel = "ELBO at convergence"),
  data.frame(k = eligible$k, value = eligible$final_rho, panel = "Ordering recovery at convergence"))
score_plot_data$panel <- factor(score_plot_data$panel,
  levels = c("ELBO after one CAVI sweep", "ELBO at convergence", "Ordering recovery at convergence"))
caption1 <- paste("Figure 1. Each point is a neighborhood on the same 300-by-12 observation matrix.",
  "ELBO is compared within each fitting stage; higher values are preferred.",
  "Recovery is absolute Spearman correlation with true positions.",
  "Orange dashed line: k=15 comparator; blue solid line: one-sweep selection.")
p1 <- ggplot(score_plot_data, aes(k, value)) +
  geom_line(color = "#555555", linewidth = .4) + geom_point(size = 1.8) +
  geom_vline(xintercept = 15, color = "#D55E00", linetype = "dashed", linewidth = .7) +
  geom_vline(xintercept = result$selected_k, color = "#0072B2", linewidth = .7) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(5, 30, 5), minor_breaks = 5:30) +
  labs(x = "Isomap number of nearest neighbors (k)", y = NULL,
    title = "Can one sweep select an Isomap start that recovers the ordering?",
    subtitle = sprintf("Selected k=%d | Highest converged ELBO: k=%d | MPCurver 0.4.0; 50 position bins",
      result$selected_k, full_choice), caption = paste(strwrap(caption1, 115), collapse = "\n")) +
  theme_minimal(base_size = 11) + theme(plot.caption = element_text(hjust = 0, size = 9))
ggsave(file.path(experiment_dir, "selection.png"), p1, width = 10, height = 9, dpi = 160, bg = "white")
ggsave(file.path(experiment_dir, "selection.pdf"), p1, width = 10, height = 9)

normalized_rank <- function(values) (rank(values, ties.method = "average") - 1) / (length(values) - 1)
aligned_rank <- function(values) {
  ranks <- normalized_rank(values)
  if (cor(truth, values, method = "spearman") < 0) 1 - ranks else ranks
}
points <- list()
for (k in unique(c(15L, result$selected_k))) {
  candidate <- result$candidates[[match(k, design$dense_neighbors)]]
  stages <- list("Raw Isomap" = candidate$raw_position,
    "After one CAVI sweep" = candidate$early$position,
    "At convergence" = candidate$final$position)
  for (stage in names(stages)) {
    points[[length(points) + 1L]] <- data.frame(truth_rank = normalized_rank(truth),
      inferred_rank = aligned_rank(stages[[stage]]), stage = stage,
      method = if (k == 15L) "Default neighborhood: k=15" else sprintf("One-sweep selected: k=%d", k),
      rho_label = sprintf("|Spearman rho| = %.4f", recovery(stages[[stage]])))
  }
}
points <- do.call(rbind, points)
points$stage <- factor(points$stage, levels = c("Raw Isomap", "After one CAVI sweep", "At convergence"))
annotations <- unique(points[c("stage", "method", "rho_label")])
caption2 <- paste("Figure 2. True and inferred normalized sample ranks for the k=15 comparator and",
  "the start selected solely by one-sweep ELBO. Columns show raw Isomap, one CAVI sweep, and convergence.",
  "Global reversal is aligned using truth for display only; the gray diagonal is perfect rank recovery.",
  "All 300 observations are retained; no truth coordinates enter fitting or ELBO selection.")
p2 <- ggplot(points, aes(truth_rank, inferred_rank)) +
  geom_abline(slope = 1, intercept = 0, color = "#999999", linewidth = .4) +
  geom_point(color = "#0072B2", alpha = .6, size = .85) +
  geom_text(data = annotations, aes(x = .04, y = .96, label = rho_label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  facet_grid(method ~ stage) + coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "True latent-position rank (normalized to [0,1])",
    y = "Inferred position rank (orientation aligned)",
    title = "One-sweep screening followed by continuation on the original failure case",
    caption = paste(strwrap(caption2, 135), collapse = "\n")) +
  theme_minimal(base_size = 11) + theme(plot.caption = element_text(hjust = 0, size = 9))
ggsave(file.path(experiment_dir, "positions.png"), p2, width = 12, height = 8, dpi = 160, bg = "white")
ggsave(file.path(experiment_dir, "positions.pdf"), p2, width = 12, height = 8)
write.csv(points, file.path(experiment_dir, "plotted_positions.csv"), row.names = FALSE)
print(candidates, row.names = FALSE)
print(timing_summary, row.names = FALSE)
print(budget_selection, row.names = FALSE)
cat(sprintf("Dense selected k=%d; sparse selected k=%d; converged-ELBO selected k=%d\n",
  result$selected_k, sparse_choice, full_choice))
