source("experiments/isomap_kmin_m1_p2_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
arguments <- commandArgs(trailingOnly = TRUE)
primary_metric <- if (length(arguments)) match.arg(arguments[1L], c("cosine", "centered_cosine")) else "cosine"
metric_title <- if (primary_metric == "cosine") "Cosine similarity, allowing global reversal" else
  "Absolute centered cosine similarity"
results <- unlist(lapply(seq_len(design$replications), function(replication)
  lapply(design$methods, function(method) readRDS(result_path(replication, method)))), recursive = FALSE)
rows <- do.call(rbind, lapply(results, function(result) {
  data.frame(replication = result$replication, method = result$method,
    k_used = result$k_used, status = result$status, converged = result$converged,
    iterations = result$iterations, graph_components = result$graph_components,
    graph_keep_count = result$graph_keep_count,
    initializer = if (is.null(result$initialization)) NA_character_ else result$initialization$method_used,
    fallback = if (is.null(result$initialization)) NA else result$initialization$fallback,
    cosine = result$metrics[["cosine"]], centered_cosine = result$metrics[["centered_cosine"]],
    spearman = result$metrics[["spearman"]], raw_cosine = result$raw_metrics[["cosine"]],
    raw_centered_cosine = result$raw_metrics[["centered_cosine"]],
    constant_cosine = result$constant_cosine,
    final_elbo = if (length(result$elbo_trace)) tail(result$elbo_trace, 1) else NA_real_,
    fit_seconds = result$elapsed_seconds, warning_count = length(result$warnings),
    warnings = paste(result$warnings, collapse = " | "), error = result$error,
    input_sha256 = result$input_sha256)
}))
rows$primary_score <- rows[[primary_metric]]
write.csv(rows, file.path(study, "metrics.csv"), row.names = FALSE)

summaries <- do.call(rbind, lapply(design$methods, function(method) {
  selected <- rows[rows$method == method, , drop = FALSE]
  finite <- is.finite(selected$primary_score)
  scores <- selected$primary_score[finite]
  data.frame(method = method, metric = primary_metric,
    total = nrow(selected), finite = sum(finite), converged = sum(selected$converged),
    failed = sum(selected$status == "failed"), fallbacks = sum(selected$fallback, na.rm = TRUE),
    median = if (length(scores)) median(scores) else NA_real_,
    mean = if (length(scores)) mean(scores) else NA_real_,
    q1 = if (length(scores)) unname(quantile(scores, .25)) else NA_real_,
    q3 = if (length(scores)) unname(quantile(scores, .75)) else NA_real_,
    min = if (length(scores)) min(scores) else NA_real_,
    max = if (length(scores)) max(scores) else NA_real_,
    median_centered_cosine = median(selected$centered_cosine, na.rm = TRUE),
    median_spearman = median(selected$spearman, na.rm = TRUE),
    min_k = min(selected$k_used, na.rm = TRUE), median_k = median(selected$k_used, na.rm = TRUE),
    max_k = max(selected$k_used, na.rm = TRUE),
    median_iterations = median(selected$iterations, na.rm = TRUE))
}))
write.csv(summaries, file.path(study, "summary.csv"), row.names = FALSE)
automatic <- rows[rows$method == "auto_kmin", ]
paired_rows <- do.call(rbind, lapply(c("fixed_k15", "fixed_k10"), function(method) {
  comparator <- rows[rows$method == method, ]
  matched <- match(automatic$replication, comparator$replication)
  data.frame(replication = automatic$replication, comparator = method,
    auto_score = automatic$primary_score, fixed_score = comparator$primary_score[matched],
    difference = automatic$primary_score - comparator$primary_score[matched])
}))
paired_summary <- do.call(rbind, lapply(split(paired_rows, paired_rows$comparator), function(values) {
  differences <- values$difference[is.finite(values$difference)]
  count <- length(differences)
  data.frame(comparator = values$comparator[1L], metric = primary_metric,
    complete_pairs = count, missing_pairs = nrow(values) - count,
    mean_difference = if (count) mean(differences) else NA_real_,
    median_difference = if (count) median(differences) else NA_real_,
    monte_carlo_se = if (count > 1L) sd(differences) / sqrt(count) else NA_real_,
    improvements = sum(differences > 1e-10), ties = sum(abs(differences) <= 1e-10),
    declines = sum(differences < -1e-10))
}))
write.csv(paired_rows, file.path(study, "paired_differences.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(study, "paired_summary.csv"), row.names = FALSE)

labels <- c(auto_kmin = "Auto kmin\n(new default)",
            fixed_k15 = "Fixed k = 15\n(previous default)",
            fixed_k10 = "Fixed k = 10\n(reference)")
colors <- c(auto_kmin = "#0072B2", fixed_k15 = "#D55E00", fixed_k10 = "#009E73")
plot_data <- rows
plot_data$method <- factor(plot_data$method, levels = design$methods)
plot_data$score <- plot_data$primary_score
write.csv(plot_data, file.path(study, "plotted_scores.csv"), row.names = FALSE)
finite_counts <- vapply(design$methods, function(method)
  sum(is.finite(plot_data$score[plot_data$method == method])), integer(1))
tick_labels <- setNames(paste0(labels[design$methods], "\n(n = ", finite_counts, ")"), design$methods)
definition <- if (primary_metric == "cosine")
  "Score = max{cos(t, q), cos(t, 1-q)} on the original [0,1] positions; no rank transform." else
  "Score = |cos(t-mean(t), q-mean(q))|, equal to absolute Pearson correlation; no rank transform."
caption <- paste("Figure 1. Twenty independent random cubic B-spline realizations with N=200, M=1, P=2",
  "and Gaussian noise SD=0.25 (average variance SNR=16).", definition,
  "Each point is a final MPCurve fit; gray lines connect identical datasets.",
  "Boxes show median and quartiles; whiskers reach the most extreme values within 1.5 IQR.",
  sprintf("Converged fits: %d/%d; failures: %d; missing scores: %d.",
    sum(rows$converged), nrow(rows), sum(rows$status == "failed"), sum(!is.finite(rows$primary_score))))
plot <- ggplot(plot_data, aes(method, score)) +
  geom_line(aes(group = replication), color = "#888888", linewidth = .35, alpha = .35, na.rm = TRUE) +
  geom_boxplot(aes(fill = method), width = .43, outlier.shape = NA, alpha = .4, na.rm = TRUE) +
  geom_point(aes(color = method), size = 2, alpha = .8, na.rm = TRUE) +
  scale_color_manual(values = colors, guide = "none") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_discrete(labels = tick_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(.01, .02))) +
  labs(x = "Isomap neighborhood setting", y = metric_title,
    title = "Position recovery with automatic and fixed Isomap neighborhoods",
    subtitle = "20 paired replicates | Independent nonmonotone curves | Fits run to convergence",
    caption = paste(strwrap(caption, 108), collapse = "\n")) +
  theme_minimal(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0, size = 9), panel.grid.minor = element_blank())
ggsave(file.path(study, "position_cosine_boxplot.png"), plot, width = 10, height = 7.5, dpi = 180, bg = "white")
ggsave(file.path(study, "position_cosine_boxplot.pdf"), plot, width = 10, height = 7.5)
save_atomic(list(primary_metric = primary_metric, definition = definition,
  summaries = summaries, paired_summary = paired_summary, plotted_count = sum(is.finite(rows$primary_score)),
  evaluated_at_utc = format(Sys.time(), tz = "UTC")), file.path(study, "evaluation.rds"))
print(summaries, digits = 5)
print(paired_summary, digits = 5)
