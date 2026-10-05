# Render paired initial/final ranks from cached full-data sensitivity fits.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(study_dir, 'exploratory_m5_ordering_b')
source(file.path(out, 'plot_ordering_comparison.R'))
paths <- file.path(out, 'sensitivity_fits', sprintf('setting_%02d.rds', 1:15))
stopifnot(all(file.exists(paths)))
results <- lapply(paths, readRDS)
stopifnot(sum(vapply(results, function(x) x$status == 'disconnected_initialization', logical(1))) == 1L)
complete <- Filter(function(x) x$status != 'disconnected_initialization', results)
stopifnot(all(vapply(complete, function(x) x$status == 'converged', logical(1))))
samples <- do.call(rbind, lapply(complete, function(x) data.frame(
 sample = seq_along(x$truth), truth = x$truth, isomap = x$isomap, final = x$final,
 noise_scale = x$setting$noise_scale, k = x$setting$k)))
comparison <- plot_ordering_comparison(samples, c('noise_scale', 'k'),
 title = 'Ordering B: Isomap initialization and final MPCurve ordering',
 subtitle = 'Rows: multiplier of the original B noise; columns: Isomap neighbors. Other feature groups and fitting controls are fixed.')
p <- comparison$plot + theme(plot.caption = element_text(hjust = 0, size = 9), plot.caption.position = 'plot') + facet_grid(noise_scale ~ k, labeller = label_both) +
 geom_label(data = data.frame(noise_scale = 0, k = 5),
  aes(x = 0.5, y = 0.53), label = 'Disconnected Isomap graph\n244/300 samples in largest component\nFull-data fit not run',
  size = 3.1, linewidth = 0.2) +
 labs(caption = paste(
 'Blue open circles: raw Isomap ranks; orange dots: final posterior-mean ranks.',
 'Gray segments connect the same sample before and after iteration.',
 'Ranks remove coordinate-spacing differences; average ranks are used for ties. Global reversal is aligned away.',
 'Each fit uses M=5 and the original feature groups for initialization; subsequent feature assignments and priors are learned.', sep = '\n'))
for (extension in c('png', 'pdf')) {
 ggsave(file.path(out, paste0('geometry_initial_final_comparison.', extension)),
  p, width = 17, height = 10, dpi = 180, bg = 'white')
}
write.csv(comparison$aligned, file.path(out, 'geometry_initial_final_positions.csv'), row.names = FALSE)
summary <- do.call(rbind, lapply(results, function(x) {
 valid <- identical(x$status, 'converged')
 data.frame(noise_scale = x$setting$noise_scale, k = x$setting$k,
  status = x$status, initial_rho = if (valid) abs(cor(x$truth,x$isomap,method='spearman')) else NA_real_,
  final_rho = if (valid) abs(cor(x$truth,x$final,method='spearman')) else NA_real_,
  initial_final_rho = if (valid) abs(cor(x$isomap,x$final,method='spearman')) else NA_real_,
  iterations = if (valid) x$iterations else NA_integer_,
  objective = if (valid) tail(x$objective,1) else NA_real_,
  ARI = if (valid) x$ARI else NA_real_, warning_count = length(x$warnings))
}))
write.csv(summary, file.path(out, 'geometry_initial_final_summary.csv'), row.names = FALSE)
# Check the observed-data controls against the earlier paired restart.
previous <- read.csv(file.path(out, 'controlled_restarts.csv'))
for (i in seq_len(nrow(previous))) {
 current <- summary[summary$noise_scale == 1 & summary$k == previous$B_neighbors[i], ]
 stopifnot(nrow(current) == 1L, abs(current$final_rho-previous$B_rho[i]) < 1e-8,
  abs(current$objective-previous$objective[i]) < 1e-5)
}
print(summary, row.names = FALSE)
