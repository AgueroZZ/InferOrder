# Plot paired Isomap and MPCurve orderings on a common rank scale.
# Each input row represents one sample; initial and final use the same orientation.
plot_ordering_comparison <- function(samples, panel_columns, title, subtitle, initial_label = 'Initial Isomap') {
  required <- c('sample', 'truth', 'isomap', 'final', panel_columns)
  stopifnot(all(required %in% names(samples)))
  groups <- interaction(samples[panel_columns], drop = TRUE)
  frames <- lapply(split(samples, groups), function(frame) {
    stopifnot(!anyDuplicated(frame$sample))
    valid <- is.finite(frame$isomap) & is.finite(frame$final)
    if (!all(valid)) stop('Comparison requires paired finite coordinates.')
    # Orient Isomap to truth, then orient the final fit to Isomap; a global
    # reversal is equivalent and should not be counted as iterative change.
    initial <- frame$isomap
    if (cor(frame$truth, initial, method = 'spearman') < 0) initial <- -initial
    final <- frame$final
    if (cor(initial, final, method = 'spearman') < 0) final <- -final
    normalized_rank <- function(x) (rank(x, ties.method = 'average') - 1) / (length(x) - 1)
    frame$initial_rank <- normalized_rank(initial)
    frame$final_rank <- normalized_rank(final)
    frame
  })
  aligned <- do.call(rbind, frames)
  initial <- transform(aligned, stage = initial_label, position = initial_rank)
  final <- transform(aligned, stage = 'Final MPCurve', position = final_rank)
  long <- rbind(initial, final)
  long$stage <- factor(long$stage, levels = c(initial_label, 'Final MPCurve'))
  annotations <- do.call(rbind, lapply(split(aligned, interaction(aligned[panel_columns], drop = TRUE)), function(x) {
    label <- sprintf('Truth |rho|: %.3f -> %.3f\nInitial-final rho: %.3f',
      abs(cor(x$truth, x$initial_rank, method = 'spearman')),
      abs(cor(x$truth, x$final_rank, method = 'spearman')),
      cor(x$initial_rank, x$final_rank, method = 'spearman'))
    cbind(x[1, panel_columns, drop = FALSE], label = label)
  }))
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = aligned,
      ggplot2::aes(x = truth, xend = truth, y = initial_rank, yend = final_rank),
      color = '#777777', alpha = 0.18, linewidth = 0.25) +
    ggplot2::geom_point(data = long,
      ggplot2::aes(truth, position, color = stage, shape = stage), alpha = 0.65, size = 1) +
    ggplot2::geom_label(data = annotations,
      ggplot2::aes(x = 0.02, y = 1.16, label = label),
      hjust = 0, vjust = 1, size = 3, linewidth = 0.15) +
    ggplot2::scale_color_manual(values = setNames(c('#0072B2', '#D55E00'), c(initial_label, 'Final MPCurve'))) +
    ggplot2::scale_shape_manual(values = setNames(c(1, 16), c(initial_label, 'Final MPCurve'))) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.19)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = 'bottom', panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = 'True latent position', y = 'Normalized sample rank',
      color = NULL, shape = NULL, title = title, subtitle = subtitle,
      caption = paste('Gray segments connect the same sample before and after iteration.',
        'Ranks remove coordinate-spacing differences; average ranks are used for ties.',
        'A global reversal is aligned away.'))
  list(plot = plot, aligned = aligned, annotations = annotations)
}
