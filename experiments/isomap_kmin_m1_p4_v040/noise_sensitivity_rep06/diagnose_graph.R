# Post hoc graph checks after the fixed-case noise comparison.
# Truth supplies diagnostic labels and an oracle path, never a fitted initializer.
helpers <- new.env(parent = globalenv())
source("experiments/isomap_kmin_m1_p4_v040/common.R", local = helpers)
suppressPackageStartupMessages(library(ggplot2))
output <- file.path(helpers$study, "noise_sensitivity_rep06")
source_input <- readRDS(helpers$input_path(6L))
record <- readRDS(file.path(output, "provenance.rds"))
stopifnot(identical(record$source_input_sha256, helpers$hash_file(helpers$input_path(6L))),
  identical(record$parent_provenance$package_source_hashes, helpers$package_source_hashes))
truth_rank <- rank(source_input$truth, ties.method = "first")
ordered_samples <- order(source_input$truth)
ordered_signal <- source_input$signal[ordered_samples, , drop = FALSE]
arc_steps <- sqrt(rowSums((ordered_signal[-1L, , drop = FALSE] -
  ordered_signal[-nrow(ordered_signal), , drop = FALSE])^2))
arc <- numeric(nrow(ordered_signal))
arc[ordered_samples] <- c(0, cumsum(arc_steps))
summaries <- edges <- components <- list()
for (sd in c(.25, .10, .05, .01, 0)) {
  X <- source_input$signal + sd * source_input$standard_noise
  dimnames(X) <- dimnames(source_input$X)
  raw <- MPCurver::isomap_ordering(X, seed = record$fitting_seed)
  distances <- as.matrix(dist(X))
  diag(distances) <- Inf
  make_graph <- function(k) {
    neighbors <- t(vapply(seq_len(nrow(X)), function(index)
      order(distances[index, ], seq_len(nrow(X)))[seq_len(k)], integer(k)))
    pairs <- cbind(rep(seq_len(nrow(X)), each = k), as.vector(t(neighbors)))
    igraph::graph_from_edgelist(pairs, directed = FALSE)
  }
  graph <- igraph::simplify(make_graph(raw$k_used))
  pairs <- igraph::as_edgelist(graph, names = FALSE)
  weights <- distances[pairs]
  igraph::E(graph)$weight <- weights
  graph_distances <- igraph::distances(graph)
  stopifnot(igraph::is_connected(graph),
    max(abs(graph_distances - raw$geodesic_to_landmark)) < 1e-12)
  if (sd > 0) {
    saved <- readRDS(file.path(output, sprintf("sd_%03d", round(sd * 1000)),
      "results", "rep06_auto_kmin.rds"))
    stopifnot(max(abs(raw$t - saved$raw_positions)) < 1e-12)
  }
  for (k in seq_len(raw$k_used)) components[[length(components) + 1L]] <- data.frame(
    noise_sd = sd, k = k, components = igraph::components(make_graph(k))$no)
  edge_data <- data.frame(noise_sd = sd, sample_1 = pairs[, 1L], sample_2 = pairs[, 2L],
    rank_1 = truth_rank[pairs[, 1L]], rank_2 = truth_rank[pairs[, 2L]],
    true_position_1 = source_input$truth[pairs[, 1L]],
    true_position_2 = source_input$truth[pairs[, 2L]], euclidean_distance = weights)
  edge_data$rank_gap <- abs(edge_data$rank_1 - edge_data$rank_2)
  edge_data$true_position_gap <- abs(edge_data$true_position_1 - edge_data$true_position_2)
  edge_data$true_order_path_distance <- abs(arc[pairs[, 1L]] - arc[pairs[, 2L]])
  edges[[length(edges) + 1L]] <- edge_data
  summaries[[length(summaries) + 1L]] <- data.frame(noise_sd = sd, k_used = raw$k_used,
    raw_spearman = abs(cor(raw$t, source_input$truth, method = "spearman")),
    edges = nrow(edge_data), max_rank_gap = max(edge_data$rank_gap),
    max_true_position_gap = max(edge_data$true_position_gap))
  if (sd == 0) noiseless_positions <- raw$t
}
summaries <- do.call(rbind, summaries)
edges <- do.call(rbind, edges)
components <- do.call(rbind, components)
write.csv(summaries, file.path(output, "graph_summary.csv"), row.names = FALSE)
write.csv(edges, file.path(output, "graph_edges.csv"), row.names = FALSE)
write.csv(components, file.path(output, "graph_components.csv"), row.names = FALSE)

# Same noiseless sampled curve with edges only between consecutive true ranks.
oracle <- as.numeric(stats::cmdscale(dist(arc), k = 1L))
oracle <- (oracle - min(oracle)) / diff(range(oracle))
oracle_spearman <- abs(cor(oracle, source_input$truth, method = "spearman"))
stopifnot(oracle_spearman > 1 - 1e-12)
if (cor(noiseless_positions, source_input$truth, method = "spearman") < 0)
  noiseless_positions <- 1 - noiseless_positions
if (cor(oracle, source_input$truth, method = "spearman") < 0) oracle <- 1 - oracle
noiseless <- data.frame(sample = seq_along(oracle), truth = source_input$truth,
  feature_1 = source_input$signal[, 1L], feature_4 = source_input$signal[, 4L],
  auto_positions = noiseless_positions, oracle_positions = oracle)
write.csv(noiseless, file.path(output, "noiseless_diagnostic_positions.csv"), row.names = FALSE)
helpers$save_atomic(list(checked_at_utc = format(Sys.time(), tz = "UTC"),
  source_input_sha256 = record$source_input_sha256,
  script_sha256 = helpers$hash_file(file.path(output, "diagnose_graph.R")),
  parent_provenance = record$parent_provenance, graph_summary = summaries,
  independent_graph_geodesics_match = TRUE, all_saved_auto_positions_match = TRUE,
  oracle_noiseless_path_spearman = oracle_spearman,
  description = "Post hoc initializer-only noiseless check; oracle uses true ordering and is not an estimated fit."),
  file.path(output, "graph_diagnostic_provenance.rds"))

# Rebuild Figure 2 from saved plot data with readable color-bar ticks.
points <- read.csv(file.path(output, "plotted_feature_points.csv"))
points$stage <- factor(points$stage, levels = c("True positions", "Auto Isomap"))
points$level <- factor(points$level, levels = unique(points$level))
caption <- function(text, width) paste(strwrap(text, width), collapse = "\n")
theme <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9),
  legend.position = "bottom")
feature_plot <- ggplot(points, aes(feature_1, feature_4, color = color_position)) +
  geom_point(size = 1.5, alpha = .85) + facet_grid(stage ~ level) + coord_equal() +
  scale_color_viridis_c(name = "Position along ordering", limits = c(0, 1),
    breaks = seq(0, 1, .25), option = "D", end = .95,
    guide = guide_colorbar(barwidth = grid::unit(7, "cm"))) +
  labs(x = "Observed feature 1", y = "Observed feature 4",
    title = "Feature 1 versus feature 4 as noise decreases",
    subtitle = "Top: colors show truth | Bottom: colors show raw automatic Isomap positions computed from all four features",
    caption = caption(paste("Figure 2. Same trajectories and noise directions, with decreasing noise SD from left to right.",
      "Axes are observed features 1 and 4; Isomap uses all four observed features. Top/bottom panels contain identical samples",
      "at each noise level. Raw Isomap positions precede binning and MPCurve fitting; only global reflection is allowed for display.",
      "Axis limits and [0,1] color limits are shared across panels."), 136)) + theme
for (extension in c("png", "pdf")) ggsave(file.path(output, paste0("feature1_feature4_noise.", extension)),
  feature_plot, width = 15, height = 9, dpi = 180, bg = "white")
print(summaries, digits = 6)
print(head(edges[order(-edges$rank_gap), ], 10L), digits = 6)
cat("Noiseless true-order oracle path Spearman:", oracle_spearman, "\n")
