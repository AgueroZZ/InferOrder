# Audit both directed neighbor lists for the noiseless cross-component edge.
helpers <- new.env(parent = globalenv())
source("experiments/isomap_kmin_m1_p4_v040/common.R", local = helpers)
suppressPackageStartupMessages(library(ggplot2))
output <- file.path(helpers$study, "noise_sensitivity_rep06")
input_path <- helpers$input_path(6L)
input <- readRDS(input_path)
record <- readRDS(file.path(output, "provenance.rds"))
stopifnot(identical(record$source_input_sha256, helpers$hash_file(input_path)),
  identical(record$parent_provenance$package_source_hashes, helpers$package_source_hashes))
X <- input$signal
nodes <- read.csv(file.path(output, "knn_graph_nodes.csv"))
nodes <- nodes[nodes$noise_sd == 0 & nodes$stage == "Selected: kmin", ]
stopifnot(identical(nodes$sample, seq_len(nrow(X))),
  max(abs(nodes$feature_1 - X[, 1L])) < 1e-13,
  max(abs(nodes$feature_4 - X[, 4L])) < 1e-13)
connecting <- read.csv(file.path(output, "knn_graph_connecting_edges.csv"))
edge <- connecting[connecting$noise_sd == 0 &
  abs(connecting$true_position_1 - connecting$true_position_2) > .5, ]
stopifnot(nrow(edge) == 1L)
queries <- c(edge$sample_1, edge$sample_2)
panels <- c("Green endpoint as query", "Blue endpoint as query")
package_neighbors <- get(".isomap_neighbors", envir = asNamespace("MPCurver"))(X, 6L)
audit <- backgrounds <- query_points <- list()
for (index in seq_along(queries)) {
  query <- queries[index]
  differences <- sweep(X, 2L, X[query, ], "-")
  distance_4d <- sqrt(rowSums(differences^2))
  distance_2d <- sqrt(rowSums(differences[, c(1L, 4L), drop = FALSE]^2))
  candidates <- setdiff(seq_len(nrow(X)), query)
  ordered <- candidates[order(distance_4d[candidates], candidates)]
  ordered_2d <- candidates[order(distance_2d[candidates], candidates)]
  stopifnot(identical(ordered[seq_len(6L)], as.integer(package_neighbors$idx[query, ])))
  rank_2d <- match(ordered, ordered_2d)
  audit[[index]] <- data.frame(panel = panels[index], query_sample = query,
    query_true_position = input$truth[query], neighbor_sample = ordered,
    neighbor_true_position = input$truth[ordered], component = nodes$prior_component[ordered],
    rank_4d = seq_along(ordered), rank_2d = rank_2d,
    distance_4d = distance_4d[ordered], distance_2d = distance_2d[ordered],
    delta_feature_1 = differences[ordered, 1L], delta_feature_2 = differences[ordered, 2L],
    delta_feature_3 = differences[ordered, 3L], delta_feature_4 = differences[ordered, 4L],
    x = X[query, 1L], y = X[query, 4L], xend = X[ordered, 1L], yend = X[ordered, 4L])
  backgrounds[[index]] <- transform(nodes, panel = panels[index])
  query_points[[index]] <- data.frame(panel = panels[index], query_sample = query,
    x = X[query, 1L], y = X[query, 4L], component = nodes$prior_component[query],
    label = sprintf("Query\nt=%.3f", input$truth[query]))
}
audit <- do.call(rbind, audit)
backgrounds <- do.call(rbind, backgrounds)
query_points <- do.call(rbind, query_points)
chosen <- audit[audit$rank_4d <= 6L, ]
crossing <- chosen[chosen$component != nodes$prior_component[chosen$query_sample], ]
stopifnot(nrow(audit) == 398L, nrow(chosen) == 12L, nrow(crossing) == 1L,
  crossing$query_sample == queries[1L], crossing$neighbor_sample == queries[2L],
  crossing$rank_4d == 6L,
  audit$rank_4d[audit$query_sample == queries[2L] & audit$neighbor_sample == queries[1L]] == 64L)
write.csv(audit, file.path(output, "shortcut_neighbor_audit.csv"), row.names = FALSE)
write.csv(chosen, file.path(output, "shortcut_six_neighbors.csv"), row.names = FALSE)
write.csv(backgrounds, file.path(output, "shortcut_neighbor_plot_points.csv"), row.names = FALSE)
for (name in c("audit", "backgrounds", "query_points", "chosen", "crossing")) {
  data <- get(name); data$panel <- factor(data$panel, levels = panels); assign(name, data)
}
query_points$label_x <- query_points$x + .75
query_points$label_y <- query_points$y + .65
crossing$label_x <- crossing$xend - .75
crossing$label_y <- crossing$yend + .45
crossing$label <- sprintf("Neighbor #6\nt=%.3f\n4D distance=%.3f", crossing$neighbor_true_position, crossing$distance_4d)
caption <- paste(strwrap(paste("Figure 4. Six nearest neighbors of each endpoint of the noiseless shortcut, computed from all four features.",
  "Black diamonds mark queries; arrows point from each query to its selected neighbors. The red arrow is selected by the green endpoint only.",
  "Background colors are the k=5 components from Figure 3. The green endpoint has five green neighbors before this blue sixth neighbor;",
  "the blue endpoint selects six blue neighbors and ranks the green endpoint 64th. Axes are features 1 and 4."), 128), collapse = "\n")
palette <- c(C1 = "#0072B2", C2 = "#E69F00", C3 = "#009E73")
plot <- ggplot() +
  geom_point(data = backgrounds, aes(feature_1, feature_4, color = prior_component), size = 1.3, alpha = .38) +
  geom_segment(data = chosen, aes(x, y, xend = xend, yend = yend), color = "#444444", linewidth = .6,
    arrow = grid::arrow(length = grid::unit(.11, "inches"))) +
  geom_point(data = chosen, aes(xend, yend, fill = component), shape = 21, color = "#222222", size = 3, stroke = .65) +
  geom_segment(data = crossing, aes(x, y, xend = xend, yend = yend), color = "#CC0000", linewidth = 1,
    arrow = grid::arrow(length = grid::unit(.14, "inches"))) +
  geom_point(data = query_points, aes(x, y, fill = component), shape = 23, color = "black", size = 4, stroke = 1.1) +
  geom_segment(data = query_points, aes(label_x, label_y, xend = x, yend = y), color = "#333333", linetype = "dotted") +
  geom_label(data = query_points, aes(label_x, label_y, label = label), size = 3.6, fill = "white") +
  geom_segment(data = crossing, aes(label_x, label_y, xend = xend, yend = yend), color = "#CC0000", linetype = "dotted") +
  geom_label(data = crossing, aes(label_x, label_y, label = label), color = "#AA0000", size = 3.5, fill = "white") +
  facet_wrap(~ panel, nrow = 1L) + coord_equal(xlim = c(-1.8, 3.3), ylim = c(-2.2, 2.9)) +
  scale_color_manual(values = palette, name = "Components at k=5") +
  scale_fill_manual(values = palette, guide = "none") +
  labs(x = "Observed feature 1 (noise SD=0)", y = "Observed feature 4 (noise SD=0)",
    title = "Which endpoint selects the shortcut?",
    subtitle = "Each panel shows only the six directed neighbor selections from its marked query",
    caption = caption) + theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9))
for (extension in c("png", "pdf")) ggsave(file.path(output, paste0("shortcut_neighbor_direction.", extension)),
  plot, width = 12, height = 7, dpi = 180, bg = "white")
helpers$save_atomic(list(verified_at_utc = format(Sys.time(), tz = "UTC"),
  source_input_sha256 = record$source_input_sha256,
  script_sha256 = helpers$hash_file(file.path(output, "inspect_shortcut_neighbors.R")),
  source_plot_nodes_sha256 = helpers$hash_file(file.path(output, "knn_graph_nodes.csv")),
  source_edge_sha256 = helpers$hash_file(file.path(output, "knn_graph_connecting_edges.csv")),
  package_source_hashes = helpers$package_source_hashes, query_samples = queries,
  independent_neighbor_lists_match_package = TRUE, ranked_candidates = nrow(audit),
  displayed_neighbors = nrow(chosen), source_artifacts_preserved =
    identical(record$source_input_sha256, helpers$hash_file(input_path))),
  file.path(output, "shortcut_neighbor_provenance.rds"))
print(chosen[, c("panel", "neighbor_sample", "neighbor_true_position", "component", "rank_4d", "distance_4d", "rank_2d", "distance_2d")],
  digits = 6, row.names = FALSE)
