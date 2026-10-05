# Exact Isomap neighbor graphs before and at the connectivity threshold.
# Layout uses observed features 1 and 4; graph distances use all four features.
helpers <- new.env(parent = globalenv())
source("experiments/isomap_kmin_m1_p4_v040/common.R", local = helpers)
suppressPackageStartupMessages(library(ggplot2))
output <- file.path(helpers$study, "noise_sensitivity_rep06")
input_path <- helpers$input_path(6L)
input <- readRDS(input_path)
record <- readRDS(file.path(output, "provenance.rds"))
stopifnot(identical(record$source_input_sha256, helpers$hash_file(input_path)),
  identical(record$parent_provenance$package_source_hashes, helpers$package_source_hashes))
summary_reference <- read.csv(file.path(output, "graph_summary.csv"))
saved_files <- c(input_path,
  list.files(output, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE),
  file.path(output, c("graph_summary.csv", "graph_edges.csv", "graph_components.csv")))
saved_files <- saved_files[!grepl("/full_fits/|/knn_graph_", saved_files)]
saved_hashes <- setNames(vapply(saved_files, helpers$hash_file, character(1)), saved_files)
package_graph <- get(".isomap_graph", envir = asNamespace("MPCurver"))
canonical_pairs <- function(graph) {
  pairs <- igraph::as_edgelist(igraph::simplify(graph), names = FALSE)
  pairs <- cbind(pmin(pairs[, 1L], pairs[, 2L]), pmax(pairs[, 1L], pairs[, 2L]))
  pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
}
pair_keys <- function(pairs) paste(pairs[, 1L], pairs[, 2L], sep = "-")
nodes <- edges <- summaries <- list()
for (sd in c(.25, .01, 0)) {
  X <- input$signal + sd * input$standard_noise
  dimnames(X) <- dimnames(input$X)
  kmin <- summary_reference$k_used[summary_reference$noise_sd == sd]
  stopifnot(length(kmin) == 1L)
  distances <- as.matrix(dist(X))
  diag(distances) <- Inf
  neighbor_ranks <- t(vapply(seq_len(nrow(X)), function(index) {
    ranks <- integer(nrow(X))
    ranks[order(distances[index, ], seq_len(nrow(X)))] <- seq_len(nrow(X))
    ranks
  }, integer(nrow(X))))
  make_graph <- function(k) {
    neighbors <- t(vapply(seq_len(nrow(X)), function(index)
      order(distances[index, ], seq_len(nrow(X)))[seq_len(k)], integer(k)))
    pairs <- cbind(rep(seq_len(nrow(X)), each = k), as.vector(t(neighbors)))
    graph <- igraph::simplify(igraph::graph_from_edgelist(pairs, directed = FALSE))
    stopifnot(identical(canonical_pairs(graph), canonical_pairs(package_graph(X, k)$graph)))
    graph
  }
  before <- make_graph(kmin - 1L)
  before_components <- igraph::components(before)
  component_order <- order(tapply(input$truth, before_components$membership, min))
  labels <- setNames(paste0("C", seq_along(component_order)), as.character(component_order))
  prior_component <- unname(labels[as.character(before_components$membership)])
  prior_sizes <- table(prior_component)
  prior_keys <- pair_keys(canonical_pairs(before))
  for (stage in c("Before: kmin - 1", "Selected: kmin")) {
    k <- if (stage == "Before: kmin - 1") kmin - 1L else kmin
    graph <- make_graph(k)
    components <- igraph::components(graph)
    pairs <- canonical_pairs(graph)
    keys <- pair_keys(pairs)
    joining <- prior_component[pairs[, 1L]] != prior_component[pairs[, 2L]]
    new_edge <- !(keys %in% prior_keys)
    rank_forward <- neighbor_ranks[pairs]
    rank_backward <- neighbor_ranks[cbind(pairs[, 2L], pairs[, 1L])]
    stopifnot(all(!joining | new_edge), all(pmin(rank_forward, rank_backward) <= k),
      all(pmin(rank_forward[new_edge], rank_backward[new_edge]) == k),
      components$no == if (stage == "Selected: kmin") 1L else before_components$no)
    column <- sprintf("Noise SD = %.2f", sd)
    nodes[[length(nodes) + 1L]] <- data.frame(noise_sd = sd, k = k, stage = stage,
      column = column, sample = seq_len(nrow(X)), feature_1 = X[, 1L], feature_4 = X[, 4L],
      true_position = input$truth, prior_component = prior_component)
    edges[[length(edges) + 1L]] <- data.frame(noise_sd = sd, k = k, stage = stage,
      column = column, sample_1 = pairs[, 1L], sample_2 = pairs[, 2L],
      x = X[pairs[, 1L], 1L], y = X[pairs[, 1L], 4L],
      xend = X[pairs[, 2L], 1L], yend = X[pairs[, 2L], 4L],
      true_position_1 = input$truth[pairs[, 1L]], true_position_2 = input$truth[pairs[, 2L]],
      component_1_before = prior_component[pairs[, 1L]],
      component_2_before = prior_component[pairs[, 2L]],
      added_at_kmin = new_edge, joins_previous_components = joining,
      neighbor_rank_1_to_2 = rank_forward, neighbor_rank_2_to_1 = rank_backward,
      euclidean_distance_4d = distances[pairs])
    summaries[[length(summaries) + 1L]] <- data.frame(noise_sd = sd, k = k,
      stage = stage, column = column, samples = nrow(X), components = components$no,
      unique_edges = nrow(pairs), added_edges = sum(new_edge),
      joining_edges = sum(joining),
      prior_component_sizes = paste(paste(names(prior_sizes), as.integer(prior_sizes), sep = "="), collapse = "; "))
  }
}
nodes <- do.call(rbind, nodes)
edges <- do.call(rbind, edges)
summaries <- do.call(rbind, summaries)
stages <- c("Before: kmin - 1", "Selected: kmin")
columns <- unique(nodes$column)
for (object in c("nodes", "edges", "summaries")) {
  data <- get(object)
  data$stage <- factor(data$stage, levels = stages)
  data$column <- factor(data$column, levels = columns)
  assign(object, data)
}
summaries$label <- sprintf("k=%d | %d component%s\n%d unique edges", summaries$k,
  summaries$components, ifelse(summaries$components == 1L, "", "s"), summaries$unique_edges)
bridges <- edges[edges$joins_previous_components, ]
bridge_endpoints <- unique(rbind(
  data.frame(stage = bridges$stage, column = bridges$column, x = bridges$x, y = bridges$y),
  data.frame(stage = bridges$stage, column = bridges$column, x = bridges$xend, y = bridges$yend)))
write.csv(nodes, file.path(output, "knn_graph_nodes.csv"), row.names = FALSE)
write.csv(edges, file.path(output, "knn_graph_edges.csv"), row.names = FALSE)
write.csv(summaries, file.path(output, "knn_graph_summary.csv"), row.names = FALSE)
write.csv(bridges, file.path(output, "knn_graph_connecting_edges.csv"), row.names = FALSE)

caption <- paste(strwrap(paste("Figure 3. Exact union-kNN graphs: an undirected edge is retained when either sample selects the other.",
  "Top: kmin-1; bottom: automatic kmin. All 200 samples and all unique edges are displayed.",
  "Nodes retain their top-panel component colors in the bottom panel; component labels restart within each column.",
  "Red edges join different top-panel components; red rings mark their endpoints. Edges use four-dimensional distances;",
  "the plotted axes are observed features 1 and 4. Projected line crossings do not create graph connections."), 132), collapse = "\n")
plot <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), color = "#6B7280", linewidth = .25, alpha = .28) +
  geom_point(data = nodes, aes(feature_1, feature_4, color = prior_component), size = 1.4, alpha = .95) +
  geom_segment(data = bridges, aes(x, y, xend = xend, yend = yend), color = "#CC0000", linewidth = 1, alpha = .9) +
  geom_point(data = bridge_endpoints, aes(x, y), shape = 21, fill = NA, color = "#CC0000", size = 2.8, stroke = .8) +
  geom_label(data = summaries, aes(x = -Inf, y = Inf, label = label), hjust = -.02, vjust = 1.02,
    size = 3.6, label.size = 0, fill = "white") +
  facet_grid(stage ~ column) + coord_equal() +
  scale_color_manual(values = c(C1 = "#0072B2", C2 = "#E69F00", C3 = "#009E73"),
    name = "Components at kmin - 1") +
  labs(x = "Observed feature 1", y = "Observed feature 4",
    title = "Isomap kNN graphs immediately before and at connectivity",
    subtitle = "Replicate 6, P=4 | Red edges connect components that are separate in the top row",
    caption = caption) + theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9))
for (extension in c("png", "pdf")) ggsave(file.path(output, paste0("knn_graph_connectivity.", extension)),
  plot, width = 13, height = 10, dpi = 180, bg = "white")
stopifnot(identical(saved_hashes, setNames(vapply(saved_files, helpers$hash_file, character(1)), saved_files)),
  nrow(nodes) == 1200L)
helpers$save_atomic(list(verified_at_utc = format(Sys.time(), tz = "UTC"),
  script_sha256 = helpers$hash_file(file.path(output, "plot_knn_graph.R")),
  source_file_sha256 = saved_hashes, package_source_hashes = helpers$package_source_hashes,
  graph_definition = "Undirected union-kNN; Euclidean distances on all four observed features; row-index tie breaking.",
  package_edge_sets_match = TRUE, source_artifacts_preserved = TRUE,
  displayed_samples = nrow(nodes), displayed_edges = nrow(edges),
  displayed_connecting_edges = nrow(bridges), summary = summaries),
  file.path(output, "knn_graph_provenance.rds"))
print(summaries, row.names = FALSE)
print(bridges[, c("noise_sd", "k", "sample_1", "sample_2", "true_position_1", "true_position_2",
  "neighbor_rank_1_to_2", "neighbor_rank_2_to_1", "euclidean_distance_4d")],
  row.names = FALSE, digits = 6)
