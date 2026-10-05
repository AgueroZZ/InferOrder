# Evaluate Samko et al.'s connectivity and degree rules on the unchanged B group.
# Each undirected union-graph edge is counted once, including reciprocal pairs.
source("experiments/isomap_elbo_screen_v040/common.R")
sample_count <- nrow(observations)
neighbor_values <- seq_len(sample_count - 1L)
nearest <- RANN::nn2(observations, observations, k = sample_count,
  treetype = "kd", eps = 0)
stopifnot(all(nearest$nn.idx[, 1L] == seq_len(sample_count)),
  all(nearest$nn.dists[, 1L] == 0))
neighbor_indices <- nearest$nn.idx[, -1L, drop = FALSE]

simple_graph <- function(indices) {
  from <- rep(seq_len(sample_count), each = ncol(indices))
  to <- as.vector(t(indices))
  stopifnot(length(from) == sample_count * ncol(indices), all(from != to))
  low <- pmin(from, to)
  high <- pmax(from, to)
  unique_pair <- !duplicated(low + sample_count * (high - 1L))
  igraph::graph_from_data_frame(
    data.frame(from = low[unique_pair], to = high[unique_pair]),
    directed = FALSE, vertices = seq_len(sample_count))
}

rows <- lapply(neighbor_values, function(k) {
  graph <- simple_graph(neighbor_indices[, seq_len(k), drop = FALSE])
  edge_count <- igraph::ecount(graph)
  degree <- igraph::degree(graph)
  components <- igraph::components(graph)$no
  stopifnot(igraph::is_simple(graph), sum(degree) == 2 * edge_count,
    min(degree) >= k)
  data.frame(k = k, components = components, connected = components == 1L,
    unique_edges = edge_count, average_degree = 2 * edge_count / sample_count,
    threshold = k + 2L, excess_degree = 2 * edge_count / sample_count - k,
    nonreciprocal_entries = 2 * edge_count - sample_count * k,
    # Compare integer numerators rather than rounding average degrees.
    degree_rule_pass = 2 * edge_count <= sample_count * (k + 2L))
})
diagnostics <- do.call(rbind, rows)
lower <- min(diagnostics$k[diagnostics$connected])
literal_upper <- max(diagnostics$k[diagnostics$degree_rule_pass])
violations <- diagnostics$k[diagnostics$k >= lower & !diagnostics$degree_rule_pass]
first_violation <- if (length(violations)) min(violations) else NA_integer_
local_upper <- if (is.na(first_violation)) literal_upper else first_violation - 1L
local_empty <- local_upper < lower

# Record every passing run; the degree condition need not be monotone in k.
run_groups <- rle(diagnostics$degree_rule_pass)
run_ends <- cumsum(run_groups$lengths)
run_starts <- c(1L, head(run_ends, -1L) + 1L)
passing_runs <- data.frame(lower = run_starts[run_groups$values],
  upper = run_ends[run_groups$values])

# Verify all-neighbor truncation against separate nearest-neighbor queries and
# the package-style multigraph followed by simplification at relevant k values.
check_values <- sort(unique(c(lower, local_upper, first_violation, 5L, 10L, 15L,
  sample_count - 1L)))
check_values <- check_values[!is.na(check_values) & check_values >= 1L]
checks <- lapply(check_values, function(k) {
  direct <- RANN::nn2(observations, observations, k = k + 1L,
    treetype = "kd", eps = 0)$nn.idx[, -1L, drop = FALSE]
  expected <- simple_graph(neighbor_indices[, seq_len(k), drop = FALSE])
  actual <- simple_graph(direct)
  pairs <- function(graph) sort(apply(igraph::as_edgelist(graph), 1L, function(edge)
    paste(sort(as.integer(edge)), collapse = ":")))
  from <- rep(seq_len(sample_count), each = k)
  multigraph <- igraph::graph_from_data_frame(
    data.frame(from = from, to = as.vector(t(direct))),
    directed = FALSE, vertices = seq_len(sample_count))
  simplified <- igraph::simplify(multigraph, remove.multiple = TRUE, remove.loops = TRUE)
  stopifnot(identical(pairs(expected), pairs(actual)),
    identical(pairs(actual), pairs(simplified)),
    igraph::ecount(multigraph) == sample_count * k)
  data.frame(k = k, independent_query_matches = TRUE,
    simplified_package_graph_matches = TRUE)
})
checks <- do.call(rbind, checks)
reference <- readRDS(file.path(experiment_dir, "results.rds"))
stopifnot(identical(reference$provenance$group_input_sha256,
  digest::digest(observations, algo = "sha256")),
  all(diff(diagnostics$components) <= 0),
  tail(diagnostics$unique_edges, 1L) == choose(sample_count, 2L),
  literal_upper == sample_count - 1L)

summary <- data.frame(interpretation = c("literal_global_maximum", "first_violation_from_connectivity"),
  k_min = lower, k_max = c(literal_upper, local_upper),
  empty_interval = c(FALSE, local_empty), first_violation = c(NA_integer_, first_violation))
write.csv(diagnostics, file.path(experiment_dir, "graph_bounds.csv"), row.names = FALSE)
write.csv(summary, file.path(experiment_dir, "graph_bounds_summary.csv"), row.names = FALSE)
write.csv(passing_runs, file.path(experiment_dir, "graph_rule_passing_runs.csv"), row.names = FALSE)
saveRDS(list(diagnostics = diagnostics, summary = summary,
  passing_runs = passing_runs, checks = checks,
  prior_selected_k_in_local_interval = !local_empty &&
    reference$selected_k >= lower && reference$selected_k <= local_upper,
  provenance = c(provenance(), list(graph_type = "simple undirected union kNN",
    nearest_neighbor_query = list(k = sample_count, treetype = "kd", eps = 0),
    graph_bounds_script_sha256 = digest::digest(
      file = file.path(experiment_dir, "graph_bounds.R"), algo = "sha256")))),
  file.path(experiment_dir, "graph_bounds.rds"))
print(summary, row.names = FALSE)
cat("All degree-rule passing runs:\n")
print(passing_runs, row.names = FALSE)
cat("First neighborhoods, including the earlier selected k=10:\n")
print(diagnostics[diagnostics$k <= max(15L, first_violation, na.rm = TRUE), ], row.names = FALSE)
cat("Independent graph checks:\n")
print(checks, row.names = FALSE)
