source("experiments/isomap_kmin_m1_p2_v040/common.R")
metadata <- readRDS(file.path(study, "provenance.rds"))
current <- provenance()
stopifnot(identical(metadata$design_hash, current$design_hash),
  identical(metadata$package_source_hashes, current$package_source_hashes),
  identical(metadata$archive_sha256, current$archive_sha256),
  identical(metadata$script_hashes, current$script_hashes))
evaluation <- readRDS(file.path(study, "evaluation.rds"))
metrics <- read.csv(file.path(study, "metrics.csv"))
plotted <- read.csv(file.path(study, "plotted_scores.csv"))
stopifnot(nrow(metrics) == design$replications * length(design$methods),
  !anyDuplicated(metrics[c("replication", "method")]),
  evaluation$primary_metric %in% c("spearman", "cosine", "centered_cosine"),
  identical(metrics$primary_score, metrics[[evaluation$primary_metric]]),
  identical(metrics$primary_score, plotted$score))
graph_checks <- list()
for (replication in seq_len(design$replications)) {
  input <- readRDS(input_path(replication))
  stopifnot(identical(input, simulate_replicate(replication)),
    identical(input$input_sha256, hash_object(input$X)),
    all(apply(input$dense_signal, 2L, is_nonmonotone)),
    abs(mean(input$dense_signal^2) - 1) < 1e-12,
    identical(dim(input$X), c(design$N, design$P)))
  for (method in design$methods) {
    result <- readRDS(result_path(replication, method))
    stopifnot(identical(result$input_sha256, input$input_sha256),
      identical(result$design_hash, design_hash),
      identical(result$package_source_hashes, package_source_hashes),
      identical(result$metrics, position_metrics(input$truth, result$positions)))
    row <- metrics[metrics$replication == replication & metrics$method == method, ]
    stopifnot(nrow(row) == 1L,
      isTRUE(all.equal(unlist(row[c("cosine", "centered_cosine", "spearman")], use.names = FALSE),
                       unname(result$metrics), tolerance = 1e-12)))
    if (result$status != "failed") stopifnot(
      abs(row$spearman - abs(cor(input$truth, result$positions, method = "spearman"))) < 1e-12)
    if (result$status != "failed") {
      stopifnot(length(result$positions) == design$N,
        all(is.finite(result$positions)),
        all(result$positions >= -1e-12 & result$positions <= 1 + 1e-12),
        result$num_bins == design$num_bins,
        length(result$elbo_trace) == result$iterations + 1L,
        min(diff(result$elbo_trace)) >= -1e-7)
      if (result$converged) stopifnot(
        abs(tail(diff(result$elbo_trace), 1L)) / length(input$X) < design$tolerance)
    }
    if (method == "auto_kmin" && !is.na(result$k_used)) {
      # Independent reference: full pairwise distances, no RANN/helper graph.
      distances <- as.matrix(dist(input$X))
      diag(distances) <- Inf
      connected <- function(k) {
        indices <- t(vapply(seq_len(design$N), function(sample)
          order(distances[sample, ], seq_len(design$N))[seq_len(k)], integer(k)))
        graph <- igraph::graph_from_edgelist(cbind(rep(seq_len(design$N), each = k),
          as.vector(t(indices))), directed = FALSE)
        igraph::is_connected(graph)
      }
      stopifnot(connected(result$k_used),
        result$k_used == 1L || !connected(result$k_used - 1L),
        result$graph_components == 1L, result$graph_keep_count == design$N)
      graph_checks[[length(graph_checks) + 1L]] <- c(replication = replication, k = result$k_used)
    }
    if (method != "auto_kmin" && !is.na(result$k_used)) stopifnot(
      result$k_used == switch(method, fixed_k15 = 15L, fixed_k10 = 10L))
  }
}

# Cosine is evaluated on coordinates; centered cosine equals Pearson up to sign.
truth <- c(0.05, 0.2, 0.6, 0.9)
stopifnot(abs(position_metrics(truth, truth)[["cosine"]] - 1) < 1e-14,
  abs(position_metrics(truth, 1 - truth)[["cosine"]] - 1) < 1e-14,
  abs(position_metrics(truth, truth^2)[["centered_cosine"]] - abs(cor(truth, truth^2))) < 1e-14,
  abs(position_metrics(truth, truth^2)[["spearman"]] - 1) < 1e-14,
  abs(position_metrics(truth, 1 - truth^2)[["spearman"]] - 1) < 1e-14)

namespace <- asNamespace("MPCurver")
source_environment <- new.env(parent = namespace)
for (path in list.files(file.path(package_source, "R"), pattern = "\\.R$", full.names = TRUE))
  sys.source(path, envir = source_environment)
functions <- ls(source_environment, all.names = TRUE)
functions <- functions[vapply(functions, function(name)
  is.function(get(name, source_environment, inherits = FALSE)), logical(1))]
stopifnot(all(vapply(functions, function(name) {
  installed <- get(name, namespace, inherits = FALSE)
  source <- get(name, source_environment, inherits = FALSE)
  identical(body(installed), body(source)) && identical(formals(installed), formals(source))
}, logical(1))))
verification <- list(verified_at_utc = format(Sys.time(), tz = "UTC"),
  replicates = design$replications, fits = nrow(metrics),
  converged = sum(metrics$converged), failures = sum(metrics$status == "failed"),
  fallback_count = sum(metrics$fallback, na.rm = TRUE),
  warning_count = sum(metrics$warning_count), primary_metric = evaluation$primary_metric,
  plotted_count = evaluation$plotted_count, graph_checks = graph_checks,
  installed_functions = length(functions), provenance = current)
save_atomic(verification, file.path(study, "verification.rds"))
writeLines(capture.output(str(verification[c("verified_at_utc", "replicates", "fits", "converged",
  "failures", "fallback_count", "warning_count", "primary_metric", "plotted_count",
  "installed_functions", "graph_checks")])), file.path(study, "verification.txt"))
cat(sprintf("Verified %d inputs, %d fits, %d kmin bounds, and %d installed functions.\n",
  verification$replicates, verification$fits, length(graph_checks), length(functions)))
