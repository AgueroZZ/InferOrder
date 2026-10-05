# Fixed trajectories and fixed errors isolate the effect of reducing noise SD.
helpers <- new.env(parent = globalenv())
source("experiments/isomap_kmin_m1_p4_v040/common.R", local = helpers)
suppressPackageStartupMessages(library(ggplot2))
parent_study <- helpers$study
output <- file.path(parent_study, "noise_sensitivity_rep06")
source_input_path <- helpers$input_path(6L)
source_input <- readRDS(source_input_path)
metadata <- readRDS(file.path(parent_study, "provenance.rds"))
current <- helpers$provenance()
stopifnot(identical(metadata$package_source_hashes, current$package_source_hashes),
  identical(metadata$script_hashes, current$script_hashes),
  identical(source_input$input_sha256, helpers$hash_object(source_input$X)))
levels <- c(.25, .10, .05, .01)
methods <- helpers$design$methods
level_root <- function(sd) file.path(output, sprintf("sd_%03d", round(sd * 1000)))
source_fits <- setNames(vapply(methods, function(method) helpers$result_path(6L, method), character(1)), methods)
record <- list(replication = 6L, noise_sd = levels,
  source_input_sha256 = helpers$hash_file(source_input_path),
  source_fits_sha256 = setNames(vapply(source_fits, helpers$hash_file, character(1)), methods),
  script_sha256 = helpers$hash_file(file.path(output, "run.R")),
  design_file_sha256 = helpers$hash_file(file.path(output, "DESIGN.md")),
  parent_provenance = metadata, fitting_seed = 202620026L)
path <- file.path(output, "provenance.rds")
if (file.exists(path)) stopifnot(identical(readRDS(path), record)) else helpers$save_atomic(record, path)

# Reuse the frozen parent fitting function, with isolated design/output bindings.
fit_environment <- function(sd) {
  environment <- new.env(parent = helpers)
  environment$study <- level_root(sd)
  environment$design <- helpers$design
  environment$design$replications <- 1L
  environment$design$noise_sd <- sd
  environment$design$source_replication <- 6L
  environment$design$fixed_noise_realization <- TRUE
  environment$design_hash <- helpers$hash_object(environment$design)
  environment$result_path <- helpers$result_path
  environment(environment$result_path) <- environment
  environment$runner <- helpers$fit_replicate
  environment(environment$runner) <- environment
  environment
}
for (sd in levels) {
  environment <- fit_environment(sd)
  for (directory in c("inputs", "results", "full_fits"))
    dir.create(file.path(environment$study, directory), recursive = TRUE, showWarnings = FALSE)
  input <- source_input
  input$X <- source_input$signal + sd * source_input$standard_noise
  dimnames(input$X) <- dimnames(source_input$X)
  input$input_sha256 <- helpers$hash_object(input$X)
  input$noise_sd <- sd
  input$design_hash <- if (sd == .25) source_input$design_hash else environment$design_hash
  if (sd == .25) stopifnot(identical(input$X, source_input$X))
  input_path <- file.path(environment$study, "inputs", "rep06.rds")
  if (file.exists(input_path)) stopifnot(identical(readRDS(input_path), input)) else
    helpers$save_atomic(input, input_path)
  for (method in methods) {
    result_path <- environment$result_path(6L, method)
    if (!file.exists(result_path)) {
      if (sd == .25) stopifnot(file.copy(source_fits[[method]], result_path)) else {
        result <- environment$runner(input, method)
        helpers$save_atomic(result, result_path)
      }
    }
    result <- readRDS(result_path)
    stopifnot(identical(result$input_sha256, input$input_sha256),
      identical(result$design_hash, input$design_hash),
      identical(result$package_source_hashes, metadata$package_source_hashes))
    cat(sprintf("SD=%.2f %-10s k=%d %s raw=%.6f final=%.6f iter=%d\n", sd, method,
      result$k_used, result$status, result$raw_metrics[["spearman"]], result$metrics[["spearman"]], result$iterations))
    flush.console()
  }
}

rows <- list()
points <- list()
graph_checks <- list()
for (sd in levels) {
  environment <- fit_environment(sd)
  input <- readRDS(file.path(environment$study, "inputs", "rep06.rds"))
  stopifnot(identical(input$truth, source_input$truth),
    identical(input$signal, source_input$signal), identical(input$dense_signal, source_input$dense_signal),
    identical(input$standard_noise, source_input$standard_noise),
    max(abs(input$X - (source_input$signal + sd * source_input$standard_noise))) == 0)
  for (method in methods) {
    result <- readRDS(environment$result_path(6L, method))
    stopifnot(identical(result$metrics, helpers$position_metrics(input$truth, result$positions)),
      identical(result$raw_metrics, helpers$position_metrics(input$truth, result$raw_positions)))
    if (result$status != "failed") {
      stopifnot(length(result$positions) == helpers$design$N,
        all(result$positions >= -1e-12 & result$positions <= 1 + 1e-12),
        min(diff(result$elbo_trace)) >= -1e-7,
        length(result$elbo_trace) == result$iterations + 1L)
      if (result$converged) stopifnot(abs(tail(diff(result$elbo_trace), 1L)) / length(input$X) < helpers$design$tolerance)
    }
    rows[[length(rows) + 1L]] <- data.frame(noise_sd = sd, method = method,
      k_used = result$k_used, raw_spearman = result$raw_metrics[["spearman"]],
      final_spearman = result$metrics[["spearman"]], status = result$status,
      converged = result$converged, iterations = result$iterations,
      warnings = paste(result$warnings, collapse = " | "), warning_count = length(result$warnings),
      graph_components = result$graph_components, retained_samples = result$graph_keep_count,
      initializer = if (is.null(result$initialization)) NA_character_ else result$initialization$method_used,
      baseline_reused = sd == .25)
    if (method == "auto_kmin") {
      # Independent connectivity reference with stable ties and no package graph helper.
      distances <- as.matrix(dist(input$X)); diag(distances) <- Inf
      connected <- function(k) {
        neighbors <- t(vapply(seq_len(nrow(input$X)), function(index)
          order(distances[index, ], seq_len(nrow(input$X)))[seq_len(k)], integer(k)))
        graph <- igraph::graph_from_edgelist(cbind(rep(seq_len(nrow(input$X)), each = k),
          as.vector(t(neighbors))), directed = FALSE)
        igraph::is_connected(graph)
      }
      stopifnot(connected(result$k_used), result$k_used == 1L || !connected(result$k_used - 1L))
      graph_checks[[length(graph_checks) + 1L]] <- c(noise_sd = sd, k_used = result$k_used)
      raw <- result$raw_positions
      if (cor(input$truth, raw, method = "spearman") < 0) raw <- 1 - raw
      for (stage in c("True positions", "Auto Isomap")) points[[length(points) + 1L]] <- data.frame(
        sample = rownames(input$X), stage = stage,
        level = sprintf("Noise SD = %.2f\nAuto k=%d; raw |Spearman|=%.3f", sd, result$k_used, result$raw_metrics[["spearman"]]),
        feature_1 = input$X[, 1L], feature_4 = input$X[, 4L],
        color_position = if (stage == "True positions") input$truth else raw)
    }
  }
}
rows <- do.call(rbind, rows)
points <- do.call(rbind, points)
write.csv(rows, file.path(output, "scores.csv"), row.names = FALSE)
points$stage <- factor(points$stage, levels = c("True positions", "Auto Isomap"))
points$level <- factor(points$level, levels = unique(points$level))
write.csv(points, file.path(output, "plotted_feature_points.csv"), row.names = FALSE)
plotted <- rbind(data.frame(noise_sd = rows$noise_sd, method = rows$method, stage = "Raw Isomap", score = rows$raw_spearman),
  data.frame(noise_sd = rows$noise_sd, method = rows$method, stage = "Final MPCurve", score = rows$final_spearman))
plotted$stage <- factor(plotted$stage, levels = c("Raw Isomap", "Final MPCurve"))
plotted$method <- factor(plotted$method, levels = methods, labels = c("Auto kmin", "Fixed k=15", "Fixed k=10"))
write.csv(plotted, file.path(output, "plotted_scores.csv"), row.names = FALSE)
theme <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9), legend.position = "bottom")
save_plot <- function(plot, stem, width, height) {
  for (extension in c("png", "pdf")) ggsave(file.path(output, paste0(stem, ".", extension)),
    plot, width = width, height = height, dpi = 180, bg = "white")
}
caption <- function(text, width) paste(strwrap(text, width), collapse = "\n")
score_plot <- ggplot(plotted, aes(noise_sd, score, color = method)) +
  geom_line(linewidth = .8) + geom_point(size = 2.5) + facet_wrap(~ stage, nrow = 1) +
  scale_x_log10(breaks = sort(levels), labels = sprintf("%.2f", sort(levels))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_color_manual(values = c("Auto kmin" = "#0072B2", "Fixed k=15" = "#D55E00", "Fixed k=10" = "#009E73"), name = NULL) +
  labs(x = "Observation-noise SD (logarithmic axis)", y = "Absolute Spearman correlation",
    title = "Noise reduction on fixed trajectories: replicate 6, P=4",
    subtitle = "Same four signals, 200 true positions, and standard-normal errors at every level",
    caption = caption(paste("Figure 1. Raw Isomap and final MPCurve recovery at SD=0.25, 0.10, 0.05, and 0.01.",
      "Only error amplitude changes. Absolute Spearman allows global reversal. The original SD=0.25 fits are reused;",
      "all reduced-noise settings use the same frozen package and fitting controls. One fixed case; no uncertainty intervals."), 116)) + theme
save_plot(score_plot, "noise_recovery", 11, 6.5)
feature_plot <- ggplot(points, aes(feature_1, feature_4, color = color_position)) +
  geom_point(size = 1.5, alpha = .85) + facet_grid(stage ~ level) + coord_equal() +
  scale_color_viridis_c(name = "Position along ordering", limits = c(0, 1), option = "D", end = .95) +
  labs(x = "Observed feature 1", y = "Observed feature 4",
    title = "Feature 1 versus feature 4 as noise decreases",
    subtitle = "Top: colors show truth | Bottom: colors show raw automatic Isomap positions computed from all four features",
    caption = caption(paste("Figure 2. Same trajectories and noise directions, with decreasing noise SD from left to right.",
      "Axes are observed features 1 and 4; Isomap uses all four observed features. Top/bottom panels contain identical samples",
      "at each noise level. Raw Isomap positions precede binning and MPCurve fitting; only global reflection is allowed for display.",
      "Axis limits and [0,1] color limits are shared across panels."), 136)) + theme
save_plot(feature_plot, "feature1_feature4_noise", 15, 9)
stopifnot(identical(record$source_input_sha256, helpers$hash_file(source_input_path)),
  identical(record$source_fits_sha256, setNames(vapply(source_fits, helpers$hash_file, character(1)), methods)))
verification <- list(verified_at_utc = format(Sys.time(), tz = "UTC"), levels = length(levels),
  new_fits = 9L, reused_fits = 3L, converged = sum(rows$converged),
  failures = sum(rows$status == "failed"), warnings = sum(rows$warning_count),
  plotted_scores = nrow(plotted), plotted_feature_points = nrow(points), graph_checks = graph_checks,
  source_input_and_fits_preserved = TRUE)
helpers$save_atomic(verification, file.path(output, "verification.rds"))
writeLines(capture.output(str(verification)), file.path(output, "verification.txt"))
print(rows, digits = 6)
