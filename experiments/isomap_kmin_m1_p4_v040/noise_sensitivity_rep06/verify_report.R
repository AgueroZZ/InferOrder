# Verify the final scientific report against saved fits and plotted data.
helpers <- new.env(parent = globalenv())
source("experiments/isomap_kmin_m1_p4_v040/common.R", local = helpers)
output <- file.path(helpers$study, "noise_sensitivity_rep06")
record <- readRDS(file.path(output, "provenance.rds"))
scores <- read.csv(file.path(output, "scores.csv"))
plotted <- read.csv(file.path(output, "plotted_scores.csv"))
points <- read.csv(file.path(output, "plotted_feature_points.csv"))
stopifnot(nrow(scores) == 12L, nrow(plotted) == 24L, nrow(points) == 1600L,
  sum(scores$baseline_reused) == 3L, all(scores$converged), all(scores$status == "converged"))
for (index in seq_len(nrow(scores))) {
  row <- scores[index, ]
  level_root <- file.path(output, sprintf("sd_%03d", round(row$noise_sd * 1000)))
  input <- readRDS(file.path(level_root, "inputs", "rep06.rds"))
  result <- readRDS(file.path(level_root, "results", paste0("rep06_", row$method, ".rds")))
  stopifnot(length(result$raw_warnings) == 0L, length(result$warnings) == 0L,
    identical(result$initialization$method_used, "isomap"), result$graph_keep_count == 200L,
    result$graph_components == 1L,
    abs(row$raw_spearman - abs(cor(input$truth, result$raw_positions, method = "spearman"))) < 1e-14,
    abs(row$final_spearman - abs(cor(input$truth, result$positions, method = "spearman"))) < 1e-14)
  selected <- plotted$noise_sd == row$noise_sd & plotted$method ==
    switch(row$method, auto_kmin = "Auto kmin", fixed_k15 = "Fixed k=15", fixed_k10 = "Fixed k=10")
  stopifnot(sum(selected) == 2L,
    abs(plotted$score[selected & plotted$stage == "Raw Isomap"] - row$raw_spearman) < 1e-14,
    abs(plotted$score[selected & plotted$stage == "Final MPCurve"] - row$final_spearman) < 1e-14)
  if (row$method == "auto_kmin") {
    observed <- points[grepl(sprintf("Noise SD = %.2f", row$noise_sd), points$level, fixed = TRUE), ]
    truth_points <- observed[observed$stage == "True positions", ]
    auto_points <- observed[observed$stage == "Auto Isomap", ]
    raw <- result$raw_positions
    if (cor(raw, input$truth, method = "spearman") < 0) raw <- 1 - raw
    stopifnot(nrow(truth_points) == 200L, nrow(auto_points) == 200L,
      identical(truth_points$sample, rownames(input$X)),
      max(abs(truth_points$feature_1 - input$X[, 1L])) < 1e-13,
      max(abs(truth_points$feature_4 - input$X[, 4L])) < 1e-13,
      max(abs(truth_points$color_position - input$truth)) < 1e-14,
      max(abs(auto_points$color_position - raw)) < 1e-14,
      identical(truth_points$feature_1, auto_points$feature_1),
      identical(truth_points$feature_4, auto_points$feature_4))
  }
}
source_fits <- setNames(vapply(helpers$design$methods, function(method)
  helpers$result_path(6L, method), character(1)), helpers$design$methods)
stopifnot(identical(record$source_input_sha256, helpers$hash_file(helpers$input_path(6L))),
  identical(record$source_fits_sha256, setNames(vapply(source_fits, helpers$hash_file, character(1)), names(source_fits))),
  identical(record$script_sha256, helpers$hash_file(file.path(output, "run.R"))),
  identical(record$design_file_sha256, helpers$hash_file(file.path(output, "DESIGN.md"))),
  identical(record$parent_provenance$archive_sha256, helpers$hash_file(helpers$archive_path)))
writeLines(c(paste("Verified at UTC:", format(Sys.time(), tz = "UTC")),
  "PASS: 12 converged fits, 9 new and 3 reused; zero raw/fitting warnings or fallback.",
  "PASS: all 24 plotted scores match source fits and absolute Spearman calculations.",
  "PASS: all 1600 feature coordinates/colors match input and original Isomap vectors.",
  "PASS: source input/fits, frozen archive, and original fit/design script hashes preserved.",
  paste("Verification script SHA256:", helpers$hash_file(file.path(output, "verify_report.R")))),
  file.path(output, "report_verification.txt"))
cat("Final report verification passed.\n")
