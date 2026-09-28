#!/usr/bin/env Rscript
# Describe correlation patterns used by the unchanged similarity initialization.
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
summary_dir <- file.path(study_dir, "main_summary")
dir.create(summary_dir, showWarnings = FALSE)
manifest <- active_manifest()

correlation_summary <- function(values) {
  stopifnot(length(values) > 0L, all(is.finite(values)))
  c(median = stats::median(values),
    q25 = unname(stats::quantile(values, 0.25, type = 7)),
    q75 = unname(stats::quantile(values, 0.75, type = 7)),
    mean = mean(values), pairs = length(values))
}

describe_input <- function(data, anchors, row, family) {
  similarity <- abs(stats::cor(data$X, method = "spearman"))
  same_group <- outer(unname(data$true_assign), unname(data$true_assign), "==")
  upper <- upper.tri(similarity, diag = FALSE)
  values <- list(within = similarity[upper & same_group],
    between = similarity[upper & !same_group],
    anchor = unlist(lapply(anchors, function(anchor) {
      others <- which(data$true_assign == data$true_assign[anchor])
      similarity[anchor, others[others != anchor]]
    }), use.names = FALSE))
  output <- row
  output$family <- family
  output$input_hash <- data$input_hash
  for (relation in names(values)) {
    statistics <- correlation_summary(values[[relation]])
    for (statistic in names(statistics)) {
      output[[paste(relation, statistic, sep = "_")]] <- unname(statistics[[statistic]])
    }
  }
  output
}

rows <- vector("list", nrow(manifest) * 2L)
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  smooth <- readRDS(file.path(study_dir, "data", paste0(row$id, ".rds")))
  monotone <- readRDS(file.path(baseline_dir, "data", paste0(row$id, ".rds")))
  stopifnot(identical(smooth$true_assign, monotone$true_assign),
    identical(dimnames(smooth$X), dimnames(monotone$X)))
  rows[[2L * i - 1L]] <- describe_input(monotone, smooth$anchor_indices, row, "monotone")
  rows[[2L * i]] <- describe_input(smooth, smooth$anchor_indices, row, "smooth")
}
per_dataset <- do.call(rbind, rows)
rownames(per_dataset) <- NULL
stopifnot(nrow(per_dataset) == 180L, all(table(per_dataset$id) == 2L))
write.csv(per_dataset, file.path(summary_dir, "initialization_correlations.csv"), row.names = FALSE)

summarize_dataset_medians <- function(data, groups) {
  combinations <- unique(data[, groups, drop = FALSE])
  output <- list()
  index <- 0L
  for (i in seq_len(nrow(combinations))) {
    included <- rep(TRUE, nrow(data))
    for (group in groups) included <- included & data[[group]] == combinations[[group]][i]
    subset <- data[included, , drop = FALSE]
    for (relation in c("within", "between", "anchor")) {
      values <- subset[[paste0(relation, "_median")]]
      index <- index + 1L
      output[[index]] <- cbind(combinations[i, , drop = FALSE],
        data.frame(relation = relation, datasets = length(values),
          median_of_dataset_medians = stats::median(values),
          q25_of_dataset_medians = unname(stats::quantile(values, 0.25, type = 7)),
          q75_of_dataset_medians = unname(stats::quantile(values, 0.75, type = 7)),
          mean_of_dataset_medians = mean(values),
          sd_of_dataset_medians = stats::sd(values),
          minimum_dataset_median = min(values), maximum_dataset_median = max(values)))
    }
  }
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}
condition <- summarize_dataset_medians(per_dataset, c("family", "true_M", "snr"))
overall <- summarize_dataset_medians(per_dataset, "family")
by_snr <- summarize_dataset_medians(per_dataset, c("family", "snr"))
write.csv(condition, file.path(summary_dir, "initialization_correlations_summary.csv"), row.names = FALSE)
write.csv(overall, file.path(summary_dir, "initialization_correlations_overall.csv"), row.names = FALSE)
write.csv(by_snr, file.path(summary_dir, "initialization_correlations_by_snr.csv"), row.names = FALSE)

paired_rows <- list()
index <- 0L
for (i in seq_len(nrow(manifest))) {
  original <- per_dataset[per_dataset$id == manifest$id[i] & per_dataset$family == "monotone", ]
  replacement <- per_dataset[per_dataset$id == manifest$id[i] & per_dataset$family == "smooth", ]
  stopifnot(nrow(original) == 1L, nrow(replacement) == 1L)
  for (relation in c("within", "between", "anchor")) {
    index <- index + 1L
    old_value <- original[[paste0(relation, "_median")]]
    new_value <- replacement[[paste0(relation, "_median")]]
    paired_rows[[index]] <- cbind(manifest[i, , drop = FALSE],
      data.frame(relation = relation, monotone_median = old_value,
        smooth_median = new_value, difference = new_value - old_value))
  }
}
paired <- do.call(rbind, paired_rows)
rownames(paired) <- NULL
paired_summary <- do.call(rbind, lapply(c("within", "between", "anchor"), function(relation) {
  values <- paired$difference[paired$relation == relation]
  data.frame(relation = relation, datasets = length(values), lower = sum(values < 0),
    higher = sum(values > 0), unchanged = sum(values == 0),
    median_difference = stats::median(values), mean_difference = mean(values))
}))
write.csv(paired, file.path(summary_dir, "initialization_correlations_paired.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(summary_dir, "initialization_correlations_paired_summary.csv"), row.names = FALSE)

jsonlite::write_json(list(quantity = "absolute Spearman correlation of noisy observed features",
  within_pairs = "all unordered feature pairs sharing a true ordering, diagonal excluded",
  between_pairs = "all unordered feature pairs assigned to different true orderings",
  anchor_pairs = "each retained monotone feature paired with the other features of its true ordering",
  summary_unit = "median correlation within each dataset; datasets weighted equally",
  quantile_type = 7L, datasets_per_family = 90L,
  design_hash = design_hash, descriptive = TRUE),
  file.path(summary_dir, "initialization_correlations_provenance.json"),
  auto_unbox = TRUE, pretty = TRUE)
cat("Descriptive initialization correlations summarized for 90 paired datasets.\n")
print(overall[, c("family", "relation", "median_of_dataset_medians")], row.names = FALSE)
print(by_snr[, c("family", "snr", "relation", "median_of_dataset_medians")], row.names = FALSE)
