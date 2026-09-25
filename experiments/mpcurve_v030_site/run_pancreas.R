#!/usr/bin/env Rscript

if (!file.exists("_workflowr.yml")) {
  stop("Run this script from the InferOrder repository root.")
}
if (!requireNamespace("MPCurver", quietly = TRUE) ||
    as.character(utils::packageVersion("MPCurver")) != "0.3.0") {
  stop("MPCurver 0.3.0 must be installed before generating site results.")
}

library(MPCurver)
source("code/plot_ordering.R")
resume <- "--resume" %in% commandArgs(trailingOnly = TRUE)

out_dir <- file.path("experiments", "mpcurve_v030_site", "results", "pancreas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

factor_file <- file.path("data", "loading_order", "pancreas_factors.RData")
metadata_file <- file.path("data", "loading_order", "pancreas.RData")
input_env <- new.env(parent = globalenv())
load(factor_file, envir = input_env)
load(metadata_file, envir = input_env)

set.seed(1)
cells <- subsample_cell_types(input_env$sample_info$celltype, n = 500)
factor_columns <- c(3L, 8L, 9L, 12L, 17L, 18L, 20L, 21L)
selected <- which(as.character(input_env$sample_info$tech[cells]) == "inDrop3")
sample_ids <- cells[selected]
X <- as.matrix(input_env$fl_snmf_ldf$L[sample_ids, factor_columns, drop = FALSE])
storage.mode(X) <- "double"
stopifnot(identical(dim(X), c(865L, 8L)), all(is.finite(X)),
          !anyDuplicated(rownames(X)))

metadata <- data.frame(
  sample_id = rownames(X),
  cell_type = as.character(input_env$sample_info$celltype[sample_ids]),
  technology = as.character(input_env$sample_info$tech[sample_ids]),
  stringsAsFactors = FALSE
)
stopifnot(identical(metadata$technology, rep("inDrop3", nrow(X))))
utils::write.csv(metadata, file.path(out_dir, "selected_samples.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  file = c(factor_file, metadata_file),
  md5 = unname(tools::md5sum(c(factor_file, metadata_file))),
  stringsAsFactors = FALSE
), file.path(out_dir, "input_hashes.csv"), row.names = FALSE)
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
           file.path(out_dir, "session_info.txt"))

settings <- data.frame(
  model = c("M1", "M2"),
  seed = c(20260925L, 20260927L),
  intrinsic_dim = c(1L, 2L),
  method = c("PCA", "fiedler"),
  K = c(50L, 50L),
  rw_q = c(2L, 2L),
  ridge = c(0, 0),
  iter = c(150L, NA_integer_),
  n_outer = c(NA_integer_, 25L),
  max_converge_iter = c(NA_integer_, 100L),
  partition_init = c(NA_character_, "similarity"),
  stringsAsFactors = FALSE
)
utils::write.csv(settings, file.path(out_dir, "fit_settings.csv"), row.names = FALSE)

cat("Selected pancreatic loading matrix:", nrow(X), "samples by", ncol(X), "factors\n")
fit1_path <- file.path(out_dir, "fit_m1.rds")
fit2_path <- file.path(out_dir, "fit_m2.rds")
if (resume) {
  cat("Loading saved M = 1 fit...\n")
  fit1 <- readRDS(fit1_path)$fit
  stopifnot(isTRUE(all.equal(fit1$data, X, check.attributes = FALSE)))
} else {
  cat("Fitting M = 1...\n")
  set.seed(20260925)
  fit1 <- fit_mpcurve(X, intrinsic_dim = 1, method = "PCA", K = 50,
                     rw_q = 2, ridge = 0, iter = 150, verbose = FALSE)
}
stopifnot(inherits(fit1, "mpcurve"), identical(fit1$algorithm, "cavi"))
position1 <- as.numeric(fit1$locations$mean$pseudotime)
stopifnot(length(position1) == nrow(X), all(is.finite(position1)))
saveRDS(list(fit = fit1, position = position1), fit1_path)

if (resume) {
  cat("Loading saved M = 2 fit...\n")
  fit2 <- readRDS(fit2_path)$fit
  stopifnot(isTRUE(all.equal(fit2$data, X, check.attributes = FALSE)))
} else {
  cat("Fitting M = 2...\n")
  set.seed(20260927)
  fit2 <- fit_mpcurve(X, intrinsic_dim = 2, method = "fiedler",
                     partition_init = "similarity", K = 50, rw_q = 2,
                     ridge = 0, n_outer = 25L, max_converge_iter = 100L,
                     verbose = FALSE)
}
stopifnot(inherits(fit2, "mpcurve"), identical(fit2$algorithm, "cavi"))
positions2 <- lapply(fit2$locations, function(item) as.numeric(item$mean$pseudotime))
stopifnot(length(positions2) == 2L,
          all(vapply(positions2, length, integer(1)) == nrow(X)),
          all(is.finite(unlist(positions2))))
weights2 <- as.matrix(fit2$partition$pi_weights)
stopifnot(all(dim(weights2) == c(ncol(X), 2L)),
          all(is.finite(weights2)), max(abs(rowSums(weights2) - 1)) < 1e-8)
saveRDS(list(fit = fit2, positions = positions2, weights = weights2),
        fit2_path)

t1_obj <- fit2$objective_history[fit2$temperature_history == 1]
summary <- data.frame(
  model = c("M1", "M2"),
  samples = nrow(X),
  factors = ncol(X),
  K = 50L,
  objective_records = c(length(fit1$elbo_trace), length(t1_obj)),
  min_within_objective_increment = c(
    if (length(fit1$elbo_trace) > 1L) min(diff(fit1$elbo_trace)) else NA_real_,
    if (length(t1_obj) > 1L) min(diff(t1_obj)) else NA_real_
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(out_dir, "fit_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  factor = paste0("Factor ", factor_columns),
  ordering_A = weights2[, 1],
  ordering_B = weights2[, 2],
  assignment = as.character(fit2$partition$assign),
  stringsAsFactors = FALSE
), file.path(out_dir, "factor_assignments.csv"), row.names = FALSE)

position_table <- cbind(metadata,
                        data.frame(M1 = position1,
                                   M2_A = positions2[[1]],
                                   M2_B = positions2[[2]]))
utils::write.csv(position_table, file.path(out_dir, "sample_positions.csv"), row.names = FALSE)

plot_cell_type_positions <- function(position, title) {
  type_levels <- names(sort(table(metadata$cell_type), decreasing = TRUE))
  y <- match(metadata$cell_type, type_levels)
  graphics::plot(position, jitter(y, amount = 0.22), pch = 16, cex = 0.55,
                 col = grDevices::adjustcolor("#2E8A99", alpha.f = 0.38),
                 xlab = "Inferred pseudotime", ylab = "Cell type",
                 yaxt = "n", ylim = c(length(type_levels) + 0.6, 0.4),
                 main = title)
  graphics::axis(2, at = seq_along(type_levels), labels = type_levels,
                 las = 2, cex.axis = 0.65)
  graphics::grid(nx = NA, ny = length(type_levels), col = "grey90")
}

grDevices::png(file.path(out_dir, "cell_type_positions.png"),
               width = 1500, height = 1500, res = 150)
op <- graphics::par(mfrow = c(3, 1), mar = c(4, 9, 3, 1))
plot_cell_type_positions(position1, "M = 1")
plot_cell_type_positions(positions2[[1]], "M = 2, ordering A")
plot_cell_type_positions(positions2[[2]], "M = 2, ordering B")
graphics::par(op)
grDevices::dev.off()

plot_loading_heatmap <- function(position, title) {
  ordered_X <- X[order(position), , drop = FALSE]
  colours <- grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)
  graphics::image(seq_len(nrow(ordered_X)), seq_len(ncol(ordered_X)),
                  ordered_X, col = colours, zlim = range(X),
                  xlab = "Cells ordered by inferred pseudotime",
                  ylab = "Semi-NMF factor column", yaxt = "n", main = title)
  graphics::axis(2, at = seq_len(ncol(X)), labels = factor_columns, las = 2)
  graphics::box()
}

grDevices::png(file.path(out_dir, "ordered_loading_heatmaps.png"),
               width = 1500, height = 1500, res = 150)
op <- graphics::par(mfrow = c(3, 1), mar = c(4, 6, 3, 1))
plot_loading_heatmap(position1, "M = 1")
plot_loading_heatmap(positions2[[1]], "M = 2, ordering A")
plot_loading_heatmap(positions2[[2]], "M = 2, ordering B")
graphics::par(op)
grDevices::dev.off()

grDevices::png(file.path(out_dir, "factor_assignments.png"),
               width = 1200, height = 680, res = 150)
graphics::barplot(t(weights2), col = c("#2E8A99", "#E76F51"), border = NA,
                  names.arg = factor_columns, xlab = "Semi-NMF factor column",
                  ylab = "Assignment probability", main = "M = 2 feature assignments")
weight_labels <- colnames(weights2)
if (is.null(weight_labels)) weight_labels <- c("A", "B")
graphics::legend("topright", legend = weight_labels,
                 fill = c("#2E8A99", "#E76F51"), bty = "n")
grDevices::dev.off()

grDevices::png(file.path(out_dir, "objective_traces.png"),
               width = 1500, height = 680, res = 150)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
graphics::plot(fit1$elbo_trace, type = "l", lwd = 2, col = "#2E8A99",
               xlab = "Recorded iteration", ylab = "ELBO", main = "M = 1")
graphics::plot(seq_along(t1_obj), t1_obj, type = "l", lwd = 2, col = "#E76F51",
               xlab = "Recorded T = 1 update", ylab = "Structural objective",
               main = "M = 2")
graphics::par(op)
grDevices::dev.off()

cat("Saved M = 1 and M = 2 results to", out_dir, "\n")
