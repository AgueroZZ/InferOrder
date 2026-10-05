# Compare all feature pairs in the previously inspected four-feature replicate.
source("experiments/isomap_kmin_m1_p4_v040/common.R")
suppressPackageStartupMessages(library(ggplot2))
replication <- 6L
input <- readRDS(input_path(replication))
stopifnot(identical(input$input_sha256, hash_object(input$X)),
  identical(input$X[, 1:2, drop = FALSE], parent_input(replication)$X),
  identical(input$truth, parent_input(replication)$truth))
feature_names <- paste("Feature", seq_len(design$P))
observed <- list()
curves <- list()
diagonal <- list()
limits <- lapply(seq_len(design$P), function(feature)
  range(c(input$X[, feature], input$dense_signal[, feature])))
for (y_feature in seq_len(design$P)) {
  for (x_feature in seq_len(design$P)) {
    if (y_feature > x_feature) {
      observed[[length(observed) + 1L]] <- data.frame(
        sample = rownames(input$X), x_index = x_feature, y_index = y_feature,
        x_feature = feature_names[x_feature], y_feature = feature_names[y_feature],
        x = input$X[, x_feature], y = input$X[, y_feature], position = input$truth)
    } else if (y_feature < x_feature) {
      curves[[length(curves) + 1L]] <- data.frame(
        grid_index = seq_along(input$grid), x_index = x_feature, y_index = y_feature,
        x_feature = feature_names[x_feature], y_feature = feature_names[y_feature],
        x = input$dense_signal[, x_feature], y = input$dense_signal[, y_feature],
        position = input$grid)
    } else {
      diagonal[[length(diagonal) + 1L]] <- data.frame(
        x_feature = feature_names[x_feature], y_feature = feature_names[y_feature],
        x = mean(limits[[x_feature]]), y = mean(limits[[y_feature]]),
        xmin = limits[[x_feature]][1L], xmax = limits[[x_feature]][2L],
        ymin = limits[[y_feature]][1L], ymax = limits[[y_feature]][2L],
        label = if (x_feature <= 2L) "Original\nfeature" else "Added\nfeature")
    }
  }
}
observed <- do.call(rbind, observed)
curves <- do.call(rbind, curves)
diagonal <- do.call(rbind, diagonal)
set_features <- function(data) {
  data$x_feature <- factor(data$x_feature, levels = feature_names)
  data$y_feature <- factor(data$y_feature, levels = feature_names)
  data
}
observed <- set_features(observed)
curves <- set_features(curves)
diagonal <- set_features(diagonal)
stopifnot(nrow(observed) == design$N * choose(design$P, 2L),
  nrow(curves) == length(input$grid) * choose(design$P, 2L),
  all(observed$y_index > observed$x_index), all(curves$y_index < curves$x_index))
for (row in seq_len(nrow(observed))) {
  index <- match(observed$sample[row], rownames(input$X))
  stopifnot(observed$x[row] == input$X[index, observed$x_index[row]],
    observed$y[row] == input$X[index, observed$y_index[row]],
    observed$position[row] == input$truth[index])
}
write.csv(observed, file.path(study, "rep06_scatter_matrix_observed.csv"), row.names = FALSE)
write.csv(curves, file.path(study, "rep06_scatter_matrix_curves.csv"), row.names = FALSE)
caption <- paste("Figure 4. Scatterplot matrix of the four features in replicate 6 (N=200, M=1).",
  "Columns specify horizontal features and rows specify vertical features. Lower triangle: observed samples with Gaussian noise SD=0.25.",
  "Upper triangle: noiseless generating curves on a 2,001-point grid. Colors show true position along the trajectory.",
  "Features 1 and 2 are preserved from P=2; features 3 and 4 were added. Each column and row shares its feature-specific axis limits.")
plot <- ggplot() +
  geom_rect(data = diagonal, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#EEEEEE", color = NA) +
  geom_point(data = observed, aes(x, y, color = position), size = 1.5, alpha = .8) +
  geom_path(data = curves, aes(x, y, color = position), linewidth = .85) +
  geom_text(data = diagonal, aes(x, y, label = label), color = "#555555", size = 4.5) +
  facet_grid(y_feature ~ x_feature, scales = "free", labeller = labeller(
    x_feature = function(labels) paste("X:", labels),
    y_feature = function(labels) paste("Y:", labels))) +
  scale_color_viridis_c(name = "True position", limits = c(0, 1),
    breaks = seq(0, 1, .25), option = "D", end = .95) +
  labs(x = "Horizontal feature (column)", y = "Vertical feature (row)",
    title = "Feature scatterplot matrix: replicate 6, P=4",
    subtitle = "Lower triangle: noisy observations | Upper triangle: noiseless generating curves",
    caption = paste(strwrap(caption, 115), collapse = "\n")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.spacing = grid::unit(.7, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom", legend.key.width = grid::unit(1.5, "cm"),
    plot.caption = element_text(hjust = 0, size = 9))
for (extension in c("png", "pdf")) ggsave(
  file.path(study, paste0("rep06_feature_scatter_matrix.", extension)), plot,
  width = 12, height = 11.5, dpi = 180, bg = "white")
save_atomic(list(replication = replication, input_file_sha256 = hash_file(input_path(replication)),
  script_sha256 = hash_file(file.path(study, "plot_scatter_matrix.R")),
  parent_provenance = readRDS(file.path(study, "provenance.rds")),
  observed_rows = nrow(observed), curve_rows = nrow(curves),
  diagonal_features = feature_names, created_at_utc = format(Sys.time(), tz = "UTC")),
  file.path(study, "rep06_scatter_matrix_provenance.rds"))
cat(sprintf("Verified %d observed points and plotted %d generating-curve points across six feature pairs.\n",
  nrow(observed), nrow(curves)))
