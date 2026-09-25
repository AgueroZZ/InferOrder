#!/usr/bin/env Rscript

if (!file.exists("_workflowr.yml")) {
  stop("Run this script from the InferOrder repository root.")
}
if (!requireNamespace("MPCurver", quietly = TRUE) ||
    as.character(utils::packageVersion("MPCurver")) != "0.3.0") {
  stop("MPCurver 0.3.0 must be installed before generating site results.")
}

library(MPCurver)

out_dir <- file.path("experiments", "mpcurve_v030_site", "results", "simulations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
           file.path(out_dir, "session_info.txt"))

sim1 <- simulate_swiss_roll_1d_2d(
  n = 1000,
  t_range = c(1.5 * pi, 4.5 * pi),
  sigma = 0.12,
  seed = 1
)
X1 <- as.matrix(sim1$obs)
set.seed(42)
fit1 <- fit_mpcurve(X1, method = "fiedler", K = 60, iter = 120)
stopifnot(inherits(fit1, "mpcurve"), identical(fit1$algorithm, "cavi"))

position1 <- as.numeric(fit1$gamma %*% seq_len(ncol(fit1$gamma)))
rho1 <- abs(stats::cor(position1, sim1$t, method = "spearman"))
elbo_diff1 <- diff(fit1$elbo_trace)
summary1 <- data.frame(
  n = nrow(X1),
  d = ncol(X1),
  K = ncol(fit1$gamma),
  abs_spearman = rho1,
  iterations = length(fit1$elbo_trace),
  min_elbo_increment = if (length(elbo_diff1)) min(elbo_diff1) else NA_real_,
  stringsAsFactors = FALSE
)
stopifnot(all(is.finite(unlist(summary1[1, 1:4]))))
saveRDS(list(simulation = sim1, fit = fit1, position = position1,
             summary = summary1), file.path(out_dir, "simulation_m1.rds"))
utils::write.csv(summary1, file.path(out_dir, "simulation_m1_summary.csv"), row.names = FALSE)

pal <- grDevices::colorRampPalette(c("#1F3A5F", "#2E8A99", "#F2C14E", "#E76F51"))(256)
to_colors <- function(x) pal[cut(x, breaks = 256, labels = FALSE)]
grDevices::png(file.path(out_dir, "simulation_m1_recovery.png"), width = 1500, height = 680, res = 150)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
graphics::plot(X1, col = to_colors(sim1$t), pch = 19, cex = 0.55,
               xlab = "Feature 1", ylab = "Feature 2", main = "Generating position")
graphics::plot(X1, col = to_colors(position1), pch = 19, cex = 0.55,
               xlab = "Feature 1", ylab = "Feature 2", main = "Inferred position")
graphics::par(op)
grDevices::dev.off()

grDevices::png(file.path(out_dir, "simulation_m1_elbo.png"), width = 1050, height = 620, res = 150)
graphics::plot(fit1$elbo_trace, type = "l", lwd = 2, col = "#2E8A99",
               xlab = "Recorded iteration", ylab = "ELBO", main = "Single-ordering CAVI objective")
grDevices::dev.off()

sim2 <- simulate_dual_trajectory(
  n = 500,
  d1 = 10,
  d2 = 10,
  d_noise = 0,
  sigma = 0.05,
  trajectory_family = c("quadratic", "monotone"),
  seed = 1
)
set.seed(42)
fit2 <- fit_mpcurve(sim2$X, intrinsic_dim = 2, method = "fiedler", K = 50,
                   rw_q = 2, ridge = 0, verbose = FALSE)
stopifnot(inherits(fit2, "mpcurve"), identical(fit2$algorithm, "cavi"))

truth2 <- as.character(sim2$true_assign)
pred2 <- as.character(fit2$partition$assign)
stopifnot(length(truth2) == ncol(sim2$X), length(pred2) == length(truth2))
direct_accuracy2 <- mean(pred2 == truth2)
swapped_accuracy2 <- mean(pred2 != truth2)
accuracy2 <- max(direct_accuracy2, swapped_accuracy2)
weights2 <- as.matrix(fit2$partition$pi_weights)
stopifnot(all(dim(weights2) == c(ncol(sim2$X), 2L)))
labels2 <- colnames(weights2)
if (is.null(labels2)) labels2 <- names(fit2$locations)
if (is.null(labels2)) labels2 <- c("A", "B")
colnames(weights2) <- labels2

truth_order_for_fit <- if (direct_accuracy2 >= swapped_accuracy2) c(1L, 2L) else c(2L, 1L)
rho2 <- vapply(seq_along(labels2), function(m) {
  p <- as.numeric(fit2$locations[[m]]$mean$pseudotime)
  abs(stats::cor(p, sim2$latent_positions[, truth_order_for_fit[m]], method = "spearman"))
}, numeric(1))

temperature2 <- fit2$temperature_history
objective2 <- fit2$objective_history
t1_obj2 <- objective2[temperature2 == 1]
summary2 <- data.frame(
  n = nrow(sim2$X),
  d = ncol(sim2$X),
  K = ncol(fit2$gamma[[1]]),
  partition_accuracy = accuracy2,
  abs_spearman_A = rho2[1],
  abs_spearman_B = rho2[2],
  t1_objective_records = length(t1_obj2),
  min_t1_increment = if (length(t1_obj2) > 1L) min(diff(t1_obj2)) else NA_real_,
  stringsAsFactors = FALSE
)
assignments2 <- data.frame(
  feature = colnames(sim2$X), truth = truth2, predicted = pred2,
  weight_A = weights2[, 1], weight_B = weights2[, 2],
  stringsAsFactors = FALSE
)
stopifnot(all(is.finite(unlist(summary2[1, 1:6]))),
          all(is.finite(weights2)), max(abs(rowSums(weights2) - 1)) < 1e-8)
saveRDS(list(simulation = sim2, fit = fit2, summary = summary2,
             assignments = assignments2,
             truth_order_for_fit = truth_order_for_fit),
        file.path(out_dir, "simulation_m2.rds"))
utils::write.csv(summary2, file.path(out_dir, "simulation_m2_summary.csv"), row.names = FALSE)
utils::write.csv(assignments2, file.path(out_dir, "simulation_m2_assignments.csv"), row.names = FALSE)

grDevices::png(file.path(out_dir, "simulation_m2_assignments.png"), width = 1350, height = 680, res = 150)
graphics::barplot(t(weights2), col = c("#2E8A99", "#E76F51"), border = NA,
                  names.arg = truth2, las = 2,
                  xlab = "Feature (label shows generating group)",
                  ylab = "Assignment probability", main = "Inferred feature assignments")
graphics::legend("topright", legend = labels2, fill = c("#2E8A99", "#E76F51"), bty = "n")
grDevices::dev.off()

grDevices::png(file.path(out_dir, "simulation_m2_objective.png"), width = 1050, height = 620, res = 150)
graphics::plot(seq_along(t1_obj2), t1_obj2, type = "l", lwd = 2, col = "#2E8A99",
               xlab = "Recorded T = 1 update", ylab = "Structural objective",
               main = "Two-ordering fixed-temperature objective")
grDevices::dev.off()

cat("M = 1 absolute Spearman:", rho1, "\n")
cat("M = 2 partition accuracy:", accuracy2, "\n")
cat("M = 2 ordering absolute Spearman:", paste(rho2, collapse = ", "), "\n")
