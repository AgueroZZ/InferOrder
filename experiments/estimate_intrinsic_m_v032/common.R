# Shared study design, deterministic data generation, and checkpoint utilities.
options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
study_dir <- normalizePath("experiments/estimate_intrinsic_m_v032")
.libPaths(c(file.path(study_dir, "library"), .libPaths()))
stopifnot(as.character(utils::packageVersion("MPCurver")) == "0.3.2")
design <- list(n = 300L, D = 60L, true_M = 3:5, snr = c(1, 4, 16),
  repetitions = 50L, max_M = 8L, K = 50L, tolerance = 1e-6, convergence = "normalized",
  initial_budget = 1500L, continuation_budget = 1500L, maximum_sweeps = 10000L,
  effective_weight_tol = 1e-12, method = "isomap", similarity_metric = "spearman",
  cluster_linkage = "single", rw_q = 2L, ridge = 0,
  T_start = 5, T_end = 1, n_outer = 25L, inner_iter = 1L,
  position_prior = "adaptive", generator = "unit_variance_monotone_v1")
design_hash <- digest::digest(design, algo = "sha256")
for (directory in c("data", "results", "candidates", "checkpoints", "logs", "full_fits", "figures")) {
  dir.create(file.path(study_dir, directory), showWarnings = FALSE)
}
atomic_save <- function(object, path, compress = FALSE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  stopifnot(file.rename(temporary, path))
}
make_manifest <- function() {
  rows <- list()
  index <- 0L
  for (M in design$true_M) for (snr in design$snr) for (replicate in seq_len(design$repetitions)) {
    index <- index + 1L
    rows[[index]] <- data.frame(id = sprintf("main_M%d_S%d_r%03d", M, snr, replicate),
      phase = "main", true_M = M, snr = snr, replicate = replicate, seed = 260925000L + index)
  }
  main <- do.call(rbind, rows)
  pilot <- main[main$replicate == 1L, ]
  pilot$phase <- "pilot"
  pilot$replicate <- 0L
  pilot$id <- sprintf("pilot_M%d_S%d", pilot$true_M, pilot$snr)
  pilot$seed <- 260926000L + seq_len(nrow(pilot))
  rbind(pilot, main)
}
load_dataset <- function(row) {
  path <- file.path(study_dir, "data", paste0(row$id, ".rds"))
  if (file.exists(path)) {
    value <- readRDS(path)
    stopifnot(identical(value$design_hash, design_hash), value$row$seed == row$seed)
    return(value)
  }
  simulation <- MPCurver::simulate_intrinsic_trajectories(n = design$n,
    d_signal = rep(design$D / row$true_M, row$true_M), d_noise = 0L,
    sigma = 0, seed = row$seed, trajectory_family = "monotone")
  signal <- scale(simulation$X, center = TRUE, scale = TRUE)
  signal <- matrix(signal, nrow = design$n, dimnames = dimnames(simulation$X))
  noise <- matrix(stats::rnorm(length(signal), sd = 1 / sqrt(row$snr)), nrow = design$n)
  X <- signal + noise
  rownames(X) <- rownames(signal) <- sprintf("sample_%03d", seq_len(design$n))
  stopifnot(identical(dim(X), c(design$n, design$D)), all(is.finite(X)),
    max(abs(apply(signal, 2, var) - 1)) < 1e-10,
    all(table(simulation$true_assign) == design$D / row$true_M))
  value <- list(row = row, design_hash = design_hash, X = X, signal = signal,
    true_assign = simulation$true_assign, latent_positions = simulation$latent_positions,
    ordering_labels = simulation$ordering_labels, noise_variance = 1 / row$snr,
    realized_noise_variance = apply(noise, 2, var),
    input_hash = digest::digest(X, algo = "sha256"))
  atomic_save(value, path)
  value
}
fit_options <- function() list(algorithm = "cavi", method = design$method, K = design$K,
  rw_q = design$rw_q, ridge = design$ridge, lambda = 1, fix_lambda = FALSE, S = NULL,
  position_prior = design$position_prior, similarity_metric = design$similarity_metric,
  cluster_linkage = design$cluster_linkage, discretization = "quantile", num_cores = 1L,
  iter = design$initial_budget, tol = design$tolerance, convergence = design$convergence,
  T_start = design$T_start, T_end = design$T_end, n_outer = design$n_outer,
  inner_iter = design$inner_iter, max_converge_iter = design$initial_budget,
  tol_outer = design$tolerance, verbose = FALSE)

# Retain the original seed registry and select the approved ten repeats.
# This version has a new design hash and separate saved artifacts.
planned_repetitions <- 10L
active_manifest <- function() {
  registry <- read.csv(file.path(study_dir, "manifest.csv"))
  registry[registry$phase == "pilot" | registry$replicate <= planned_repetitions, , drop = FALSE]
}
