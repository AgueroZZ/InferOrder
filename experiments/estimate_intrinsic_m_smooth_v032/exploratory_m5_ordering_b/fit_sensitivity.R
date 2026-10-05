# Full-data controlled fits for Isomap or PCA initialization/noise comparisons.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir, 'exploratory_m5_ordering_b')
arguments <- commandArgs(trailingOnly = TRUE)
use_pca <- '--pca' %in% arguments
arguments <- arguments[arguments != '--pca']
result_dir <- file.path(out, if (use_pca) 'pca_fits' else 'sensitivity_fits')
dir.create(result_dir, showWarnings = FALSE)
id <- 'main_M5_S4_r001'
d <- load_fixed_dataset(manifest[manifest$id == id, ])
settings <- expand.grid(noise_scale = c(0, 0.5, 1), k = c(5L, 10L, 15L, 20L, 30L))
if (use_pca) settings <- data.frame(noise_scale = c(0, 0.5, 1), k = NA_integer_)
indices <- as.integer(arguments)
stopifnot(length(indices) > 0, all(indices %in% seq_len(nrow(settings))))
seed <- d$row$seed + design$seed_offset
set.seed(seed)
original <- MPCurver:::init_m_trajectories_cavi(X = d$X, M = 5L,
 methods = rep('isomap', 5), K = 50L, rw_q = 2L, ridge = 0,
 discretization = 'quantile', partition_init = 'similarity',
 similarity_metric = 'spline_r2', spline_r2_df = 5L,
 cluster_linkage = 'single', num_iter = 2L)
b <- which(d$true_assign == 'B')
slot <- which(vapply(original$init_info, function(x) setequal(x$feature_idx, b), logical(1)))
stopifnot(length(slot) == 1)
for (i in indices) {
 z <- settings[i, ]; path <- file.path(result_dir, sprintf('setting_%02d.rds', i))
 if (file.exists(path)) next
 X <- d$X
 X[, b] <- d$signal[, b] + z$noise_scale * (d$X[, b] - d$signal[, b])
 set.seed(seed)
 warnings <- character()
 result <- withCallingHandlers({
  iso <- if (use_pca) MPCurver:::PCA_ordering(X[, b], component = 1L) else
   MPCurver:::isomap_ordering(X[, b], k = z$k)
  record <- list(setting = z, input_hash = d$input_hash,
   modified_input_hash = digest::digest(X, algo = 'sha256'), seed = seed,
   package_version = as.character(packageVersion('MPCurver')),
   package_commit = design$package_source_commit,
   initial_method = if (use_pca) 'PCA' else 'isomap',
   initial_coordinates = iso$t,
   isomap = if (use_pca) NULL else iso$t, components = iso$n_components,
   truth = d$latent_positions[, 2], slot = slot)
  if (any(!is.finite(iso$t))) {
   record$status <- 'disconnected_initialization'
   record
  } else {
   sub <- MPCurver:::.cavi_build_from_ordering(X = X[, b], ordering_vec = iso$t,
    K = 50L, rw_q = 2L, ridge = 0, max_iter = 2L,
    tol = 1e-6, discretization = 'quantile', strict_K = TRUE)
   # Preserve each group's subset initialization; expand it to the modified
   # full matrix without iterations, exactly as similarity initialization does.
   fits <- lapply(seq_len(5), function(m) {
    gamma <- if (m == slot) sub$gamma else original$fits[[m]]$gamma
    MPCurver:::cavi(X = X, K = 50L, responsibilities_init = gamma,
     position_prior_init = colMeans(gamma), rw_q = 2L, ridge = 0,
     max_iter = 0L, convergence = 'relative', verbose = FALSE)
   })
   set.seed(seed)
   f <- MPCurver::fit_mpcurve(X = X, intrinsic_dim = 5L, algorithm = 'cavi',
    method = 'isomap', K = 50L, rw_q = 2L, ridge = 0, fits_init = fits,
    partition_prior = 'adaptive', position_prior = 'adaptive',
    lambda = 1, fix_lambda = FALSE, discretization = 'quantile', num_cores = 1L,
    iter = 1500L, tol = 1e-6, convergence = 'normalized',
    T_start = 5, T_end = 1, n_outer = 25L, inner_iter = 1L,
    max_converge_iter = 1500L, tol_outer = 1e-6, verbose = FALSE)
   while (!isTRUE(f$fit$converged) && f$fit$iter < 9999L) {
    f <- MPCurver::do_mpcurve(f, iter = min(1500L, 9999L - f$fit$iter),
     tol = 1e-6, tol_outer = 1e-6, convergence = 'normalized', verbose = FALSE)
   }
   record$status <- if (isTRUE(f$fit$converged)) 'converged' else 'iteration_limit'
   record$final <- as.numeric(f$fit$fits[[slot]]$gamma %*% seq(0, 1, length.out = 50))
   record$objective <- f$fit$objective_history
   record$iterations <- f$fit$iter
   record$ARI <- mclust::adjustedRandIndex(d$true_assign, f$fit$assign)
   record$assignments <- f$fit$assign
   record
  }
 }, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
 })
 result$warnings <- warnings
 atomic_save(result, path, compress = TRUE)
 cat(sprintf('Setting %02d: noise=%s k=%d %s; iterations=%s\n',
  i, z$noise_scale, z$k, result$status,
  if (is.null(result$iterations)) 'NA' else result$iterations))
 flush.console()
}
