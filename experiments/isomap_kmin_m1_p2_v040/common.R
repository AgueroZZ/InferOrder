# Shared design and helpers for the paired M=1, P=2 internal experiment.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(MPCurver))
study <- "experiments/isomap_kmin_m1_p2_v040"
package_source <- normalizePath("../MPCurver")
archive_path <- file.path(study, "source", "MPCurver_0.4.0.9000.tar.gz")
stopifnot(as.character(packageVersion("MPCurver")) == "0.4.0.9000",
  identical(normalizePath(find.package("MPCurver")),
            normalizePath(file.path(study, "library", "MPCurver"))),
  is.null(formals(isomap_ordering)$num_neighbors))

design <- list(replications = 20L, N = 200L, P = 2L, M = 1L,
  spline_df = 8L, degree = 3L, dense_grid_size = 2001L,
  monotonicity_increment_tolerance = 1e-8, noise_sd = 0.25,
  generating_seed_base = 202610020L, fit_seed_base = 202620020L,
  fit_order_seed = 202630020L, num_bins = 50L,
  tolerance = 1e-6, block_budget = 2000L, total_budget = 10000L,
  methods = c("auto_kmin", "fixed_k15", "fixed_k10"),
  control = mpcurve_control(rw_order = 2L, ridge = 0, lambda_init = 1,
                            convergence = "normalized"),
  discretization = "quantile", on_failure = "pca", keep = "giant",
  data_orientation = "samples by features", observed_standardization = FALSE)

hash_object <- function(value) digest::digest(value, algo = "sha256")
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
design_hash <- hash_object(design)
source_paths <- c(file.path(package_source, "DESCRIPTION"),
  file.path(package_source, "NAMESPACE"),
  list.files(file.path(package_source, "R"), pattern = "\\.R$", full.names = TRUE))
package_source_hashes <- setNames(vapply(source_paths, hash_file, character(1)),
  substring(source_paths, nchar(package_source) + 2L))
script_paths <- file.path(study,
  c("common.R", "prepare.R", "run.R", "report.R", "verify.R", "run_r.sh", "DESIGN.md"))

provenance <- function() list(design = design, design_hash = design_hash,
  package_version = as.character(packageVersion("MPCurver")),
  package_base_commit = system2("git", c("-C", package_source, "rev-parse", "HEAD"), stdout = TRUE),
  inferorder_base_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  package_source_hashes = package_source_hashes,
  archive_sha256 = hash_file(archive_path),
  script_hashes = setNames(vapply(script_paths, hash_file, character(1)), basename(script_paths)),
  session = sessionInfo())

input_path <- function(replication) file.path(study, "inputs", sprintf("rep%02d.rds", replication))
result_path <- function(replication, method, full = FALSE) file.path(study,
  if (full) "full_fits" else "results", sprintf("rep%02d_%s.rds", replication, method))
save_atomic <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop("Could not save ", path)
}

is_nonmonotone <- function(curve) {
  increments <- diff(curve)
  any(increments > design$monotonicity_increment_tolerance) &&
    any(increments < -design$monotonicity_increment_tolerance)
}

simulate_replicate <- function(replication) {
  seed <- design$generating_seed_base + as.integer(replication)
  set.seed(seed)
  truth <- runif(design$N)
  grid <- seq(0, 1, length.out = design$dense_grid_size)
  basis_grid <- splines::bs(grid, df = design$spline_df, degree = design$degree,
    intercept = TRUE, Boundary.knots = c(0, 1))
  basis <- predict(basis_grid, truth)
  attempts <- 0L
  repeat {
    attempts <- attempts + 1L
    if (attempts > 10000L) stop("Nonmonotone generator exhausted its budget.")
    coefficients <- matrix(rnorm(design$spline_df * design$P), design$spline_df)
    dense <- basis_grid %*% coefficients
    centers <- colMeans(dense)
    centered <- sweep(dense, 2L, centers)
    signal_scale <- sqrt(mean(centered^2))
    dense_signal <- centered / signal_scale
    if (all(apply(dense_signal, 2L, is_nonmonotone))) break
  }
  signal <- sweep(basis %*% coefficients, 2L, centers) / signal_scale
  standard_noise <- matrix(rnorm(design$N * design$P), design$N)
  X <- signal + design$noise_sd * standard_noise
  dimnames(X) <- list(sprintf("sample_%03d", seq_len(design$N)),
                      sprintf("feature_%d", seq_len(design$P)))
  list(replication = as.integer(replication), seed = seed,
    design_hash = design_hash, truth = truth, coefficients = coefficients,
    grid = grid, centers = centers, signal_scale = signal_scale,
    dense_signal = dense_signal, signal = signal, standard_noise = standard_noise,
    shape_attempts = attempts, X = X, input_sha256 = hash_object(X))
}

cosine <- function(a, b) {
  denominator <- sqrt(sum(a^2) * sum(b^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(a * b) / denominator
}
position_metrics <- function(truth, estimated) {
  if (length(truth) != length(estimated) || any(!is.finite(estimated))) {
    return(c(cosine = NA_real_, centered_cosine = NA_real_, spearman = NA_real_))
  }
  direct <- cosine(truth, estimated)
  reflected <- cosine(truth, 1 - estimated)
  ordinary <- if (all(is.na(c(direct, reflected)))) NA_real_ else max(c(direct, reflected), na.rm = TRUE)
  centered <- abs(cosine(truth - mean(truth), estimated - mean(estimated)))
  spearman <- if (length(unique(estimated)) < 2L) NA_real_ else
    abs(cor(truth, estimated, method = "spearman"))
  metrics <- c(cosine = ordinary, centered_cosine = centered, spearman = spearman)
  pmin(metrics, 1)
}

method_arguments <- function(method, seed) {
  arguments <- list(seed = seed)
  if (method != "auto_kmin") arguments$num_neighbors <-
    switch(method, fixed_k15 = 15L, fixed_k10 = 10L, stop("Unknown method."))
  arguments
}

fit_replicate <- function(input, method) {
  seed <- design$fit_seed_base + input$replication
  arguments <- method_arguments(method, seed)
  warnings <- character()
  raw_warnings <- character()
  raw <- withCallingHandlers(tryCatch(do.call(isomap_ordering,
    c(list(X = input$X), arguments)), error = function(error) NULL),
    warning = function(warning) {
      raw_warnings <<- c(raw_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(tryCatch({
    state <- fit_mpcurve(input$X, num_bins = design$num_bins,
      intrinsic_dim = design$M, initial_method = "isomap", position_prior = "adaptive",
      max_iter = design$block_budget, tol = design$tolerance,
      init_control = mpcurve_init_control(discretization = design$discretization,
        method_args = arguments), control = design$control)
    while (!isTRUE(state$converged) && state$fit$iter < design$total_budget) {
      state <- do_mpcurve(state,
        max_iter = min(design$block_budget, design$total_budget - state$fit$iter),
        tol = design$tolerance)
    }
    state
  }, error = identity), warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  elapsed <- proc.time()[["elapsed"]] - started
  failed <- inherits(fit, "error")
  estimated <- if (failed) rep(NA_real_, design$N) else as.numeric(fitted_positions(fit))
  result <- list(replication = input$replication, method = method, seed = seed,
    design_hash = design_hash, input_sha256 = input$input_sha256,
    package_source_hashes = package_source_hashes,
    status = if (failed) "failed" else if (isTRUE(fit$converged)) "converged" else "budget_limited",
    error = if (failed) conditionMessage(fit) else "",
    warnings = warnings, raw_warnings = raw_warnings,
    k_used = if (failed) if (is.null(raw)) NA_integer_ else raw$k_used else
      fit$fit$init_info$k_used %||% NA_integer_,
    graph_components = if (is.null(raw)) NA_integer_ else as.integer(raw$n_components),
    graph_keep_count = if (is.null(raw)) NA_integer_ else length(raw$keep_idx),
    raw_positions = if (is.null(raw)) rep(NA_real_, design$N) else as.numeric(raw$t),
    positions = estimated, metrics = position_metrics(input$truth, estimated),
    raw_metrics = position_metrics(input$truth,
      if (is.null(raw)) rep(NA_real_, design$N) else as.numeric(raw$t)),
    constant_cosine = cosine(input$truth, rep(0.5, design$N)),
    elbo_trace = if (failed) numeric() else fit$elbo_trace,
    iterations = if (failed) NA_integer_ else fit$fit$iter,
    converged = if (failed) FALSE else isTRUE(fit$converged),
    num_bins = if (failed) NA_integer_ else fit$K,
    initialization = if (failed) NULL else fit$fit$init_info,
    sigma2 = if (failed) NULL else fit$params$sigma2,
    precision = if (failed) NULL else fit$fit$lambda_vec,
    elapsed_seconds = elapsed)
  if (!failed) save_atomic(fit, result_path(input$replication, method, full = TRUE))
  result
}

`%||%` <- function(value, alternative) if (is.null(value)) alternative else value
