source("experiments/isomap_kmin_m1_p2_v040/common.R")
metadata <- readRDS(file.path(study, "provenance.rds"))
stopifnot(identical(metadata$design_hash, design_hash),
  identical(metadata$package_source_hashes, package_source_hashes),
  identical(metadata$script_hashes, provenance()$script_hashes))
set.seed(design$fit_order_seed)
fit_order <- do.call(rbind, lapply(seq_len(design$replications), function(replication)
  data.frame(replication = replication, method = sample(design$methods))))
write.csv(fit_order, file.path(study, "fit_order.csv"), row.names = FALSE)
for (index in seq_len(nrow(fit_order))) {
  replication <- fit_order$replication[index]
  method <- fit_order$method[index]
  input <- readRDS(input_path(replication))
  stopifnot(identical(input$input_sha256, hash_object(input$X)))
  path <- result_path(replication, method)
  if (file.exists(path)) {
    result <- readRDS(path)
    stopifnot(identical(result$input_sha256, input$input_sha256),
      identical(result$design_hash, design_hash),
      identical(result$package_source_hashes, package_source_hashes))
  } else {
    result <- fit_replicate(input, method)
    save_atomic(result, path)
  }
  cat(sprintf("%02d/%02d: rep=%02d, %-10s k=%2d, %-14s iter=%4d, cosine=%.6f, ELBO=%.6f\n",
    index, nrow(fit_order), replication, method, result$k_used, result$status,
    result$iterations, result$metrics[["cosine"]],
    if (length(result$elbo_trace)) tail(result$elbo_trace, 1) else NA_real_))
  flush.console()
}
