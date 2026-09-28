#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (length(value)) sub(paste0("^--", name, "="), "", value[1]) else default
}
phase <- option("phase", "pilot")
stopifnot(phase %in% c("pilot", "main"))
manifest <- active_manifest()
rows <- manifest[manifest$phase == phase, , drop = FALSE]
first <- as.integer(option("first", "1"))
last <- as.integer(option("last", as.character(nrow(rows))))
stopifnot(first >= 1L, last <= nrow(rows), last >= first)
rows <- rows[seq.int(first, last), , drop = FALSE]
worker <- as.integer(option("worker", "1"))
workers <- as.integer(option("workers", "1"))
stopifnot(workers %in% 1:2, worker >= 1L, worker <= workers)
rows <- rows[(seq_len(nrow(rows)) - 1L) %% workers == worker - 1L, , drop = FALSE]
if (phase == "main") {
  gate <- file.path(study_dir, "execution_authorization.json")
  stopifnot(file.exists(gate), isTRUE(jsonlite::read_json(gate)$approved_for_main))
}
for (id in rows$id) {
  completed <- vapply(c("adaptive", "forward"), function(method)
    file.exists(file.path(study_dir, "results", paste0(id, "_", method, ".rds"))), logical(1))
  if (all(completed)) {
    cat("Reusing", id, "\n")
    next
  }
  cat(format(Sys.time()), "Dispatching", id, "\n")
  flush.console()
  log <- file.path(study_dir, "logs", paste0(id, ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("experiments/estimate_intrinsic_m_v032/run_dataset.R", id), stdout = log, stderr = log)
  if (status != 0L) stop("Dataset process failed: ", id, "; inspect ", log)
  cat(format(Sys.time()), "Completed", id, "\n")
  flush.console()
}
