#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
status <- jsonlite::read_json(file.path(study_dir, "main_summary", "status.json"))
stopifnot(isTRUE(status$complete), status$completed_method_runs == 180L)
source(file.path(study_dir, "compare_baseline.R"))
source(file.path(study_dir, "render_design.R"))
source(file.path(study_dir, "render_examples.R"))
rmarkdown::render(file.path(study_dir, "analysis_report.Rmd"),
  knit_root_dir = getwd(), quiet = TRUE)
