#!/usr/bin/env Rscript
source("experiments/estimate_intrinsic_m_v032/common.R")
source(file.path(study_dir, "render_design.R"))
rmarkdown::render(file.path(study_dir, "analysis_report.Rmd"),
  knit_root_dir = getwd(), quiet = TRUE)
