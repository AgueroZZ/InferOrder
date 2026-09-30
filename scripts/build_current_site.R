#!/usr/bin/env Rscript

local({
  pages <- c(
    "index.Rmd",
    "method.Rmd",
    "simulation_m1.Rmd",
    "simulation_m1_comparison.Rmd",
    "simulation_summary.Rmd",
    "simulation_m2.Rmd",
    "estimate_intrinsic_m.Rmd",
    "estimate_intrinsic_m_smooth.Rmd",
    "fitness.Rmd",
    "pancreas.Rmd"
  )

  stopifnot(file.exists("_workflowr.yml"), dir.exists("analysis"), dir.exists("docs"))
  stopifnot(requireNamespace("MPCurver", quietly = TRUE))
  stopifnot(requireNamespace("workflowr", quietly = TRUE))
  stopifnot(requireNamespace("rmarkdown", quietly = TRUE))

  for (page in pages) {
    input <- file.path("analysis", page)
    stopifnot(file.exists(input))
    message("Rendering ", input)
    rmarkdown::render(
      input,
      output_dir = "docs",
      knit_root_dir = getwd(),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
    output <- file.path("docs", sub("\\.Rmd$", ".html", page))
    lines <- readLines(output, warn = FALSE)
    writeLines(sub("[[:blank:]]+$", "", lines), output, useBytes = TRUE)
  }

  # Keep the published exploratory playground without rerunning its model fits.
  # It is maintained separately from these saved-result pages.
  expected <- c(sub("\\.Rmd$", ".html", pages), "ordering_playground.html")
  actual <- list.files("docs", pattern = "\\.html$", full.names = FALSE)
  if (!setequal(actual, expected)) {
    stop(
      "The public docs directory must contain only the current HTML pages. ",
      "Expected: ", paste(expected, collapse = ", "), "; found: ",
      paste(actual, collapse = ", ")
    )
  }
  message("Rendered and checked ", length(pages), " current pages.")
})
