#!/usr/bin/env Rscript

local({
  pages <- c(
    "index.Rmd",
    "method.Rmd",
    "simulation_m1.Rmd",
    "simulation_m2.Rmd",
    "fitness.Rmd",
    "pancreas.Rmd"
  )

  stopifnot(file.exists("_workflowr.yml"), dir.exists("analysis"), dir.exists("docs"))
  stopifnot(utils::packageVersion("MPCurver") == "0.3.0")
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

  expected <- sub("\\.Rmd$", ".html", pages)
  actual <- list.files("docs", pattern = "\\.html$", full.names = FALSE)
  if (!setequal(actual, expected)) {
    stop(
      "The public docs directory must contain only the six current HTML pages. ",
      "Expected: ", paste(expected, collapse = ", "), "; found: ",
      paste(actual, collapse = ", ")
    )
  }
  message("Rendered and checked ", length(pages), " current pages.")
})
