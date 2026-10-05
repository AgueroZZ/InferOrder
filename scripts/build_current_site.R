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

  # Refresh navigation on retained exploratory pages without rerunning fits.
  retained <- c("ordering_playground.html", "explore_initialization.html")
  navbar_pattern <- '(?s)<div class="navbar navbar-default.*?</div><!--/\\.navbar -->'
  index_html <- paste(readLines("docs/index.html", warn = FALSE), collapse = "\n")
  navbar <- regmatches(index_html, regexpr(navbar_pattern, index_html, perl = TRUE))
  stopifnot(length(navbar) == 1L, nzchar(navbar))
  for (page in retained) {
    path <- file.path("docs", page)
    html <- paste(readLines(path, warn = FALSE), collapse = "\n")
    match <- regexpr(navbar_pattern, html, perl = TRUE)
    stopifnot(match[1] > 0)
    regmatches(html, match) <- navbar
    writeLines(html, path, useBytes = TRUE)
  }
  expected <- c(sub("\\.Rmd$", ".html", pages), retained)
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
