#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
suppressWarnings(suppressMessages({
  if (!requireNamespace("pkgdown", quietly = TRUE)) {
    install.packages("pkgdown", repos = "https://cloud.r-project.org")
  }
  library(pkgdown)
}))

# Keep things in-process and avoid example execution to sidestep native crashes
try(pkgdown::build_reference(examples = FALSE, run_dont_run = FALSE, devel = FALSE), silent = TRUE)
try(pkgdown::build_articles(), silent = TRUE)
pkgdown::build_site(preview = FALSE, new_process = FALSE)

cat("\nDone. Site in ./docs.\n")

