test_that("native helper types cannot leak external linkage", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  lint <- file.path(root, "tools", "check-native-helper-types.R")
  testthat::skip_if_not(file.exists(lint), "source-tree lint is run by CI")

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(lint, root),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  expect_equal(status, 0L, info = paste(output, collapse = "\n"))
  expect_true(any(grepl("native helper linkage lint: ok", output, fixed = TRUE)))
})

test_that("tiny SLIC completes in an isolated bounded subprocess", {
  dll <- getLoadedDLLs()[["neurocluster"]][["path"]]
  script <- tempfile("neurocluster-tiny-slic-", fileext = ".R")
  on.exit(unlink(script), add = TRUE)

  writeLines(c(
    "library(Rcpp)",
    "library(RcppParallel)",
    sprintf("dyn.load(%s)", encodeString(normalizePath(dll, mustWork = TRUE), quote = "'")),
    "feats <- matrix(seq_len(16), nrow = 8, ncol = 2)",
    "coords <- cbind(rep(0:1, 4), rep(rep(0:1, each = 2), 2), rep(0:1, each = 4))",
    paste0(
      "ans <- .Call('_neurocluster_slic4d_core', PACKAGE = 'neurocluster', ",
      "feats, coords, 0:7, c(2L, 2L, 2L), c(1, 1, 1), 2L, 10, 1L, 0, 1L, ",
      "'mask_grid', FALSE, 0L, 26L, TRUE, FALSE, 2L, numeric(), 1L, FALSE)"
    ),
    "stopifnot(length(ans$labels) == 8L, all(is.finite(ans$center_feats)))",
    "cat('tiny-slic-ok\\n')"
  ), script)

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    script,
    stdout = TRUE,
    stderr = TRUE,
    timeout = 60
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  expect_equal(status, 0L, info = paste(output, collapse = "\n"))
  expect_true(any(grepl("tiny-slic-ok", output, fixed = TRUE)))
})
