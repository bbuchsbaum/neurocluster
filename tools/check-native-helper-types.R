#!/usr/bin/env Rscript

# Fail when a C++ helper type is declared at translation-unit scope. A type in
# the global namespace gives its inline members external linkage; repeating a
# common name such as GridIndex, UnionFind, or Edge in another translation unit
# can therefore violate the one-definition rule even though neither type is
# part of the package API.

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)

if (length(args) >= 1L) {
  root <- normalizePath(args[[1L]], mustWork = TRUE)
} else if (length(file_arg) == 1L) {
  script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
} else {
  root <- normalizePath(".", mustWork = TRUE)
}

src_dir <- file.path(root, "src")
files <- sort(list.files(src_dir, pattern = "\\.(cc|cpp|cxx)$", full.names = TRUE))
if (!length(files)) {
  stop("No C++ translation units found under ", src_dir, call. = FALSE)
}

strip_comments <- function(lines) {
  out <- character(length(lines))
  in_block <- FALSE

  for (i in seq_along(lines)) {
    chars <- strsplit(lines[[i]], "", fixed = TRUE)[[1L]]
    kept <- character()
    j <- 1L

    while (j <= length(chars)) {
      next_char <- if (j < length(chars)) chars[[j + 1L]] else ""

      if (in_block) {
        if (chars[[j]] == "*" && next_char == "/") {
          in_block <- FALSE
          j <- j + 2L
        } else {
          j <- j + 1L
        }
      } else if (chars[[j]] == "/" && next_char == "*") {
        in_block <- TRUE
        j <- j + 2L
      } else if (chars[[j]] == "/" && next_char == "/") {
        break
      } else {
        kept <- c(kept, chars[[j]])
        j <- j + 1L
      }
    }

    out[[i]] <- paste0(kept, collapse = "")
  }

  out
}

global_types <- list()

for (file in files) {
  lines <- strip_comments(readLines(file, warn = FALSE))
  depth <- 0L

  for (line_no in seq_along(lines)) {
    line <- lines[[line_no]]
    code <- gsub('"([^"\\\\]|\\\\.)*"', '""', line, perl = TRUE)
    code <- gsub("'([^'\\\\]|\\\\.)*'", "''", code, perl = TRUE)

    match <- regexec(
      "^[[:space:]]*(?:template[[:space:]]*<[^;]+>[[:space:]]*)?(struct|class)[[:space:]]+([A-Za-z_][A-Za-z0-9_]*)",
      code,
      perl = TRUE
    )
    fields <- regmatches(code, match)[[1L]]
    if (length(fields) == 3L && depth == 0L) {
      global_types[[length(global_types) + 1L]] <- data.frame(
        file = sub(paste0("^", root, "/?"), "", file),
        line = line_no,
        kind = fields[[2L]],
        name = fields[[3L]],
        stringsAsFactors = FALSE
      )
    }

    opens <- lengths(regmatches(code, gregexpr("\\{", code, perl = TRUE)))
    closes <- lengths(regmatches(code, gregexpr("\\}", code, perl = TRUE)))
    if (identical(opens, 1L) && !grepl("\\{", code, perl = TRUE)) opens <- 0L
    if (identical(closes, 1L) && !grepl("\\}", code, perl = TRUE)) closes <- 0L
    depth <- depth + opens - closes

    if (depth < 0L) {
      stop("Unbalanced braces while scanning ", file, ":", line_no, call. = FALSE)
    }
  }
}

if (length(global_types)) {
  exposed <- do.call(rbind, global_types)
  duplicate_names <- unique(exposed$name[duplicated(exposed$name) |
                                           duplicated(exposed$name, fromLast = TRUE)])
  lines <- sprintf("%s:%d: externally linked %s %s",
                   exposed$file, exposed$line, exposed$kind, exposed$name)
  if (length(duplicate_names)) {
    lines <- c(lines, paste0("duplicate helper names: ",
                             paste(sort(duplicate_names), collapse = ", ")))
  }
  cat(paste(lines, collapse = "\n"), "\n", file = stderr())
  quit(save = "no", status = 1L)
}

cat(sprintf("native helper linkage lint: ok (%d translation units)\n", length(files)))
