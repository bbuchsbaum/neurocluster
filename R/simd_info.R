#' SIMD / CPU Feature Report (Development Helper)
#'
#' Reports build-time SIMD feature macros and (when available) a small set of
#' runtime CPU feature flags, plus the currently-selected int8 dot-product
#' kernel used by neurocluster's correlation-SLIC code paths.
#'
#' This is primarily a development/benchmarking aid. For definitive evidence
#' that an ISA-specific kernel is in use, disassemble the package shared
#' library and look for the corresponding instructions.
#'
#' @return A named list with fields `compiler`, `arch`, `compile`, `runtime`,
#'   `dot_i8_kernel`, `shared_library`, and `notes`.
#' @export
neurocluster_simd_info <- function() {
  simd_info_cpp <- getFromNamespace("simd_info_cpp", "neurocluster")
  info <- simd_info_cpp()
  info$shared_library <- .neurocluster_shared_library_path()
  info
}

.neurocluster_shared_library_path <- function() {
  dlls <- getLoadedDLLs()
  if ("neurocluster" %in% names(dlls)) {
    # DLLInfo has a custom `$` method; use `[[` to access fields.
    return(normalizePath(dlls[["neurocluster"]][["path"]], winslash = "/", mustWork = FALSE))
  }

  # Installed package layout
  libdir <- system.file("libs", package = "neurocluster")
  if (nzchar(libdir)) {
    p <- file.path(libdir, paste0("neurocluster", .Platform$dynlib.ext))
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }

  # Source checkout layout (pkgload / devtools)
  p <- file.path(getwd(), "src", paste0("neurocluster", .Platform$dynlib.ext))
  if (file.exists(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }

  NA_character_
}
