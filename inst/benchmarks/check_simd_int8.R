# Disassemble neurocluster shared library and grep for int8 SIMD dot-product
# instructions.
#
# Usage:
#   Rscript inst/benchmarks/check_simd_int8.R
#
# Notes:
# - This script is best-effort and platform-dependent.
# - "Real SIMD int8 dot kernels" means the corresponding ISA instructions are
#   present in the compiled binary.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  } else {
    library(neurocluster)
  }
})

info <- neurocluster_simd_info()
cat("== neurocluster_simd_info() ==\n")
print(info)

so <- info$shared_library
if (is.na(so) || !nzchar(so) || !file.exists(so)) {
  stop("Could not locate neurocluster shared library (info$shared_library).")
}

cat("\n== Shared library ==\n")
cat(so, "\n", sep = "")

os <- Sys.info()[["sysname"]]
is_windows <- identical(.Platform$OS.type, "windows")

pick_tool <- function() {
  if (!is_windows && identical(os, "Darwin") && nzchar(Sys.which("otool"))) {
    return(list(cmd = "otool", args = c("-tvV", so), label = "otool -tvV"))
  }
  if (nzchar(Sys.which("llvm-objdump"))) {
    return(list(cmd = "llvm-objdump", args = c("-d", so), label = "llvm-objdump -d"))
  }
  if (nzchar(Sys.which("objdump"))) {
    return(list(cmd = "objdump", args = c("-d", so), label = "objdump -d"))
  }
  NULL
}

tool <- pick_tool()
if (is_windows || is.null(tool)) {
  cat("\nNo supported disassembler found on this platform.\n")
  cat("Try one of: otool (macOS), llvm-objdump, objdump.\n")
  quit(status = 0)
}

cat("\n== Disassembly tool ==\n")
cat(tool$label, "\n")

cat("\nDisassembling... (this can take a few seconds)\n")
out <- suppressWarnings(system2(tool$cmd, tool$args, stdout = TRUE, stderr = TRUE))
if (!length(out)) {
  stop("Disassembler returned no output.")
}

strip_comments <- function(lines) {
  # otool uses ';' for comments; objdump often uses '#'.
  lines <- sub(";.*$", "", lines)
  lines <- sub("#.*$", "", lines)
  lines
}

out_clean <- strip_comments(out)

patterns <- c(
  # x86 AVX-512 VNNI int8 dot
  "vpdpbusd",
  # x86 AVX2 int8 dot patterns (common lowering)
  "vpmaddub",   # matches vpmaddubs/vpmaddubsw variants
  "vpmaddwd",
  # ARM NEON dotprod
  "sdot",
  "udot",
  # ARM matrix multiply accumulate (sometimes used in dot-style kernels)
  "smmla",
  "usmmla"
)

counts <- vapply(
  patterns,
  function(p) sum(grepl(p, out_clean, ignore.case = TRUE, fixed = FALSE)),
  integer(1)
)
hits <- counts[counts > 0]

cat("\n== Instruction hits ==\n")
if (!length(hits)) {
  cat("No matches found for:\n")
  cat(paste0("  - ", patterns, collapse = "\n"), "\n")
  cat("\nThis usually means the binary does not contain dedicated int8 SIMD dot instructions.\n")
  cat("If you expect them, ensure you built with ISA-specific kernels and the relevant flags/dispatch.\n")
} else {
  for (nm in names(hits)) {
    cat(sprintf("  %s: %d\n", nm, hits[[nm]]))
  }
  cat("\nIf you see vpdpbusd/vpmaddub*/sdot/udot, the compiler emitted corresponding SIMD instructions.\n")
}
