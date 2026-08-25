#!/usr/bin/env Rscript

# Fixed-fixture comparison of the pinned pre-remediation Rcpp::List queue and
# seed-after-pop implementation (`snic_main`) with the ownership-first,
# lightweight candidate (`snic_main_optimized`). The legacy implementation is
# retained only as the executable before-state for this benchmark.

suppressPackageStartupMessages(
  devtools::load_all(".", quiet = TRUE, recompile = FALSE)
)

make_fixture <- function() {
  dims <- c(24L, 24L, 4L)
  coords_1 <- arrayInd(seq_len(prod(dims)), dims)
  coords_0 <- coords_1 - 1L
  n <- nrow(coords_0)
  features <- vapply(seq_len(n), function(index) {
    xyz <- coords_1[index, ]
    sin(seq_len(12L) * 0.17 + xyz[[1]] * 0.09) +
      cos(seq_len(12L) * 0.11 + xyz[[2]] * 0.07) +
      xyz[[3]] * seq_len(12L) * 0.002
  }, numeric(12L))
  features <- neurocluster:::.snic_normalize_features(features)
  seed_grid <- expand.grid(
    x = seq(2L, dims[[1]] - 1L, length.out = 6L),
    y = seq(2L, dims[[2]] - 1L, length.out = 6L),
    z = seq_len(dims[[3]])
  )
  seed_grid <- unique(round(as.matrix(seed_grid)))
  seeds <- match(
    paste(seed_grid[, 1], seed_grid[, 2], seed_grid[, 3]),
    paste(coords_1[, 1], coords_1[, 2], coords_1[, 3])
  ) - 1L
  seeds <- as.integer(seeds[!is.na(seeds)])
  mask <- array(1L, dims)
  lookup <- array(seq_len(n) - 1L, dims)
  norm_coords <- matrix(as.numeric(coords_0), ncol = 3L)
  extents <- apply(norm_coords, 2L, function(x) max(diff(range(x)), 1e-6))
  spacing <- (prod(extents) / length(seeds))^(1 / 3) * 1.1
  list(
    dims = dims,
    n = n,
    K = length(seeds),
    seeds = seeds,
    features = features,
    mask = mask,
    lookup = lookup,
    coords = matrix(as.integer(coords_0), ncol = 3L),
    norm_coords = norm_coords,
    centroids = norm_coords[seeds + 1L, , drop = FALSE],
    s = spacing^2
  )
}

fixture <- make_fixture()

run_method <- function(method) {
  labels <- array(0L, fixture$dims)
  if (method == "baseline") {
    # compactness=1 makes the legacy priority monotone-equivalent to a 0.5/0.5
    # feature-spatial mixture before centroid accounting diverges.
    snic_main(
      labels, fixture$mask, fixture$centroids, fixture$seeds,
      fixture$coords, fixture$norm_coords, fixture$features,
      fixture$K, fixture$s, 1, fixture$lookup
    )
  } else {
    snic_main_optimized(
      labels, fixture$mask, fixture$centroids, fixture$seeds,
      fixture$coords, fixture$norm_coords, fixture$features,
      fixture$K, fixture$s, 0.5, fixture$lookup
    )
  }
}

visible_allocated_bytes <- function(method) {
  profile <- tempfile("snic-rprofmem-")
  on.exit(unlink(profile), add = TRUE)
  Rprofmem(profile)
  on.exit(Rprofmem(NULL), add = TRUE)
  invisible(run_method(method))
  Rprofmem(NULL)
  lines <- readLines(profile, warn = FALSE)
  sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  sum(sizes, na.rm = TRUE)
}

benchmark_method <- function(method, repetitions = 7L) {
  invisible(run_method(method))
  gc()
  elapsed <- numeric(repetitions)
  result <- NULL
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- system.time({
      result <- run_method(method)
    })[["elapsed"]]
  }
  counts <- attr(result, "snic_centroid_counts")
  labels <- as.integer(result)
  data.frame(
    method = method,
    n_voxels = fixture$n,
    n_features = nrow(fixture$features),
    requested_k = fixture$K,
    actual_k = length(unique(labels)),
    seeds_owned = identical(labels[fixture$seeds + 1L], seq_len(fixture$K)),
    centroid_counts_match = identical(
      as.integer(counts), tabulate(labels, nbins = fixture$K)
    ),
    queue_pushes = as.numeric(attr(result, "snic_queue_pushes")),
    queue_pops = as.numeric(attr(result, "snic_queue_pops")),
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    r_visible_alloc_bytes = visible_allocated_bytes(method),
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  benchmark_method("baseline"),
  benchmark_method("candidate")
)
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
