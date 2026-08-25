#!/usr/bin/env Rscript

# Exact masked-grid construction and Dijkstra propagation versus the pinned
# pre-remediation Euclidean kNN graph at approximately matched directed edges.

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE, recompile = FALSE))

reps <- as.integer(Sys.getenv("RENA_G3S_BENCH_REPS", "5"))
if (is.na(reps) || reps < 1L) stop("RENA_G3S_BENCH_REPS must be positive")

median_time <- function(expr) {
  timings <- numeric(reps)
  result <- NULL
  for (i in seq_len(reps)) {
    timing <- system.time(result <- force(expr()))
    timings[i] <- unname(timing[["elapsed"]])
  }
  list(seconds = stats::median(timings), result = result)
}

dims <- c(24L, 24L, 4L)
mask_values <- array(1, dims)
mask_values[12L, 7:18, 2:3] <- 0
space <- neuroim2::NeuroSpace(dims, spacing = c(1.2, 1.2, 2.5))
mask <- neuroim2::NeuroVol(mask_values, space)
mask_idx <- which(mask_values > 0)
coords <- neurocluster:::.cluster4d_index_to_coord(mask, mask_idx)
grid_coords <- arrayInd(mask_idx, dims)
n <- length(mask_idx)

grid_build <- median_time(function() {
  neurocluster:::build_grid_neighbors_g3s(mask, mask_idx, 26L)
})
grid <- grid_build$result
directed_grid_edges <- sum(grid$nn.index > 0L)
matched_k <- max(1L, min(n - 1L, as.integer(round(directed_grid_edges / n))))
knn_build <- median_time(function() FNN::get.knn(coords, k = matched_k))
knn <- knn_build$result

set.seed(825)
feature_rows <- cbind(
  sin(grid_coords[, 1L] * 0.19),
  cos(grid_coords[, 2L] * 0.17),
  sin(grid_coords[, 3L] * 0.53),
  matrix(rnorm(n * 5L, sd = 0.25), n, 5L)
)
feature_rows <- feature_rows / sqrt(rowSums(feature_rows^2))
feature_mat <- t(feature_rows)
seed_indices <- unique(as.integer(round(seq(1, n, length.out = 32L))))
spatial_scale <- neurocluster:::.g3s_spatial_scale(
  mask, mask_idx, length(seed_indices)
)$scale

grid_propagation <- median_time(function() {
  g3s_propagate_cpp(
    feature_mat, seed_indices, grid$nn.index, grid$nn.dist,
    alpha = 0.5, compactness = spatial_scale
  )
})
knn_propagation <- median_time(function() {
  g3s_propagate_cpp(
    feature_mat, seed_indices, knn$nn.index, knn$nn.dist,
    alpha = 0.5, compactness = spatial_scale
  )
})

knn_left <- rep(seq_len(n), each = matched_k)
knn_right <- as.integer(t(knn$nn.index))
delta <- abs(grid_coords[knn_left, , drop = FALSE] -
               grid_coords[knn_right, , drop = FALSE])
knn_non_grid_edges <- sum(apply(delta, 1L, max) > 1L)

connected_piece_count <- function(labels) {
  normalized <- neurocluster:::.exact_k_connected_labels(
    labels, grid$graph, grid$edges
  )
  length(unique(normalized))
}

grid_labels <- as.integer(grid_propagation$result)
knn_labels <- as.integer(knn_propagation$result)
metrics <- data.frame(
  metric = c(
    "masked_voxels",
    "grid_directed_edges",
    "knn_directed_edges",
    "matched_knn_k",
    "grid_build_seconds_median",
    "knn_build_seconds_median",
    "grid_propagation_seconds_median",
    "knn_propagation_seconds_median",
    "knn_non_grid_edges",
    "grid_connected_pieces",
    "knn_connected_pieces",
    "requested_seed_labels"
  ),
  value = c(
    n,
    directed_grid_edges,
    length(knn$nn.index),
    matched_k,
    grid_build$seconds,
    knn_build$seconds,
    grid_propagation$seconds,
    knn_propagation$seconds,
    knn_non_grid_edges,
    connected_piece_count(grid_labels),
    connected_piece_count(knn_labels),
    length(seed_indices)
  )
)
print(metrics, row.names = FALSE)

if (any(grid_labels <= 0L) || connected_piece_count(grid_labels) != length(seed_indices)) {
  stop("exact-grid propagation connectivity contract failed")
}
if (knn_non_grid_edges < 1L) {
  stop("pinned kNN comparison did not reproduce a non-grid edge")
}
