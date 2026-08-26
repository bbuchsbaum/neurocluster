library(testthat)
library(neuroim2)

acsc_contract_neuro <- function(features, dims) {
  features <- as.matrix(features)
  stopifnot(nrow(features) == prod(dims))
  list(
    vec = NeuroVec(
      array(features, dim = c(dims, ncol(features))),
      NeuroSpace(c(dims, ncol(features)))
    ),
    mask = NeuroVol(array(1, dim = dims), NeuroSpace(dims))
  )
}

direct_biweight_midcor <- function(x, y) {
  transform_one <- function(z) {
    residual <- z - stats::median(z)
    scale <- stats::median(abs(residual))
    if (scale <= sqrt(.Machine$double.eps)) return(rep(0, length(z)))
    u <- residual / (9 * scale)
    weight <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
    residual * weight
  }
  xw <- transform_one(x)
  yw <- transform_one(y)
  denominator <- sqrt(sum(xw^2) * sum(yw^2))
  if (denominator <= sqrt(.Machine$double.eps)) return(0)
  sum(xw * yw) / denominator
}

sorted_acsc_edges <- function(graph) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  if (!nrow(edges)) return(edges)
  edge_id <- t(apply(cbind(as.integer(edges$from), as.integer(edges$to)), 1L, sort))
  edges$from <- edge_id[, 1L]
  edges$to <- edge_id[, 2L]
  edges[order(edges$from, edges$to), , drop = FALSE]
}

test_that("robust ACSC preprocessing preserves N by T shape", {
  set.seed(2101)
  fixture <- acsc_contract_neuro(matrix(rnorm(18L * 10L), 18L, 10L), c(18L, 1L, 1L))
  observed <- preprocess_time_series(fixture$vec, fixture$mask, "robust")
  expect_identical(dim(observed), c(18L, 10L))
  expect_true(all(is.finite(observed)))
})

test_that("Pearson, Spearman, and robust similarities match independent oracles", {
  x <- rbind(
    c(1, 2, 3, 4, 5, 6, 7),
    c(2, 1, 4, 3, 7, 5, 6),
    c(1, 4, 9, 16, 25, 36, 49)
  )
  pearson <- acsc_similarity_matrix(x, "pearson")
  spearman <- acsc_similarity_matrix(x, "spearman")
  robust <- acsc_similarity_matrix(x, "robust")

  expect_equal(pearson, stats::cor(t(x), method = "pearson"), tolerance = 1e-12)
  expect_equal(spearman, stats::cor(t(x), method = "spearman"), tolerance = 1e-12)
  robust_oracle <- outer(seq_len(nrow(x)), seq_len(nrow(x)), Vectorize(function(i, j) {
    direct_biweight_midcor(x[i, ], x[j, ])
  }))
  expect_equal(robust, robust_oracle, tolerance = 1e-12)
  expect_equal(spearman[1, 3], 1, tolerance = 1e-12)
  expect_lt(pearson[1, 3], 0.98)
})

test_that("the three ACSC metric paths differ under ranks and outliers", {
  clean <- seq_len(15L)
  monotone <- clean^3
  contaminated <- clean
  contaminated[15L] <- 5000
  data <- rbind(clean, monotone, contaminated)

  pearson <- acsc_similarity_matrix(data, "pearson")
  spearman <- acsc_similarity_matrix(data, "spearman")
  robust <- acsc_similarity_matrix(data, "robust")
  expect_equal(spearman[1, 2], 1, tolerance = 1e-12)
  expect_lt(pearson[1, 2], 0.95)
  expect_gt(robust[1, 3], pearson[1, 3] + 0.4)
  expect_gt(robust[1, 3], 0.9)
  expect_false(isTRUE(all.equal(pearson, spearman)))
  expect_false(isTRUE(all.equal(pearson, robust)))
})

test_that("alpha zero selects only spatial adjacency before feature weighting", {
  grid <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1, 0))
  spatial <- grid + 1
  features_a <- rbind(
    c(1, 2, 3, 4, 5), c(5, 4, 3, 2, 1),
    c(1, 3, 5, 2, 4), c(2, 5, 1, 4, 3)
  )
  features_b <- features_a[c(4, 2, 1, 3), rev(seq_len(ncol(features_a)))] * 17
  summary_a <- list(summaries = features_a, spatial = spatial, grid = grid)
  summary_b <- list(summaries = features_b, spatial = spatial, grid = grid)

  graph_a <- build_acsc_graph(summary_a, 2L, 0, "binary", 1, "pearson")
  graph_b <- build_acsc_graph(summary_b, 2L, 0, "binary", 1, "robust")
  edges_a <- sorted_acsc_edges(graph_a)
  edges_b <- sorted_acsc_edges(graph_b)
  expect_identical(edges_a[, c("from", "to", "weight", "edge_source")],
                   edges_b[, c("from", "to", "weight", "edge_source")])
  expect_identical(nrow(edges_a), 4L)
  expect_true(all(edges_a$weight == 1))
  expect_true(all(edges_a$edge_source == "spatial"))
})

test_that("ACSC graph weights expose the declared correlation, not Euclidean proxy", {
  summaries <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(1, 4, 9, 16, 25, 36),
    c(6, 2, 5, 1, 4, 3)
  )
  block_summary <- list(
    summaries = summaries,
    spatial = cbind(seq_len(3L), 1, 1),
    grid = cbind(seq_len(3L) - 1L, 0, 0)
  )

  for (metric in c("pearson", "spearman", "robust")) {
    graph <- build_acsc_graph(
      block_summary, ann_k = 2L, alpha = 1,
      spatial_weighting = "binary", block_size = 1,
      correlation_metric = metric
    )
    edges <- sorted_acsc_edges(graph)
    oracle <- switch(
      metric,
      pearson = stats::cor(t(summaries), method = "pearson"),
      spearman = stats::cor(t(summaries), method = "spearman"),
      robust = outer(seq_len(3L), seq_len(3L), Vectorize(function(i, j) {
        direct_biweight_midcor(summaries[i, ], summaries[j, ])
      }))
    )
    expected <- mapply(function(i, j) oracle[i, j], edges$from, edges$to)
    expect_equal(edges$correlation, expected, tolerance = 1e-12)
    expect_equal(edges$weight, (expected + 1) / 2, tolerance = 1e-12)
    expect_identical(igraph::graph_attr(graph, "correlation_metric"), metric)
  }
})

test_that("one voxel and one block return standardized exact-K results", {
  set.seed(2102)
  one <- acsc_contract_neuro(matrix(rnorm(10L), 1L, 10L), c(1L, 1L, 1L))
  direct <- expect_silent(acsc(one$vec, one$mask, K = 1L))
  unified <- expect_silent(cluster4d(one$vec, one$mask, n_clusters = 1L, method = "acsc"))
  for (result in list(direct, unified)) {
    expect_s3_class(result, "cluster4d_result")
    expect_identical(result$labels, 1L)
    expect_identical(dim(result$centers), c(1L, 10L))
    expect_identical(dim(result$coord_centers), c(1L, 3L))
    expect_identical(result$actual_k, 1L)
  }

  one_block <- acsc_contract_neuro(matrix(rnorm(8L * 12L), 8L, 12L), c(8L, 1L, 1L))
  split <- acsc(
    one_block$vec, one_block$mask, block_size = 20,
    K = 3L, refine = FALSE
  )
  expect_identical(split$actual_k, 3L)
  expect_identical(sort(unique(split$labels)), 1:3)
})

test_that("ACSC flood-fill normalization and exact-K preserve contiguity", {
  set.seed(2103)
  dims <- c(6L, 4L, 1L)
  time <- seq(0, 2 * pi, length.out = 16L)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  region <- 1L + (coords[, 1L] > 2L) + (coords[, 1L] > 4L)
  prototypes <- rbind(sin(time), cos(time), sin(2 * time))
  features <- prototypes[region, , drop = FALSE] +
    matrix(rnorm(prod(dims) * length(time), sd = 0.05), prod(dims), length(time))
  fixture <- acsc_contract_neuro(features, dims)
  result <- acsc(
    fixture$vec, fixture$mask, K = 3L, block_size = 2L,
    ann_k = 3L, correlation_metric = "pearson", refine = TRUE,
    max_refine_iter = 2L
  )

  graph_info <- .exact_k_graph(fixture$mask, 6L)
  flood_filled <- .exact_k_connected_labels(
    result$labels, graph_info$graph, graph_info$edges
  )
  expect_identical(result$actual_k, 3L)
  expect_identical(flood_filled, result$labels)
  expect_identical(sort(unique(result$labels)), 1:3)
})

test_that("ACSC fails closed when exact K violates mask components", {
  set.seed(2104)
  dims <- c(5L, 1L, 1L)
  features <- matrix(rnorm(prod(dims) * 8L), prod(dims), 8L)
  fixture <- acsc_contract_neuro(features, dims)
  mask_values <- array(c(1, 1, 0, 1, 1), dims)
  mask <- NeuroVol(mask_values, NeuroSpace(dims))
  expect_error(
    acsc(fixture$vec, mask, K = 1L, refine = FALSE),
    class = "cluster4d_exact_k_infeasible"
  )
})

test_that("ACSC centroids retain the sorted label mapping", {
  labels <- c(2L, 2L, 1L, 1L)
  features <- rbind(c(10, 0), c(10, 0), c(0, 10), c(0, 10))
  centroids <- compute_cluster_centroids(labels, features)

  expect_identical(names(centroids), c("1", "2"))
  expect_equal(centroids[["1"]], c(0, 1), tolerance = 1e-12)
  expect_equal(centroids[["2"]], c(1, 0), tolerance = 1e-12)
})

test_that("ACSC Louvain is deterministic and preserves caller RNG state", {
  graph <- igraph::make_ring(24L)
  graph <- igraph::set_edge_attr(graph, "weight", value = rep(1, 24L))

  set.seed(2105)
  before <- .Random.seed
  first <- run_louvain_clustering(graph)
  expect_true(identical(.Random.seed, before))

  set.seed(9999)
  second_before <- .Random.seed
  second <- run_louvain_clustering(graph)
  expect_true(identical(.Random.seed, second_before))
  expect_identical(igraph::membership(first), igraph::membership(second))
})

test_that("full ACSC results do not depend on incoming RNG state", {
  set.seed(2106)
  dims <- c(6L, 4L, 1L)
  features <- matrix(rnorm(prod(dims) * 16L), prod(dims), 16L)
  fixture <- acsc_contract_neuro(features, dims)

  set.seed(1)
  before_first <- .Random.seed
  first <- suppressWarnings(acsc(
    fixture$vec, fixture$mask, K = 4L, block_size = 2L,
    ann_k = 3L, refine = FALSE
  ))
  expect_true(identical(.Random.seed, before_first))

  set.seed(999)
  before_second <- .Random.seed
  second <- suppressWarnings(acsc(
    fixture$vec, fixture$mask, K = 4L, block_size = 2L,
    ann_k = 3L, refine = FALSE
  ))
  expect_true(identical(.Random.seed, before_second))
  expect_identical(first$labels, second$labels)
})
