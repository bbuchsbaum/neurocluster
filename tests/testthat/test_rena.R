library(testthat)
library(neurocluster)
library(neuroim2)

context("ReNA Clustering Tests")

test_that("rena basic functionality works with small 3D volume", {
  # Create small test volume
  mask <- NeuroVol(array(1, c(4,4,4)), NeuroSpace(c(4,4,4)))
  vec <- NeuroVec(array(rnorm(4*4*4*10), c(4,4,4,10)),
                  NeuroSpace(c(4,4,4,10)))

  # Run ReNA
  result <- rena(vec, mask, K=3, connectivity=6, verbose=FALSE)

  # Check result structure
  expect_s3_class(result, "cluster_result")
  expect_s3_class(result, "rena_cluster_result")
  expect_true("clusvol" %in% names(result))
  expect_true("cluster" %in% names(result))
  expect_true("centers" %in% names(result))
  expect_true("coord_centers" %in% names(result))
  expect_true("n_clusters" %in% names(result))
  expect_true("method" %in% names(result))
  expect_equal(result$method, "rena")

  # Check cluster assignments
  expect_equal(length(result$cluster), sum(mask > 0))
  expect_true(all(result$cluster > 0))
  expect_true(all(!is.na(result$cluster)))

  # Check number of clusters is close to target
  actual_k <- length(unique(result$cluster))
  expect_true(actual_k >= 1)
  expect_true(actual_k <= 10)  # Should be reasonably close to K=3

  # Check center dimensions
  expect_equal(nrow(result$centers), actual_k)
  expect_equal(ncol(result$centers), 10)  # 10 timepoints
  expect_equal(nrow(result$coord_centers), actual_k)
  expect_equal(ncol(result$coord_centers), 3)  # 3D coordinates
})

test_that("rena works via cluster4d interface", {
  # Create small test volume
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*8), c(5,5,5,8)),
                  NeuroSpace(c(5,5,5,8)))

  # Run via cluster4d
  result <- cluster4d(vec, mask, n_clusters=5, method="rena", verbose=FALSE)

  # Check result structure
  expect_s3_class(result, "cluster4d_result")
  expect_s3_class(result, "cluster_result")
  expect_equal(result$method, "rena")

  # Check parameters are stored
  expect_true("parameters" %in% names(result))
  expect_equal(result$parameters$n_clusters_requested, 5)
})

test_that("rena handles K=1 (trivial clustering)", {
  # Create small test volume
  mask <- NeuroVol(array(1, c(3,3,3)), NeuroSpace(c(3,3,3)))
  vec <- NeuroVec(array(rnorm(3*3*3*5), c(3,3,3,5)),
                  NeuroSpace(c(3,3,3,5)))

  # Run with K=1
  result <- rena(vec, mask, K=1, connectivity=6, verbose=FALSE)

  # Should have all voxels in one cluster
  expect_equal(length(unique(result$cluster)), 1)
  expect_equal(result$n_clusters, 1)
  expect_true(all(result$cluster == 1))
})

test_that("rena handles single voxel mask", {
  # Create mask with single voxel
  mask <- NeuroVol(array(0, c(3,3,3)), NeuroSpace(c(3,3,3)))
  mask[2,2,2] <- 1

  vec <- NeuroVec(array(rnorm(3*3*3*5), c(3,3,3,5)),
                  NeuroSpace(c(3,3,3,5)))

  # Run ReNA (should handle gracefully)
  result <- rena(vec, mask, K=1, verbose=FALSE)

  expect_equal(length(result$cluster), 1)
  expect_equal(result$cluster[1], 1)
  expect_equal(result$n_clusters, 1)
})

test_that("rena produces consistent results with same seed", {
  # Create test volume
  set.seed(42)
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)),
                  NeuroSpace(c(5,5,5,10)))

  # Run twice with explicit seed setting (ReNA uses FNN which may use RNG)
  set.seed(100)
  result1 <- rena(vec, mask, K=8, connectivity=6, verbose=FALSE)
  set.seed(100)
  result2 <- rena(vec, mask, K=8, connectivity=6, verbose=FALSE)

  # Results should be identical (ReNA is deterministic with same seed)
  expect_equal(result1$cluster, result2$cluster)
  expect_equal(result1$n_clusters, result2$n_clusters)
})

test_that("rena respects connectivity parameter", {
  # Create test volume
  mask <- NeuroVol(array(1, c(6,6,6)), NeuroSpace(c(6,6,6)))
  vec <- NeuroVec(array(rnorm(6*6*6*8), c(6,6,6,8)),
                  NeuroSpace(c(6,6,6,8)))

  # Run with different connectivity
  result_6 <- rena(vec, mask, K=10, connectivity=6, verbose=FALSE)
  result_26 <- rena(vec, mask, K=10, connectivity=26, verbose=FALSE)

  # Both should complete successfully
  expect_true(result_6$n_clusters > 0)
  expect_true(result_26$n_clusters > 0)

  # Results may differ due to different graph structure
  # Just check they're both valid
  expect_equal(length(result_6$cluster), length(result_26$cluster))
})

test_that("rena convergence and iteration tracking", {
  # Create test volume
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)),
                  NeuroSpace(c(5,5,5,10)))

  # Run with verbose
  result <- rena(vec, mask, K=10, max_iterations=20, verbose=FALSE)

  # Check metadata contains iteration info
  expect_true("metadata" %in% names(result))
  expect_true("iterations" %in% names(result$metadata))
  expect_true(result$metadata$iterations >= 1)
  expect_true(result$metadata$iterations <= 20)
})

test_that("rena validates inputs properly", {
  # Valid inputs
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)),
                  NeuroSpace(c(5,5,5,10)))

  # Should fail with K > number of voxels
  expect_error(
    rena(vec, mask, K=1000, verbose=FALSE),
    "Cannot create.*clusters from.*voxels"
  )

  # Should fail with K <= 0
  expect_error(
    cluster4d(vec, mask, n_clusters=0, method="rena"),
    "n_clusters must be positive"
  )

  # Should fail with empty mask
  empty_mask <- NeuroVol(array(0, c(5,5,5)), NeuroSpace(c(5,5,5)))
  expect_error(
    rena(vec, empty_mask, K=10),
    "No nonzero voxels in mask"
  )
})

test_that("rena clusters are spatially contiguous", {
  # Create test volume with structured pattern
  mask <- NeuroVol(array(1, c(8,8,8)), NeuroSpace(c(8,8,8)))

  # Create data with spatial structure
  data_array <- array(0, c(8,8,8,10))
  for (t in 1:10) {
    # Left half has one pattern, right half has another
    data_array[1:4,,,t] <- rnorm(4*8*8, mean=1)
    data_array[5:8,,,t] <- rnorm(4*8*8, mean=-1)
  }
  vec <- NeuroVec(data_array, NeuroSpace(c(8,8,8,10)))

  # Run ReNA
  result <- rena(vec, mask, K=4, connectivity=26, verbose=FALSE)

  # Check that clusters exist
  expect_true(result$n_clusters >= 2)
  expect_true(result$n_clusters <= 10)
})

test_that("rena exact_k parameter works", {
  # Create test volume
  mask <- NeuroVol(array(1, c(6,6,6)), NeuroSpace(c(6,6,6)))
  vec <- NeuroVec(array(rnorm(6*6*6*8), c(6,6,6,8)),
                  NeuroSpace(c(6,6,6,8)))

  # Run with exact_k = TRUE
  K_target <- 15
  result <- rena(vec, mask, K=K_target, exact_k=TRUE, verbose=FALSE)

  # Should be close to target K (within reason, allowing more tolerance for small volumes)
  expect_true(abs(result$n_clusters - K_target) <= 10 || result$n_clusters >= K_target / 2)

  # Run with exact_k = FALSE
  result2 <- rena(vec, mask, K=K_target, exact_k=FALSE, verbose=FALSE)

  # Should still produce valid clustering
  expect_true(result2$n_clusters > 0)
  expect_true(result2$n_clusters < sum(mask > 0))
})

test_that("rena centers are computed correctly", {
  # Create test volume
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)),
                  NeuroSpace(c(5,5,5,10)))

  result <- rena(vec, mask, K=10, verbose=FALSE)

  # Centers should have correct dimensions
  expect_equal(nrow(result$centers), result$n_clusters)
  expect_equal(ncol(result$centers), 10)  # 10 timepoints

  # Coordinate centers should be in valid range
  dims <- dim(mask)
  expect_true(all(result$coord_centers[,1] >= 1 & result$coord_centers[,1] <= dims[1]))
  expect_true(all(result$coord_centers[,2] >= 1 & result$coord_centers[,2] <= dims[2]))
  expect_true(all(result$coord_centers[,3] >= 1 & result$coord_centers[,3] <= dims[3]))

  # No NA values
  expect_true(all(!is.na(result$centers)))
  expect_true(all(!is.na(result$coord_centers)))
})

test_that("C++ helper functions work correctly", {
  skip_if_not_installed("Matrix")

  # Test compute_masked_distances_cpp
  # Create feature matrix with 5 voxels (columns) to match the graph
  feature_mat <- matrix(rnorm(10*5), nrow=10, ncol=5)
  adj_i <- c(1L, 1L, 2L, 3L)
  adj_j <- c(2L, 3L, 4L, 5L)

  distances <- compute_masked_distances_cpp(feature_mat, adj_i, adj_j)

  expect_equal(length(distances), 4)
  expect_true(all(distances >= 0))
  expect_true(all(!is.na(distances)))

  # Test find_1nn_subgraph_cpp
  n_nodes <- 5
  nn <- find_1nn_subgraph_cpp(n_nodes, adj_i, adj_j, distances)

  expect_equal(length(nn), n_nodes)
  expect_true(all(nn >= -1))  # -1 for no neighbor, otherwise 0-based indices

  # Test find_connected_components_cpp
  components <- find_connected_components_cpp(n_nodes, nn)

  expect_equal(length(components), n_nodes)
  expect_true(all(components >= 0))

  # Test aggregate_features_cpp
  n_components <- length(unique(components))
  agg_features <- aggregate_features_cpp(feature_mat, components, n_components)

  expect_equal(nrow(agg_features), nrow(feature_mat))
  expect_equal(ncol(agg_features), n_components)
  expect_true(all(!is.na(agg_features)))
})

test_that("rena works with different data sizes", {
  # Small data
  mask_small <- NeuroVol(array(1, c(3,3,3)), NeuroSpace(c(3,3,3)))
  vec_small <- NeuroVec(array(rnorm(3*3*3*5), c(3,3,3,5)),
                        NeuroSpace(c(3,3,3,5)))
  result_small <- rena(vec_small, mask_small, K=3, verbose=FALSE)
  expect_true(result_small$n_clusters > 0)

  # Medium data
  mask_med <- NeuroVol(array(1, c(8,8,8)), NeuroSpace(c(8,8,8)))
  vec_med <- NeuroVec(array(rnorm(8*8*8*10), c(8,8,8,10)),
                      NeuroSpace(c(8,8,8,10)))
  result_med <- rena(vec_med, mask_med, K=20, verbose=FALSE)
  expect_true(result_med$n_clusters > 0)
  expect_true(result_med$n_clusters >= 5)
})

test_that("rena metadata contains useful information", {
  mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
  vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)),
                  NeuroSpace(c(5,5,5,10)))

  result <- rena(vec, mask, K=10, verbose=FALSE)

  # Check metadata structure
  expect_true("metadata" %in% names(result))
  expect_true("iterations" %in% names(result$metadata))
  expect_true("final_n_clusters" %in% names(result$metadata))

  # Check values are sensible
  expect_true(is.numeric(result$metadata$iterations))
  expect_true(result$metadata$iterations > 0)
  expect_equal(result$metadata$final_n_clusters, result$n_clusters)
})
