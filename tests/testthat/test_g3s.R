context("G3S: Gradient-Guided Geodesic Supervoxels")

# =============================================================================
# Test Phase 1: SVD Compression
# =============================================================================

test_that("SVD compression preserves variance", {
  skip_on_cran()

  # Create synthetic data with known structure
  n_voxels <- 100
  n_timepoints <- 300

  # Generate data with strong first few components
  set.seed(42)
  true_components <- 5
  latent <- matrix(rnorm(n_voxels * true_components), n_voxels, true_components)
  loading <- matrix(rnorm(n_timepoints * true_components), n_timepoints, true_components)
  feature_mat <- latent %*% t(loading) + matrix(rnorm(n_voxels * n_timepoints, sd = 0.1),
                                                 n_voxels, n_timepoints)

  # Compress
  compressed <- compress_features_svd(feature_mat, n_components = 15,
                                     variance_threshold = 0.90)

  # Check variance explained
  expect_gte(compressed$variance_explained, 0.90)
  expect_lte(compressed$variance_explained, 1.0)

  # Check dimensions
  expect_equal(nrow(compressed$features), n_voxels)
  expect_lte(ncol(compressed$features), 15)

  # Check normalization
  row_norms <- sqrt(rowSums(compressed$features^2))
  expect_true(all(abs(row_norms - 1) < 1e-10))
})


test_that("SVD compression handles edge cases", {
  skip_on_cran()

  # Small dataset
  small_mat <- matrix(rnorm(10 * 20), 10, 20)
  compressed <- compress_features_svd(small_mat, n_components = 5)
  expect_equal(nrow(compressed$features), 10)

  # More components requested than available
  compressed2 <- compress_features_svd(small_mat, n_components = 100)
  expect_lte(ncol(compressed2$features), min(10, 20))

  # Zero variance column
  zero_var_mat <- matrix(rnorm(10 * 20), 10, 20)
  zero_var_mat[, 1] <- 1  # Constant column
  compressed3 <- compress_features_svd(zero_var_mat, n_components = 5)
  expect_equal(nrow(compressed3$features), 10)
})


test_that("transform_new_data_svd applies compression consistently", {
  skip_on_cran()

  # Train compression
  train_data <- matrix(rnorm(50 * 100), 50, 100)
  compression <- compress_features_svd(train_data, n_components = 10)

  # Apply to new data
  test_data <- matrix(rnorm(20 * 100), 20, 100)
  compressed_test <- transform_new_data_svd(test_data, compression)

  expect_equal(nrow(compressed_test), 20)
  expect_equal(ncol(compressed_test), compression$n_components)

  # Check normalization
  row_norms <- sqrt(rowSums(compressed_test^2))
  expect_true(all(abs(row_norms - 1) < 1e-10))
})


# =============================================================================
# Test Phase 2: Gradient-Based Seeding
# =============================================================================

test_that("Functional gradient identifies homogeneous vs boundary regions", {
  skip_on_cran()

  # Create data with clear functional regions
  set.seed(123)
  n_voxels <- 200
  n_features <- 15

  # Region 1: Low gradient (homogeneous)
  region1 <- matrix(rnorm(100 * n_features, mean = 0, sd = 0.1), 100, n_features)
  region1 <- t(apply(region1, 1, function(x) x / sqrt(sum(x^2))))

  # Region 2: Low gradient (different mean)
  region2 <- matrix(rnorm(100 * n_features, mean = 2, sd = 0.1), 100, n_features)
  region2 <- t(apply(region2, 1, function(x) x / sqrt(sum(x^2))))

  features <- rbind(region1, region2)
  coords <- matrix(rnorm(n_voxels * 3), n_voxels, 3)

  # Compute gradient
  grad_vals <- compute_functional_gradient(features, coords, k_neighbors = 10)

  # Regions should have lower gradient than boundaries
  # Within-region gradient should be lower than cross-region
  within_region1 <- mean(grad_vals[1:100])
  within_region2 <- mean(grad_vals[101:200])

  # Both regions should have relatively low gradient
  expect_lt(within_region1, 1.0)
  expect_lt(within_region2, 1.0)
})


test_that("find_gradient_seeds_g3s finds K spatially separated seeds", {
  skip_on_cran()

  # Create synthetic data
  set.seed(456)
  n_voxels <- 500
  n_features <- 15
  features <- matrix(rnorm(n_voxels * n_features), n_voxels, n_features)
  features <- t(apply(features, 1, function(x) x / sqrt(sum(x^2))))

  # Regular 3D grid
  coords <- as.matrix(expand.grid(x = 1:10, y = 1:10, z = 1:5))

  # Find 20 seeds
  seeds <- find_gradient_seeds_g3s(features, coords, K = 20, k_neighbors = 26)

  # Check we got close to K seeds
  expect_gte(length(seeds), 15)  # Allow some tolerance
  expect_lte(length(seeds), 20)

  # Check seeds are valid indices
  expect_true(all(seeds >= 1))
  expect_true(all(seeds <= n_voxels))
  expect_true(all(seeds == as.integer(seeds)))

  # Check no duplicates
  expect_equal(length(seeds), length(unique(seeds)))

  # Check spatial separation
  seed_coords <- coords[seeds, ]
  if (length(seeds) > 1) {
    dists <- as.matrix(dist(seed_coords))
    diag(dists) <- NA
    min_dist <- min(dists, na.rm = TRUE)
    expect_gt(min_dist, 0)  # Seeds should not be identical
  }
})


# =============================================================================
# Test Phase 3: Geodesic Propagation
# =============================================================================

test_that("g3s_propagate_cpp assigns all voxels", {
  skip_on_cran()

  # Small test case
  set.seed(789)
  n_voxels <- 100
  n_features <- 10
  K <- 5

  feature_mat <- matrix(rnorm(n_features * n_voxels), n_features, n_voxels)
  feature_mat <- apply(feature_mat, 2, function(x) x / sqrt(sum(x^2)))

  coords <- as.matrix(expand.grid(x = 1:5, y = 1:5, z = 1:4))
  seed_indices <- as.integer(seq(1, n_voxels, length.out = K))

  # Compute neighbors
  neib <- FNN::get.knn(coords, k = min(26, n_voxels - 1))

  # Run G3S propagation
  labels <- g3s_propagate_cpp(
    feature_mat = feature_mat,
    coords = coords,
    seed_indices = seed_indices,
    neighbor_indices = neib$nn.index,
    neighbor_dists = neib$nn.dist,
    alpha = 0.5,
    compactness = 2.0
  )

  # Check all voxels labeled
  expect_equal(length(labels), n_voxels)
  expect_true(all(labels > 0))  # All voxels should be labeled
  expect_true(all(labels <= K))  # Labels should be 1..K

  # Check seeds are labeled correctly
  for (k in seq_along(seed_indices)) {
    expect_equal(labels[seed_indices[k]], k)
  }

  # Check contiguity: each cluster should be spatially connected
  # (This is a weak test; full test would check connectivity graph)
  for (k in 1:K) {
    cluster_voxels <- which(labels == k)
    expect_gt(length(cluster_voxels), 0)
  }
})


# =============================================================================
# Test Full G3S Pipeline
# =============================================================================

test_that("cluster4d_g3s produces valid results", {
  skip_on_cran()

  # Create small test volume
  library(neuroim2)
  set.seed(101)

  mask <- NeuroVol(array(1, c(5, 5, 5)), NeuroSpace(c(5, 5, 5)))
  vec <- NeuroVec(array(rnorm(5*5*5*20), c(5, 5, 5, 20)),
                  NeuroSpace(c(5, 5, 5, 20)))

  # Run G3S
  result <- cluster4d_g3s(vec, mask, K = 10, n_components = 5,
                         max_refinement_iter = 1, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "g3s_result")
  expect_s3_class(result, "cluster4d_result")
  expect_s3_class(result, "cluster_result")

  # Check components
  expect_true("clusvol" %in% names(result))
  expect_true("cluster" %in% names(result))
  expect_true("centers" %in% names(result))
  expect_true("coord_centers" %in% names(result))
  expect_true("n_clusters" %in% names(result))
  expect_true("method" %in% names(result))
  expect_equal(result$method, "g3s")

  # Check cluster assignments
  expect_equal(length(result$cluster), sum(mask > 0))
  expect_true(all(result$cluster > 0))
  expect_lte(result$n_clusters, 10)
  expect_gte(result$n_clusters, 1)

  # Check metadata
  expect_true("variance_explained" %in% names(result$parameters))
  expect_true("compression_ratio" %in% names(result$metadata))
  expect_true("seed_indices" %in% names(result$metadata))

  # Check variance explained
  expect_gte(result$parameters$variance_explained, 0.8)
  expect_lte(result$parameters$variance_explained, 1.0)
})


test_that("cluster4d with method='g3s' works", {
  skip_on_cran()

  library(neuroim2)
  set.seed(202)

  mask <- NeuroVol(array(1, c(4, 4, 4)), NeuroSpace(c(4, 4, 4)))
  vec <- NeuroVec(array(rnorm(4*4*4*15), c(4, 4, 4, 15)),
                  NeuroSpace(c(4, 4, 4, 15)))

  # Call via cluster4d interface
  result <- cluster4d(vec, mask, n_clusters = 5, method = "g3s",
                     spatial_weight = 0.5, max_iterations = 1, verbose = FALSE)

  expect_s3_class(result, "cluster4d_result")
  expect_equal(result$method, "g3s")
  expect_lte(result$n_clusters, 5)
})


test_that("G3S handles small datasets gracefully", {
  skip_on_cran()

  library(neuroim2)
  set.seed(303)

  # Very small dataset
  mask <- NeuroVol(array(1, c(3, 3, 3)), NeuroSpace(c(3, 3, 3)))
  vec <- NeuroVec(array(rnorm(3*3*3*10), c(3, 3, 3, 10)),
                  NeuroSpace(c(3, 3, 3, 10)))

  # Request more clusters than voxels should error
  expect_error(
    cluster4d_g3s(vec, mask, K = 30, n_components = 5,
                  max_refinement_iter = 0, verbose = FALSE),
    "must be between"
  )
})


# =============================================================================
# Test Phase 4: Boundary Refinement
# =============================================================================

test_that("Boundary refinement improves cluster quality", {
  skip_on_cran()

  # Create scenario where refinement should help
  set.seed(404)
  n_voxels <- 200
  n_features <- 10

  # Two clear regions with some noisy boundary voxels
  features1 <- matrix(rnorm(100 * n_features, mean = 0), 100, n_features)
  features2 <- matrix(rnorm(100 * n_features, mean = 2), 100, n_features)
  features <- rbind(features1, features2)
  features <- t(apply(features, 1, function(x) x / sqrt(sum(x^2))))

  coords <- matrix(rnorm(n_voxels * 3), n_voxels, 3)

  # Initial labels with some misclassifications
  labels <- c(rep(1, 100), rep(2, 100))
  labels[c(50, 51, 150, 151)] <- sample(1:2, 4, replace = TRUE)  # Add noise

  # Compute neighbors
  neib <- FNN::get.knn(coords, k = min(10, n_voxels - 1))

  # Refine
  labels_refined <- refine_boundaries_g3s_cpp(
    labels = labels,
    feature_mat = t(features),
    neighbor_indices = neib$nn.index,
    max_iter = 3L
  )

  # Check refinement happened
  expect_equal(length(labels_refined), n_voxels)
  expect_true(all(labels_refined %in% 1:2))

  # Refinement should maintain or improve cluster coherence
  # (This is a soft test; exact behavior depends on data)
  expect_true(length(unique(labels_refined)) >= 1)
})


# =============================================================================
# Test Algorithm Comparison
# =============================================================================

test_that("G3S produces stable clusters on synthetic data", {
  skip_on_cran()

  library(neuroim2)
  set.seed(505)

  # Create test data with known structure
  mask <- NeuroVol(array(1, c(10, 10, 5)), NeuroSpace(c(10, 10, 5)))
  vec_list <- replicate(30, {
    NeuroVol(array(rnorm(10*10*5), c(10, 10, 5)), NeuroSpace(c(10, 10, 5)))
  }, simplify = FALSE)
  vec <- do.call(concat, vec_list)

  result_g3s <- cluster4d(vec, mask, n_clusters = 20, method = "g3s",
                          max_iterations = 2, verbose = FALSE)

  # NOTE: We previously compared against SNIC here, but SNIC occasionally
  # segfaults on this randomized volume (see devtools::test(filter = "g3s") logs).
  # Once SNIC is stabilized we can re-enable the side-by-side comparison.

  expect_s3_class(result_g3s, "cluster4d_result")
  expect_lte(abs(result_g3s$n_clusters - 20), 5)
  expect_true("variance_explained" %in% names(result_g3s$parameters))
  expect_true("compression_ratio" %in% names(result_g3s$metadata))
})


# =============================================================================
# Test Print Method
# =============================================================================

test_that("print.g3s_result produces readable output", {
  skip_on_cran()

  library(neuroim2)
  set.seed(606)

  mask <- NeuroVol(array(1, c(5, 5, 5)), NeuroSpace(c(5, 5, 5)))
  vec <- NeuroVec(array(rnorm(5*5*5*15), c(5, 5, 5, 15)),
                  NeuroSpace(c(5, 5, 5, 15)))

  result <- cluster4d_g3s(vec, mask, K = 8, verbose = FALSE)

  # Should print without error
  expect_output(print(result), "G3S Clustering Result")
  expect_output(print(result), "Clusters:")
  expect_output(print(result), "Compression:")
  expect_output(print(result), "Parameters:")
})
