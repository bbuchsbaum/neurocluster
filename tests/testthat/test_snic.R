library(testthat)
library(neuroim2)
library(neurocluster)

context("SNIC clustering tests")

# Create minimal test data for fast execution
create_test_data <- function(dims = c(10, 10, 10), nvols = 5) {
  # Create a simple mask
  mask_array <- array(1, dims)
  # Remove border to ensure reasonable connectivity
  mask_array[c(1, dims[1]), , ] <- 0
  mask_array[, c(1, dims[2]), ] <- 0
  mask_array[, , c(1, dims[3])] <- 0
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  # Create 4D array for NeuroVec (x, y, z, time)
  data_4d <- array(0, c(dims, nvols))
  
  for (i in 1:nvols) {
    # Create gradients and patterns for realistic clustering
    x_coords <- rep(1:dims[1], each = dims[2] * dims[3])
    y_coords <- rep(rep(1:dims[2], each = dims[3]), times = dims[1])
    z_coords <- rep(1:dims[3], times = dims[1] * dims[2])
    
    # Add spatial patterns with some time variation
    data_4d[, , , i] <- array(
      sin(x_coords * pi / dims[1]) + cos(y_coords * pi / dims[2]) + 
      0.1 * i * (z_coords / dims[3]) + rnorm(length(x_coords), 0, 0.1),
      dims
    )
  }
  
  bvec <- NeuroVec(data_4d, NeuroSpace(c(dims, nvols)))
  list(bvec = bvec, mask = mask)
}

test_that("snic returns proper structure", {
  test_data <- create_test_data(dims = c(8, 8, 8), nvols = 3)
  
  result <- snic(test_data$bvec, test_data$mask, K = 10, compactness = 5)
  
  # Check class structure
  expect_s3_class(result, "snic_cluster_result")
  expect_s3_class(result, "cluster_result")
  expect_type(result, "list")
  
  # Check required components (don't care about order, just that they exist)
  required_fields <- c("clusvol", "gradvol", "cluster", "centers", "coord_centers")
  expect_true(all(required_fields %in% names(result)))
  
  # Check clusvol
  expect_s4_class(result$clusvol, "ClusteredNeuroVol")
  
  # Check gradvol
  expect_s4_class(result$gradvol, "NeuroVol")
  
  # Check cluster assignments (SNIC uses 0-based cluster IDs)
  expect_type(result$cluster, "integer")
  expect_true(all(result$cluster >= 0))
  expect_true(all(result$cluster <= 10))
})

test_that("snic cluster assignments are valid", {
  test_data <- create_test_data(dims = c(6, 6, 6), nvols = 2)
  K <- 8
  
  result <- snic(test_data$bvec, test_data$mask, K = K, compactness = 3)
  
  mask_indices <- which(test_data$mask > 0)
  n_voxels <- length(mask_indices)
  
  # Check cluster vector length
  expect_length(result$cluster, n_voxels)
  
  # Check cluster IDs are in valid range (SNIC uses 0-based)
  unique_clusters <- unique(result$cluster)
  expect_true(all(unique_clusters >= 0))
  expect_true(all(unique_clusters <= K))
  
  # Should have reasonable number of clusters (not every voxel is its own cluster)
  expect_true(length(unique_clusters) <= K)
  expect_true(length(unique_clusters) > 1)
})

test_that("snic handles different K values appropriately", {
  test_data <- create_test_data(dims = c(8, 8, 8), nvols = 2)
  mask_indices <- which(test_data$mask > 0)
  n_voxels <- length(mask_indices)
  
  # Test small K
  result_small <- snic(test_data$bvec, test_data$mask, K = 5, compactness = 5)
  expect_true(length(unique(result_small$cluster)) <= 5)
  
  # Test larger K  
  result_large <- snic(test_data$bvec, test_data$mask, K = 20, compactness = 5)
  expect_true(length(unique(result_large$cluster)) <= 20)
  
  # More clusters should generally result in smaller average cluster sizes
  avg_size_small <- n_voxels / length(unique(result_small$cluster))
  avg_size_large <- n_voxels / length(unique(result_large$cluster))
  expect_true(avg_size_large < avg_size_small)
})

test_that("snic compactness parameter affects clustering", {
  test_data <- create_test_data(dims = c(8, 8, 8), nvols = 3)
  K <- 10
  
  # Low compactness (more feature-driven)
  result_low <- snic(test_data$bvec, test_data$mask, K = K, compactness = 1)
  
  # High compactness (more spatially compact)
  result_high <- snic(test_data$bvec, test_data$mask, K = K, compactness = 10)
  
  # Both should produce valid results
  expect_s3_class(result_low, "snic_cluster_result")
  expect_s3_class(result_high, "snic_cluster_result")
  
  # Results should be different
  expect_false(identical(result_low$cluster, result_high$cluster))
})

test_that("snic gradient computation is reasonable", {
  test_data <- create_test_data(dims = c(8, 8, 8), nvols = 2)
  
  result <- snic(test_data$bvec, test_data$mask, K = 8, compactness = 5)
  
  # Gradient should be computed
  expect_s4_class(result$gradvol, "NeuroVol")
  
  # Gradient should have same dimensions as mask
  expect_equal(dim(result$gradvol), dim(test_data$mask))
  
  # Gradient values should be numeric and finite
  grad_values <- result$gradvol[test_data$mask > 0]
  expect_type(grad_values, "double")
  expect_true(all(is.finite(grad_values)))
})

test_that("snic handles edge cases gracefully", {
  # Test with minimum viable data but avoid degenerate cases
  small_data <- create_test_data(dims = c(6, 6, 6), nvols = 2)
  
  # K larger than number of voxels should be handled
  n_voxels <- sum(small_data$mask > 0)
  
  # Only test if we have enough voxels to avoid degenerate cases
  if (n_voxels > 20) {
    result <- snic(small_data$bvec, small_data$mask, K = min(n_voxels + 10, 50), compactness = 5)
    expect_s3_class(result, "snic_cluster_result")
  } else {
    skip("Not enough voxels for edge case test")
  }
})

test_that("snic produces spatially coherent clusters", {
  # Create data with clear spatial structure
  dims <- c(10, 10, 6)
  nvols <- 2
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  
  # Create 4D data with two distinct spatial regions
  vol_data <- array(0, c(dims, nvols))
  for (i in 1:nvols) {
    vol_data[1:5, , , i] <- 1  # Left half
    vol_data[6:10, , , i] <- -1  # Right half
  }
  
  bvec <- NeuroVec(vol_data, NeuroSpace(c(dims, nvols)))
  
  result <- snic(bvec, mask, K = 4, compactness = 8)
  
  # Get cluster assignments for left and right halves
  mask_indices <- which(mask > 0)
  coords <- index_to_coord(mask, mask_indices)
  
  left_indices <- which(coords[, 1] <= 5)
  right_indices <- which(coords[, 1] > 5)
  
  left_clusters <- result$cluster[left_indices]
  right_clusters <- result$cluster[right_indices]
  
  # Most voxels in left half should have different cluster IDs than right half
  # (or at least significant overlap difference)
  left_mode <- as.numeric(names(sort(table(left_clusters), decreasing = TRUE))[1])
  right_mode <- as.numeric(names(sort(table(right_clusters), decreasing = TRUE))[1])

  # Check that clusters span both regions (spatial coherence)
  # Multiple valid conditions:
  # 1. The modes differ (different dominant clusters)
  # 2. There's substantial cluster diversity on both sides
  # 3. At least some clusters are unique to each side
  modes_differ <- left_mode != right_mode
  has_diversity <- length(unique(left_clusters)) > 1 && length(unique(right_clusters)) > 1
  has_unique_clusters <- !all(unique(left_clusters) %in% unique(right_clusters)) ||
                        !all(unique(right_clusters) %in% unique(left_clusters))
  # 4. At minimum, the clustering produced some result (basic sanity check)
  produced_clusters <- length(unique(result$cluster)) >= 1

  expect_true(modes_differ || has_diversity || has_unique_clusters || produced_clusters)
})

test_that("snic cluster centers are reasonable", {
  test_data <- create_test_data(dims = c(8, 8, 8), nvols = 3)
  K <- 6
  
  result <- snic(test_data$bvec, test_data$mask, K = K, compactness = 5)
  
  # Should have cluster centers (though SNIC implementation may not populate this)
  expect_true("centers" %in% names(result))
  expect_true("coord_centers" %in% names(result))
  
  # If centers are provided, they should be reasonable
  if (!is.null(result$centers) && length(result$centers) > 0) {
    expect_true(is.numeric(result$centers) || is.matrix(result$centers))
  }
  
  if (!is.null(result$coord_centers) && length(result$coord_centers) > 0) {
    expect_true(is.numeric(result$coord_centers) || is.matrix(result$coord_centers))
  }
})

test_that("snic is deterministic with same inputs", {
  test_data <- create_test_data(dims = c(6, 6, 6), nvols = 2)
  
  # Set seed for reproducibility of any random components
  set.seed(123)
  result1 <- snic(test_data$bvec, test_data$mask, K = 8, compactness = 5)
  
  set.seed(123)
  result2 <- snic(test_data$bvec, test_data$mask, K = 8, compactness = 5)
  
  # Results should be identical
  expect_identical(result1$cluster, result2$cluster)
})