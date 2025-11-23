library(testthat)
library(neurocluster)
library(neuroim2)

# Test helper function
create_test_data <- function(dims = c(10, 10, 10), n_time = 50, seed = 123) {
  set.seed(seed)
  
  # Create a spherical mask
  mask_array <- array(FALSE, dims)
  center <- dims / 2
  radius <- min(dims) / 3
  
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        if (sum(((c(i,j,k) - center) / radius)^2) <= 1) {
          mask_array[i,j,k] <- TRUE
        }
      }
    }
  }
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))
  
  # Use neuroim2::simulate_fmri for proper NeuroVec creation
  vec <- neuroim2::simulate_fmri(mask, n_time = n_time, seed = seed)
  
  nvoxels <- sum(mask_array)
  
  list(vec = vec, mask = mask, nvoxels = nvoxels)
}

test_that("supervoxels_flash3d works with basic inputs", {
  data <- create_test_data(c(10, 10, 10), n_time = 30)
  
  result <- supervoxels_flash3d(data$vec, data$mask, K = 10, verbose = FALSE)
  
  expect_s3_class(result, "cluster_result")
  expect_s3_class(result, "flash3d_result")
  expect_s4_class(result$clusvol, "ClusteredNeuroVol")
  expect_equal(length(result$cluster), data$nvoxels)
  expect_lte(result$K, 10)
  expect_equal(result$method, "flash3d")
  expect_equal(nrow(result$centers), result$K)
  expect_equal(ncol(result$centers), 30)
  expect_equal(nrow(result$coord_centers), result$K)
  expect_equal(ncol(result$coord_centers), 3)
})

test_that("supervoxels_flash3d handles different bit widths", {
  data <- create_test_data(c(8, 8, 8), n_time = 20)
  
  # Test 64-bit hashing
  result_64 <- supervoxels_flash3d(data$vec, data$mask, K = 5, bits = 64, verbose = FALSE)
  expect_s3_class(result_64, "cluster_result")
  expect_lte(result_64$K, 5)
  
  # Test 128-bit hashing
  result_128 <- supervoxels_flash3d(data$vec, data$mask, K = 5, bits = 128, verbose = FALSE)
  expect_s3_class(result_128, "cluster_result")
  expect_lte(result_128$K, 5)
})

test_that("supervoxels_flash3d handles different DCT coefficients", {
  data <- create_test_data(c(8, 8, 8), n_time = 25)
  
  # Test minimum DCT coefficients
  result_min <- supervoxels_flash3d(data$vec, data$mask, K = 4, dctM = 4, verbose = FALSE)
  expect_s3_class(result_min, "cluster_result")
  
  # Test maximum DCT coefficients
  result_max <- supervoxels_flash3d(data$vec, data$mask, K = 4, dctM = 32, verbose = FALSE)
  expect_s3_class(result_max, "cluster_result")
})

test_that("supervoxels_flash3d handles lambda parameters", {
  data <- create_test_data(c(8, 8, 8), n_time = 20)
  
  # High spatial weight (more compact clusters)
  result_spatial <- supervoxels_flash3d(data$vec, data$mask, K = 5, 
                                        lambda_s = 2.0, lambda_t = 0.5, 
                                        verbose = FALSE)
  expect_s3_class(result_spatial, "cluster_result")
  
  # High temporal weight (more feature-similar clusters)
  result_temporal <- supervoxels_flash3d(data$vec, data$mask, K = 5,
                                         lambda_s = 0.1, lambda_t = 2.0,
                                         verbose = FALSE)
  expect_s3_class(result_temporal, "cluster_result")
})

test_that("supervoxels_flash3d handles multiple rounds", {
  data <- create_test_data(c(8, 8, 8), n_time = 20)
  
  # Single round
  result_1 <- supervoxels_flash3d(data$vec, data$mask, K = 4, rounds = 1, verbose = FALSE)
  expect_s3_class(result_1, "cluster_result")
  
  # Multiple rounds
  result_3 <- supervoxels_flash3d(data$vec, data$mask, K = 4, rounds = 3, verbose = FALSE)
  expect_s3_class(result_3, "cluster_result")
})

test_that("supervoxels_flash3d handles voxel scaling", {
  data <- create_test_data(c(8, 8, 8), n_time = 20)
  
  # Anisotropic voxels
  result <- supervoxels_flash3d(data$vec, data$mask, K = 5,
                                vox_scale = c(1, 1, 3), # Z dimension has 3x spacing
                                verbose = FALSE)
  expect_s3_class(result, "cluster_result")
  expect_lte(result$K, 5)
})

test_that("supervoxels_flash3d handles barrier volume", {
  data <- create_test_data(c(8, 8, 8), n_time = 20)
  
  # Create a barrier that discourages crossing the middle
  barrier <- array(0, dim = c(8, 8, 8))
  barrier[4:5, , ] <- 1.0  # High barrier in middle X slices
  
  result <- supervoxels_flash3d(data$vec, data$mask, K = 4,
                                lambda_g = 0.5, barrier = barrier,
                                verbose = FALSE)
  expect_s3_class(result, "cluster_result")
  expect_lte(result$K, 4)
})

test_that("supervoxels_flash3d handles edge cases", {
  # Very small data
  data <- create_test_data(c(5, 5, 5), n_time = 10)
  
  # K = 1 (single cluster)
  result_1 <- supervoxels_flash3d(data$vec, data$mask, K = 1, verbose = FALSE)
  expect_equal(result_1$K, 1)
  expect_true(all(result_1$cluster == 1))
  
  # K close to nvoxels
  nvox <- data$nvoxels
  if (nvox > 10) {
    result_many <- supervoxels_flash3d(data$vec, data$mask, K = nvox - 2, verbose = FALSE)
    expect_lte(result_many$K, nvox - 2)
  }
})

test_that("supervoxels_flash3d produces valid cluster assignments", {
  data <- create_test_data(c(8, 8, 8), n_time = 25)
  
  result <- supervoxels_flash3d(data$vec, data$mask, K = 6, verbose = FALSE)
  
  # Check all voxels are assigned
  expect_false(any(is.na(result$cluster)))
  
  # Check cluster IDs are valid
  expect_true(all(result$cluster %in% 1:result$K))
  
  # Check each cluster has at least one voxel
  cluster_counts <- table(result$cluster)
  expect_true(all(cluster_counts > 0))
  
  # Check centers are finite
  expect_true(all(is.finite(result$centers)))
  expect_true(all(is.finite(result$coord_centers)))
})

test_that("supervoxels_flash3d input validation works", {
  data <- create_test_data(c(6, 6, 6), n_time = 15)
  
  # Invalid K
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 0))
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = -5))
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = data$nvoxels + 10))
  
  # Invalid bits
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 5, bits = 32))
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 5, bits = 256))
  
  # Invalid dctM
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 5, dctM = 3))
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 5, dctM = 33))
  
  # Invalid rounds
  expect_error(supervoxels_flash3d(data$vec, data$mask, K = 5, rounds = 0))
  
  # Empty mask
  empty_mask <- NeuroVol(array(0, c(6, 6, 6)), NeuroSpace(c(6, 6, 6)))
  expect_error(supervoxels_flash3d(data$vec, empty_mask, K = 5))
})

test_that("supervoxels_flash3d handles K = nvoxels correctly", {
  data <- create_test_data(c(5, 5, 5), n_time = 10)
  
  expect_warning(
    result <- supervoxels_flash3d(data$vec, data$mask, K = data$nvoxels, verbose = FALSE),
    "K equals number of voxels"
  )
  
  expect_equal(result$K, data$nvoxels)
})

test_that("supervoxels_flash3d is deterministic", {
  data <- create_test_data(c(8, 8, 8), n_time = 20, seed = 42)
  
  # Run twice with same parameters
  result1 <- supervoxels_flash3d(data$vec, data$mask, K = 5, 
                                 lambda_s = 0.6, lambda_t = 1.0,
                                 bits = 64, dctM = 12, rounds = 2,
                                 verbose = FALSE)
  
  result2 <- supervoxels_flash3d(data$vec, data$mask, K = 5,
                                 lambda_s = 0.6, lambda_t = 1.0,
                                 bits = 64, dctM = 12, rounds = 2,
                                 verbose = FALSE)
  
  # Results should be identical
  expect_equal(result1$cluster, result2$cluster)
  expect_equal(result1$K, result2$K)
  expect_equal(result1$centers, result2$centers)
  expect_equal(result1$coord_centers, result2$coord_centers)
})