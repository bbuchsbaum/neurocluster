library(testthat)
library(neurocluster)
library(neuroim2)

# Helper function to create simulated fMRI data
create_test_data <- function(dims, n_time = 50, mask_prop = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create a spherical mask
  mask_array <- array(FALSE, dims)
  center <- dims / 2
  radius <- min(dims) / 3
  
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        if (sum(((c(i,j,k) - center) / radius)^2) <= mask_prop) {
          mask_array[i,j,k] <- TRUE
        }
      }
    }
  }
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))
  
  # Use neuroim2::simulate_fmri for realistic data
  bvec <- neuroim2::simulate_fmri(mask, n_time = n_time, seed = seed)
  
  list(mask = mask, bvec = bvec, nvoxels = sum(mask_array))
}

test_that("supervoxels handles K > nvoxels correctly", {
  # Create very small data
  data <- create_test_data(c(8,8,8), n_time = 20, mask_prop = 0.3, seed = 123)
  
  # Should fail when K > nvoxels
  expect_error(
    supervoxels(data$bvec, data$mask, K = data$nvoxels + 10),
    "Cannot create .* clusters from .* masked voxels"
  )
})

test_that("supervoxels handles K = nvoxels correctly", {
  # Create very small data
  data <- create_test_data(c(8,8,8), n_time = 20, mask_prop = 0.3, seed = 124)
  
  # Should return warning and create one cluster per voxel
  expect_warning(
    res <- supervoxels(data$bvec, data$mask, K = data$nvoxels, iterations = 2),
    "K equals number of voxels"
  )
  
  # Each voxel should have unique cluster
  expect_equal(length(unique(res$cluster)), data$nvoxels)
})

test_that("supervoxels works with small data and various K values", {
  # Small data: 10x10x10
  data <- create_test_data(c(10,10,10), n_time = 30, mask_prop = 0.6, seed = 125)
  
  # Test various K values
  k_values <- c(5, 10, 20, 50)
  
  for (k in k_values) {
    if (k < data$nvoxels) {
      res <- supervoxels(data$bvec, data$mask, K = k, 
                        iterations = 10, verbose = FALSE)
      expect_s3_class(res, "cluster_result")
      expect_true(length(unique(res$cluster)) <= k)
      expect_equal(nrow(res$centers), length(unique(res$cluster)))
    }
  }
})

test_that("supervoxels works with medium data and various K values", {
  # Medium data: 32x32x20
  data <- create_test_data(c(32,32,20), n_time = 50, mask_prop = 0.7, seed = 126)
  
  # Test various K values
  k_values <- c(50, 100, 200)
  
  for (k in k_values) {
    if (k < data$nvoxels) {
      res <- supervoxels(data$bvec, data$mask, K = k, 
                        iterations = 15, verbose = FALSE,
                        converge_thresh = 0.01)  # Less strict for faster testing
      expect_s3_class(res, "cluster_result")
      expect_true(length(unique(res$cluster)) <= k)
    }
  }
})

test_that("supervoxels converges properly with different parameters", {
  # Create moderate size data
  data <- create_test_data(c(20,20,15), n_time = 40, mask_prop = 0.7, seed = 127)
  
  # Test with different convergence thresholds
  res1 <- supervoxels(data$bvec, data$mask, K = 30,
                     iterations = 50, converge_thresh = 0.1,  # Loose threshold
                     verbose = FALSE)
  
  res2 <- supervoxels(data$bvec, data$mask, K = 30,
                     iterations = 50, converge_thresh = 0.001,  # Strict threshold
                     verbose = FALSE)
  
  # Both should produce valid results
  expect_s3_class(res1, "cluster_result")
  expect_s3_class(res2, "cluster_result")
})

test_that("supervoxels handles different sigma and alpha parameters", {
  data <- create_test_data(c(15,15,15), n_time = 30, seed = 128)
  
  param_grid <- expand.grid(
    sigma1 = c(0.5, 1, 2),
    sigma2 = c(2, 4, 6),
    alpha = c(0.2, 0.5, 0.8)
  )
  
  # Test a subset of parameter combinations
  test_rows <- sample(nrow(param_grid), min(6, nrow(param_grid)))
  
  for (i in test_rows) {
    params <- param_grid[i, ]
    
    res <- supervoxels(data$bvec, data$mask, K = 20,
                      sigma1 = params$sigma1,
                      sigma2 = params$sigma2,
                      alpha = params$alpha,
                      iterations = 10,
                      verbose = FALSE)
    
    expect_s3_class(res, "cluster_result")
  }
})

test_that("supervoxels parallel and sequential produce similar results", {
  data <- create_test_data(c(25,25,15), n_time = 30, seed = 129)
  
  # Ensure we have enough voxels for parallel to activate
  if (data$nvoxels > 1000) {
    set.seed(456)
    res_par <- supervoxels(data$bvec, data$mask, K = 50,
                          iterations = 10, parallel = TRUE,
                          verbose = FALSE)
    
    set.seed(456)
    res_seq <- supervoxels(data$bvec, data$mask, K = 50,
                          iterations = 10, parallel = FALSE,
                          verbose = FALSE)
    
    # Check that both produce valid results (exact match not guaranteed due to numerical differences)
    expect_s3_class(res_par, "cluster_result")
    expect_s3_class(res_seq, "cluster_result")
    expect_equal(length(unique(res_par$cluster)), length(unique(res_seq$cluster)), tolerance = 5)
  }
})

test_that("supervoxels handles edge cases gracefully", {
  # Test with single time point - create proper NeuroVec
  dims <- c(10,10,10)
  mask <- NeuroVol(array(TRUE, dims), NeuroSpace(dims, c(1,1,1)))
  
  # Create a single time point as a list and concatenate
  vol_list <- list(NeuroVol(array(rnorm(prod(dims)), dims), NeuroSpace(dims, c(1,1,1))))
  bvec <- do.call(concat, vol_list)
  
  # Should work even with single time point
  res <- supervoxels(bvec, mask, K = 10, iterations = 5, verbose = FALSE)
  expect_s3_class(res, "cluster_result")
  
  # Test with very sparse mask
  sparse_mask_array <- array(FALSE, c(10,10,10))
  sparse_mask_array[5,5,1:5] <- TRUE  # Only 5 voxels
  sparse_mask <- NeuroVol(sparse_mask_array, NeuroSpace(c(10,10,10), c(1,1,1)))
  
  res <- supervoxels(bvec, sparse_mask, K = 2, iterations = 5, verbose = FALSE)
  expect_true(length(unique(res$cluster)) <= 2)
})

test_that("supervoxels doesn't hang with large K on small data", {
  # This is the key test for the hanging issue
  data <- create_test_data(c(20,20,10), n_time = 50, mask_prop = 0.5, seed = 130)
  
  # Calculate reasonable K values
  max_k <- min(500, floor(data$nvoxels * 0.8))
  
  # This should complete without hanging
  time_taken <- system.time({
    res <- supervoxels(data$bvec, data$mask, K = max_k,
                      iterations = 20,
                      verbose = FALSE,
                      converge_thresh = 0.01)
  })
  
  # Should complete in reasonable time (< 60 seconds for small data)
  expect_true(time_taken[3] < 60, 
             info = paste("Took", time_taken[3], "seconds"))
  
  expect_s3_class(res, "cluster_result")
  expect_true(length(unique(res$cluster)) <= max_k)
})

test_that("supervoxels verbose mode provides useful information", {
  data <- create_test_data(c(10,10,10), n_time = 20, seed = 131)
  
  # Capture messages
  messages <- capture_messages({
    res <- supervoxels(data$bvec, data$mask, K = 10,
                      iterations = 10, verbose = TRUE,
                      converge_thresh = 0.05)
  })
  
  # Should include convergence information
  expect_true(any(grepl("% of voxels", messages)) || 
              any(grepl("converged", messages)) ||
              any(grepl("no improvement", messages)))
})

test_that("supervoxels handles constant time series", {
  data <- create_test_data(c(10,10,10), n_time = 20, seed = 132)
  
  # Create constant time series (all values the same)
  const_array <- array(1, c(10,10,10,20))
  vec_list <- lapply(1:20, function(i) {
    NeuroVol(const_array[,,,i], NeuroSpace(c(10,10,10), c(3,3,3)))
  })
  bvec_const <- do.call(concat, vec_list)
  
  # Should handle constant values gracefully (with warning about scaling)
  expect_warning({
    res <- supervoxels(bvec_const, data$mask, K = 10,
                      iterations = 5, verbose = FALSE)
  }, "NA|scaling")
  
  expect_s3_class(res, "cluster_result")
})

test_that("supervoxels accepts 3D NeuroVol input", {
  # Test that supervoxels can handle 3D NeuroVol (structural images)
  # This is useful for T1-weighted MRI segmentation
  dims <- c(12, 12, 8)
  vol_data <- array(rnorm(prod(dims)), dims)
  vol <- NeuroVol(vol_data, NeuroSpace(dims, c(2, 2, 2)))
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims, c(2, 2, 2)))

  # Should work with 3D NeuroVol input (automatically converted internally)
  res <- supervoxels(vol, mask, K = 8, iterations = 5, verbose = FALSE)

  expect_s3_class(res, "cluster_result")
  expect_true(length(unique(res$cluster)) <= 8)
  expect_true(length(res$cluster) == prod(dims))
})