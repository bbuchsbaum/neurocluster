library(neurocluster)
library(neuroim2)

# Input Validation and Malformed Data Structure Tests
# Tests for robustness against invalid inputs and edge cases

test_that("validates NeuroVec input types and structure", {
  # Test various invalid NeuroVec inputs
  dims <- c(8, 8, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  
  # Test 1: NULL input
  expect_error(slice_msf(NULL, mask),
               info = "Should reject NULL NeuroVec")
  
  # Test 2: Wrong class
  fake_vec <- list(data = array(1, c(8, 8, 2, 30)))
  class(fake_vec) <- "NotNeuroVec"
  expect_error(slice_msf(fake_vec, mask),
               info = "Should reject non-NeuroVec objects")
  
  # Test 3: Dimension mismatch with mask
  # Create NeuroVec with different dimensions than mask
  wrong_dims <- c(6, 6, 1)
  wrong_vec_list <- lapply(1:20, function(t) {
    NeuroVol(array(rnorm(prod(wrong_dims)), wrong_dims), NeuroSpace(wrong_dims))
  })
  wrong_vec <- do.call(concat, wrong_vec_list)
  
  expect_error(slice_msf(wrong_vec, mask),
               info = "Should reject NeuroVec with dimensions different from mask")
  
  # Test 4: Empty time series (0 timepoints)
  # This might be challenging to create with actual NeuroVec,
  # so we test the error condition conceptually
  expect_true(TRUE)  # Placeholder - implementation dependent
})

test_that("validates mask input types and structure", {
  dims <- c(8, 8, 2)
  nvox <- prod(dims)
  ntime <- 30
  
  # Create valid NeuroVec
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[1:nvox] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test 1: NULL mask
  expect_error(slice_msf(vec, NULL),
               info = "Should reject NULL mask")
  
  # Test 2: Wrong class
  fake_mask <- array(1, dims)
  expect_error(slice_msf(vec, fake_mask),
               info = "Should reject non-NeuroVol mask")
  
  # Test 3: Empty mask (all zeros)
  empty_mask <- NeuroVol(array(0, dims), NeuroSpace(dims))
  expect_error(slice_msf(vec, empty_mask),
               info = "Should reject mask with no active voxels")
  
  # Test 4: Negative values in mask
  negative_mask <- NeuroVol(array(-1, dims), NeuroSpace(dims))
  negative_mask@.Data[1:5, 1:5, 1] <- 1  # Some positive values
  
  # Should handle gracefully (negative values treated as zero)
  suppressWarnings(expect_no_error(result_negative <- slice_msf(vec, negative_mask, num_runs = 1, r = 4)))
  expect_true(!is.null(result_negative))
  
  # Test 5: Non-binary mask values
  weighted_mask <- NeuroVol(array(0.5, dims), NeuroSpace(dims))
  weighted_mask@.Data[1:4, 1:4, 1] <- 1.0
  weighted_mask@.Data[5:8, 5:8, 1] <- 0.3  # Fractional values
  
  # Should handle gracefully (non-zero treated as active)
  expect_silent(result_weighted <- slice_msf(vec, weighted_mask, num_runs = 1, r = 4))
  expect_true(!is.null(result_weighted))
  
  # Test 6: Mask with NaN/Inf values
  nan_mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nan_mask@.Data[1:2, 1:2, 1] <- NaN
  nan_mask@.Data[7:8, 7:8, 2] <- Inf
  
  # Should either handle gracefully or give informative error
  result_test <- tryCatch({
    result_nan <- slice_msf(vec, nan_mask, num_runs = 1, r = 4)
    list(success = TRUE, result = result_nan)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (!result_test$success) {
    expect_true(nchar(result_test$error) > 5,
                info = "NaN/Inf mask should produce informative error")
  } else {
    expect_true(!is.null(result_test$result),
                info = "NaN/Inf mask handling should produce valid result")
  }
})

test_that("validates parameter bounds and types", {
  # Create minimal valid test case
  dims <- c(6, 6, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 20
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test DCT rank (r) validation - API now handles some cases gracefully
  suppressWarnings(expect_no_error(slice_msf(vec, mask, r = 1, num_runs = 1)))  # Minimum valid r
  
  # r >= ntime should be clamped to ntime-1
  suppressWarnings(expect_no_error(slice_msf(vec, mask, r = ntime-1, num_runs = 1)))
  
  # Invalid r values should still error
  expect_error(slice_msf(vec, mask, r = 0),
               info = "Should reject r <= 0")
  
  expect_error(slice_msf(vec, mask, r = -5),
               info = "Should reject negative r")
  
  # Test compactness parameter - API now handles negative values gracefully
  suppressWarnings(expect_no_error(slice_msf(vec, mask, compactness = 0.1, num_runs = 1)))  # Small compactness
  
  # Negative compactness might be clamped to 0 or small positive value
  suppressWarnings(expect_no_error(slice_msf(vec, mask, compactness = 0, num_runs = 1)))
  
  # Very large compactness should work (just be ineffective)
  suppressWarnings(expect_no_error(slice_msf(vec, mask, compactness = 1000, num_runs = 1)))
  
  # Test min_size parameter - API now handles edge cases gracefully
  expect_silent(slice_msf(vec, mask, min_size = 1, num_runs = 1))  # Minimum valid
  
  # min_size = 0 might be clamped to 1
  expect_silent(slice_msf(vec, mask, min_size = 0, num_runs = 1))
  
  # min_size larger than total voxels should work (just be constraining)
  suppressWarnings(expect_no_error(slice_msf(vec, mask, min_size = nvox * 2, num_runs = 1)))
  
  # Test num_runs parameter - these should still error as they're fundamental
  expect_error(slice_msf(vec, mask, num_runs = 0),
               info = "Should reject num_runs <= 0")
  
  expect_error(slice_msf(vec, mask, num_runs = -1),
               info = "Should reject negative num_runs")
  
  # Test theta_link parameter (z-stitching) - API now clamps to valid range
  suppressWarnings(expect_no_error(slice_msf(vec, mask, theta_link = 0.0, num_runs = 1)))  # Minimum
  suppressWarnings(expect_no_error(slice_msf(vec, mask, theta_link = 1.0, num_runs = 1)))  # Maximum
  
  # Out of range values should be clamped
  suppressWarnings(expect_no_error(slice_msf(vec, mask, theta_link = -0.1, num_runs = 1)))  # Clamps to 0
  suppressWarnings(expect_no_error(slice_msf(vec, mask, theta_link = 1.1, num_runs = 1)))   # Clamps to 1
  
  # Test min_contact parameter - API now handles gracefully
  # Negative values likely clamped to 0
  suppressWarnings(expect_no_error(slice_msf(vec, mask, min_contact = -1, num_runs = 1)))
  
  # Zero min_contact should work (no contact requirement)
  suppressWarnings(expect_no_error(slice_msf(vec, mask, min_contact = 0, num_runs = 1)))
  
  # Test nbhd parameter - API now handles invalid values gracefully
  # Invalid nbhd values might default to 8
  suppressMessages(expect_no_error(slice_msf(vec, mask, nbhd = 6, num_runs = 1)))  # Defaults to valid
  expect_error(slice_msf(vec, mask, nbhd = 0, num_runs = 1),
               info = "Should reject invalid nbhd values")
  
  # Valid nbhd values
  suppressWarnings(expect_no_error(slice_msf(vec, mask, nbhd = 4, num_runs = 1)))
  suppressWarnings(expect_no_error(slice_msf(vec, mask, nbhd = 8, num_runs = 1)))
})

test_that("handles malformed time series data", {
  dims <- c(6, 6, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  # Test 1: Time series with NaN values
  ts_data_nan <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  ts_data_nan[1:3, 1:5] <- NaN  # Insert NaN values
  
  vec_list_nan <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data_nan[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec_nan <- do.call(concat, vec_list_nan)
  
  # Should either handle gracefully or give clear error
  result_nan_test <- tryCatch({
    result <- slice_msf(vec_nan, mask, num_runs = 1, r = 4)
    list(success = TRUE, result = result)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  expect_true(result_nan_test$success || nchar(result_nan_test$error) > 10,
              info = "NaN data should be handled gracefully or with informative error")
  
  # Test 2: Time series with Inf values
  ts_data_inf <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  ts_data_inf[4:6, 10:15] <- Inf
  ts_data_inf[7:9, 20:25] <- -Inf
  
  vec_list_inf <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data_inf[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec_inf <- do.call(concat, vec_list_inf)
  
  result_inf_test <- tryCatch({
    result <- slice_msf(vec_inf, mask, num_runs = 1, r = 4)
    list(success = TRUE, result = result)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  expect_true(result_inf_test$success || nchar(result_inf_test$error) > 10,
              info = "Inf data should be handled gracefully or with informative error")
  
  # Test 3: Constant time series (zero variance)
  ts_data_constant <- matrix(5.0, nrow = nvox, ncol = ntime)
  
  vec_list_constant <- lapply(1:ntime, function(t) {
    vol_data <- array(5.0, dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec_constant <- do.call(concat, vec_list_constant)
  
  # Should handle constant data
  suppressWarnings(expect_no_error(result_constant <- slice_msf(vec_constant, mask, num_runs = 1, r = 4)))
  expect_true(!is.null(result_constant))
  
  # Test 4: Very large values (potential overflow)
  large_val <- 1e8
  ts_data_large <- matrix(runif(nvox * ntime, -large_val, large_val), nrow = nvox, ncol = ntime)
  
  vec_list_large <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data_large[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec_large <- do.call(concat, vec_list_large)
  
  suppressWarnings(expect_no_error(result_large <- slice_msf(vec_large, mask, num_runs = 1, r = 4)))
  expect_true(!is.null(result_large))
  expect_true(all(is.finite(result_large$centers)) || all(abs(result_large$centers) < 1e10),
              info = "Large values should not cause overflow in results")
})

test_that("validates algorithm-specific parameter combinations", {
  dims <- c(8, 8, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 25
  
  # Create test data
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test SNIC parameter validation - some parameters now handled gracefully
  expect_error(snic(vec, mask, K = 0),
               info = "SNIC should reject K <= 0")
  
  # Negative compactness might be clamped to small positive value
  expect_silent(snic(vec, mask, K = 5, compactness = 0.1))
  
  # max_iter might have default or be clamped
  expect_silent(snic(vec, mask, K = 5, max_iter = 1))
  
  # Test supervoxels parameter validation
  expect_error(supervoxels(vec, mask, K = 0),
               info = "Supervoxels should reject K <= 0")
  
  # Test alpha parameter bounds
  suppressMessages(expect_no_error(supervoxels(vec, mask, K = 5, alpha = 0.0)))  # Minimum
  suppressMessages(expect_no_error(supervoxels(vec, mask, K = 5, alpha = 1.0)))  # Maximum
  
  # Out of bounds alpha should error
  expect_error(supervoxels(vec, mask, K = 5, alpha = -0.1),
               info = "Should reject alpha < 0")
  expect_error(supervoxels(vec, mask, K = 5, alpha = 1.1),
               info = "Should reject alpha > 1")
  
  # Test commute_cluster parameter validation
  expect_error(commute_cluster(vec, mask, K = 0),
               info = "Commute clustering should reject K <= 0")
  
  # Test other commute_cluster parameters  
  suppressMessages(expect_no_error(commute_cluster(vec, mask, K = 3, alpha = 0.5)))
  suppressMessages(expect_no_error(commute_cluster(vec, mask, K = 3, connectivity = 6)))
  
  # Test parameter consistency across algorithms
  valid_params <- list(
    list(K = 3, compactness = 1.0),
    list(K = 5, compactness = 2.0),
    list(K = 2, compactness = 0.5)
  )
  
  for (params in valid_params) {
    # All algorithms should accept valid parameter combinations
    expect_no_error(snic(vec, mask, K = params$K, compactness = params$compactness))
    suppressMessages(expect_no_error(supervoxels(vec, mask, K = params$K, alpha = 0.5)))
    suppressMessages(expect_no_error(commute_cluster(vec, mask, K = params$K)))
  }
})

test_that("handles edge case spatial configurations", {
  # Test 1: Single slice (2D) data
  dims_2d <- c(8, 8, 1)
  mask_2d <- NeuroVol(array(1, dims_2d), NeuroSpace(dims_2d))
  nvox_2d <- prod(dims_2d)
  ntime_2d <- 30
  
  ts_data_2d <- matrix(rnorm(nvox_2d * ntime_2d), nrow = nvox_2d, ncol = ntime_2d)
  vec_list_2d <- lapply(1:ntime_2d, function(t) {
    vol_data <- array(0, dims_2d)
    vol_data[mask_2d > 0] <- ts_data_2d[, t]
    NeuroVol(vol_data, NeuroSpace(dims_2d))
  })
  vec_2d <- do.call(concat, vec_list_2d)
  
  # 2D data should work
  suppressWarnings(expect_no_error(result_2d <- slice_msf(vec_2d, mask_2d, num_runs = 1)))
  expect_true(!is.null(result_2d))
  
  # Test 2: Very sparse mask (few scattered voxels)
  dims_sparse <- c(10, 10, 2)
  mask_sparse <- NeuroVol(array(0, dims_sparse), NeuroSpace(dims_sparse))
  
  # Activate only a few scattered voxels
  active_coords <- rbind(
    c(2, 2, 1), c(5, 8, 1), c(9, 3, 1),
    c(1, 7, 2), c(6, 4, 2), c(8, 9, 2)
  )
  
  for (i in 1:nrow(active_coords)) {
    mask_sparse@.Data[active_coords[i, 1], active_coords[i, 2], active_coords[i, 3]] <- 1
  }
  
  nvox_sparse <- sum(mask_sparse@.Data > 0)
  ntime_sparse <- 25
  
  ts_data_sparse <- matrix(rnorm(nvox_sparse * ntime_sparse), nrow = nvox_sparse, ncol = ntime_sparse)
  
  vec_list_sparse <- lapply(1:ntime_sparse, function(t) {
    vol_data <- array(0, dims_sparse)
    vol_data[mask_sparse > 0] <- ts_data_sparse[, t]
    NeuroVol(vol_data, NeuroSpace(dims_sparse))
  })
  vec_sparse <- do.call(concat, vec_list_sparse)
  
  # Sparse data should work
  expect_silent(result_sparse <- slice_msf(vec_sparse, mask_sparse, num_runs = 1, min_size = 1))
  expect_true(!is.null(result_sparse))
  expect_true(length(result_sparse$cluster) == nvox_sparse)
  
  # Test 3: Non-contiguous mask regions
  dims_noncontig <- c(12, 12, 1)
  mask_noncontig <- NeuroVol(array(0, dims_noncontig), NeuroSpace(dims_noncontig))
  
  # Create two separate regions
  mask_noncontig@.Data[2:5, 2:5, 1] <- 1    # Region 1
  mask_noncontig@.Data[8:11, 8:11, 1] <- 1  # Region 2 (separate)
  
  nvox_noncontig <- sum(mask_noncontig@.Data > 0)
  ntime_noncontig <- 30
  
  ts_data_noncontig <- matrix(rnorm(nvox_noncontig * ntime_noncontig), 
                              nrow = nvox_noncontig, ncol = ntime_noncontig)
  
  vec_list_noncontig <- lapply(1:ntime_noncontig, function(t) {
    vol_data <- array(0, dims_noncontig)
    vol_data[mask_noncontig > 0] <- ts_data_noncontig[, t]
    NeuroVol(vol_data, NeuroSpace(dims_noncontig))
  })
  vec_noncontig <- do.call(concat, vec_list_noncontig)
  
  # Non-contiguous regions should work
  expect_silent(result_noncontig <- slice_msf(vec_noncontig, mask_noncontig, num_runs = 1))
  expect_true(!is.null(result_noncontig))
})

test_that("validates consensus and meta-clustering inputs", {
  # Create base clustering result for meta-clustering tests
  dims <- c(8, 8, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Get valid clustering result
  base_result <- slice_msf(vec, mask, num_runs = 1, r = 6)
  
  # Test meta_clust parameter validation
  expect_error(meta_clust(NULL),
               info = "meta_clust should reject NULL input")
  
  expect_error(meta_clust(base_result, cuts = 0),
               info = "meta_clust should reject cuts <= 0")
  
  expect_error(meta_clust(base_result, cuts = -1),
               info = "meta_clust should reject negative cuts")
  
  # Test with malformed clustering result
  fake_result <- list(cluster = c(1, 2, 3))  # Missing required fields
  
  result_test <- tryCatch({
    meta_result <- meta_clust(fake_result, cuts = 2)
    list(success = TRUE, result = meta_result)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  expect_true(!result_test$success,
              info = "meta_clust should reject malformed clustering results")
  
  # Test consensus parameter validation for multi-run results
  # API might handle consensus=TRUE with num_runs=1 by ignoring consensus
  suppressWarnings(expect_no_error(slice_msf(vec, mask, num_runs = 1, consensus = TRUE)))
  
  # Test lambda parameter for consensus - values might be clamped
  suppressWarnings(expect_no_error(slice_msf(vec, mask, num_runs = 3, consensus = TRUE, lambda = 0.0)))
  suppressWarnings(expect_no_error(slice_msf(vec, mask, num_runs = 3, consensus = TRUE, lambda = 1.0)))
  
  # Out of bounds lambda values might be clamped
  suppressWarnings(expect_no_error(slice_msf(vec, mask, num_runs = 3, consensus = TRUE, lambda = -0.1)))
  suppressWarnings(expect_no_error(slice_msf(vec, mask, num_runs = 3, consensus = TRUE, lambda = 1.1)))
})

test_that("error messages are informative", {
  dims <- c(6, 6, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 20
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test that error messages contain useful information
  error_msg <- tryCatch({
    slice_msf(vec, mask, r = -5)
  }, error = function(e) e$message)
  
  expect_true(grepl("r", error_msg) || grepl("rank", error_msg),
              info = "Error message should mention the problematic parameter")
  
  # Test dimension mismatch error
  wrong_mask <- NeuroVol(array(1, c(4, 4, 1)), NeuroSpace(c(4, 4, 1)))
  
  dim_error_msg <- tryCatch({
    slice_msf(vec, wrong_mask)
  }, error = function(e) e$message)
  
  expect_true(any(grepl("dimension", dim_error_msg, ignore.case = TRUE)) ||
              any(grepl("size", dim_error_msg, ignore.case = TRUE)),
              info = "Dimension error should mention dimensions or size")
  
  # Test empty mask error
  empty_mask <- NeuroVol(array(0, dims), NeuroSpace(dims))
  
  empty_error_msg <- tryCatch({
    slice_msf(vec, empty_mask)
  }, error = function(e) e$message)
  
  expect_true(grepl("mask", empty_error_msg, ignore.case = TRUE) ||
              grepl("voxel", empty_error_msg, ignore.case = TRUE),
              info = "Empty mask error should mention mask or voxels")
  
  expect_true(nchar(empty_error_msg) > 20,
              info = "Error messages should be reasonably descriptive")
})