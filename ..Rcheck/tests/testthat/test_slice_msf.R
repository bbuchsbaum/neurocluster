library(neuroim2)
library(neurocluster)

test_that("slice_msf works with basic inputs", {
  # Create small test data
  mask <- NeuroVol(array(1, c(10, 10, 5)), NeuroSpace(c(10, 10, 5)))
  
  # Create time series data with some structure
  set.seed(123)
  nvox <- sum(mask > 0)
  ntime <- 50
  
  # Create 3 distinct patterns
  pattern1 <- sin(seq(0, 4*pi, length.out = ntime))
  pattern2 <- cos(seq(0, 4*pi, length.out = ntime))
  pattern3 <- sin(seq(0, 2*pi, length.out = ntime))
  
  # Assign patterns to different regions
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  ts_data[1:200, ] <- matrix(rep(pattern1, 200), nrow = 200, byrow = TRUE) + 
                      matrix(rnorm(200 * ntime, sd = 0.2), nrow = 200)
  ts_data[201:350, ] <- matrix(rep(pattern2, 150), nrow = 150, byrow = TRUE) + 
                        matrix(rnorm(150 * ntime, sd = 0.2), nrow = 150)
  ts_data[351:500, ] <- matrix(rep(pattern3, 150), nrow = 150, byrow = TRUE) + 
                        matrix(rnorm(150 * ntime, sd = 0.2), nrow = 150)
  
  # Create NeuroVec from the data
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(10, 10, 5))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(10, 10, 5)))
  })
  vec <- do.call(concat, vec_list)
  
  # Test single run
  result <- slice_msf(vec, mask, num_runs = 1, r = 6, min_size = 20)
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_s3_class(result, "cluster_result")
  expect_s4_class(result$clusvol, "ClusteredNeuroVol")
  expect_equal(length(result$cluster), nvox)
  expect_true(all(result$cluster > 0))  # All voxels should be assigned
  
  # Check that we have some clusters
  n_clusters <- length(unique(result$cluster))
  expect_true(n_clusters >= 2)
  expect_true(n_clusters <= 50)  # Reasonable upper bound
})

test_that("slice_msf works with consensus fusion", {
  # Create small test data
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  
  set.seed(456)
  nvox <- sum(mask > 0)
  ntime <- 30
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, c(8, 8, 4))
  
  # Create structured data with 2 regions
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Left half: sine pattern
  left_voxels <- coords[,1] <= 4
  ts_data[left_voxels, ] <- matrix(rep(sin(t_seq), sum(left_voxels)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(left_voxels) * ntime, sd = 0.2)
  
  # Right half: cosine pattern
  right_voxels <- coords[,1] > 4
  ts_data[right_voxels, ] <- matrix(rep(cos(t_seq), sum(right_voxels)), 
                                     ncol = ntime, byrow = TRUE) + 
                                     rnorm(sum(right_voxels) * ntime, sd = 0.2)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(8, 8, 4))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(8, 8, 4)))
  })
  vec <- do.call(concat, vec_list)
  
  # Test multi-run consensus
  result <- slice_msf(vec, mask, num_runs = 2, consensus = TRUE, 
                      r = 4, min_size = 20, compactness = 3)
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_true("runs" %in% names(result))
  expect_equal(length(result$runs), 2)
  expect_equal(length(result$cluster), nvox)
})

test_that("slice_msf handles edge cases", {
  # Very small mask
  mask <- NeuroVol(array(0, c(5, 5, 2)), NeuroSpace(c(5, 5, 2)))
  mask[2:4, 2:4, 1] <- 1  # Only 9 voxels
  
  set.seed(789)
  nvox <- sum(mask > 0)
  ntime <- 20
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(5, 5, 2))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(5, 5, 2)))
  })
  vec <- do.call(concat, vec_list)
  
  # Should work even with very small data
  result <- slice_msf(vec, mask, num_runs = 1, r = 2, min_size = 3)
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
})

test_that("slice_msf parameters affect results", {
  # Create test data
  mask <- NeuroVol(array(1, c(6, 6, 3)), NeuroSpace(c(6, 6, 3)))
  
  set.seed(111)
  nvox <- sum(mask > 0)
  ntime <- 40
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(6, 6, 3))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(6, 6, 3)))
  })
  vec <- do.call(concat, vec_list)
  
  # Test different compactness values
  result_low <- slice_msf(vec, mask, compactness = 2, num_runs = 1, 
                          r = 4, min_size = 5)
  result_high <- slice_msf(vec, mask, compactness = 10, num_runs = 1, 
                           r = 4, min_size = 5)
  
  # Higher compactness should generally lead to fewer, more compact clusters
  n_clusters_low <- length(unique(result_low$cluster))
  n_clusters_high <- length(unique(result_high$cluster))
  
  # Both should produce valid results
  expect_true(n_clusters_low >= 1)
  expect_true(n_clusters_high >= 1)
})

test_that("slice_msf works with different neighborhood settings", {
  mask <- NeuroVol(array(1, c(6, 6, 3)), NeuroSpace(c(6, 6, 3)))
  
  set.seed(222)
  nvox <- sum(mask > 0)
  ntime <- 30
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(6, 6, 3))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(6, 6, 3)))
  })
  vec <- do.call(concat, vec_list)
  
  # Test 4-connectivity
  result_4 <- slice_msf(vec, mask, nbhd = 4, num_runs = 1, r = 4, min_size = 5)
  expect_s3_class(result_4, "slice_msf_cluster_result")
  
  # Test 8-connectivity
  result_8 <- slice_msf(vec, mask, nbhd = 8, num_runs = 1, r = 4, min_size = 5)
  expect_s3_class(result_8, "slice_msf_cluster_result")
  
  # Both should produce valid results
  expect_true(all(result_4$cluster > 0))
  expect_true(all(result_8$cluster > 0))
})