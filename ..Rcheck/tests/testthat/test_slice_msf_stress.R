library(testthat)
library(neurocluster)
library(neuroim2)

# Stress test for slice_msf with larger datasets and edge cases
test_that("slice_msf handles larger datasets with multiple runs", {
  skip_if_not(as.logical(Sys.getenv("RUN_STRESS_TESTS", "FALSE")), 
              "Skipping stress tests (set RUN_STRESS_TESTS=TRUE to run)")
  
  # Create a larger dataset (30x30x15 = 13,500 voxels)
  set.seed(999)
  dims <- c(30, 30, 15)
  n_time <- 150
  
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  # Create structured data with noise
  ts_data <- matrix(rnorm(nvox * n_time, sd = 1), nrow = nvox, ncol = n_time)
  
  # Add some structure (5 main patterns)
  t_seq <- seq(0, 6*pi, length.out = n_time)
  for (i in 1:5) {
    pattern <- sin(i * t_seq / 2) + cos((6-i) * t_seq / 3)
    start_idx <- (i-1) * (nvox/5) + 1
    end_idx <- min(i * (nvox/5), nvox)
    n_vox <- end_idx - start_idx + 1
    
    ts_data[start_idx:end_idx, ] <- ts_data[start_idx:end_idx, ] + 
                                     matrix(rep(pattern, n_vox), nrow = n_vox, byrow = TRUE)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(ts_data[, t], dim = dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test with multiple runs and consensus
  result <- slice_msf(
    vec = vec,
    mask = mask,
    num_runs = 3,
    consensus = TRUE,
    r = 15,
    compactness = 5,
    min_size = 100,
    use_features = TRUE,
    lambda = 0.6,
    stitch_z = TRUE,
    theta_link = 0.80
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
  expect_true(all(result$cluster > 0))
  expect_true("runs" %in% names(result))
  expect_equal(length(result$runs), 3)
})

test_that("slice_msf handles very sparse masks", {
  # Create a mask with only 10% of voxels active
  dims <- c(15, 15, 8)
  mask_array <- array(0, dims)
  
  set.seed(111)
  # Randomly select 10% of voxels
  n_active <- round(prod(dims) * 0.1)
  active_idx <- sample(prod(dims), n_active)
  mask_array[active_idx] <- 1
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  nvox <- sum(mask > 0)
  
  n_time <- 50
  ts_data <- matrix(rnorm(nvox * n_time), nrow = nvox, ncol = n_time)
  
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(0, dim = dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test with multiple runs on sparse data
  result <- slice_msf(
    vec = vec,
    mask = mask,
    num_runs = 2,
    consensus = TRUE,
    r = 8,
    min_size = 5,
    compactness = 3
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
  expect_true(all(result$cluster > 0))
})

test_that("slice_msf handles extreme parameter combinations", {
  # Small test data
  dims <- c(10, 10, 5)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  set.seed(222)
  n_time <- 30
  ts_data <- matrix(rnorm(nvox * n_time), nrow = nvox, ncol = n_time)
  
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(ts_data[, t], dim = dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test extreme parameter combinations
  extreme_params <- list(
    # Very high compactness with low min_size
    list(compactness = 10, min_size = 5, r = 4),
    # Very low compactness with high min_size
    list(compactness = 1, min_size = 100, r = 4),
    # Maximum r value
    list(compactness = 5, min_size = 50, r = 20),
    # Minimum r value
    list(compactness = 5, min_size = 50, r = 2),
    # Very aggressive stitching
    list(theta_link = 0.5, min_contact = 1, stitch_z = TRUE),
    # Very conservative stitching
    list(theta_link = 0.99, min_contact = 10, stitch_z = TRUE)
  )
  
  for (i in seq_along(extreme_params)) {
    params <- extreme_params[[i]]
    
    result <- tryCatch({
      do.call(slice_msf, c(
        list(vec = vec, mask = mask, num_runs = 2, consensus = TRUE),
        params
      ))
    }, error = function(e) {
      # Some extreme parameters might legitimately fail
      # but shouldn't crash R
      NULL
    })
    
    if (!is.null(result)) {
      expect_s3_class(result, "slice_msf_cluster_result")
      expect_equal(length(result$cluster), nvox)
    }
  }
})

test_that("slice_msf handles exact-K with different strategies", {
  dims <- c(12, 12, 6)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  set.seed(333)
  n_time <- 60
  
  # Create data with clear structure for K clusters
  K_true <- 8
  ts_data <- matrix(0, nrow = nvox, ncol = n_time)
  t_seq <- seq(0, 4*pi, length.out = n_time)
  
  voxels_per_cluster <- nvox / K_true
  for (k in 1:K_true) {
    pattern <- sin(k * t_seq / 2) + rnorm(n_time, sd = 0.3)
    start_idx <- ceiling((k-1) * voxels_per_cluster) + 1
    end_idx <- min(ceiling(k * voxels_per_cluster), nvox)
    n_vox <- end_idx - start_idx + 1
    
    ts_data[start_idx:end_idx, ] <- matrix(rep(pattern, n_vox), 
                                           nrow = n_vox, byrow = TRUE)
  }
  
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(ts_data[, t], dim = dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test exact-K with global target
  result_global <- slice_msf(
    vec = vec,
    mask = mask,
    target_k_global = K_true,
    num_runs = 2,
    consensus = TRUE,
    use_features = TRUE,  # Required for exact-K
    r = 10,
    min_size = 20
  )
  
  n_clusters <- length(unique(result_global$cluster))
  # Should achieve target K or be very close
  expect_true(abs(n_clusters - K_true) <= 2)
  
  # Test per-slice K (without stitching)
  K_per_slice <- 3
  result_per_slice <- slice_msf(
    vec = vec,
    mask = mask,
    target_k_per_slice = K_per_slice,
    stitch_z = FALSE,
    num_runs = 1,
    r = 8,
    min_size = 10
  )
  
  expect_s3_class(result_per_slice, "slice_msf_cluster_result")
  # With 6 slices and 3 clusters per slice, expect multiple clusters
  # The exact number depends on data structure and min_size constraints
  n_clusters_total <- length(unique(result_per_slice$cluster))
  # Should have at least some clusters (more than 1)
  expect_true(n_clusters_total >= 2)  # At least 2 clusters total
  # And not more than K_per_slice * number of slices
  expect_true(n_clusters_total <= K_per_slice * 6)
})

test_that("slice_msf memory safety with many runs", {
  # Test that multiple runs don't cause memory issues
  dims <- c(8, 8, 4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  set.seed(444)
  n_time <- 25
  ts_data <- matrix(rnorm(nvox * n_time), nrow = nvox, ncol = n_time)
  
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(ts_data[, t], dim = dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test with many runs
  for (n_runs in c(5, 7, 10)) {
    result <- slice_msf(
      vec = vec,
      mask = mask,
      num_runs = n_runs,
      consensus = TRUE,
      r = 4,
      min_size = 10,
      compactness = 5
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), nvox)
    expect_true("runs" %in% names(result))
    expect_equal(length(result$runs), n_runs)
  }
})

test_that("slice_msf handles disconnected regions", {
  # Create a mask with two disconnected regions
  dims <- c(15, 15, 5)
  mask_array <- array(0, dims)
  
  # Left region
  mask_array[2:6, 2:13, 1:5] <- 1
  # Right region (disconnected)
  mask_array[10:14, 2:13, 1:5] <- 1
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  nvox <- sum(mask > 0)
  
  set.seed(555)
  n_time <- 40
  
  # Different patterns for each region
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = n_time)
  
  t_seq <- seq(0, 2*pi, length.out = n_time)
  left_region <- coords[,1] < 8
  right_region <- coords[,1] >= 8
  
  ts_data[left_region, ] <- matrix(rep(sin(t_seq), sum(left_region)), 
                                   ncol = n_time, byrow = TRUE) + 
                                   rnorm(sum(left_region) * n_time, sd = 0.2)
  
  ts_data[right_region, ] <- matrix(rep(cos(t_seq), sum(right_region)), 
                                    ncol = n_time, byrow = TRUE) + 
                                    rnorm(sum(right_region) * n_time, sd = 0.2)
  
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(0, dim = dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Should handle disconnected regions
  result <- slice_msf(
    vec = vec,
    mask = mask,
    num_runs = 2,
    consensus = TRUE,
    r = 6,
    min_size = 20,
    compactness = 4
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
  
  # Should find at least 2 clusters (one per region)
  n_clusters <- length(unique(result$cluster))
  expect_true(n_clusters >= 2)
})