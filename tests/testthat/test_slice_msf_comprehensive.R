library(testthat)
library(neurocluster)
library(neuroim2)

# Helper function to create test data with controllable structure
create_slice_test_data <- function(dims = c(10, 10, 5), 
                                  n_time = 50, 
                                  n_patterns = 3,
                                  noise_sd = 0.2,
                                  seed = 123) {
  set.seed(seed)
  
  # Create mask
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  # Generate different temporal patterns
  t_seq <- seq(0, 4*pi, length.out = n_time)
  patterns <- list(
    sin(t_seq),
    cos(t_seq),
    sin(2 * t_seq),
    cos(2 * t_seq),
    sin(0.5 * t_seq)
  )
  
  # Take only requested number of patterns
  patterns <- patterns[seq_len(min(n_patterns, length(patterns)))]
  
  # Assign patterns to spatial regions
  ts_data <- matrix(0, nrow = nvox, ncol = n_time)
  voxels_per_pattern <- nvox / n_patterns
  
  for (i in seq_along(patterns)) {
    start_idx <- ceiling((i - 1) * voxels_per_pattern) + 1
    end_idx <- min(ceiling(i * voxels_per_pattern), nvox)
    n_vox <- end_idx - start_idx + 1
    
    ts_data[start_idx:end_idx, ] <- matrix(rep(patterns[[i]], n_vox), 
                                           nrow = n_vox, byrow = TRUE) + 
                                    matrix(rnorm(n_vox * n_time, sd = noise_sd), 
                                          nrow = n_vox)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(ts_data[, t], dim = dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  list(vec = vec, mask = mask, nvox = nvox, n_time = n_time)
}

test_that("slice_msf handles multiple runs with various parameter combinations", {
  # Small dataset for fast testing
  data <- create_slice_test_data(dims = c(8, 8, 4), n_time = 30, seed = 42)
  
  # Test matrix of parameter combinations
  test_params <- expand.grid(
    num_runs = c(2, 3),
    consensus = c(TRUE, FALSE),
    r = c(4, 8),
    compactness = c(2, 5),
    stitch_z = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  
  # Test each combination
  for (i in 1:nrow(test_params)) {
    params <- test_params[i, ]
    
    # Skip non-sensical combinations
    if (params$num_runs == 1 && params$consensus) next
    
    result <- tryCatch({
      slice_msf(
        vec = data$vec,
        mask = data$mask,
        num_runs = params$num_runs,
        consensus = params$consensus,
        r = params$r,
        compactness = params$compactness,
        stitch_z = params$stitch_z,
        min_size = 10,
        theta_link = 0.85,
        min_contact = 1
      )
    }, error = function(e) {
      fail(sprintf("slice_msf failed with params: num_runs=%d, consensus=%s, r=%d, compactness=%d, stitch_z=%s\nError: %s",
                  params$num_runs, params$consensus, params$r, params$compactness, 
                  params$stitch_z, e$message))
      NULL
    })
    
    if (!is.null(result)) {
      expect_s3_class(result, "slice_msf_cluster_result")
      expect_equal(length(result$cluster), data$nvox)
      expect_true(all(result$cluster > 0), 
                  info = sprintf("Failed for params row %d", i))
      
      # Check runs are stored when appropriate
      if (params$num_runs > 1 && params$consensus) {
        expect_true("runs" %in% names(result))
        expect_equal(length(result$runs), params$num_runs)
      }
    }
  }
})

test_that("slice_msf works with medium-sized datasets", {
  # Medium dataset
  data <- create_slice_test_data(dims = c(20, 20, 10), n_time = 100, seed = 123)
  
  # Test with multiple runs
  result <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    num_runs = 3,
    consensus = TRUE,
    r = 12,
    compactness = 5,
    min_size = 50,
    use_features = TRUE,
    lambda = 0.7
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), data$nvox)
  expect_true(all(result$cluster > 0))
  expect_true("runs" %in% names(result))
  expect_equal(length(result$runs), 3)
  
  # Check cluster properties
  n_clusters <- length(unique(result$cluster))
  expect_true(n_clusters >= 2)
  expect_true(n_clusters <= 100)
})

test_that("slice_msf handles target_k parameters correctly", {
  data <- create_slice_test_data(dims = c(10, 10, 5), n_time = 40, seed = 456)
  
  # Test target_k_global (requires use_features=TRUE when using consensus)
  result_global <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    target_k_global = 10,
    num_runs = 2,
    consensus = TRUE,
    use_features = TRUE,  # Required for exact-K with consensus
    r = 8,
    min_size = 20
  )
  
  n_clusters_global <- length(unique(result_global$cluster))
  # Should be close to target (allowing some flexibility due to algorithm constraints)
  expect_true(abs(n_clusters_global - 10) <= 5)
  
  # Test target_k_per_slice (requires stitch_z = FALSE)
  result_per_slice <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    target_k_per_slice = 5,
    stitch_z = FALSE,
    num_runs = 1,
    r = 6,
    min_size = 10
  )
  
  expect_s3_class(result_per_slice, "slice_msf_cluster_result")
  expect_equal(length(result_per_slice$cluster), data$nvox)
})

test_that("slice_msf consensus fusion parameters work correctly", {
  data <- create_slice_test_data(dims = c(8, 8, 4), n_time = 30, seed = 789)
  
  # Test use_features parameter
  result_no_features <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    num_runs = 3,
    consensus = TRUE,
    use_features = FALSE,
    lambda = 0.7,
    r = 6,
    min_size = 15
  )
  
  result_with_features <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    num_runs = 3,
    consensus = TRUE,
    use_features = TRUE,
    lambda = 0.5,
    r = 6,
    min_size = 15
  )
  
  # Both should produce valid results
  expect_s3_class(result_no_features, "slice_msf_cluster_result")
  expect_s3_class(result_with_features, "slice_msf_cluster_result")
  
  # Test different lambda values
  lambdas <- c(0.3, 0.5, 0.7, 0.9)
  for (lam in lambdas) {
    result <- slice_msf(
      vec = data$vec,
      mask = data$mask,
      num_runs = 2,
      consensus = TRUE,
      use_features = TRUE,
      lambda = lam,
      r = 4,
      min_size = 10
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), data$nvox)
  }
})

test_that("slice_msf handles edge cases gracefully", {
  # Very small mask with limited voxels
  small_mask <- NeuroVol(array(0, c(5, 5, 3)), NeuroSpace(c(5, 5, 3)))
  small_mask[2:3, 2:3, 1:2] <- 1  # Only 8 voxels
  
  nvox <- sum(small_mask > 0)
  ntime <- 20
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(5, 5, 3))
    vol_data[small_mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(5, 5, 3)))
  })
  vec <- do.call(concat, vec_list)
  
  # Should handle very small data
  result <- slice_msf(
    vec = vec,
    mask = small_mask,
    num_runs = 2,
    consensus = TRUE,
    r = 2,
    min_size = 2
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
  
  # Test with single slice
  single_slice_mask <- NeuroVol(array(0, c(10, 10, 1)), NeuroSpace(c(10, 10, 1)))
  single_slice_mask[3:8, 3:8, 1] <- 1
  
  nvox_single <- sum(single_slice_mask > 0)
  ntime <- 25
  ts_data_single <- matrix(rnorm(nvox_single * ntime), nrow = nvox_single, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(10, 10, 1))
    vol_data[single_slice_mask > 0] <- ts_data_single[, t]
    NeuroVol(vol_data, NeuroSpace(c(10, 10, 1)))
  })
  vec_single <- do.call(concat, vec_list)
  
  result_single <- slice_msf(
    vec = vec_single,
    mask = single_slice_mask,
    num_runs = 2,
    consensus = TRUE,
    stitch_z = FALSE,  # No z-stitching for single slice
    r = 4,
    min_size = 5
  )
  
  expect_s3_class(result_single, "slice_msf_cluster_result")
  expect_equal(length(result_single$cluster), nvox_single)
})

test_that("slice_msf reliability weighting works", {
  data <- create_slice_test_data(dims = c(8, 8, 4), n_time = 40, seed = 111)
  
  # Test different gamma values
  gammas <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  
  for (g in gammas) {
    result <- slice_msf(
      vec = data$vec,
      mask = data$mask,
      gamma = g,
      num_runs = 2,
      consensus = TRUE,
      r = 6,
      min_size = 15
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), data$nvox)
  }
})

test_that("slice_msf z-stitching parameters work correctly", {
  data <- create_slice_test_data(dims = c(10, 10, 6), n_time = 40, seed = 222)
  
  # Test different theta_link values
  theta_links <- c(0.7, 0.75, 0.85, 0.9, 0.95)
  
  for (theta in theta_links) {
    result <- slice_msf(
      vec = data$vec,
      mask = data$mask,
      stitch_z = TRUE,
      theta_link = theta,
      min_contact = 1,
      num_runs = 2,
      consensus = TRUE,
      r = 6,
      min_size = 20
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), data$nvox)
  }
  
  # Test different min_contact values
  min_contacts <- c(1, 2, 3, 5)
  
  for (mc in min_contacts) {
    result <- slice_msf(
      vec = data$vec,
      mask = data$mask,
      stitch_z = TRUE,
      theta_link = 0.85,
      min_contact = mc,
      num_runs = 1,
      r = 6,
      min_size = 20
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), data$nvox)
  }
})

test_that("slice_msf produces consistent results across runs", {
  data <- create_slice_test_data(dims = c(10, 10, 5), n_time = 50, seed = 333)
  
  # Run twice with same parameters
  result1 <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    num_runs = 3,
    consensus = TRUE,
    r = 8,
    compactness = 5,
    min_size = 25,
    use_features = TRUE,
    lambda = 0.6
  )
  
  result2 <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    num_runs = 3,
    consensus = TRUE,
    r = 8,
    compactness = 5,
    min_size = 25,
    use_features = TRUE,
    lambda = 0.6
  )
  
  # Check both produce valid results
  expect_s3_class(result1, "slice_msf_cluster_result")
  expect_s3_class(result2, "slice_msf_cluster_result")
  
  # Number of clusters should be similar (allowing for some stochasticity)
  n_clusters1 <- length(unique(result1$cluster))
  n_clusters2 <- length(unique(result2$cluster))
  expect_true(abs(n_clusters1 - n_clusters2) <= 5)
})

test_that("slice_msf neighborhood connectivity options work", {
  data <- create_slice_test_data(dims = c(8, 8, 4), n_time = 30, seed = 444)
  
  # Test different neighborhood options
  nbhds <- c(4, 6, 8)
  
  for (nb in nbhds) {
    result <- slice_msf(
      vec = data$vec,
      mask = data$mask,
      nbhd = nb,
      num_runs = 2,
      consensus = TRUE,
      r = 5,
      min_size = 15
    )
    
    expect_s3_class(result, "slice_msf_cluster_result")
    expect_equal(length(result$cluster), data$nvox)
  }
})

test_that("slice_msf handles sparse masks correctly", {
  # Create sparse mask (checkerboard pattern)
  sparse_mask <- NeuroVol(array(0, c(10, 10, 4)), NeuroSpace(c(10, 10, 4)))
  for (i in 1:10) {
    for (j in 1:10) {
      for (k in 1:4) {
        if ((i + j + k) %% 2 == 0) {
          sparse_mask[i, j, k] <- 1
        }
      }
    }
  }
  
  nvox <- sum(sparse_mask > 0)
  ntime <- 30
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(10, 10, 4))
    vol_data[sparse_mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(10, 10, 4)))
  })
  vec <- do.call(concat, vec_list)
  
  result <- slice_msf(
    vec = vec,
    mask = sparse_mask,
    num_runs = 2,
    consensus = TRUE,
    r = 6,
    min_size = 10
  )
  
  expect_s3_class(result, "slice_msf_cluster_result")
  expect_equal(length(result$cluster), nvox)
  expect_true(all(result$cluster > 0))
})

test_that("z-axis smoothing increases vertical feature agreement", {
  dims <- c(2, 2, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  T <- 12
  t_seq <- seq(0, 2 * pi, length.out = T)
  base_pattern <- sin(t_seq)
  alt_pattern <- cos(t_seq)

  TS <- matrix(0, nrow = T, ncol = nvox)
  idx <- 1
  for (z in 1:dims[3]) {
    for (y in 1:dims[2]) {
      for (x in 1:dims[1]) {
        TS[, idx] <- if (z == 2) alt_pattern else base_pattern
        idx <- idx + 1
      }
    }
  }

  mask_flat <- as.integer(mask@.Data)
  vol_dim <- dims
  voxel_dim <- spacing(mask)

  base <- slice_msf_runwise(
    TS = TS,
    mask = mask_flat,
    vol_dim = vol_dim,
    r = 4,
    fh_scale = 0.3,
    min_size = 1,
    nbhd = 4,
    stitch_z = TRUE,
    theta_link = 0.9,
    min_contact = 1,
    rows_are_time = TRUE,
    gamma = 1.0,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0,
    target_k_global = -1,
    target_k_per_slice = -1,
    z_mult = 0.0
  )

  smoothed <- slice_msf_runwise(
    TS = TS,
    mask = mask_flat,
    vol_dim = vol_dim,
    r = 4,
    fh_scale = 0.3,
    min_size = 1,
    nbhd = 4,
    stitch_z = TRUE,
    theta_link = 0.9,
    min_contact = 1,
    rows_are_time = TRUE,
    gamma = 1.0,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0,
    target_k_global = -1,
    target_k_per_slice = -1,
    z_mult = 0.5
  )

  idx3d <- function(x, y, z) {
    (x - 1) + dims[1] * ((y - 1) + dims[2] * (z - 1)) + 1
  }

  vertical_mean <- function(U) {
    dots <- numeric()
    for (x in 1:dims[1]) {
      for (y in 1:dims[2]) {
        for (z in 1:(dims[3] - 1)) {
          i <- idx3d(x, y, z)
          j <- idx3d(x, y, z + 1)
          dots <- c(dots, sum(U[, i] * U[, j]))
        }
      }
    }
    mean(dots)
  }

  base_mean <- vertical_mean(base$sketch)
  smooth_mean <- vertical_mean(smoothed$sketch)

  expect_gt(smooth_mean, base_mean + 0.05)
})
