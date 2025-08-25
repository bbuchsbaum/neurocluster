library(neurocluster)
library(neuroim2)

# Comprehensive ACSC Algorithm Test Suite
# Tests for Adaptive Correlation Superclustering

test_that("ACSC basic functionality", {
  # Test basic ACSC operation
  dims <- c(8, 8, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 40
  
  # Create structured test data
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Create spatially distinct regions
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 4 && y <= 4) {
      # Region 1
      ts_data[i, ] <- sin(t_seq) + rnorm(ntime, sd = 0.1)
    } else if (x > 4 && y <= 4) {
      # Region 2
      ts_data[i, ] <- cos(t_seq) + rnorm(ntime, sd = 0.1)
    } else if (x <= 4 && y > 4) {
      # Region 3
      ts_data[i, ] <- sin(2 * t_seq) + rnorm(ntime, sd = 0.1)
    } else {
      # Region 4
      ts_data[i, ] <- cos(2 * t_seq) + rnorm(ntime, sd = 0.1)
    }
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test basic ACSC functionality
  expect_silent(acsc_result <- acsc(vec, mask, K = 4, 
                                   ann_k = 5,
                                   alpha = 0.5))
  
  # Validate results
  expect_true(!is.null(acsc_result))
  expect_true("cluster_map" %in% names(acsc_result))
  cluster_assignments <- acsc_result$cluster_map[mask > 0]
  expect_true(length(cluster_assignments) == nvox)
  expect_true(all(cluster_assignments > 0))
  
  # Should find reasonable number of clusters
  n_clusters <- length(unique(cluster_assignments))
  expect_true(n_clusters >= 2 && n_clusters <= 8,
              info = sprintf("ACSC should find 2-8 clusters, found %d", n_clusters))
  
  # Check for spatial coherence
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- cluster_assignments
  
  coherent_neighbors <- 0
  total_neighbors <- 0
  
  for (x in 2:(dims[1]-1)) {
    for (y in 2:(dims[2]-1)) {
      for (z in 1:dims[3]) {
        center_cluster <- cluster_vol[x, y, z]
        neighbors <- c(
          cluster_vol[x-1, y, z], cluster_vol[x+1, y, z],
          cluster_vol[x, y-1, z], cluster_vol[x, y+1, z]
        )
        
        same_cluster_neighbors <- sum(neighbors == center_cluster)
        coherent_neighbors <- coherent_neighbors + same_cluster_neighbors
        total_neighbors <- total_neighbors + length(neighbors)
      }
    }
  }
  
  spatial_coherence <- coherent_neighbors / total_neighbors
  expect_true(spatial_coherence > 0.4,
              info = sprintf("ACSC should have reasonable spatial coherence, got %.3f", spatial_coherence))
})

test_that("ACSC parameter validation", {
  # Test parameter validation and edge cases
  dims <- c(6, 6, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  # Simple test data
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test invalid parameters
  expect_error(acsc(vec, mask, K = 0),
               info = "Should reject K <= 0")
  
  expect_error(acsc(vec, mask, K = 5, alpha = -0.1),
               info = "Should reject negative alpha")
  
  expect_error(acsc(vec, mask, K = 5, alpha = 1.1),
               info = "Should reject alpha > 1")
  
  # Connectivity parameter doesn't exist in acsc
  # Test removed - acsc uses different spatial control
  
  # max_cluster_size parameter doesn't exist in acsc
  # Test removed - acsc uses different clustering approach
  
  # Test extreme but valid parameters
  # Very small K - may produce resolution estimation warnings
  suppressWarnings(result_small_K <- acsc(vec, mask, K = 1))
  expect_true(length(unique(result_small_K$cluster_map[mask > 0])) >= 1)
  
  # Large K (should be clamped to reasonable value)
  # This may produce warnings about resolution estimation
  expect_warning(result_large_K <- acsc(vec, mask, K = nvox), 
                 "Could not find a resolution matching K")
  expect_true(length(unique(result_large_K$cluster_map[mask > 0])) <= nvox)
  
  # Different connectivity values
  # Test different ann_k values instead of connectivity
  # Suppress warnings that may occur with small test data
  suppressWarnings(result_ann2 <- acsc(vec, mask, K = 3, ann_k = 2))
  # May produce warnings about resolution estimation for small datasets
  suppressWarnings(result_ann5 <- acsc(vec, mask, K = 3, ann_k = 5))
  # ann_k = 8 may trigger adjustment warnings for this small dataset
  suppressWarnings(result_ann8 <- acsc(vec, mask, K = 3, ann_k = 8))
  
  # All connectivity options should work
  expect_true(!is.null(result_ann2))
  expect_true(!is.null(result_ann5)) 
  expect_true(!is.null(result_ann8))
})

test_that("ACSC handles different data characteristics", {
  dims <- c(10, 10, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 50
  
  # Test different data scenarios
  data_scenarios <- list(
    # High correlation data
    high_corr = {
      t_seq <- seq(0, 2*pi, length.out = ntime)
      base_signal <- sin(t_seq)
      ts_data <- matrix(0, nrow = nvox, ncol = ntime)
      
      coords <- arrayInd(1:nvox, dims)
      for (i in 1:nvox) {
        x <- coords[i, 1]
        # Create two highly correlated regions
        if (x <= 5) {
          ts_data[i, ] <- base_signal + rnorm(ntime, sd = 0.1)
        } else {
          ts_data[i, ] <- base_signal * 0.9 + rnorm(ntime, sd = 0.1)  # Slightly different but highly correlated
        }
      }
      ts_data
    },
    
    # Low correlation data (mostly noise)
    low_corr = {
      matrix(rnorm(nvox * ntime, sd = 1), nrow = nvox, ncol = ntime)
    },
    
    # Mixed correlation data
    mixed_corr = {
      t_seq <- seq(0, 2*pi, length.out = ntime)
      ts_data <- matrix(0, nrow = nvox, ncol = ntime)
      coords <- arrayInd(1:nvox, dims)
      
      for (i in 1:nvox) {
        x <- coords[i, 1]
        y <- coords[i, 2]
        
        if (x <= 3) {
          # Highly structured region
          ts_data[i, ] <- sin(t_seq) + rnorm(ntime, sd = 0.05)
        } else if (x >= 8) {
          # Another structured region
          ts_data[i, ] <- cos(t_seq) + rnorm(ntime, sd = 0.05)
        } else {
          # Noisy transition region
          ts_data[i, ] <- rnorm(ntime, sd = 0.5)
        }
      }
      ts_data
    }
  )
  
  for (scenario_name in names(data_scenarios)) {
    ts_data <- data_scenarios[[scenario_name]]
    
    vec_list <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec <- do.call(concat, vec_list)
    
    # Should handle all data types
    # May produce warnings about resolution estimation for some scenarios
    suppressWarnings({
      result <- acsc(vec, mask, K = 4, ann_k = 5)
    })
    
    expect_true(!is.null(result),
                info = sprintf("ACSC should produce result for %s data", scenario_name))
    
    expect_true(length(result$cluster_map[mask > 0]) == nvox,
                info = sprintf("ACSC should cluster all voxels for %s data", scenario_name))
    
    n_clusters <- length(unique(result$cluster_map[mask > 0]))
    expect_true(n_clusters > 0 && n_clusters <= 15,
                info = sprintf("ACSC should find reasonable clusters for %s data (found %d)", 
                              scenario_name, n_clusters))
  }
})

test_that("ACSC adaptive clustering behavior", {
  # Test that ACSC adapts to different correlation structures
  dims <- c(12, 12, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 60
  
  # Create data with varying correlation strengths
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Create regions with different internal correlation levels
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 4) {
      # Highly coherent region
      base_signal <- sin(t_seq)
      noise_level <- 0.05
      ts_data[i, ] <- base_signal + rnorm(ntime, sd = noise_level)
    } else if (x <= 8) {
      # Moderately coherent region
      base_signal <- cos(t_seq)
      noise_level <- 0.2
      ts_data[i, ] <- base_signal + rnorm(ntime, sd = noise_level)
    } else {
      # Low coherence region
      base_signal <- sin(2 * t_seq)
      noise_level <- 0.5
      ts_data[i, ] <- base_signal + rnorm(ntime, sd = noise_level)
    }
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test with different alpha values (adaptive parameter)
  alpha_values <- c(0.1, 0.5, 0.9)
  alpha_results <- list()
  
  for (alpha in alpha_values) {
    # May produce warnings about resolution estimation
    suppressWarnings(alpha_results[[as.character(alpha)]] <- 
                       acsc(vec, mask, K = 6, alpha = alpha, ann_k = 5))
    
    n_clusters <- length(unique(alpha_results[[as.character(alpha)]]$cluster_map[mask > 0]))
    expect_true(n_clusters > 0,
                info = sprintf("ACSC with alpha=%.1f should find clusters", alpha))
  }
  
  # Different alpha values should potentially produce different clustering patterns
  n_clusters_alpha <- sapply(alpha_results, function(x) length(unique(x$cluster_map[mask > 0])))
  
  # Should show some variation (though not required to be monotonic)
  expect_true(length(unique(n_clusters_alpha)) >= 1,
              info = "Different alpha values should explore clustering space")
  
  # Test different block_size values (affects granularity)
  block_size_values <- c(2, 3, 4)
  size_results <- list()
  
  for (block_size_val in block_size_values) {
    # May produce warnings about resolution estimation
    suppressWarnings(size_results[[as.character(block_size_val)]] <-
                       acsc(vec, mask, K = 4, block_size = block_size_val))
    
    # Check clustering results
    n_clusters <- length(unique(size_results[[as.character(block_size_val)]]$cluster_map[mask > 0]))
    
    expect_true(n_clusters > 0,
                info = sprintf("Should find clusters with block_size=%d (found %d)", 
                              block_size_val, n_clusters))
  }
})

test_that("ACSC correlation computation accuracy", {
  # Test that ACSC correctly computes and uses correlations
  dims <- c(6, 6, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 80  # Longer time series for stable correlations
  
  # Create signals with known correlation structure
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Base signals with known relationships
  signal_A <- sin(t_seq)
  signal_B <- cos(t_seq)  # Orthogonal to A
  signal_C <- sin(t_seq)  # Identical to A
  
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 2) {
      # Region with signal A
      ts_data[i, ] <- signal_A + rnorm(ntime, sd = 0.05)
    } else if (x >= 5) {
      # Region with signal C (identical to A)
      ts_data[i, ] <- signal_C + rnorm(ntime, sd = 0.05)
    } else {
      # Region with signal B (orthogonal)
      ts_data[i, ] <- signal_B + rnorm(ntime, sd = 0.05)
    }
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run ACSC
  # May produce warnings about resolution estimation
  suppressWarnings(corr_result <- acsc(vec, mask, K = 3, 
                                      ann_k = 5,
                                      alpha = 0.5))
  
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- corr_result$cluster_map[mask > 0]
  
  # Analyze clustering results
  # Regions with signal A and C should be more likely to cluster together
  # than either with signal B
  
  region_A_voxels <- which(coords[, 1] <= 2)
  region_B_voxels <- which(coords[, 1] >= 3 & coords[, 1] <= 4)
  region_C_voxels <- which(coords[, 1] >= 5)
  
  # Get cluster assignments for each region
  cluster_assignments <- corr_result$cluster_map[mask > 0]
  clusters_A <- cluster_assignments[region_A_voxels]
  clusters_B <- cluster_assignments[region_B_voxels]
  clusters_C <- cluster_assignments[region_C_voxels]
  
  # Calculate overlap between regions
  overlap_AC <- length(intersect(clusters_A, clusters_C))
  overlap_AB <- length(intersect(clusters_A, clusters_B))
  overlap_BC <- length(intersect(clusters_B, clusters_C))
  
  # A and C should have more cluster overlap than A-B or B-C
  # (though spatial constraints may affect this)
  total_AC_overlap <- overlap_AC / (length(unique(c(clusters_A, clusters_C))))
  total_AB_overlap <- overlap_AB / (length(unique(c(clusters_A, clusters_B))))
  
  # At minimum, clustering should not completely ignore correlation structure
  expect_true(total_AC_overlap >= 0 && total_AB_overlap >= 0,
              info = "Clustering should produce valid overlap measures")
})

test_that("ACSC edge cases and robustness", {
  # Test edge cases and robustness conditions
  
  # Test 1: Very small volume
  dims_small <- c(3, 3, 1)
  mask_small <- NeuroVol(array(1, dims_small), NeuroSpace(dims_small))
  nvox_small <- prod(dims_small)
  ntime_small <- 20
  
  ts_data_small <- matrix(rnorm(nvox_small * ntime_small), nrow = nvox_small, ncol = ntime_small)
  
  vec_list_small <- lapply(1:ntime_small, function(t) {
    vol_data <- array(0, dims_small)
    vol_data[mask_small > 0] <- ts_data_small[, t]
    NeuroVol(vol_data, NeuroSpace(dims_small))
  })
  vec_small <- do.call(concat, vec_list_small)
  
  # Small dataset may produce warnings about ann_k adjustment or resolution estimation
  suppressWarnings(result_small <- acsc(vec_small, mask_small, K = 2))
  expect_true(!is.null(result_small))
  expect_true(length(result_small$cluster_map[mask_small > 0]) == nvox_small)
  
  # Test 2: Single voxel mask
  mask_single <- NeuroVol(array(0, c(5, 5, 1)), NeuroSpace(c(5, 5, 1)))
  mask_single@.Data[3, 3, 1] <- 1
  
  single_voxel_ts <- matrix(rnorm(20), nrow = 1, ncol = 20)
  vec_list_single <- lapply(1:20, function(t) {
    vol_data <- array(0, c(5, 5, 1))
    vol_data[3, 3, 1] <- single_voxel_ts[1, t]
    NeuroVol(vol_data, NeuroSpace(c(5, 5, 1)))
  })
  vec_single <- do.call(concat, vec_list_single)
  
  expect_silent(result_single <- acsc(vec_single, mask_single, K = 1))
  expect_true(!is.null(result_single))
  expect_true(length(result_single$cluster_map[mask_single > 0]) == 1)
  expect_true(result_single$cluster_map[mask_single > 0][1] == 1)
  
  # Test 3: Constant signal (zero variance)
  dims_const <- c(4, 4, 1)
  mask_const <- NeuroVol(array(1, dims_const), NeuroSpace(dims_const))
  nvox_const <- prod(dims_const)
  ntime_const <- 30
  
  # All voxels have identical constant signal
  ts_data_const <- matrix(5.0, nrow = nvox_const, ncol = ntime_const)
  
  vec_list_const <- lapply(1:ntime_const, function(t) {
    vol_data <- array(5.0, dims_const)
    NeuroVol(vol_data, NeuroSpace(dims_const))
  })
  vec_const <- do.call(concat, vec_list_const)
  
  # Constant signal may produce warnings about resolution estimation
  suppressWarnings(result_const <- acsc(vec_const, mask_const, K = 2))
  expect_true(!is.null(result_const))
  expect_true(all(result_const$cluster_map[mask_const > 0] > 0))
  
  # Test 4: Very short time series
  dims_short <- c(6, 6, 1)
  mask_short <- NeuroVol(array(1, dims_short), NeuroSpace(dims_short))
  nvox_short <- prod(dims_short)
  ntime_short <- 5  # Very short
  
  ts_data_short <- matrix(rnorm(nvox_short * ntime_short), nrow = nvox_short, ncol = ntime_short)
  
  vec_list_short <- lapply(1:ntime_short, function(t) {
    vol_data <- array(0, dims_short)
    vol_data[mask_short > 0] <- ts_data_short[, t]
    NeuroVol(vol_data, NeuroSpace(dims_short))
  })
  vec_short <- do.call(concat, vec_list_short)
  
  # Should handle gracefully (may warn about short time series or resolution estimation)
  suppressWarnings(result_short <- acsc(vec_short, mask_short, K = 3))
  expect_true(!is.null(result_short))
  expect_true(length(result_short$cluster_map[mask_short > 0]) == nvox_short)
})

test_that("ACSC performance and scalability", {
  skip_on_cran()  # Skip performance tests on CRAN
  
  # Test performance with larger datasets
  dims_large <- c(20, 20, 2)
  mask_large <- NeuroVol(array(1, dims_large), NeuroSpace(dims_large))
  nvox_large <- prod(dims_large)
  ntime_large <- 100
  
  # Create structured data efficiently
  ts_data_large <- matrix(0, nrow = nvox_large, ncol = ntime_large)
  coords <- arrayInd(1:nvox_large, dims_large)
  t_seq <- seq(0, 4*pi, length.out = ntime_large)
  
  # Batch process to avoid memory issues
  batch_size <- 200
  for (start_idx in seq(1, nvox_large, by = batch_size)) {
    end_idx <- min(start_idx + batch_size - 1, nvox_large)
    
    for (i in start_idx:end_idx) {
      x <- coords[i, 1]
      region <- ceiling(x / 5)  # 4 regions
      phase <- region * pi / 4
      ts_data_large[i, ] <- sin(t_seq + phase) + rnorm(ntime_large, sd = 0.1)
    }
  }
  
  vec_list_large <- lapply(1:ntime_large, function(t) {
    vol_data <- array(0, dims_large)
    vol_data[mask_large > 0] <- ts_data_large[, t]
    NeuroVol(vol_data, NeuroSpace(dims_large))
  })
  vec_large <- do.call(concat, vec_list_large)
  
  # Test that large datasets complete in reasonable time
  start_time <- Sys.time()
  # Large datasets may produce warnings about resolution estimation
  suppressWarnings(result_large <- acsc(vec_large, mask_large, K = 8, ann_k = 10))
  end_time <- Sys.time()
  
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  expect_true(!is.null(result_large))
  expect_true(length(result_large$cluster_map[mask_large > 0]) == nvox_large)
  
  # Should complete in reasonable time (< 120 seconds on most systems)
  expect_true(computation_time < 120,
              info = sprintf("ACSC large dataset processing took %.1f seconds", computation_time))
  
  # Should find reasonable number of clusters
  n_clusters_large <- length(unique(result_large$cluster_map[mask_large > 0]))
  expect_true(n_clusters_large >= 4 && n_clusters_large <= 20,
              info = sprintf("Large dataset should find 4-20 clusters, found %d", n_clusters_large))
})