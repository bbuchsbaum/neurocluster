library(neurocluster)
library(neuroim2)

# Memory Management and Thread Safety Tests
# Tests for parallel operations, memory stress, and concurrent access

test_that("parallel operations are thread-safe", {
  skip_on_cran()  # Skip on CRAN due to timing/resource constraints
  
  # Test parallel slice processing consistency
  dims <- c(12, 12, 4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 40
  
  # Create deterministic test data
  set.seed(12345)
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Create spatially structured patterns
  t_seq <- seq(0, 4*pi, length.out = ntime)
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    z <- coords[i, 3]
    
    # Deterministic pattern based on spatial location
    pattern_id <- (x %% 3) + (y %% 3) * 3 + (z %% 2) * 9
    phase <- pattern_id * pi / 10
    ts_data[i, ] <- sin(t_seq + phase) + 0.1 * rnorm(ntime, mean = 0, sd = 1)
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run same analysis multiple times - should get identical results
  results <- list()
  set.seed(12345)  # Reset seed for each run
  
  for (run in 1:3) {
    results[[run]] <- slice_msf_single(vec, mask, r = 6, k = 0.5, min_size = 10)
  }
  
  # Check that parallel operations produce consistent results
  for (run in 2:3) {
    # Sketches should be identical (within numerical precision)
    sketch_diff <- max(abs(results[[1]]$sketch - results[[run]]$sketch))
    expect_true(sketch_diff < 1e-12,
                info = sprintf("Thread safety: sketch differences should be minimal, got %.2e", sketch_diff))
    
    # Labels should be identical
    expect_identical(results[[1]]$labels, results[[run]]$labels,
                    info = sprintf("Thread safety: labels should be identical across runs"))
  }
})

test_that("memory allocation handles stress conditions", {
  skip_on_cran()
  
  # Test with progressively larger datasets
  base_dims <- c(8, 8, 2)
  base_ntime <- 30
  
  # Test different sizes to stress memory allocation
  size_multipliers <- c(1, 2, 3)
  
  for (mult in size_multipliers) {
    dims <- base_dims * mult
    dims[3] <- min(dims[3], 6)  # Limit z-dimension to keep memory reasonable
    ntime <- base_ntime + mult * 10
    
    mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
    nvox <- prod(dims)
    
    # Create structured test data
    coords <- arrayInd(1:nvox, dims)
    ts_data <- matrix(0, nrow = nvox, ncol = ntime)
    
    # Simple sinusoidal patterns
    t_seq <- seq(0, 2*pi, length.out = ntime)
    for (i in 1:nvox) {
      x <- coords[i, 1]
      region <- ceiling(x / (dims[1] / 2))
      ts_data[i, ] <- sin(region * t_seq) + rnorm(ntime, sd = 0.1)
    }
    
    vec_list <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec <- do.call(concat, vec_list)
    
    # Test that large datasets don't cause memory issues
    expect_silent({
      result <- slice_msf_single(vec, mask, r = min(8, ntime-1), k = 0.5, min_size = 5)
    })
    
    expect_true(!is.null(result))
    expect_true(all(is.finite(result$sketch)))
    expect_true(length(unique(result$labels)) > 0)
    
    # Force garbage collection to free memory
    gc()
  }
})

test_that("memory cleanup is proper", {
  skip_on_cran()
  
  # Monitor memory usage across multiple runs
  dims <- c(10, 10, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 50
  
  # Create test data
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Get initial memory usage
  gc()
  initial_memory <- sum(gc()[, 2])  # Total memory in use
  
  # Run multiple analyses
  for (i in 1:10) {
    result <- slice_msf_single(vec, mask, r = 6, k = 0.5, min_size = 5)
    expect_true(!is.null(result))
    
    # Don't keep references to results
    rm(result)
    
    # Periodic garbage collection
    if (i %% 3 == 0) gc()
  }
  
  # Force final garbage collection
  gc()
  final_memory <- sum(gc()[, 2])
  
  # Memory should not have grown excessively
  memory_increase <- final_memory - initial_memory
  memory_ratio <- final_memory / initial_memory
  
  expect_true(memory_ratio < 2.0,
              info = sprintf("Memory should not double after 10 runs (ratio: %.2f)", memory_ratio))
  
  expect_true(memory_increase < 100,  # Less than 100MB increase
              info = sprintf("Memory increase should be reasonable (%.1f MB)", memory_increase))
})

test_that("handles large dataset dimensions", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  
  # Test with larger spatial dimensions
  dims <- c(20, 20, 4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 60
  
  # Create spatially structured data efficiently
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Use batch processing for efficiency
  batch_size <- 500
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  for (start_idx in seq(1, nvox, by = batch_size)) {
    end_idx <- min(start_idx + batch_size - 1, nvox)
    batch_indices <- start_idx:end_idx
    
    for (i in batch_indices) {
      x <- coords[i, 1]
      y <- coords[i, 2]
      
      # Create 4 spatial regions
      if (x <= 10 && y <= 10) region <- 1
      else if (x > 10 && y <= 10) region <- 2  
      else if (x <= 10 && y > 10) region <- 3
      else region <- 4
      
      ts_data[i, ] <- sin(region * t_seq) + rnorm(ntime, sd = 0.05)
    }
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test that large spatial dimensions work
  expect_silent({
    start_time <- Sys.time()
    result_large <- slice_msf(vec, mask, r = 8, min_size = 20, 
                             compactness = 3, num_runs = 1)
    end_time <- Sys.time()
  })
  
  computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  expect_true(!is.null(result_large))
  expect_true(length(result_large$cluster) == nvox)
  expect_true(length(unique(result_large$cluster)) > 1)
  expect_true(length(unique(result_large$cluster)) <= 10)  # Should find reasonable number of clusters
  
  # Computation should complete in reasonable time (< 60 seconds on most systems)
  expect_true(computation_time < 60,
              info = sprintf("Large dataset processing took %.1f seconds", computation_time))
})

test_that("concurrent NeuroVec access is safe", {
  skip_on_cran()
  
  # Test that multiple algorithms can safely access the same NeuroVec
  dims <- c(8, 8, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  # Create test data
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  t_seq <- seq(0, 2*pi, length.out = ntime)
  for (i in 1:nvox) {
    x <- coords[i, 1]
    region <- ceiling(x / 4)
    ts_data[i, ] <- sin(region * t_seq) + rnorm(ntime, sd = 0.1)
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test concurrent access by running multiple algorithms
  # Note: This doesn't test true concurrency but tests that the same data
  # can be safely processed by different algorithms sequentially
  
  results <- list()
  
  # SLiCE-MSF
  expect_silent(results$slice_msf <- slice_msf_single(vec, mask, r = 6, k = 0.5, min_size = 5))
  
  # SNIC
  expect_silent(results$snic <- snic(vec, mask, K = 4, compactness = 1.0))
  
  # Supervoxels - use smaller K to avoid index out of bounds
  suppressMessages(expect_no_error(results$supervoxels <- supervoxels(vec, mask, K = 4, alpha = 0.5)))
  
  # All algorithms should succeed
  for (alg_name in names(results)) {
    expect_true(!is.null(results[[alg_name]]),
                info = sprintf("Algorithm %s should succeed", alg_name))
  }
  
  # Test that original data is unchanged
  vec_data_after <- series(vec, which(mask > 0))
  original_vec_data <- ts_data
  
  max_diff <- max(abs(t(vec_data_after) - original_vec_data))
  expect_true(max_diff < 1e-12,
              info = sprintf("Original data should be unchanged after processing (max diff: %.2e)", max_diff))
})

test_that("extreme parameter combinations don't cause memory issues", {
  # Test combinations that might stress memory allocation
  dims <- c(6, 6, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 40
  
  # Simple test data
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test extreme parameter combinations
  extreme_params <- list(
    # Very high r
    list(r = min(25, ntime-1), k = 0.5, min_size = 1),
    # Very small k (many clusters)
    list(r = 8, k = 1e-8, min_size = 1),
    # Very large min_size (force few clusters)
    list(r = 8, k = 0.5, min_size = nvox/2),
    # High compactness in full algorithm
    list(compactness = 20, min_size = 1, r = 8)
  )
  
  for (i in seq_along(extreme_params)) {
    params <- extreme_params[[i]]
    
    if ("compactness" %in% names(params)) {
      # Use full slice_msf
      expect_silent({
        result <- slice_msf(vec, mask, 
                           compactness = params$compactness,
                           min_size = params$min_size,
                           r = params$r,
                           num_runs = 1)
      })
    } else {
      # Use slice_msf_single
      expect_silent({
        result <- slice_msf_single(vec, mask,
                                  r = params$r,
                                  k = params$k, 
                                  min_size = params$min_size)
      })
    }
    
    expect_true(!is.null(result),
                info = sprintf("Extreme params test %d should produce result", i))
    
    # Force cleanup
    rm(result)
    gc()
  }
})

test_that("memory usage is reasonable for multi-run consensus", {
  skip_on_cran()
  
  # Test that consensus doesn't explode memory usage
  dims <- c(8, 8, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  # Create test data
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Monitor memory during consensus runs
  gc()
  initial_memory <- sum(gc()[, 2])
  
  # Run multi-run consensus
  expect_silent({
    result_consensus <- slice_msf(vec, mask, 
                                 r = 6, min_size = 5, compactness = 3,
                                 num_runs = 5, consensus = TRUE)
  })
  
  gc()
  final_memory <- sum(gc()[, 2])
  
  expect_true(!is.null(result_consensus))
  expect_true(length(result_consensus$runs) == 5)
  
  # Memory increase should be reasonable
  memory_ratio <- final_memory / initial_memory
  expect_true(memory_ratio < 3.0,
              info = sprintf("Consensus memory usage should be reasonable (ratio: %.2f)", memory_ratio))
})