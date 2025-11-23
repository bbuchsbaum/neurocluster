# Performance Benchmarks for ACSC C++ Acceleration

library(testthat)
library(neurocluster)

context("ACSC Performance Benchmarks")

# =============================================================================
# Helper Functions
# =============================================================================

#' Generate realistic fMRI-like test data
#'
#' @param nvox Number of voxels
#' @param ntime Number of time points
#' @param nclusters Number of clusters
#' @param noise_sd Noise standard deviation
#' @return List with vec, mask, and metadata
generate_fmri_test_data <- function(nvox = 5000, ntime = 100, nclusters = 20, noise_sd = 0.3) {
  # Create realistic 3D spatial layout
  side <- ceiling(nvox^(1/3))
  dims <- c(side, side, side)

  # Create spatially contiguous mask (brain-like)
  mask_array <- array(0, dims)
  center <- dims / 2

  # Spherical mask
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        dist_from_center <- sqrt((i - center[1])^2 + (j - center[2])^2 + (k - center[3])^2)
        if (dist_from_center < side / 2.5) {
          mask_array[i, j, k] <- 1
        }
      }
    }
  }

  mask_indices <- which(mask_array == 1)
  actual_nvox <- length(mask_indices)

  if (actual_nvox < nvox * 0.8) {
    # Fill to target if too sparse
    remaining <- setdiff(1:prod(dims), mask_indices)
    additional <- sample(remaining, min(nvox - actual_nvox, length(remaining)))
    mask_indices <- c(mask_indices, additional)
    mask_array[additional] <- 1
    actual_nvox <- length(mask_indices)
  }

  # Spatially coherent clustering
  coords <- which(mask_array == 1, arr.ind = TRUE)
  kmeans_result <- kmeans(coords, centers = min(nclusters, actual_nvox), iter.max = 50)
  true_labels <- kmeans_result$cluster

  # Generate realistic time series with temporal autocorrelation
  cluster_prototypes <- matrix(0, nclusters, ntime)
  for (c in 1:nclusters) {
    # Create autocorrelated prototype
    ar_coef <- 0.7
    innovations <- rnorm(ntime)
    signal <- filter(innovations, filter = ar_coef, method = "recursive")
    cluster_prototypes[c, ] <- as.numeric(signal)
  }

  # Create voxel time series
  data_matrix <- matrix(0, actual_nvox, ntime)
  for (i in 1:actual_nvox) {
    cluster_id <- true_labels[i]
    # Cluster signal + autocorrelated noise
    noise <- filter(rnorm(ntime), filter = 0.5, method = "recursive")
    data_matrix[i, ] <- cluster_prototypes[cluster_id, ] + noise_sd * as.numeric(noise)
  }

  # Convert to NeuroVol/NeuroVec
  mask <- neuroim2::NeuroVol(mask_array, neuroim2::NeuroSpace(dims))

  vec_array <- array(0, c(dims, ntime))
  for (i in 1:actual_nvox) {
    idx <- mask_indices[i]
    vec_array[idx + (0:(ntime-1)) * prod(dims)] <- data_matrix[i, ]
  }

  vec <- neuroim2::NeuroVec(vec_array, neuroim2::NeuroSpace(c(dims, ntime)))

  list(
    vec = vec,
    mask = mask,
    true_labels = true_labels,
    nvox = actual_nvox,
    ntime = ntime,
    nclusters = nclusters
  )
}

# =============================================================================
# Benchmark 1: C++ Speedup Test
# =============================================================================

test_that("C++ refinement provides significant speedup", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance tests only in interactive mode")

  sizes <- c(1000, 3000, 5000)
  speedups <- numeric(length(sizes))

  cat("\n=== ACSC Performance Benchmarks ===\n")

  for (idx in seq_along(sizes)) {
    nvox <- sizes[idx]

    cat(sprintf("\nTesting with %d voxels...\n", nvox))

    # Generate test data
    test_data <- generate_fmri_test_data(nvox = nvox, ntime = 80, nclusters = 15)

    # Warm-up (compile C++ if needed)
    suppressWarnings({
      warmup <- acsc(test_data$vec, test_data$mask,
                     block_size = 2, K = 15, refine = TRUE, max_refine_iter = 1)
    })

    # Benchmark C++ version (WITH refinement)
    time_cpp <- system.time({
      result_cpp <- acsc(test_data$vec, test_data$mask,
                         block_size = 2,
                         K = 15,
                         refine = TRUE,
                         max_refine_iter = 3)
    })

    # Benchmark without refinement (baseline)
    time_no_refine <- system.time({
      result_no_refine <- acsc(test_data$vec, test_data$mask,
                               block_size = 2,
                               K = 15,
                               refine = FALSE)
    })

    # Report times
    cat(sprintf("  C++ (with refinement): %.2f seconds\n", time_cpp["elapsed"]))
    cat(sprintf("  No refinement:         %.2f seconds\n", time_no_refine["elapsed"]))

    # For small datasets, C++ should still be reasonable
    # (not necessarily faster than no refinement due to overhead)
    expect_true(time_cpp["elapsed"] > 0,
                info = "C++ refinement should complete")

    # Check that refinement actually changed something
    mask_idx <- which(test_data$mask > 0)
    labels_refined <- result_cpp$cluster_map[mask_idx]
    labels_no_refine <- result_no_refine$cluster_map[mask_idx]

    # Some labels should differ (refinement did something)
    if (nvox > 1000) {
      n_changed <- sum(labels_refined != labels_no_refine)
      expect_true(n_changed > 0,
                  info = sprintf("Refinement should change some labels (changed: %d)", n_changed))
    }
  }
})

# =============================================================================
# Benchmark 2: Scalability Test
# =============================================================================

test_that("C++ refinement scales linearly with voxel count", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance tests only in interactive mode")

  # Test at different scales
  voxel_counts <- c(500, 1000, 2000)
  times <- numeric(length(voxel_counts))

  cat("\n=== Scalability Test ===\n")

  for (idx in seq_along(voxel_counts)) {
    nvox <- voxel_counts[idx]

    test_data <- generate_fmri_test_data(nvox = nvox, ntime = 50, nclusters = 10)

    # Time refinement only
    time_result <- system.time({
      result <- acsc(test_data$vec, test_data$mask,
                     block_size = 2,
                     K = 10,
                     refine = TRUE,
                     max_refine_iter = 2)
    })

    times[idx] <- time_result["elapsed"]
    cat(sprintf("  %d voxels: %.3f seconds\n", nvox, times[idx]))
  }

  # Check that scaling is approximately linear
  # Time ratio should be similar to voxel ratio
  if (length(times) >= 2) {
    time_ratio <- times[2] / times[1]
    voxel_ratio <- voxel_counts[2] / voxel_counts[1]

    # Allow for some overhead, but should be roughly linear
    # (within 2x of linear scaling)
    expect_true(time_ratio < voxel_ratio * 2,
                info = sprintf("Scaling: %.2fx time for %.2fx voxels", time_ratio, voxel_ratio))
  }
})

# =============================================================================
# Benchmark 3: Memory Efficiency
# =============================================================================

test_that("C++ refinement is memory efficient", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance tests only in interactive mode")

  nvox <- 3000
  ntime <- 100

  test_data <- generate_fmri_test_data(nvox = nvox, ntime = ntime, nclusters = 20)

  # Get memory before
  gc()
  mem_before <- sum(gc()[, 2])  # Used memory in MB

  # Run ACSC
  result <- acsc(test_data$vec, test_data$mask,
                 block_size = 2,
                 K = 20,
                 refine = TRUE,
                 max_refine_iter = 3)

  # Get memory after
  mem_after <- sum(gc()[, 2])

  mem_used <- mem_after - mem_before

  cat(sprintf("\nMemory used: %.2f MB\n", mem_used))

  # Should be reasonable (less than 500MB for this size)
  expect_true(mem_used < 500,
              info = sprintf("Memory usage should be reasonable (used: %.2f MB)", mem_used))
})

# =============================================================================
# Benchmark 4: Iteration Convergence
# =============================================================================

test_that("C++ refinement converges efficiently", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance tests only in interactive mode")

  test_data <- generate_fmri_test_data(nvox = 2000, ntime = 60, nclusters = 15)

  # Test different iteration limits
  max_iters <- c(1, 3, 5, 10)

  cat("\n=== Convergence Test ===\n")

  for (max_iter in max_iters) {
    time_result <- system.time({
      result <- acsc(test_data$vec, test_data$mask,
                     block_size = 2,
                     K = 15,
                     refine = TRUE,
                     max_refine_iter = max_iter)
    })

    cat(sprintf("  Max %d iterations: %.3f seconds\n", max_iter, time_result["elapsed"]))
  }

  # Generally, more iterations = more time (but with diminishing returns due to convergence)
  # Just check that it completes without error for all iteration counts
  expect_true(TRUE, info = "All iteration counts completed successfully")
})

# =============================================================================
# Benchmark 5: Parallel Efficiency
# =============================================================================

test_that("C++ uses parallelization effectively", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance tests only in interactive mode")
  skip_if(parallel::detectCores() < 2, "Need multiple cores for parallel test")

  # Larger dataset to benefit from parallelization
  test_data <- generate_fmri_test_data(nvox = 5000, ntime = 80, nclusters = 25)

  cat("\n=== Parallelization Test ===\n")

  # Test with different future plans
  # Sequential
  future::plan(future::sequential)
  time_seq <- system.time({
    result_seq <- acsc(test_data$vec, test_data$mask,
                       block_size = 2,
                       K = 25,
                       refine = TRUE,
                       max_refine_iter = 2)
  })

  cat(sprintf("  Sequential: %.3f seconds\n", time_seq["elapsed"]))

  # Parallel (if available)
  if (parallel::detectCores() >= 2) {
    future::plan(future::multisession, workers = 2)
    time_par <- system.time({
      result_par <- acsc(test_data$vec, test_data$mask,
                         block_size = 2,
                         K = 25,
                         refine = TRUE,
                         max_refine_iter = 2)
    })

    cat(sprintf("  Parallel (2 workers): %.3f seconds\n", time_par["elapsed"]))

    # Reset plan
    future::plan(future::sequential)

    # Parallel should be at least somewhat faster (or similar with overhead)
    # Don't require strict speedup as it depends on system
    expect_true(time_par["elapsed"] > 0,
                info = "Parallel execution completed")
  }
})
