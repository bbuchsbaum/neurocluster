library(testthat)
library(neurocluster)
library(neuroim2)

# Benchmark test comparing FLASH-3D with other algorithms
# Note: These tests focus on correctness and relative performance
# For full benchmarking, use a dedicated benchmark script

create_benchmark_data <- function(dims, n_time, seed = 123) {
  set.seed(seed)
  
  # Create spherical mask
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
  
  # Create structured data with 4 distinct regions
  nvoxels <- sum(mask_array)
  mask_idx <- which(mask_array)
  coords <- arrayInd(mask_idx, dims)
  ts_data <- matrix(0, nrow = nvoxels, ncol = n_time)
  
  t_seq <- seq(0, 4*pi, length.out = n_time)
  
  for (i in 1:nvoxels) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    z <- coords[i, 3]
    
    # Create 4 regions with different temporal patterns
    if (x < center[1] && y < center[2]) {
      # Region 1: Low frequency sine
      ts_data[i, ] <- sin(t_seq) + rnorm(n_time, sd = 0.1)
    } else if (x >= center[1] && y < center[2]) {
      # Region 2: High frequency sine
      ts_data[i, ] <- sin(3 * t_seq) + rnorm(n_time, sd = 0.1)
    } else if (x < center[1] && y >= center[2]) {
      # Region 3: Cosine
      ts_data[i, ] <- cos(t_seq) + rnorm(n_time, sd = 0.1)
    } else {
      # Region 4: Linear trend
      ts_data[i, ] <- seq(-1, 1, length.out = n_time) + rnorm(n_time, sd = 0.1)
    }
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:n_time, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask_array] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  list(vec = vec, mask = mask, nvoxels = nvoxels, ground_truth = coords)
}

# Helper to calculate cluster compactness
calculate_compactness <- function(result, mask) {
  mask_idx <- which(mask > 0)
  coords <- index_to_coord(mask, mask_idx)
  
  # Get number of clusters
  if (!is.null(result$K)) {
    n_clusters <- result$K
    cluster_ids <- 1:n_clusters
  } else {
    cluster_ids <- sort(unique(result$cluster))
    n_clusters <- length(cluster_ids)
  }
  
  compactness <- numeric(n_clusters)
  for (i in seq_along(cluster_ids)) {
    k <- cluster_ids[i]
    cluster_voxels <- which(result$cluster == k)
    if (length(cluster_voxels) > 0) {
      cluster_coords <- coords[cluster_voxels, , drop = FALSE]
      centroid <- colMeans(cluster_coords)
      distances <- rowSums((sweep(cluster_coords, 2, centroid))^2)
      compactness[i] <- mean(distances)
    }
  }
  mean(compactness, na.rm = TRUE)
}

# Helper to calculate temporal coherence
calculate_temporal_coherence <- function(result, vec, mask) {
  mask_idx <- which(mask > 0)
  ts_matrix <- series(vec, mask_idx)  # Returns T x Nmask
  
  # Get number of clusters
  if (!is.null(result$K)) {
    n_clusters <- result$K
    cluster_ids <- 1:n_clusters
  } else {
    cluster_ids <- sort(unique(result$cluster))
    n_clusters <- length(cluster_ids)
  }
  
  coherence <- numeric(n_clusters)
  for (i in seq_along(cluster_ids)) {
    k <- cluster_ids[i]
    cluster_voxels <- which(result$cluster == k)
    if (length(cluster_voxels) > 1) {
      # ts_matrix is T x Nmask, so select columns for cluster voxels
      cluster_ts <- ts_matrix[, cluster_voxels, drop = FALSE]
      # Correlate across voxels (columns)
      cor_mat <- cor(cluster_ts)
      coherence[i] <- mean(cor_mat[upper.tri(cor_mat)])
    } else if (length(cluster_voxels) == 1) {
      coherence[i] <- 1
    }
  }
  mean(coherence, na.rm = TRUE)
}

test_that("Algorithm comparison on small data", {
  skip_if_not(requireNamespace("future", quietly = TRUE))
  
  # Small dataset for quick testing
  data <- create_benchmark_data(c(12, 12, 12), n_time = 30, seed = 42)
  K <- 8
  
  # FLASH-3D
  time_flash3d <- system.time({
    result_flash3d <- supervoxels_flash3d(data$vec, data$mask, K = K, 
                                          verbose = FALSE)
  })
  
  # Supervoxels (traditional iterative method)
  time_supervoxels <- system.time({
    result_supervoxels <- supervoxels(data$vec, data$mask, K = K,
                                      iterations = 5, verbose = FALSE)
  })
  
  # SNIC
  time_snic <- system.time({
    result_snic <- snic(data$vec, data$mask, K = K, compactness = 5)
  })
  
  # SLIC4D
  time_slic4d <- system.time({
    result_slic4d <- slic4d_supervoxels(data$vec, data$mask, K = K,
                                        compactness = 10, max_iter = 5,
                                        verbose = FALSE)
  })
  
  # All methods should produce valid results
  expect_s3_class(result_flash3d, "cluster_result")
  expect_s3_class(result_supervoxels, "cluster_result")
  expect_s3_class(result_snic, "cluster_result")
  expect_s3_class(result_slic4d, "cluster_result")
  
  # All should have reasonable number of clusters
  expect_true(result_flash3d$K <= K * 1.5)
  
  # Check supervoxels
  if (!is.null(result_supervoxels$K)) {
    expect_true(result_supervoxels$K <= K * 1.5)
  } else {
    expect_true(length(unique(result_supervoxels$cluster)) <= K * 1.5)
  }
  
  # Check snic
  if (!is.null(result_snic$K)) {
    expect_true(result_snic$K <= K * 1.5)
  } else {
    expect_true(length(unique(result_snic$cluster)) <= K * 1.5)
  }
  
  # Check slic4d
  if (!is.null(result_slic4d$K)) {
    expect_true(result_slic4d$K <= K * 1.5)
  } else {
    expect_true(length(unique(result_slic4d$cluster)) <= K * 1.5)
  }
  
  # Calculate quality metrics
  compact_flash3d <- calculate_compactness(result_flash3d, data$mask)
  compact_supervoxels <- calculate_compactness(result_supervoxels, data$mask)
  compact_snic <- calculate_compactness(result_snic, data$mask)
  compact_slic4d <- calculate_compactness(result_slic4d, data$mask)
  
  coherence_flash3d <- calculate_temporal_coherence(result_flash3d, data$vec, data$mask)
  coherence_supervoxels <- calculate_temporal_coherence(result_supervoxels, data$vec, data$mask)
  coherence_snic <- calculate_temporal_coherence(result_snic, data$vec, data$mask)
  coherence_slic4d <- calculate_temporal_coherence(result_slic4d, data$vec, data$mask)
  
  # Print benchmark results (for manual inspection)
  cat("\n=== Small Data Benchmark (12x12x12, 30 timepoints, K=", K, ") ===\n")
  cat("Execution times (seconds):\n")
  cat("  FLASH-3D:    ", time_flash3d["elapsed"], "\n")
  cat("  Supervoxels: ", time_supervoxels["elapsed"], "\n")
  cat("  SNIC:        ", time_snic["elapsed"], "\n")
  cat("  SLIC4D:      ", time_slic4d["elapsed"], "\n")
  
  cat("\nCluster counts:\n")
  cat("  FLASH-3D:    ", result_flash3d$K, "\n")
  cat("  Supervoxels: ", result_supervoxels$K, "\n")
  cat("  SNIC:        ", result_snic$K, "\n")
  cat("  SLIC4D:      ", result_slic4d$K, "\n")
  
  cat("\nSpatial compactness (lower is better):\n")
  cat("  FLASH-3D:    ", round(compact_flash3d, 3), "\n")
  cat("  Supervoxels: ", round(compact_supervoxels, 3), "\n")
  cat("  SNIC:        ", round(compact_snic, 3), "\n")
  cat("  SLIC4D:      ", round(compact_slic4d, 3), "\n")
  
  cat("\nTemporal coherence (higher is better):\n")
  cat("  FLASH-3D:    ", round(coherence_flash3d, 3), "\n")
  cat("  Supervoxels: ", round(coherence_supervoxels, 3), "\n")
  cat("  SNIC:        ", round(coherence_snic, 3), "\n")
  cat("  SLIC4D:      ", round(coherence_slic4d, 3), "\n")
  
  # Basic expectations - all should produce reasonable results
  expect_true(all(c(compact_flash3d, compact_supervoxels, compact_snic, compact_slic4d) > 0))
  expect_true(all(c(coherence_flash3d, coherence_supervoxels, coherence_snic, coherence_slic4d) > 0))
  expect_true(all(c(coherence_flash3d, coherence_supervoxels, coherence_snic, coherence_slic4d) <= 1))
})

test_that("Algorithm scaling with data size", {
  skip_if_not(requireNamespace("future", quietly = TRUE))
  skip_on_cran()  # Skip on CRAN due to time constraints
  
  # Test scaling with different data sizes
  sizes <- list(
    small = c(10, 10, 10),
    medium = c(15, 15, 15)
  )
  
  times <- list()
  
  for (size_name in names(sizes)) {
    dims <- sizes[[size_name]]
    data <- create_benchmark_data(dims, n_time = 25, seed = 123)
    K <- 10
    
    times[[size_name]] <- list(
      flash3d = system.time({
        supervoxels_flash3d(data$vec, data$mask, K = K, verbose = FALSE)
      })["elapsed"],
      
      supervoxels = system.time({
        supervoxels(data$vec, data$mask, K = K, iterations = 3, verbose = FALSE)
      })["elapsed"],
      
      snic = system.time({
        snic(data$vec, data$mask, K = K, compactness = 5)
      })["elapsed"],
      
      slic4d = system.time({
        slic4d_supervoxels(data$vec, data$mask, K = K, 
                          compactness = 10, max_iter = 3, verbose = FALSE)
      })["elapsed"]
    )
  }
  
  # Calculate scaling factors
  scaling <- list()
  for (alg in c("flash3d", "supervoxels", "snic", "slic4d")) {
    scaling[[alg]] <- times$medium[[alg]] / times$small[[alg]]
  }
  
  cat("\n=== Scaling Analysis (10³ to 15³ voxels) ===\n")
  cat("Time scaling factors (lower is better):\n")
  cat("  FLASH-3D:    ", round(scaling$flash3d, 2), "x\n")
  cat("  Supervoxels: ", round(scaling$supervoxels, 2), "x\n")
  cat("  SNIC:        ", round(scaling$snic, 2), "x\n")
  cat("  SLIC4D:      ", round(scaling$slic4d, 2), "x\n")
  
  # All algorithms should scale reasonably (not worse than cubic)
  # Volume increases by (15/10)³ = 3.375x
  expect_true(all(unlist(scaling) < 10))  # Should scale better than 10x
})

test_that("FLASH-3D parameter sensitivity", {
  data <- create_benchmark_data(c(10, 10, 10), n_time = 20, seed = 456)
  K <- 6
  
  # Test different bit widths
  result_64 <- supervoxels_flash3d(data$vec, data$mask, K = K, 
                                   bits = 64, verbose = FALSE)
  result_128 <- supervoxels_flash3d(data$vec, data$mask, K = K,
                                    bits = 128, verbose = FALSE)
  
  coherence_64 <- calculate_temporal_coherence(result_64, data$vec, data$mask)
  coherence_128 <- calculate_temporal_coherence(result_128, data$vec, data$mask)
  
  # 128-bit should provide similar or better temporal coherence
  expect_true(coherence_128 >= coherence_64 * 0.9)  # Allow 10% tolerance
  
  # Test different DCT coefficients
  result_dct_low <- supervoxels_flash3d(data$vec, data$mask, K = K,
                                        dctM = 6, verbose = FALSE)
  result_dct_high <- supervoxels_flash3d(data$vec, data$mask, K = K,
                                         dctM = 16, verbose = FALSE)
  
  # Both should produce valid results
  expect_s3_class(result_dct_low, "cluster_result")
  expect_s3_class(result_dct_high, "cluster_result")
  
  # Test lambda balance
  result_spatial <- supervoxels_flash3d(data$vec, data$mask, K = K,
                                        lambda_s = 2.0, lambda_t = 0.5,
                                        verbose = FALSE)
  result_temporal <- supervoxels_flash3d(data$vec, data$mask, K = K,
                                         lambda_s = 0.2, lambda_t = 2.0,
                                         verbose = FALSE)
  
  compact_spatial <- calculate_compactness(result_spatial, data$mask)
  compact_temporal <- calculate_compactness(result_temporal, data$mask)
  coherence_spatial <- calculate_temporal_coherence(result_spatial, data$vec, data$mask)
  coherence_temporal <- calculate_temporal_coherence(result_temporal, data$vec, data$mask)
  
  # Spatial-weighted should be more compact
  expect_true(compact_spatial <= compact_temporal * 1.2)  # Allow some tolerance
  
  # Temporal-weighted should have better coherence
  expect_true(coherence_temporal >= coherence_spatial * 0.95)  # Allow some tolerance
  
  cat("\n=== FLASH-3D Parameter Sensitivity ===\n")
  cat("Bit width comparison:\n")
  cat("  64-bit coherence:  ", round(coherence_64, 3), "\n")
  cat("  128-bit coherence: ", round(coherence_128, 3), "\n")
  cat("\nLambda weighting:\n")
  cat("  Spatial-weighted compactness:  ", round(compact_spatial, 3), "\n")
  cat("  Temporal-weighted compactness: ", round(compact_temporal, 3), "\n")
  cat("  Spatial-weighted coherence:    ", round(coherence_spatial, 3), "\n")
  cat("  Temporal-weighted coherence:   ", round(coherence_temporal, 3), "\n")
})

test_that("Algorithms handle structured data correctly", {
  # Create data with very clear structure
  dims <- c(12, 12, 6)
  set.seed(789)
  
  mask_array <- array(TRUE, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  # Create 4 quadrants with distinct signals
  nvoxels <- prod(dims)
  mask_idx <- 1:nvoxels
  coords <- arrayInd(mask_idx, dims)
  ts_data <- matrix(0, nrow = nvoxels, ncol = 40)
  
  true_labels <- integer(nvoxels)
  
  for (i in 1:nvoxels) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 6 && y <= 6) {
      ts_data[i, ] <- rep(c(1, -1), 20) + rnorm(40, sd = 0.05)
      true_labels[i] <- 1
    } else if (x > 6 && y <= 6) {
      ts_data[i, ] <- rep(c(-1, 1), 20) + rnorm(40, sd = 0.05)
      true_labels[i] <- 2
    } else if (x <= 6 && y > 6) {
      ts_data[i, ] <- rep(c(1, 1, -1, -1), 10) + rnorm(40, sd = 0.05)
      true_labels[i] <- 3
    } else {
      ts_data[i, ] <- rep(c(-1, -1, 1, 1), 10) + rnorm(40, sd = 0.05)
      true_labels[i] <- 4
    }
  }
  
  vec_list <- lapply(1:40, function(t) {
    vol_data <- array(ts_data[, t], dims)
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test all algorithms with K=4 (matching true structure)
  result_flash3d <- supervoxels_flash3d(vec, mask, K = 4, verbose = FALSE)
  result_supervoxels <- supervoxels(vec, mask, K = 4, iterations = 5, verbose = FALSE)
  result_snic <- snic(vec, mask, K = 4, compactness = 5)
  result_slic4d <- slic4d_supervoxels(vec, mask, K = 4, 
                                      compactness = 10, max_iter = 5,
                                      verbose = FALSE)
  
  # Calculate how well each algorithm recovers the structure
  # Using normalized mutual information or adjusted Rand index would be ideal
  # For simplicity, we'll check temporal coherence which should be very high
  
  coherence_flash3d <- calculate_temporal_coherence(result_flash3d, vec, mask)
  coherence_supervoxels <- calculate_temporal_coherence(result_supervoxels, vec, mask)
  coherence_snic <- calculate_temporal_coherence(result_snic, vec, mask)
  coherence_slic4d <- calculate_temporal_coherence(result_slic4d, vec, mask)
  
  cat("\n=== Structured Data Recovery ===\n")
  cat("Temporal coherence (should be high for structured data):\n")
  cat("  FLASH-3D:    ", round(coherence_flash3d, 3), "\n")
  cat("  Supervoxels: ", round(coherence_supervoxels, 3), "\n")
  cat("  SNIC:        ", round(coherence_snic, 3), "\n")
  cat("  SLIC4D:      ", round(coherence_slic4d, 3), "\n")
  
  # All algorithms should achieve good coherence on structured data
  expect_true(all(c(coherence_flash3d, coherence_supervoxels, 
                   coherence_snic, coherence_slic4d) > 0.7))
})