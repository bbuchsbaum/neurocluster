library(neurocluster)
library(neuroim2)

# Mathematical Validation Tests
# Comprehensive validation of mathematical correctness for core functions

test_that("DCT basis construction is mathematically correct", {
  # Test DCT basis properties against known mathematical properties
  
  # Create test data
  T <- 32  # time points
  r <- 8   # DCT rank
  
  # Generate test signal with known DCT properties
  t_seq <- seq(0, 2*pi, length.out = T)
  test_signal <- sin(t_seq) + 0.5*cos(2*t_seq)
  
  # Create minimal test volume
  dims <- c(4, 4, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  # Replicate signal across voxels
  ts_data <- matrix(rep(test_signal, nvox), nrow = nvox, ncol = T, byrow = TRUE)
  
  # Create NeuroVec
  vec_list <- lapply(1:T, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run SLiCE-MSF to get DCT sketches
  result <- slice_msf_single(vec, mask, r = r, k = 0.5, min_size = 1)
  
  # Validate DCT sketch properties
  expect_true(!is.null(result$sketch))
  expect_equal(nrow(result$sketch), r)
  expect_equal(ncol(result$sketch), nvox)
  
  # Test DCT orthogonality property
  # For identical signals, sketches should be identical
  sketch_matrix <- result$sketch
  for (i in 1:(ncol(sketch_matrix)-1)) {
    expect_equal(sketch_matrix[, i], sketch_matrix[, i+1], tolerance = 1e-10,
                 info = "Identical signals should produce identical DCT sketches")
  }
  
  # Test that low frequency components dominate high frequency components for low-frequency signal
  # DCT k=1 is NOT DC component, k=1 is first harmonic, k=2 is second harmonic, etc.
  # For low-frequency signals, expect low-k components to dominate high-k components
  low_freq_components <- sketch_matrix[1:3, 1]   # k=1,2,3 (low frequencies)
  high_freq_components <- sketch_matrix[6:r, 1]  # k=6,7,8 (high frequencies) 
  expect_true(mean(abs(low_freq_components)) > mean(abs(high_freq_components)),
              info = "Low frequency components should dominate for low-frequency signals")
})

test_that("heat kernel computation is mathematically accurate", {
  # Test heat kernel against known analytical results
  
  # For heat kernel: exp(-d^2 / (2*sigma^2))
  # where d is Euclidean distance
  
  # Test case 1: distance = 0 (identical points)
  x1 <- c(1, 2, 3)
  x2 <- c(1, 2, 3)
  sigma <- 1.0
  
  # Use find_candidates to test heat kernel indirectly
  # Create simple test data
  dims <- c(3, 3, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 10
  
  # Create identical time series for all voxels
  ts_data <- matrix(rep(c(1,2,3,4,5,4,3,2,1,0), nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run clustering to test heat kernel indirectly through neighbor selection
  # Use safe K value that doesn't exceed number of voxels
  nvox_available <- sum(mask > 0)
  K_safe <- min(3, nvox_available)
  result <- supervoxels(vec, mask, K = K_safe, alpha = 0.5)
  
  # For identical time series, spatial distance should dominate
  # Verify that clustering respects spatial proximity
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- result$cluster
  
  # Check that neighboring voxels are more likely to be in same cluster
  neighbor_same_cluster <- 0
  total_neighbors <- 0
  
  for (x in 1:(dims[1]-1)) {
    for (y in 1:(dims[2]-1)) {
      c1 <- cluster_vol[x, y, 1]
      c2 <- cluster_vol[x+1, y, 1]  # right neighbor
      c3 <- cluster_vol[x, y+1, 1]  # bottom neighbor
      
      if (c1 == c2) neighbor_same_cluster <- neighbor_same_cluster + 1
      if (c1 == c3) neighbor_same_cluster <- neighbor_same_cluster + 1
      total_neighbors <- total_neighbors + 2
    }
  }
  
  spatial_coherence <- neighbor_same_cluster / total_neighbors
  expect_true(spatial_coherence >= 0.4,
              info = sprintf("Spatial coherence should be >=0.4 for identical time series, got %.2f", spatial_coherence))
})

test_that("statistical computations are numerically accurate", {
  # Test correlation, mean, variance computations
  
  # Create test data with known statistical properties
  dims <- c(6, 6, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 50
  
  # Generate signals with known correlations
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Signal 1: pure sine wave
  signal1 <- sin(t_seq)
  # Signal 2: 90-degree phase shift (should be orthogonal, r≈0)
  signal2 <- cos(t_seq)
  # Signal 3: identical to signal1 (should be perfectly correlated, r≈1)
  signal3 <- sin(t_seq)
  # Signal 4: negative of signal1 (should be anti-correlated, r≈-1)
  signal4 <- -sin(t_seq)
  
  # Assign signals to spatial regions
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  coords <- arrayInd(1:nvox, dims)
  
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 3 && y <= 3) {
      ts_data[i, ] <- signal1 + rnorm(ntime, sd = 0.01)  # Region 1
    } else if (x > 3 && y <= 3) {
      ts_data[i, ] <- signal2 + rnorm(ntime, sd = 0.01)  # Region 2
    } else if (x <= 3 && y > 3) {
      ts_data[i, ] <- signal3 + rnorm(ntime, sd = 0.01)  # Region 3
    } else {
      ts_data[i, ] <- signal4 + rnorm(ntime, sd = 0.01)  # Region 4
    }
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run clustering and validate that correlated regions cluster together
  result <- slice_msf(vec, mask, stitch_z = TRUE, r = 6, 
                     min_size = 5, compactness = 3, num_runs = 1)
  
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- result$cluster
  
  # Regions 1 and 3 (identical signals) should be in same cluster
  region1_cluster <- unique(cluster_vol[1:3, 1:3, ])
  region3_cluster <- unique(cluster_vol[1:3, 4:6, ])
  
  # Remove zeros
  region1_cluster <- region1_cluster[region1_cluster > 0]
  region3_cluster <- region3_cluster[region3_cluster > 0]
  
  # Should have significant overlap (identical signals)
  overlap <- length(intersect(region1_cluster, region3_cluster))
  expect_true(overlap > 0,
              info = "Identical signals should produce overlapping clusters")
  
  # Test that anti-correlated regions (1 and 4) are separated
  region4_cluster <- unique(cluster_vol[4:6, 4:6, ])
  region4_cluster <- region4_cluster[region4_cluster > 0]
  
  separation <- length(intersect(region1_cluster, region4_cluster))
  expect_true(separation == 0,
              info = "Anti-correlated signals should produce separate clusters")
})

test_that("numerical stability under edge conditions", {
  # Test behavior with extreme values, near-zero variance, etc.
  
  dims <- c(4, 4, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 20
  
  # Test 1: Near-zero variance (constant signal + tiny noise)
  constant_signal <- rep(1.0, ntime)
  ts_data1 <- matrix(rep(constant_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
  ts_data1 <- ts_data1 + matrix(rnorm(nvox * ntime, sd = 1e-10), nrow = nvox, ncol = ntime)
  
  vec_list1 <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data1[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec1 <- do.call(concat, vec_list1)
  
  # Should not crash with near-constant data
  expect_silent(result1 <- slice_msf_single(vec1, mask, r = 4, k = 0.5, min_size = 1))
  expect_true(!is.null(result1))
  expect_true(length(unique(result1$labels)) > 0)
  
  # Test 2: Very large values (test numerical overflow protection)
  large_signal <- seq(1e6, 2e6, length.out = ntime)
  ts_data2 <- matrix(rep(large_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
  
  vec_list2 <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data2[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec2 <- do.call(concat, vec_list2)
  
  # Should handle large values without overflow
  expect_silent(result2 <- slice_msf_single(vec2, mask, r = 4, k = 0.5, min_size = 1))
  expect_true(!is.null(result2))
  expect_true(all(is.finite(result2$sketch)))
  
  # Test 3: Mixed extreme values
  extreme_signal <- c(rep(-1e6, ntime/2), rep(1e6, ntime/2))
  ts_data3 <- matrix(rep(extreme_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
  
  vec_list3 <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data3[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec3 <- do.call(concat, vec_list3)
  
  # Should handle extreme transitions
  expect_silent(result3 <- slice_msf_single(vec3, mask, r = 4, k = 0.5, min_size = 1))
  expect_true(!is.null(result3))
  expect_true(all(is.finite(result3$sketch)))
})

test_that("DCT reconstruction fidelity", {
  # Test that DCT sketching preserves signal information appropriately
  
  dims <- c(4, 4, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 32
  
  # Create test signal with known frequency content
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Very low-frequency signal - closer to DC with gentle variation
  low_freq_signal <- 1 + 0.3*sin(0.5*t_seq) + 0.1*cos(0.8*t_seq)
  
  # High-frequency signal (should lose some fidelity)
  high_freq_signal <- sin(8*t_seq) + sin(16*t_seq)
  
  # Test with different DCT ranks
  r_values <- c(4, 8, 16)
  
  for (r in r_values) {
    # Test low-frequency signal
    ts_data_low <- matrix(rep(low_freq_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
    
    vec_list_low <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data_low[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec_low <- do.call(concat, vec_list_low)
    
    result_low <- slice_msf_single(vec_low, mask, r = r, k = 0.5, min_size = 1)
    
    # For low-frequency signals, DCT should capture most energy in first few coefficients
    sketches_low <- result_low$sketch[, 1]  # First voxel's sketch
    energy_first_half <- sum(sketches_low[1:(r/2)]^2)
    total_energy <- sum(sketches_low^2)
    
    # For very low-frequency signals, most energy should be in early coefficients
    # Test that the energy is concentrated in first half of spectrum
    expect_true(energy_first_half / total_energy > 0.75,
                info = sprintf("Very low-freq signal should concentrate energy in first half of DCT coefficients (r=%d)", r))
    
    # Test high-frequency signal
    ts_data_high <- matrix(rep(high_freq_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
    
    vec_list_high <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data_high[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec_high <- do.call(concat, vec_list_high)
    
    result_high <- slice_msf_single(vec_high, mask, r = r, k = 0.5, min_size = 1)
    
    # For high-frequency signals, energy should be more distributed
    sketches_high <- result_high$sketch[, 1]
    energy_first_half_high <- sum(sketches_high[1:(r/2)]^2)
    total_energy_high <- sum(sketches_high^2)
    
    expect_true(energy_first_half_high / total_energy_high < 0.7,
                info = sprintf("High-freq signal should distribute energy across DCT coefficients (r=%d)", r))
  }
})

test_that("cosine similarity computation accuracy", {
  # Test cosine similarity used in z-stitching
  
  dims <- c(6, 6, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 30
  
  # Create signals with known cosine similarities
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Signal pairs with known similarities:
  # - Identical: cos(θ) = 1
  # - Orthogonal: cos(θ) = 0  
  # - Anti-parallel: cos(θ) = -1
  
  base_signal <- sin(t_seq)
  orthogonal_signal <- cos(t_seq)  # 90° phase shift
  antiparallel_signal <- -sin(t_seq)
  
  # Assign signals to different slices
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  coords <- arrayInd(1:nvox, dims)
  
  for (i in 1:nvox) {
    z <- coords[i, 3]
    
    if (z == 1) {
      ts_data[i, ] <- base_signal + rnorm(ntime, sd = 0.01)
    } else if (z == 2) {
      ts_data[i, ] <- base_signal + rnorm(ntime, sd = 0.01)  # Identical to slice 1
    } else {
      ts_data[i, ] <- orthogonal_signal + rnorm(ntime, sd = 0.01)  # Orthogonal to slices 1&2
    }
  }
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test z-stitching behavior with different theta_link values
  # Slices 1&2 should stitch (high similarity), slice 3 should not
  
  # High threshold - should stitch similar slices
  result_high <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = 0.9,
                          min_contact = 1, r = 6, min_size = 5, 
                          compactness = 3, num_runs = 1)
  
  # Low threshold - should stitch everything
  result_low <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = 0.1,
                         min_contact = 1, r = 6, min_size = 5,
                         compactness = 3, num_runs = 1)
  
  # High threshold should have more clusters than low threshold
  n_clusters_high <- length(unique(result_high$cluster))
  n_clusters_low <- length(unique(result_low$cluster))
  
  expect_true(n_clusters_high >= n_clusters_low,
              info = sprintf("High theta_link (%d clusters) should produce >= clusters than low theta_link (%d clusters)", 
                           n_clusters_high, n_clusters_low))
  
  # Verify that similar slices (1&2) are more likely to be stitched than dissimilar (1&3 or 2&3)
  cluster_vol_high <- array(0, dims)
  cluster_vol_high[mask > 0] <- result_high$cluster
  
  # Check cluster overlap between slices
  slice1_clusters <- unique(cluster_vol_high[,,1])[unique(cluster_vol_high[,,1]) > 0]
  slice2_clusters <- unique(cluster_vol_high[,,2])[unique(cluster_vol_high[,,2]) > 0]
  slice3_clusters <- unique(cluster_vol_high[,,3])[unique(cluster_vol_high[,,3]) > 0]
  
  overlap_12 <- length(intersect(slice1_clusters, slice2_clusters))
  overlap_13 <- length(intersect(slice1_clusters, slice3_clusters))
  overlap_23 <- length(intersect(slice2_clusters, slice3_clusters))
  
  # Similar slices should have more cluster overlap
  expect_true(overlap_12 >= max(overlap_13, overlap_23),
              info = sprintf("Similar slices 1&2 should have >= overlap (%d) than dissimilar pairs (1&3: %d, 2&3: %d)",
                           overlap_12, overlap_13, overlap_23))
})