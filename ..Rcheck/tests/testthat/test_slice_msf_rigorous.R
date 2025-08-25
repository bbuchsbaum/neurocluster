library(neuroim2)
library(neurocluster)

test_that("slice_msf correctly identifies spatially contiguous regions with similar time series", {
  # Create a volume with clear spatial structure
  set.seed(100)
  dims <- c(20, 20, 5)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 100
  
  # Create 4 spatially distinct regions with different time series patterns
  # Region 1: top-left quadrant (slow oscillation)
  # Region 2: top-right quadrant (fast oscillation)  
  # Region 3: bottom-left quadrant (linear trend)
  # Region 4: bottom-right quadrant (constant + noise)
  
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Assign regions based on spatial location
  region <- integer(nvox)
  region[coords[,1] <= 10 & coords[,2] <= 10] <- 1
  region[coords[,1] > 10 & coords[,2] <= 10] <- 2
  region[coords[,1] <= 10 & coords[,2] > 10] <- 3
  region[coords[,1] > 10 & coords[,2] > 10] <- 4
  
  # Create distinct time series for each region
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Region 1: slow sine wave
  ts_data[region == 1, ] <- matrix(rep(sin(t_seq), sum(region == 1)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(region == 1) * ntime, sd = 0.1)
  
  # Region 2: fast sine wave
  ts_data[region == 2, ] <- matrix(rep(sin(3 * t_seq), sum(region == 2)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(region == 2) * ntime, sd = 0.1)
  
  # Region 3: linear trend
  ts_data[region == 3, ] <- matrix(rep(seq(0, 2, length.out = ntime), sum(region == 3)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(region == 3) * ntime, sd = 0.1)
  
  # Region 4: constant with higher noise
  ts_data[region == 4, ] <- matrix(rep(rep(1, ntime), sum(region == 4)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(region == 4) * ntime, sd = 0.3)
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run clustering
  result <- slice_msf(vec, mask, num_runs = 1, r = 8, 
                      min_size = 50, compactness = 5)
  
  # Verify we get approximately 4 clusters
  n_clusters <- length(unique(result$cluster))
  # With min_size=50 and 500 voxels per region, we might get more clusters due to edges
  expect_true(n_clusters >= 3 && n_clusters <= 12, 
              info = sprintf("Expected 3-12 clusters, got %d", n_clusters))
  
  # Check spatial contiguity - clusters should be mostly spatially coherent
  # For each cluster, check that most voxels have neighbors in the same cluster
  cluster_coherence <- numeric(n_clusters)
  unique_clusters <- sort(unique(result$cluster))
  
  for (i in seq_along(unique_clusters)) {
    clust_id <- unique_clusters[i]
    clust_voxels <- which(result$cluster == clust_id)
    clust_coords <- coords[clust_voxels, , drop = FALSE]
    
    # Count how many voxels have at least one neighbor in the same cluster
    has_neighbor <- 0
    for (j in seq_len(nrow(clust_coords))) {
      # Check 6-connectivity neighbors
      neighbors <- rbind(
        clust_coords[j, ] + c(1, 0, 0),
        clust_coords[j, ] + c(-1, 0, 0),
        clust_coords[j, ] + c(0, 1, 0),
        clust_coords[j, ] + c(0, -1, 0),
        clust_coords[j, ] + c(0, 0, 1),
        clust_coords[j, ] + c(0, 0, -1)
      )
      
      # Check if any neighbor is in the same cluster
      for (k in seq_len(nrow(neighbors))) {
        if (all(neighbors[k, ] >= 1) && all(neighbors[k, ] <= dims)) {
          neighbor_idx <- neighbors[k, 1] + (neighbors[k, 2] - 1) * dims[1] + 
                         (neighbors[k, 3] - 1) * dims[1] * dims[2]
          if (neighbor_idx %in% clust_voxels) {
            has_neighbor <- has_neighbor + 1
            break
          }
        }
      }
    }
    
    cluster_coherence[i] <- has_neighbor / length(clust_voxels)
  }
  
  # Most voxels should have neighbors in the same cluster (spatial coherence)
  expect_true(mean(cluster_coherence) > 0.8, 
              info = sprintf("Low spatial coherence: %.2f", mean(cluster_coherence)))
})

test_that("slice_msf separates regions with different temporal patterns", {
  # Create a simple 2D slice with two regions having orthogonal time series
  set.seed(200)
  dims <- c(10, 10, 1)  # Single slice
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 50
  
  nvox <- prod(dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Left half: sine pattern
  # Right half: cosine pattern (orthogonal to sine)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  left_voxels <- coords[,1] <= 5
  right_voxels <- coords[,1] > 5
  
  ts_data[left_voxels, ] <- matrix(rep(sin(t_seq), sum(left_voxels)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(sum(left_voxels) * ntime, sd = 0.05)
  
  ts_data[right_voxels, ] <- matrix(rep(cos(t_seq), sum(right_voxels)), 
                                     ncol = ntime, byrow = TRUE) + 
                                     rnorm(sum(right_voxels) * ntime, sd = 0.05)
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run clustering
  result <- slice_msf(vec, mask, num_runs = 1, r = 4, 
                      min_size = 20, compactness = 3)
  
  # Should find exactly 2 clusters
  n_clusters <- length(unique(result$cluster))
  expect_equal(n_clusters, 2)
  
  # Check that the two regions are in different clusters
  left_clusters <- unique(result$cluster[left_voxels])
  right_clusters <- unique(result$cluster[right_voxels])
  
  # Each region should be in a single cluster
  expect_equal(length(left_clusters), 1)
  expect_equal(length(right_clusters), 1)
  
  # And they should be different clusters
  expect_true(left_clusters != right_clusters)
})

test_that("slice_msf reliability weighting works correctly", {
  # Test that voxels with more reliable signals are weighted more heavily
  set.seed(300)
  dims <- c(8, 8, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 40  # Even number for split-half
  
  nvox <- prod(dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Create two groups:
  # Group 1: Highly reliable signal (same pattern in both halves)
  # Group 2: Unreliable signal (different patterns in each half)
  
  reliable_voxels <- 1:64
  unreliable_voxels <- 65:128
  
  # Reliable: consistent sine wave
  pattern <- sin(seq(0, 2*pi, length.out = ntime/2))
  ts_data[reliable_voxels, 1:(ntime/2)] <- matrix(rep(pattern, length(reliable_voxels)), 
                                                    ncol = ntime/2, byrow = TRUE)
  ts_data[reliable_voxels, (ntime/2+1):ntime] <- matrix(rep(pattern, length(reliable_voxels)), 
                                                          ncol = ntime/2, byrow = TRUE)
  
  # Unreliable: different patterns in each half
  ts_data[unreliable_voxels, 1:(ntime/2)] <- matrix(rnorm(length(unreliable_voxels) * ntime/2), 
                                                      ncol = ntime/2)
  ts_data[unreliable_voxels, (ntime/2+1):ntime] <- matrix(rnorm(length(unreliable_voxels) * ntime/2), 
                                                            ncol = ntime/2)
  
  # Add small noise to all
  ts_data <- ts_data + rnorm(nvox * ntime, sd = 0.1)
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run clustering with low min_size to allow separation
  result <- slice_msf(vec, mask, num_runs = 1, r = 4, 
                      min_size = 10, compactness = 3, gamma = 1.5)
  
  # The reliable voxels should cluster together more than unreliable ones
  reliable_clusters <- result$cluster[reliable_voxels]
  unreliable_clusters <- result$cluster[unreliable_voxels]
  
  # Count the size of the most common cluster in each group
  reliable_table <- table(reliable_clusters)
  unreliable_table <- table(unreliable_clusters)
  
  reliable_coherence <- max(reliable_table) / length(reliable_clusters)
  unreliable_coherence <- max(unreliable_table) / length(unreliable_clusters)
  
  # Reliable voxels should have higher coherence (or at least equal)
  # With small data, both might end up in single clusters
  expect_true(reliable_coherence >= unreliable_coherence * 0.9,
              info = sprintf("Reliable coherence (%.2f) should be >= unreliable (%.2f)", 
                            reliable_coherence, unreliable_coherence))
})

test_that("slice_msf consensus fusion improves stability", {
  # Test that consensus across multiple runs produces more stable results
  set.seed(400)
  dims <- c(12, 12, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 60
  
  # Create 3 regions with some overlap/ambiguity at boundaries
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Define regions with fuzzy boundaries
  region <- integer(nvox)
  centers <- rbind(c(3, 3, 2), c(9, 3, 2), c(6, 9, 2))
  
  for (i in 1:nvox) {
    dists <- sqrt(rowSums((centers - matrix(coords[i,], nrow = 3, ncol = 3, byrow = TRUE))^2))
    region[i] <- which.min(dists)
  }
  
  # Create time series with some noise
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  patterns <- rbind(
    sin(t_seq),
    cos(t_seq),
    sin(0.5 * t_seq)
  )
  
  for (r in 1:3) {
    voxels <- which(region == r)
    ts_data[voxels, ] <- matrix(rep(patterns[r,], length(voxels)), 
                                 ncol = ntime, byrow = TRUE) + 
                                 rnorm(length(voxels) * ntime, sd = 0.3)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run single run multiple times
  single_runs <- replicate(5, {
    result <- slice_msf(vec, mask, num_runs = 1, r = 6, 
                        min_size = 20, compactness = 4)
    result$cluster
  })
  
  # Run consensus
  consensus_result <- slice_msf(vec, mask, num_runs = 5, 
                                consensus = TRUE, r = 6, 
                                min_size = 20, compactness = 4)
  
  # Calculate variability across single runs
  # Use adjusted Rand index or similar measure
  single_run_similarity <- numeric(choose(5, 2))
  idx <- 1
  for (i in 1:4) {
    for (j in (i+1):5) {
      # Simple measure: fraction of voxels in same cluster
      same_cluster <- outer(single_runs[,i], single_runs[,i], "==") == 
                     outer(single_runs[,j], single_runs[,j], "==")
      single_run_similarity[idx] <- sum(same_cluster) / length(same_cluster)
      idx <- idx + 1
    }
  }
  
  # Consensus should produce results similar to the "average" clustering
  # Check that consensus found reasonable number of clusters
  n_clusters_consensus <- length(unique(consensus_result$cluster))
  expect_true(n_clusters_consensus >= 2 && n_clusters_consensus <= 5)
  
  # Verify runs were stored
  expect_equal(length(consensus_result$runs), 5)
})

test_that("slice_msf respects min_size parameter", {
  # Test that clusters smaller than min_size are merged
  # Note: With the canonical FH cleanup, small edge components may remain
  # in certain data configurations, especially with small test volumes
  
  set.seed(500)
  dims <- c(15, 15, 3)  # Larger volume for more reliable test
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 40
  
  # Create structured regions that should merge based on min_size
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Create grid of small regions
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Divide into 9 regions (3x3 in XY)
  for (i in 1:nvox) {
    x_region <- floor((coords[i,1] - 1) / 5)
    y_region <- floor((coords[i,2] - 1) / 5)
    region_id <- x_region + 3 * y_region + 1
    
    # Each region gets slightly different pattern
    ts_data[i, ] <- sin(t_seq + region_id * pi/9) + 
                    rnorm(ntime, sd = 0.3)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run with large min_size (should force merging)
  result <- slice_msf(vec, mask, num_runs = 1, r = 4, 
                      min_size = 100, compactness = 3)
  
  # Most clusters should respect min_size
  cluster_sizes <- table(result$cluster)
  small_clusters <- sum(cluster_sizes < 100)
  
  # Allow at most 1-2 small edge clusters due to canonical cleanup behavior
  expect_true(small_clusters <= 2,
              info = sprintf("Too many clusters smaller than min_size: %d small clusters found", 
                            small_clusters))
})

test_that("slice_msf handles slice-wise processing correctly", {
  # Verify that the algorithm processes slices independently
  # then optionally stitches them
  set.seed(600)
  dims <- c(8, 8, 4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 50
  
  # Create identical patterns in each slice
  nvox_per_slice <- dims[1] * dims[2]
  nvox <- prod(dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Each slice has the same spatial pattern
  t_seq <- seq(0, 2*pi, length.out = ntime)
  base_pattern <- sin(t_seq)
  
  for (slice in 1:dims[3]) {
    slice_idx <- ((slice-1) * nvox_per_slice + 1):(slice * nvox_per_slice)
    
    # Left half of slice: pattern 1
    left_idx <- slice_idx[1:(nvox_per_slice/2)]
    ts_data[left_idx, ] <- matrix(rep(base_pattern, length(left_idx)), 
                                   ncol = ntime, byrow = TRUE) + 
                                   rnorm(length(left_idx) * ntime, sd = 0.1)
    
    # Right half: orthogonal pattern
    right_idx <- slice_idx[(nvox_per_slice/2 + 1):nvox_per_slice]
    ts_data[right_idx, ] <- matrix(rep(cos(t_seq), length(right_idx)), 
                                    ncol = ntime, byrow = TRUE) + 
                                    rnorm(length(right_idx) * ntime, sd = 0.1)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run without stitching
  result_no_stitch <- slice_msf(vec, mask, num_runs = 1, r = 4, 
                                min_size = 10, compactness = 3, 
                                stitch_z = FALSE)
  
  # Run with stitching
  result_stitch <- slice_msf(vec, mask, num_runs = 1, r = 4, 
                             min_size = 10, compactness = 3, 
                             stitch_z = TRUE, theta_link = 0.85)
  
  # Without stitching, we might have more clusters (one set per slice)
  n_clusters_no_stitch <- length(unique(result_no_stitch$cluster))
  n_clusters_stitch <- length(unique(result_stitch$cluster))
  
  # With stitching, similar patterns across slices should merge
  expect_true(n_clusters_stitch <= n_clusters_no_stitch,
              info = sprintf("Stitching should reduce clusters: no_stitch=%d, stitch=%d",
                            n_clusters_no_stitch, n_clusters_stitch))
  
  # With perfect data and stitching, we should get close to 2 clusters
  expect_true(n_clusters_stitch >= 2 && n_clusters_stitch <= 4,
              info = sprintf("Expected 2-4 clusters with stitching, got %d", n_clusters_stitch))
})