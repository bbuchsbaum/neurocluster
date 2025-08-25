library(neuroim2)
library(neurocluster)

test_that("slice_msf achieves exact K with target_k_global", {
  # Create test data with clear structure
  set.seed(100)
  dims <- c(12, 12, 4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 50
  
  # Create 6 distinct regions with different patterns
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Divide into 6 regions (3x2 grid in XY)
  region <- integer(nvox)
  region[coords[,1] <= 4 & coords[,2] <= 6] <- 1
  region[coords[,1] > 4 & coords[,1] <= 8 & coords[,2] <= 6] <- 2
  region[coords[,1] > 8 & coords[,2] <= 6] <- 3
  region[coords[,1] <= 4 & coords[,2] > 6] <- 4
  region[coords[,1] > 4 & coords[,1] <= 8 & coords[,2] > 6] <- 5
  region[coords[,1] > 8 & coords[,2] > 6] <- 6
  
  # Create distinct time series for each region
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Different patterns for each region
  patterns <- rbind(
    sin(t_seq),                    # Region 1
    cos(t_seq),                    # Region 2
    sin(0.5 * t_seq),              # Region 3
    cos(0.5 * t_seq),              # Region 4
    sin(2 * t_seq),                # Region 5
    0.5 * sin(t_seq) + 0.5 * cos(t_seq)  # Region 6
  )
  
  for (r in 1:6) {
    voxels <- which(region == r)
    ts_data[voxels, ] <- matrix(rep(patterns[r,], length(voxels)), 
                                 ncol = ntime, byrow = TRUE) + 
                                 rnorm(length(voxels) * ntime, sd = 0.1)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test exact K = 4 (should merge some regions)
  result <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1, 
                      r = 8, min_size = 20, compactness = 4)
  
  n_clusters <- length(unique(result$cluster))
  expect_equal(n_clusters, 4, 
               info = "Should have exactly 4 clusters with target_k_global=4")
  
  # Test exact K = 2 (should merge more regions)
  result2 <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                       r = 8, min_size = 20, compactness = 4)
  
  n_clusters2 <- length(unique(result2$cluster))
  expect_equal(n_clusters2, 2,
               info = "Should have exactly 2 clusters with target_k_global=2")
})

test_that("slice_msf achieves exact K per slice with target_k_per_slice", {
  # Create test data with structure in each slice
  set.seed(200)
  dims <- c(10, 10, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 40
  
  nvox <- prod(dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Each slice has left/right regions with different patterns
  for (z in 1:dims[3]) {
    for (y in 1:dims[2]) {
      for (x in 1:dims[1]) {
        idx <- x + (y-1)*dims[1] + (z-1)*dims[1]*dims[2]
        if (x <= 5) {
          # Left: sine wave
          ts_data[idx, ] <- sin(t_seq + z*pi/4) + rnorm(ntime, sd = 0.15)
        } else {
          # Right: cosine wave
          ts_data[idx, ] <- cos(t_seq + z*pi/4) + rnorm(ntime, sd = 0.15)
        }
      }
    }
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test exact K = 2 per slice (should preserve left/right in each slice)
  result <- slice_msf(vec, mask, target_k_per_slice = 2, 
                      stitch_z = FALSE,  # Important: no z-stitching
                      num_runs = 1, r = 6, min_size = 10, compactness = 3)
  
  # Check each slice has exactly 2 clusters
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- result$cluster
  
  for (z in 1:dims[3]) {
    slice_data <- cluster_vol[,,z]
    slice_clusters <- unique(slice_data[slice_data > 0])
    expect_equal(length(slice_clusters), 2,
                 info = sprintf("Slice %d should have exactly 2 clusters", z))
  }
})

test_that("slice_msf exact K preserves spatial contiguity", {
  # Create data that would naturally form many small clusters
  set.seed(300)
  dims <- c(15, 15, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 30
  
  nvox <- prod(dims)
  # Add structured noise that creates many natural clusters
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  # Add some local spatial structure
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  for (i in 1:nvox) {
    # Add spatial smoothing based on coordinates
    spatial_factor <- sin(coords[i,1] * 0.5) * cos(coords[i,2] * 0.5)
    ts_data[i, ] <- ts_data[i, ] + spatial_factor
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run with target K = 3
  result <- slice_msf(vec, mask, target_k_global = 3, num_runs = 1,
                      r = 4, min_size = 20, compactness = 5)
  
  n_clusters <- length(unique(result$cluster))
  expect_equal(n_clusters, 3)
  
  # Check spatial contiguity - each cluster should be connected
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- result$cluster
  
  for (clust_id in unique(result$cluster)) {
    # Get all voxels in this cluster
    clust_mask <- (cluster_vol == clust_id)
    
    # Use flood fill to check connectivity
    # Find first voxel
    first_idx <- which(clust_mask, arr.ind = TRUE)[1,]
    visited <- array(FALSE, dims)
    
    # Flood fill from first voxel
    to_visit <- list(first_idx)
    n_visited <- 0
    
    while (length(to_visit) > 0) {
      current <- to_visit[[1]]
      to_visit <- to_visit[-1]
      
      if (visited[current[1], current[2], current[3]]) next
      if (!clust_mask[current[1], current[2], current[3]]) next
      
      visited[current[1], current[2], current[3]] <- TRUE
      n_visited <- n_visited + 1
      
      # Add 6-connected neighbors (and 8-connected in-plane)
      neighbors <- list(
        c(-1, 0, 0), c(1, 0, 0),  # X neighbors
        c(0, -1, 0), c(0, 1, 0),  # Y neighbors
        c(0, 0, -1), c(0, 0, 1),  # Z neighbors
        c(-1, -1, 0), c(1, 1, 0), # Diagonal in XY
        c(-1, 1, 0), c(1, -1, 0)  # Diagonal in XY
      )
      
      for (offset in neighbors) {
        nx <- current[1] + offset[1]
        ny <- current[2] + offset[2]
        nz <- current[3] + offset[3]
        if (nx >= 1 && nx <= dims[1] && 
            ny >= 1 && ny <= dims[2] && 
            nz >= 1 && nz <= dims[3]) {
          to_visit <- append(to_visit, list(c(nx, ny, nz)))
        }
      }
    }
    
    # All voxels in cluster should be visited
    expect_equal(n_visited, sum(clust_mask),
                 info = sprintf("Cluster %d should be spatially contiguous", clust_id))
  }
})

test_that("slice_msf consensus with exact K requires use_features", {
  # Small test data
  mask <- NeuroVol(array(1, c(6, 6, 2)), NeuroSpace(c(6, 6, 2)))
  nvox <- sum(mask > 0)
  ntime <- 20
  
  ts_data <- matrix(rnorm(nvox * ntime), nrow = nvox, ncol = ntime)
  
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dim = c(6, 6, 2))
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(c(6, 6, 2)))
  })
  vec <- do.call(concat, vec_list)
  
  # Should work without exact K
  result1 <- slice_msf(vec, mask, num_runs = 2, consensus = TRUE,
                       r = 4, min_size = 5, use_features = FALSE)
  expect_s3_class(result1, "slice_msf_cluster_result")
  
  # Should fail with exact K but no features
  expect_error(
    slice_msf(vec, mask, target_k_global = 3, num_runs = 2, 
              consensus = TRUE, r = 4, min_size = 5, use_features = FALSE),
    "Exact-K in consensus requires use_features=TRUE"
  )
  
  # Should work with exact K and features
  result2 <- slice_msf(vec, mask, target_k_global = 3, num_runs = 2,
                       consensus = TRUE, r = 4, min_size = 5, use_features = TRUE)
  expect_equal(length(unique(result2$cluster)), 3)
})