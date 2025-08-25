library(neuroim2)
library(neurocluster)

test_that("slice_msf perfectly recovers structured clusters with exact K", {
  # Create a 30x30x10 volume with 9 distinct regions (3x3 grid pattern)
  set.seed(1000)
  dims <- c(30, 30, 10)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 100
  
  # Define 9 regions in a 3x3 grid pattern
  # Each region is 10x10 voxels in XY, extending through all Z slices
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Assign each voxel to one of 9 regions based on its XY coordinates
  region <- integer(nvox)
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    # Determine grid position (0-2 for both x and y)
    grid_x <- min(floor((x - 1) / 10), 2)
    grid_y <- min(floor((y - 1) / 10), 2)
    
    # Convert to region ID (1-9)
    region[i] <- grid_x + 3 * grid_y + 1
  }
  
  # Create 9 distinct temporal patterns
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Define unique patterns for each region - make them MORE distinct
  # Use orthogonal patterns to ensure natural clustering finds all 9
  patterns <- rbind(
    sin(t_seq),                           # Region 1
    cos(t_seq),                           # Region 2
    sin(2 * t_seq),                       # Region 3
    cos(2 * t_seq),                       # Region 4
    sin(3 * t_seq),                       # Region 5
    cos(3 * t_seq),                       # Region 6
    sin(0.5 * t_seq),                     # Region 7
    cos(0.5 * t_seq),                     # Region 8
    sin(t_seq) * cos(2 * t_seq)           # Region 9: complex pattern
  )
  
  # Assign patterns to regions with minimal noise for near-perfect separation
  for (r in 1:9) {
    voxels <- which(region == r)
    # Very low noise to ensure clear separation
    ts_data[voxels, ] <- matrix(rep(patterns[r,], length(voxels)), 
                                 ncol = ntime, byrow = TRUE) + 
                                 rnorm(length(voxels) * ntime, sd = 0.01)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test 1: Natural clustering should find approximately 9 clusters
  result_natural <- slice_msf(vec, mask, num_runs = 1, r = 12, 
                              min_size = 50, compactness = 5)
  n_natural <- length(unique(result_natural$cluster))
  cat(sprintf("Natural clustering found %d clusters\n", n_natural))
  
  # Test 2: Exact K = 9 should perfectly recover the 9 regions
  result_k9 <- slice_msf(vec, mask, target_k_global = 9, num_runs = 1, 
                         r = 12, min_size = 50, compactness = 5)
  
  expect_equal(length(unique(result_k9$cluster)), 9,
               info = "Should have exactly 9 clusters")
  
  # Check that each true region maps predominantly to one cluster
  # Calculate adjusted Rand index or similar measure
  cluster_mapping <- matrix(0, nrow = 9, ncol = 9)
  
  for (true_region in 1:9) {
    true_voxels <- which(region == true_region)
    for (found_cluster in 1:9) {
      found_voxels <- which(result_k9$cluster == found_cluster)
      cluster_mapping[true_region, found_cluster] <- length(intersect(true_voxels, found_voxels))
    }
  }
  
  # Each true region should map predominantly to one cluster
  # Check that the mapping is nearly one-to-one
  for (true_region in 1:9) {
    region_voxels <- sum(cluster_mapping[true_region, ])
    max_overlap <- max(cluster_mapping[true_region, ])
    purity <- max_overlap / region_voxels
    
    expect_true(purity > 0.95,
                info = sprintf("Region %d should map predominantly to one cluster (purity: %.2f)", 
                              true_region, purity))
  }
  
  # Test 3: Different K values should merge regions sensibly
  result_k5 <- slice_msf(vec, mask, target_k_global = 5, num_runs = 1,
                         r = 12, min_size = 50, compactness = 5)
  expect_equal(length(unique(result_k5$cluster)), 5)
  
  result_k3 <- slice_msf(vec, mask, target_k_global = 3, num_runs = 1,
                         r = 12, min_size = 50, compactness = 5)
  expect_equal(length(unique(result_k3$cluster)), 3)
  
  # Test 4: Per-slice exact K
  # With 10 slices and K=2 per slice, we might get up to 20 clusters
  result_per_slice <- slice_msf(vec, mask, target_k_per_slice = 2,
                                stitch_z = FALSE,  # No stitching
                                num_runs = 1, r = 12, min_size = 20, 
                                compactness = 5)
  
  # Check each slice has exactly 2 clusters
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- result_per_slice$cluster
  
  for (z in 1:dims[3]) {
    slice_data <- cluster_vol[,,z]
    slice_clusters <- unique(slice_data[slice_data > 0])
    expect_equal(length(slice_clusters), 2,
                 info = sprintf("Slice %d should have exactly 2 clusters", z))
  }
})

test_that("slice_msf exact K preserves spatial structure", {
  # Create a simpler 2x2 grid pattern for clearer testing
  set.seed(2000)
  dims <- c(20, 20, 5)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 80
  
  # Define 4 quadrants
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  region <- integer(nvox)
  region[coords[,1] <= 10 & coords[,2] <= 10] <- 1  # Top-left
  region[coords[,1] > 10 & coords[,2] <= 10] <- 2   # Top-right
  region[coords[,1] <= 10 & coords[,2] > 10] <- 3   # Bottom-left
  region[coords[,1] > 10 & coords[,2] > 10] <- 4    # Bottom-right
  
  # Create very distinct patterns
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  patterns <- rbind(
    sin(t_seq),          # Region 1
    -sin(t_seq),         # Region 2 (opposite of 1)
    cos(t_seq),          # Region 3
    -cos(t_seq)          # Region 4 (opposite of 3)
  )
  
  for (r in 1:4) {
    voxels <- which(region == r)
    ts_data[voxels, ] <- matrix(rep(patterns[r,], length(voxels)), 
                                 ncol = ntime, byrow = TRUE) + 
                                 rnorm(length(voxels) * ntime, sd = 0.02)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test exact K = 4
  result <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1,
                      r = 8, min_size = 30, compactness = 5)
  
  expect_equal(length(unique(result$cluster)), 4)
  
  # Verify spatial coherence - each cluster should be mostly in one quadrant
  cluster_ids <- sort(unique(result$cluster))
  
  for (cid in cluster_ids) {
    cluster_voxels <- which(result$cluster == cid)
    cluster_regions <- region[cluster_voxels]
    
    # Most voxels in this cluster should belong to the same true region
    region_table <- table(cluster_regions)
    dominant_region <- as.numeric(names(region_table)[which.max(region_table)])
    coherence <- max(region_table) / sum(region_table)
    
    expect_true(coherence > 0.90,
                info = sprintf("Cluster %d should be spatially coherent (%.2f)", cid, coherence))
  }
  
  # Test K = 2 - should merge opposite corners
  result_k2 <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                         r = 8, min_size = 30, compactness = 5)
  
  expect_equal(length(unique(result_k2$cluster)), 2)
  
  # The merging should be sensible based on pattern similarity
  # Regions 1&2 have opposite patterns (sin vs -sin)
  # Regions 3&4 have opposite patterns (cos vs -cos)
  # So we might see (1,3) vs (2,4) or (1,4) vs (2,3) groupings
})

test_that("slice_msf exact K works with consensus fusion", {
  # Small test for consensus with exact K
  set.seed(3000)
  dims <- c(12, 12, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  ntime <- 40
  
  # Create 3 clear regions
  nvox <- prod(dims)
  mask_idx <- which(mask > 0)
  coords <- arrayInd(mask_idx, dims)
  
  # Three horizontal bands
  region <- integer(nvox)
  region[coords[,2] <= 4] <- 1
  region[coords[,2] > 4 & coords[,2] <= 8] <- 2
  region[coords[,2] > 8] <- 3
  
  # Distinct patterns
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  for (r in 1:3) {
    voxels <- which(region == r)
    if (r == 1) {
      pattern <- sin(t_seq)
    } else if (r == 2) {
      pattern <- cos(t_seq)
    } else {
      pattern <- sin(2 * t_seq)
    }
    
    ts_data[voxels, ] <- matrix(rep(pattern, length(voxels)), 
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
  
  # Test consensus with exact K = 3
  result <- slice_msf(vec, mask, target_k_global = 3, num_runs = 3,
                      consensus = TRUE, use_features = TRUE,
                      r = 6, min_size = 20, compactness = 4)
  
  expect_equal(length(unique(result$cluster)), 3,
               info = "Consensus should produce exactly 3 clusters")
  
  # Verify the clusters correspond well to true regions
  for (true_region in 1:3) {
    true_voxels <- which(region == true_region)
    
    # Find which cluster best matches this region
    overlaps <- numeric(3)
    for (found_cluster in sort(unique(result$cluster))) {
      found_voxels <- which(result$cluster == found_cluster)
      overlaps[found_cluster] <- length(intersect(true_voxels, found_voxels))
    }
    
    # The best matching cluster should have high overlap
    best_overlap <- max(overlaps)
    total_region_voxels <- length(true_voxels)
    recovery_rate <- best_overlap / total_region_voxels
    
    expect_true(recovery_rate > 0.85,
                info = sprintf("Region %d should be well recovered (%.2f)", 
                              true_region, recovery_rate))
  }
})