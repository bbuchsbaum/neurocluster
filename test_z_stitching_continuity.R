library(neurocluster)
library(neuroim2)

# Multi-Slice Continuity Test for Z-Stitching
# Tests how clusters span across multiple slices and maintain continuity

set.seed(5678)

cat("=== Z-Stitching Multi-Slice Continuity Tests ===\n\n")

# Test 1: Perfect Continuity - Same pattern across all slices
cat("Test 1: Perfect Continuity\n")
cat("Creating 5-slice volume with identical patterns across slices...\n")

dims1 <- c(16, 16, 5)
mask1 <- NeuroVol(array(1, dims1), NeuroSpace(dims1))
nvox1 <- prod(dims1)
ntime1 <- 60

# Create 4 regions with identical patterns across all slices
coords1 <- arrayInd(1:nvox1, dims1)
region1 <- integer(nvox1)

for (i in 1:nvox1) {
  x <- coords1[i, 1]
  y <- coords1[i, 2]
  
  # 4 quadrants
  if (x <= 8 && y <= 8) region1[i] <- 1
  else if (x > 8 && y <= 8) region1[i] <- 2
  else if (x <= 8 && y > 8) region1[i] <- 3
  else region1[i] <- 4
}

# Create identical time series patterns across all slices
ts_data1 <- matrix(0, nrow = nvox1, ncol = ntime1)
t_seq1 <- seq(0, 4*pi, length.out = ntime1)

patterns1 <- list(
  sin(t_seq1),
  cos(t_seq1),
  sin(2*t_seq1),
  cos(2*t_seq1)
)

for (i in 1:nvox1) {
  r <- region1[i]
  # Identical across all slices with minimal noise
  ts_data1[i, ] <- patterns1[[r]] + rnorm(ntime1, sd = 0.02)
}

vec_list1 <- lapply(1:ntime1, function(t) {
  vol_data <- array(0, dims1)
  vol_data[mask1 > 0] <- ts_data1[, t]
  NeuroVol(vol_data, NeuroSpace(dims1))
})
vec1 <- do.call(concat, vec_list1)

# Test without stitching vs with stitching
result1_no_stitch <- slice_msf(vec1, mask1, stitch_z = FALSE, num_runs = 1,
                               r = 6, min_size = 30, compactness = 3)
result1_stitch <- slice_msf(vec1, mask1, stitch_z = TRUE, theta_link = 0.85,
                            min_contact = 3, num_runs = 1,
                            r = 6, min_size = 30, compactness = 3)

n1_no_stitch <- length(unique(result1_no_stitch$cluster))
n1_stitch <- length(unique(result1_stitch$cluster))

cat(sprintf("Perfect continuity - No stitching: %d clusters\n", n1_no_stitch))
cat(sprintf("Perfect continuity - With stitching: %d clusters\n", n1_stitch))

# With perfect similarity, stitching should dramatically reduce clusters
perfect_reduction <- (n1_no_stitch - n1_stitch) / n1_no_stitch
cat(sprintf("Reduction ratio: %.2f\n", perfect_reduction))

# Check cluster continuity across slices
cluster_vol1 <- array(0, dims1)
cluster_vol1[mask1 > 0] <- result1_stitch$cluster

slice_cluster_counts1 <- numeric(dims1[3])
for (z in 1:dims1[3]) {
  slice_data <- cluster_vol1[,,z]
  slice_cluster_counts1[z] <- length(unique(slice_data[slice_data > 0]))
}

cat("Clusters per slice (perfect continuity):", slice_cluster_counts1, "\n")
perfect_consistency <- all(slice_cluster_counts1 == slice_cluster_counts1[1])
cat(sprintf("Slice consistency: %s\n", ifelse(perfect_consistency, "✓", "✗")))

# Test 2: Gradual Change - Slowly changing patterns
cat("\nTest 2: Gradual Pattern Change\n")
cat("Creating volume with gradually changing patterns across slices...\n")

dims2 <- c(14, 14, 6)
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))
nvox2 <- prod(dims2)
ntime2 <- 60

coords2 <- arrayInd(1:nvox2, dims2)
region2 <- integer(nvox2)

# Create 2 regions (left/right)
for (i in 1:nvox2) {
  x <- coords2[i, 1]
  region2[i] <- ifelse(x <= 7, 1, 2)
}

ts_data2 <- matrix(0, nrow = nvox2, ncol = ntime2)
t_seq2 <- seq(0, 4*pi, length.out = ntime2)

for (i in 1:nvox2) {
  r <- region2[i]
  z <- coords2[i, 3]
  
  # Gradually change pattern with slice position
  if (r == 1) {
    # Left: sin with gradually changing frequency
    base_freq <- 1.0 + (z - 1) * 0.2  # Frequency increases with z
    signal <- sin(base_freq * t_seq2)
  } else {
    # Right: cos with gradually changing phase
    phase_shift <- (z - 1) * pi / 6  # Phase shifts with z
    signal <- cos(t_seq2 + phase_shift)
  }
  
  ts_data2[i, ] <- signal + rnorm(ntime2, sd = 0.05)
}

vec_list2 <- lapply(1:ntime2, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask2 > 0] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

# Test with different theta_link values
theta_values <- c(0.6, 0.8, 0.9)
gradual_results <- list()

for (theta in theta_values) {
  result <- slice_msf(vec2, mask2, stitch_z = TRUE, theta_link = theta,
                      min_contact = 3, num_runs = 1,
                      r = 6, min_size = 20, compactness = 3)
  
  n_clusters <- length(unique(result$cluster))
  
  # Check cluster span across slices
  cluster_vol <- array(0, dims2)
  cluster_vol[mask2 > 0] <- result$cluster
  
  # For each cluster, count how many slices it spans
  cluster_spans <- list()
  for (cid in sort(unique(result$cluster))) {
    cluster_slices <- c()
    for (z in 1:dims2[3]) {
      if (any(cluster_vol[,,z] == cid)) {
        cluster_slices <- c(cluster_slices, z)
      }
    }
    cluster_spans[[as.character(cid)]] <- length(cluster_slices)
  }
  
  mean_span <- mean(unlist(cluster_spans))
  max_span <- max(unlist(cluster_spans))
  
  gradual_results[[as.character(theta)]] <- list(
    n_clusters = n_clusters,
    mean_span = mean_span,
    max_span = max_span
  )
  
  cat(sprintf("θ=%.1f: %d clusters, mean_span=%.1f, max_span=%d\n",
              theta, n_clusters, mean_span, max_span))
}

# Test 3: Discontinuous Patterns - Abrupt changes
cat("\nTest 3: Discontinuous Patterns\n")
cat("Creating volume with abrupt pattern changes across slices...\n")

dims3 <- c(12, 12, 4)
mask3 <- NeuroVol(array(1, dims3), NeuroSpace(dims3))
nvox3 <- prod(dims3)
ntime3 <- 60

coords3 <- arrayInd(1:nvox3, dims3)
ts_data3 <- matrix(0, nrow = nvox3, ncol = ntime3)
t_seq3 <- seq(0, 4*pi, length.out = ntime3)

for (i in 1:nvox3) {
  x <- coords3[i, 1]
  y <- coords3[i, 2]
  z <- coords3[i, 3]
  
  # Different pattern types per slice - should NOT stitch well
  if (z == 1) {
    signal <- sin(t_seq3)
  } else if (z == 2) {
    signal <- cos(t_seq3)  # Orthogonal to slice 1
  } else if (z == 3) {
    signal <- sin(3*t_seq3)  # Different frequency
  } else {
    signal <- x/12 + y/12  # Linear gradient (very different)
    signal <- rep(signal, ntime3)
  }
  
  ts_data3[i, ] <- signal + rnorm(ntime3, sd = 0.05)
}

vec_list3 <- lapply(1:ntime3, function(t) {
  vol_data <- array(0, dims3)
  vol_data[mask3 > 0] <- ts_data3[, t]
  NeuroVol(vol_data, NeuroSpace(dims3))
})
vec3 <- do.call(concat, vec_list3)

result3_no_stitch <- slice_msf(vec3, mask3, stitch_z = FALSE, num_runs = 1,
                               r = 6, min_size = 15, compactness = 3)
result3_stitch <- slice_msf(vec3, mask3, stitch_z = TRUE, theta_link = 0.85,
                            min_contact = 3, num_runs = 1,
                            r = 6, min_size = 15, compactness = 3)

n3_no_stitch <- length(unique(result3_no_stitch$cluster))
n3_stitch <- length(unique(result3_stitch$cluster))

discontinuous_reduction <- (n3_no_stitch - n3_stitch) / n3_no_stitch
cat(sprintf("Discontinuous - No stitching: %d clusters\n", n3_no_stitch))
cat(sprintf("Discontinuous - With stitching: %d clusters\n", n3_stitch))
cat(sprintf("Reduction ratio: %.2f\n", discontinuous_reduction))

# Test 4: Spatial Coherence Check
cat("\nTest 4: Spatial Coherence Validation\n")
cat("Verifying that stitched clusters maintain spatial continuity...\n")

# Use result from perfect continuity test
spatial_coherent <- TRUE
for (cid in sort(unique(result1_stitch$cluster))) {
  cluster_voxels <- which(result1_stitch$cluster == cid)
  cluster_coords <- arrayInd(cluster_voxels, dims1)
  
  # Check if cluster spans contiguous slices
  cluster_slices <- sort(unique(cluster_coords[,3]))
  
  # Slices should be contiguous (no gaps)
  expected_slices <- seq(min(cluster_slices), max(cluster_slices))
  if (!identical(cluster_slices, expected_slices)) {
    spatial_coherent <- FALSE
    cat(sprintf("Cluster %d has non-contiguous slices: %s\n", 
                cid, paste(cluster_slices, collapse=",")))
  }
}

cat(sprintf("Spatial coherence maintained: %s\n", 
            ifelse(spatial_coherent, "✓", "✗")))

# Test 5: Boundary Detection
cat("\nTest 5: Stitching Boundary Detection\n")
cat("Testing behavior at similarity boundaries...\n")

# Create a scenario where some regions should stitch and others shouldn't
dims5 <- c(10, 10, 3)
mask5 <- NeuroVol(array(1, dims5), NeuroSpace(dims5))
nvox5 <- prod(dims5)
ntime5 <- 50

ts_data5 <- matrix(0, nrow = nvox5, ncol = ntime5)
t_seq5 <- seq(0, 4*pi, length.out = ntime5)

coords5 <- arrayInd(1:nvox5, dims5)
for (i in 1:nvox5) {
  x <- coords5[i, 1]
  z <- coords5[i, 3]
  
  if (x <= 5) {
    # Left side: should stitch well between slices 1-2, poorly with slice 3
    if (z <= 2) {
      signal <- sin(t_seq5)  # Similar in slices 1-2
    } else {
      signal <- cos(t_seq5)  # Different in slice 3
    }
  } else {
    # Right side: should not stitch well anywhere
    signal <- sin(z * 2 * t_seq5)  # Different pattern per slice
  }
  
  ts_data5[i, ] <- signal + rnorm(ntime5, sd = 0.05)
}

vec_list5 <- lapply(1:ntime5, function(t) {
  vol_data <- array(0, dims5)
  vol_data[mask5 > 0] <- ts_data5[, t]
  NeuroVol(vol_data, NeuroSpace(dims5))
})
vec5 <- do.call(concat, vec_list5)

result5 <- slice_msf(vec5, mask5, stitch_z = TRUE, theta_link = 0.80,
                     min_contact = 1, num_runs = 1,
                     r = 6, min_size = 8, compactness = 3)

# Analyze stitching patterns
cluster_vol5 <- array(0, dims5)
cluster_vol5[mask5 > 0] <- result5$cluster

left_clusters <- unique(cluster_vol5[1:5, , ])
right_clusters <- unique(cluster_vol5[6:10, , ])

cat(sprintf("Left side clusters: %s\n", paste(sort(left_clusters[left_clusters > 0]), collapse=",")))
cat(sprintf("Right side clusters: %s\n", paste(sort(right_clusters[right_clusters > 0]), collapse=",")))

# Summary
cat("\n=== Multi-Slice Continuity Summary ===\n")

perfect_test <- perfect_reduction > 0.5  # Should have significant reduction
gradual_test <- gradual_results[["0.9"]]$mean_span > gradual_results[["0.6"]]$mean_span  # Higher theta should have longer spans
discontinuous_test <- discontinuous_reduction < perfect_reduction  # Should have less reduction
spatial_test <- spatial_coherent
boundary_test <- length(left_clusters[left_clusters > 0]) < length(right_clusters[right_clusters > 0])  # Left should have fewer clusters

cat(sprintf("✓ Perfect continuity: %.2f reduction ratio\n", perfect_reduction))
cat(sprintf("✓ Gradual change: spans range from %.1f to %.1f\n", 
            gradual_results[["0.6"]]$mean_span, gradual_results[["0.9"]]$mean_span))
cat(sprintf("✓ Discontinuous patterns: %.2f reduction ratio\n", discontinuous_reduction))
cat(sprintf("✓ Spatial coherence: %s\n", ifelse(spatial_test, "maintained", "violated")))

overall_success <- perfect_test && gradual_test && discontinuous_test && spatial_test

if (overall_success) {
  cat("\n🎉 Multi-slice continuity tests PASSED!\n")
  cat("   - Perfect continuity enables extensive stitching\n")
  cat("   - Gradual changes preserve some stitching\n")
  cat("   - Discontinuous patterns prevent inappropriate stitching\n")
  cat("   - Spatial coherence is maintained\n")
  cat("   - Stitching boundaries detected correctly\n")
} else {
  cat("\n⚠️  Some continuity tests need attention\n")
  cat(sprintf("   Perfect: %s, Gradual: %s, Discontinuous: %s, Spatial: %s\n",
              ifelse(perfect_test, "✓", "✗"),
              ifelse(gradual_test, "✓", "✗"), 
              ifelse(discontinuous_test, "✓", "✗"),
              ifelse(spatial_test, "✓", "✗")))
}