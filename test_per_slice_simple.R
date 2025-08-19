library(neurocluster)
library(neuroim2)

# Simple test for per-slice K
set.seed(100)
dims <- c(8, 8, 3)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 30

# Create simple left/right pattern in each slice
nvox <- prod(dims)
coords <- arrayInd(1:nvox, dims)

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 2*pi, length.out = ntime)

# Each slice has left/right with slightly different patterns
for (i in 1:nvox) {
  z <- coords[i,3]
  x <- coords[i,1]
  
  if (x <= 4) {
    # Left side - different frequency per slice
    ts_data[i, ] <- sin(t_seq * z) + rnorm(ntime, sd = 0.05)
  } else {
    # Right side - different phase per slice
    ts_data[i, ] <- cos(t_seq + z*pi/4) + rnorm(ntime, sd = 0.05)
  }
}

# Create NeuroVec
vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

cat("Test 1: Natural clustering (should find ~2 per slice):\n")
result_natural <- slice_msf(vec, mask, stitch_z = FALSE, num_runs = 1, 
                           r = 4, min_size = 5, compactness = 3)

cluster_vol <- array(0, dims)
cluster_vol[mask > 0] <- result_natural$cluster

for (z in 1:dims[3]) {
  slice_data <- cluster_vol[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("  Slice %d: %d clusters (natural)\n", z, length(slice_clusters)))
}

cat("\nTest 2: Per-slice K = 2:\n")
result_k2 <- slice_msf(vec, mask, target_k_per_slice = 2,
                       stitch_z = FALSE, num_runs = 1, 
                       r = 4, min_size = 5, compactness = 3)

cluster_vol2 <- array(0, dims)
cluster_vol2[mask > 0] <- result_k2$cluster

for (z in 1:dims[3]) {
  slice_data <- cluster_vol2[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("  Slice %d: %d clusters (exact K=2)\n", z, length(slice_clusters)))
}

# Test with more complex pattern
cat("\nTest 3: 4 regions per slice, merge to K=2:\n")
dims2 <- c(10, 10, 2)
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))
nvox2 <- prod(dims2)
coords2 <- arrayInd(1:nvox2, dims2)

# Create 4 quadrants
ts_data2 <- matrix(0, nrow = nvox2, ncol = ntime)

for (i in 1:nvox2) {
  x <- coords2[i,1]
  y <- coords2[i,2]
  
  if (x <= 5 && y <= 5) {
    # Top-left
    ts_data2[i, ] <- sin(t_seq) + rnorm(ntime, sd = 0.05)
  } else if (x > 5 && y <= 5) {
    # Top-right
    ts_data2[i, ] <- cos(t_seq) + rnorm(ntime, sd = 0.05)
  } else if (x <= 5 && y > 5) {
    # Bottom-left
    ts_data2[i, ] <- -sin(t_seq) + rnorm(ntime, sd = 0.05)
  } else {
    # Bottom-right
    ts_data2[i, ] <- -cos(t_seq) + rnorm(ntime, sd = 0.05)
  }
}

vec_list2 <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask2 > 0] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

# Natural clustering
result_nat2 <- slice_msf(vec2, mask2, stitch_z = FALSE, num_runs = 1,
                         r = 4, min_size = 10, compactness = 3)

cat("  Natural clustering: ", length(unique(result_nat2$cluster)), " total clusters\n")

# Force K=2 per slice
result_k2_2 <- slice_msf(vec2, mask2, target_k_per_slice = 2,
                         stitch_z = FALSE, num_runs = 1,
                         r = 4, min_size = 10, compactness = 3)

cluster_vol3 <- array(0, dims2)
cluster_vol3[mask2 > 0] <- result_k2_2$cluster

for (z in 1:dims2[3]) {
  slice_data <- cluster_vol3[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("  Slice %d: %d clusters\n", z, length(slice_clusters)))
}