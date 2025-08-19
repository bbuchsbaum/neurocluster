library(neurocluster)
library(neuroim2)

# Test case where we need to merge clusters (exact K < natural clusters)
set.seed(100)
dims <- c(15, 15, 3)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 50

# Create 6 distinct regions that should naturally cluster separately
nvox <- prod(dims)
mask_idx <- which(mask > 0)
coords <- arrayInd(mask_idx, dims)

# Create a 3x2 grid of regions (6 total)
region <- integer(nvox)
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  grid_x <- min(floor((x - 1) / 5), 2)  # 3 columns
  grid_y <- min(floor((y - 1) / 7.5), 1)  # 2 rows
  
  region[i] <- grid_x + 3 * grid_y + 1
}

cat("Region distribution:\n")
print(table(region))

# Create very distinct patterns to ensure they cluster separately
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)

# Make patterns orthogonal and distinct
patterns <- rbind(
  sin(t_seq),                    # Region 1
  cos(t_seq),                    # Region 2
  sin(2 * t_seq),                # Region 3
  cos(2 * t_seq),                # Region 4
  sin(3 * t_seq),                # Region 5
  cos(3 * t_seq)                 # Region 6
)

# Very low noise to ensure separation
for (r in 1:6) {
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

# Test natural clustering with low compactness to get more clusters
cat("\n1. Natural clustering (low compactness for more clusters):\n")
result_natural <- slice_msf(vec, mask, num_runs = 1, r = 6, 
                           min_size = 20, compactness = 0.5)
n_natural <- length(unique(result_natural$cluster))
cat(sprintf("  Found %d clusters\n", n_natural))

# Test exact K = 3 (should merge to 3 clusters)
cat("\n2. Exact K = 3 (should merge 6 regions to 3 clusters):\n")
result_k3 <- slice_msf(vec, mask, target_k_global = 3, num_runs = 1,
                       r = 6, min_size = 20, compactness = 0.5)
n_k3 <- length(unique(result_k3$cluster))
cat(sprintf("  Found %d clusters\n", n_k3))

if (n_k3 == 3) {
  cat("  SUCCESS: Exact K = 3 working correctly!\n")
  
  # Check which regions were merged
  cat("\n  Region to cluster mapping:\n")
  for (reg in 1:6) {
    reg_voxels <- which(region == reg)
    reg_clusters <- unique(result_k3$cluster[reg_voxels])
    cat(sprintf("    Region %d -> Cluster %s\n", reg, 
                paste(reg_clusters, collapse=", ")))
  }
} else {
  cat("  FAILED: Did not achieve target K\n")
}

# Test exact K = 2
cat("\n3. Exact K = 2 (should merge to 2 clusters):\n")
result_k2 <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                       r = 6, min_size = 20, compactness = 0.5)
n_k2 <- length(unique(result_k2$cluster))
cat(sprintf("  Found %d clusters\n", n_k2))

# Test with per-slice K on multi-slice data
cat("\n4. Per-slice K = 2 (each slice should have 2 clusters):\n")
result_per_slice <- slice_msf(vec, mask, target_k_per_slice = 2,
                              stitch_z = FALSE,
                              num_runs = 1, r = 6, min_size = 10, 
                              compactness = 0.5)

cluster_vol <- array(0, dims)
cluster_vol[mask > 0] <- result_per_slice$cluster

for (z in 1:dims[3]) {
  slice_data <- cluster_vol[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("  Slice %d: %d clusters\n", z, length(slice_clusters)))
}

# Test rows_are_time = FALSE
cat("\n5. Testing rows_are_time = FALSE:\n")
# Transpose the time series matrix
TS_transposed <- t(ts_data)  # Now N x T

# Create result directly using slice_msf_runwise
result_transpose <- slice_msf_runwise(
  TS = TS_transposed,
  mask = as.integer(mask@.Data),
  vol_dim = dim(mask),
  r = 6,
  fh_scale = 2.0 / (0.5 + 1.0),  # compactness=0.5 mapped
  min_size = 20,
  nbhd = 8,
  rows_are_time = FALSE,
  target_k_global = 3
)

n_transpose <- length(unique(result_transpose$labels[result_transpose$labels > 0]))
cat(sprintf("  With transposed matrix (rows_are_time=FALSE): %d clusters\n", n_transpose))

if (n_transpose == 3) {
  cat("  SUCCESS: rows_are_time=FALSE working correctly!\n")
} else {
  cat("  FAILED: rows_are_time=FALSE not working as expected\n")
}