library(neurocluster)
library(neuroim2)

# Debug exact K functionality
set.seed(1000)
dims <- c(30, 30, 10)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 100

# Define 9 regions in a 3x3 grid pattern
nvox <- prod(dims)
mask_idx <- which(mask > 0)
coords <- arrayInd(mask_idx, dims)

# Assign each voxel to one of 9 regions
region <- integer(nvox)
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  grid_x <- min(floor((x - 1) / 10), 2)
  grid_y <- min(floor((y - 1) / 10), 2)
  
  region[i] <- grid_x + 3 * grid_y + 1
}

# Check region distribution
cat("Region distribution:\n")
print(table(region))

# Create 9 distinct temporal patterns
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)

patterns <- rbind(
  sin(t_seq),                           # Region 1
  cos(t_seq),                           # Region 2
  sin(2 * t_seq),                       # Region 3
  cos(2 * t_seq),                       # Region 4
  sin(0.5 * t_seq),                     # Region 5
  cos(0.5 * t_seq),                     # Region 6
  0.5 * sin(t_seq) + 0.5 * cos(t_seq), # Region 7
  sin(t_seq + pi/4),                    # Region 8
  cos(t_seq + pi/4)                     # Region 9
)

for (r in 1:9) {
  voxels <- which(region == r)
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

# Test natural clustering
cat("\n1. Natural clustering:\n")
result_natural <- slice_msf(vec, mask, num_runs = 1, r = 12, 
                            min_size = 50, compactness = 5)
cat(sprintf("Found %d clusters\n", length(unique(result_natural$cluster))))

# Test exact K = 9
cat("\n2. Exact K = 9:\n")
result_k9 <- slice_msf(vec, mask, target_k_global = 9, num_runs = 1, 
                       r = 12, min_size = 50, compactness = 5)
cat(sprintf("Found %d clusters\n", length(unique(result_k9$cluster))))

# Check cluster sizes
cat("\nCluster sizes:\n")
print(table(result_k9$cluster))

# Test per-slice with smaller example
cat("\n3. Testing per-slice K = 2:\n")
# Create smaller test case
dims2 <- c(10, 10, 3)
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))
nvox2 <- prod(dims2)
ntime2 <- 40

# Simple left/right pattern
ts_data2 <- matrix(0, nrow = nvox2, ncol = ntime2)
coords2 <- arrayInd(1:nvox2, dims2)
t_seq2 <- seq(0, 2*pi, length.out = ntime2)

# Left half: sine, Right half: cosine
for (i in 1:nvox2) {
  if (coords2[i,1] <= 5) {
    ts_data2[i, ] <- sin(t_seq2) + rnorm(ntime2, sd = 0.1)
  } else {
    ts_data2[i, ] <- cos(t_seq2) + rnorm(ntime2, sd = 0.1)
  }
}

vec_list2 <- lapply(1:ntime2, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask2 > 0] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

result_per_slice <- slice_msf(vec2, mask2, target_k_per_slice = 2,
                              stitch_z = FALSE,
                              num_runs = 1, r = 6, min_size = 10, 
                              compactness = 3)

# Check clusters per slice
cluster_vol <- array(0, dims2)
cluster_vol[mask2 > 0] <- result_per_slice$cluster

cat("\nClusters per slice:\n")
for (z in 1:dims2[3]) {
  slice_data <- cluster_vol[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("Slice %d: %d clusters (%s)\n", z, length(slice_clusters), 
              paste(slice_clusters, collapse=", ")))
}

cat("\nTotal unique clusters:", length(unique(result_per_slice$cluster)), "\n")