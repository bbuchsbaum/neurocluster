library(neurocluster)
library(neuroim2)

# Test case that should have 9 regions but might naturally cluster to fewer
set.seed(1000)
dims <- c(27, 27, 3)  # 9x9 per region, 3 slices
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 80

# Create perfect 3x3 grid
nvox <- prod(dims)
mask_idx <- which(mask > 0)
coords <- arrayInd(mask_idx, dims)

region <- integer(nvox)
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  grid_x <- min(floor((x - 1) / 9), 2)
  grid_y <- min(floor((y - 1) / 9), 2)
  
  region[i] <- grid_x + 3 * grid_y + 1
}

cat("Region distribution:\n")
print(table(region))

# Create patterns that are similar within "rows" to encourage merging
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)

# Make regions 1,2,3 similar (top row)
# Make regions 4,5,6 similar (middle row)  
# Make regions 7,8,9 similar (bottom row)
patterns <- rbind(
  sin(t_seq),                      # Region 1
  sin(t_seq + 0.1),                # Region 2 (very similar to 1)
  sin(t_seq + 0.2),                # Region 3 (very similar to 1,2)
  cos(t_seq),                      # Region 4
  cos(t_seq + 0.1),                # Region 5 (very similar to 4)
  cos(t_seq + 0.2),                # Region 6 (very similar to 4,5)
  sin(2*t_seq),                    # Region 7
  sin(2*t_seq + 0.1),              # Region 8 (very similar to 7)
  sin(2*t_seq + 0.2)               # Region 9 (very similar to 7,8)
)

# Add small noise
for (r in 1:9) {
  voxels <- which(region == r)
  ts_data[voxels, ] <- matrix(rep(patterns[r,], length(voxels)), 
                               ncol = ntime, byrow = TRUE) + 
                               rnorm(length(voxels) * ntime, sd = 0.05)
}

# Create NeuroVec
vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test natural clustering
cat("\n1. Natural clustering (should merge similar regions):\n")
result_natural <- slice_msf(vec, mask, num_runs = 1, r = 8, 
                           min_size = 20, compactness = 3)
n_natural <- length(unique(result_natural$cluster))
cat(sprintf("  Found %d clusters\n", n_natural))

# Print cluster sizes
cat("  Cluster sizes: ")
cat(paste(sort(table(result_natural$cluster)), collapse=", "))
cat("\n")

# Test exact K = 9 (should separate all regions)
cat("\n2. Exact K = 9 (force all regions separate):\n")
result_k9 <- slice_msf(vec, mask, target_k_global = 9, num_runs = 1, 
                       r = 8, min_size = 20, compactness = 3)
n_k9 <- length(unique(result_k9$cluster))
cat(sprintf("  Found %d clusters\n", n_k9))

if (n_k9 < 9) {
  cat("\n  ISSUE: Not reaching target K=9\n")
  cat("  Cluster sizes: ")
  cat(paste(sort(table(result_k9$cluster)), collapse=", "))
  cat("\n")
  
  # Try with even lower compactness to get more initial clusters
  cat("\n3. Try with very low compactness:\n")
  result_k9_low <- slice_msf(vec, mask, target_k_global = 9, num_runs = 1, 
                             r = 8, min_size = 20, compactness = 0.1)
  cat(sprintf("  Found %d clusters\n", length(unique(result_k9_low$cluster))))
}

# Test exact K = 3 (should merge to 3 rows)
cat("\n4. Exact K = 3 (merge to rows):\n")
result_k3 <- slice_msf(vec, mask, target_k_global = 3, num_runs = 1,
                       r = 8, min_size = 20, compactness = 3)
cat(sprintf("  Found %d clusters\n", length(unique(result_k3$cluster))))

# Check how regions map to clusters
cluster_mapping <- table(region, result_k3$cluster)
cat("\n  Region -> Cluster mapping:\n")
print(cluster_mapping)