library(neurocluster)
library(neuroim2)

# Simple 2x2 test case
set.seed(100)
dims <- c(10, 10, 2)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 50

# Create 4 perfect quadrants
nvox <- prod(dims)
mask_idx <- which(mask > 0)
coords <- arrayInd(mask_idx, dims)

region <- integer(nvox)
region[coords[,1] <= 5 & coords[,2] <= 5] <- 1  # Top-left
region[coords[,1] > 5 & coords[,2] <= 5] <- 2   # Top-right  
region[coords[,1] <= 5 & coords[,2] > 5] <- 3   # Bottom-left
region[coords[,1] > 5 & coords[,2] > 5] <- 4    # Bottom-right

cat("Region distribution:\n")
print(table(region))

# Create very distinct patterns
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 2*pi, length.out = ntime)

# Orthogonal patterns
patterns <- rbind(
  sin(t_seq),     # Region 1
  cos(t_seq),     # Region 2  
  -sin(t_seq),    # Region 3
  -cos(t_seq)     # Region 4
)

# Assign with minimal noise
for (r in 1:4) {
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

# Test different parameters
cat("\n1. Testing with different compactness values:\n")
for (comp in c(1, 3, 5, 10)) {
  result <- slice_msf(vec, mask, num_runs = 1, r = 6, 
                      min_size = 10, compactness = comp)
  cat(sprintf("  Compactness=%d -> %d clusters\n", 
              comp, length(unique(result$cluster))))
}

cat("\n2. Testing exact K:\n")
# Try exact K = 4
result_k4 <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1, 
                       r = 6, min_size = 10, compactness = 3)
cat(sprintf("  target_k_global=4 -> %d clusters\n", 
            length(unique(result_k4$cluster))))

# Check if it's a min_size issue
cat("\n3. Testing with smaller min_size:\n")
result_k4_small <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1, 
                             r = 6, min_size = 5, compactness = 3)
cat(sprintf("  target_k_global=4, min_size=5 -> %d clusters\n", 
            length(unique(result_k4_small$cluster))))

# Try forcing more initial clusters with higher compactness
cat("\n4. Testing with lower compactness (more initial clusters):\n")
result_k4_low <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1, 
                           r = 6, min_size = 5, compactness = 0.5)
cat(sprintf("  target_k_global=4, compactness=0.5 -> %d clusters\n", 
            length(unique(result_k4_low$cluster))))

# Check initial clustering without exact K
cat("\n5. Natural clustering with very low compactness:\n")
result_natural <- slice_msf(vec, mask, num_runs = 1, r = 6, 
                           min_size = 5, compactness = 0.1)
cat(sprintf("  compactness=0.1 -> %d clusters (before exact-K)\n", 
            length(unique(result_natural$cluster))))