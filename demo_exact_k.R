library(neurocluster)
library(neuroim2)

# Demonstration of exact-K functionality in SLiCE-MSF

# Create test data with clear structure
set.seed(123)
dims <- c(20, 20, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 100

# Create 4 distinct regions with different temporal patterns
nvox <- prod(dims)
mask_idx <- which(mask > 0)
coords <- arrayInd(mask_idx, dims)

# Assign regions based on spatial location
region <- integer(nvox)
region[coords[,1] <= 10 & coords[,2] <= 10] <- 1  # Top-left
region[coords[,1] > 10 & coords[,2] <= 10] <- 2   # Top-right
region[coords[,1] <= 10 & coords[,2] > 10] <- 3   # Bottom-left
region[coords[,1] > 10 & coords[,2] > 10] <- 4    # Bottom-right

# Create distinct time series for each region
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)

# Different temporal patterns
ts_data[region == 1, ] <- matrix(rep(sin(t_seq), sum(region == 1)), 
                                  ncol = ntime, byrow = TRUE) + 
                                  rnorm(sum(region == 1) * ntime, sd = 0.1)

ts_data[region == 2, ] <- matrix(rep(sin(3 * t_seq), sum(region == 2)), 
                                  ncol = ntime, byrow = TRUE) + 
                                  rnorm(sum(region == 2) * ntime, sd = 0.1)

ts_data[region == 3, ] <- matrix(rep(seq(0, 2, length.out = ntime), sum(region == 3)), 
                                  ncol = ntime, byrow = TRUE) + 
                                  rnorm(sum(region == 3) * ntime, sd = 0.1)

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

cat("Created test data with 4 distinct regions\n\n")

# 1. Natural clustering (no target K)
cat("1. Natural FH clustering (no target K):\n")
result_natural <- slice_msf(vec, mask, num_runs = 1, r = 8, 
                           min_size = 50, compactness = 5)
cat(sprintf("   Found %d clusters\n\n", length(unique(result_natural$cluster))))

# 2. Exact K = 4 (should preserve the 4 regions)
cat("2. Exact K = 4 (target_k_global = 4):\n")
result_k4 <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1, 
                       r = 8, min_size = 50, compactness = 5)
cat(sprintf("   Found exactly %d clusters\n\n", length(unique(result_k4$cluster))))

# 3. Exact K = 2 (should merge regions)
cat("3. Exact K = 2 (target_k_global = 2):\n")
result_k2 <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                       r = 8, min_size = 50, compactness = 5)
cat(sprintf("   Found exactly %d clusters\n\n", length(unique(result_k2$cluster))))

# 4. Per-slice clustering
cat("4. Exact K = 2 per slice (target_k_per_slice = 2):\n")
result_per_slice <- slice_msf(vec, mask, target_k_per_slice = 2, 
                              stitch_z = FALSE,  # Important!
                              num_runs = 1, r = 8, min_size = 20, compactness = 5)

# Check clusters per slice
cluster_vol <- array(0, dims)
cluster_vol[mask > 0] <- result_per_slice$cluster

for (z in 1:dims[3]) {
  slice_data <- cluster_vol[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("   Slice %d: %d clusters\n", z, length(slice_clusters)))
}

cat("\n5. Consensus with exact K (requires use_features=TRUE):\n")
result_consensus <- slice_msf(vec, mask, target_k_global = 3, 
                              num_runs = 3, consensus = TRUE,
                              use_features = TRUE,  # Required for exact-K in consensus
                              r = 8, min_size = 50, compactness = 5)
cat(sprintf("   Consensus clustering found exactly %d clusters\n", 
            length(unique(result_consensus$cluster))))

cat("\nDemonstration complete!\n")
cat("The exact-K functionality allows you to:\n")
cat("- Get exactly K clusters globally (target_k_global)\n")
cat("- Get exactly K clusters per slice (target_k_per_slice)\n")
cat("- Maintain spatial contiguity while achieving the target K\n")
cat("- Use with consensus fusion (requires use_features=TRUE)\n")