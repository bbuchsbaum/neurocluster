library(neuroim2)
library(neurocluster)

# Test SLIC4D with simulated data
cat("Testing SLIC4D implementation...\n")

# Create small test data  
dims <- c(20, 20, 10)  # 3D spatial dimensions
n_time <- 50  # Time points
cat("Creating test data with dimensions:", dims, "x", n_time, "\n")

# Create a simple mask first
mask_arr <- array(TRUE, dim = dims)
mask_arr[1:5, 1:5, ] <- FALSE  # Exclude corner

# Create a NeuroSpace for the 3D mask
sp3d <- NeuroSpace(dims, spacing = c(2, 2, 2))

# Create mask as NeuroVol
mask <- NeuroVol(mask_arr, space = sp3d)

cat("Mask has", sum(mask_arr), "voxels\n")

# Create 4D data directly
data_arr <- array(rnorm(prod(dims) * n_time), dim = c(dims, n_time))

# Add some structure to make clusters meaningful
for (i in 1:5) {
  x_range <- ((i-1)*4 + 1):min(i*4, dims[1])
  y_range <- 1:dims[2]
  z_range <- 1:dims[3]
  data_arr[x_range, y_range, z_range, ] <- 
    data_arr[x_range, y_range, z_range, ] + rnorm(1, mean = i, sd = 0.5)
}

# Create NeuroVec with 3D space (time is 4th dim without spacing)
bvec <- NeuroVec(data_arr, space = sp3d)

# Test 1: Basic usage with default parameters
cat("\nTest 1: Basic usage with K=50\n")
result1 <- slic4d_supervoxels(bvec, mask, K = 50, verbose = TRUE)

cat("Result structure:\n")
cat("  - Number of clusters:", length(unique(result1$cluster)), "\n")
cat("  - Centers shape:", dim(result1$centers), "\n")
cat("  - Coord centers shape:", dim(result1$coord_centers), "\n")
cat("  - ClusteredNeuroVol dims:", dim(result1$clusvol), "\n")

# Test 2: With feature normalization
cat("\nTest 2: With L2 normalization\n")
result2 <- slic4d_supervoxels(bvec, mask, K = 30, 
                              feature_norm = "l2",
                              compactness = 20)

cat("  - Number of clusters:", length(unique(result2$cluster)), "\n")

# Test 3: With random projection
cat("\nTest 3: With random projection (n_components=10)\n")
result3 <- slic4d_supervoxels(bvec, mask, K = 25,
                              n_components = 10,
                              max_iter = 5)

cat("  - Number of clusters:", length(unique(result3$cluster)), "\n")
cat("  - Centers shape:", dim(result3$centers), "\n")

# Test 4: Check consistency with other clustering functions
cat("\nTest 4: Checking return structure consistency\n")
cat("  - Has 'cluster_result' class:", "cluster_result" %in% class(result1), "\n")
cat("  - Has required fields:", 
    all(c("clusvol", "cluster", "centers", "coord_centers") %in% names(result1)), "\n")

# Test 5: Performance comparison
cat("\nTest 5: Timing comparison\n")
t1 <- system.time({
  result_fast <- slic4d_supervoxels(bvec, mask, K = 100, max_iter = 3)
})
cat("  - SLIC4D time for K=100:", t1["elapsed"], "seconds\n")

# Test 6: Compare with different seed methods
cat("\nTest 6: Testing farthest seed method\n")
result_farthest <- slic4d_supervoxels(bvec, mask, K = 20,
                                      seed_method = "farthest",
                                      max_iter = 3)
cat("  - Number of clusters (farthest):", length(unique(result_farthest$cluster)), "\n")

cat("\nAll tests completed successfully!\n")