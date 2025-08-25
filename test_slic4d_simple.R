library(neuroim2)
library(neurocluster)

# Test SLIC4D with existing test data
cat("Testing SLIC4D implementation with test data...\n")

# Load existing test data
mask <- read_vol("testdata/mask.nii")
scan <- read_vec("testdata/rscan01.nii.gz")

cat("Data dimensions:", dim(scan), "\n")
cat("Mask dimensions:", dim(mask), "\n")
cat("Number of masked voxels:", sum(as.logical(mask)), "\n")

# Test 1: Basic usage
cat("\nTest 1: Basic SLIC4D with K=100\n")
result <- slic4d_supervoxels(scan, mask, K = 100, 
                             compactness = 15,
                             max_iter = 5,
                             verbose = TRUE)

cat("\nResults:\n")
cat("  - Number of clusters:", length(unique(result$cluster)), "\n")
cat("  - Centers shape:", dim(result$centers), "\n")
cat("  - Coord centers shape:", dim(result$coord_centers), "\n")
cat("  - Class:", class(result), "\n")
cat("  - Fields:", names(result), "\n")

# Check structure consistency
is_consistent <- all(c("clusvol", "cluster", "centers", "coord_centers") %in% names(result))
cat("\nStructure is consistent with other clustering functions:", is_consistent, "\n")

# Test 2: With random projection for speed
cat("\nTest 2: With random projection (n_components=20)\n")
t1 <- system.time({
  result2 <- slic4d_supervoxels(scan, mask, K = 50,
                                n_components = 20,
                                max_iter = 3)
})
cat("  - Time:", t1["elapsed"], "seconds\n")
cat("  - Number of clusters:", length(unique(result2$cluster)), "\n")

cat("\nSLIC4D tests completed successfully!\n")