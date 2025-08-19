library(neurocluster)
library(neuroim2)

# Use actual test data
test_data_dir <- file.path(getwd(), "testdata")
if (!dir.exists(test_data_dir)) {
  # If not found, try relative to test directory (for interactive testing)
  test_data_dir <- file.path(dirname(dirname(getwd())), "testdata")
}

bvec <- read_vec(file.path(test_data_dir, "rscan01.nii.gz"))
mask <- read_vol(file.path(test_data_dir, "mask.nii"))

# Check which supervoxel function is being called
print("Checking which supervoxel function would be called:")
print(methods("supervoxels"))

# Try to call the specific method
result <- neurocluster:::supervoxels.NeuroVec(bvec, mask, K=50, iterations=5)
print("Direct method call result:")
print(names(result))
print(paste("Has centers:", !is.null(result$centers)))
if (!is.null(result$centers)) {
  print(paste("Centers dim:", paste(dim(result$centers), collapse="x")))
} else {
  print("Centers is NULL")
}