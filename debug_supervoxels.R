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

# Try supervoxels
tryCatch({
  cres <- supervoxels(bvec, mask, K=50, iterations=5)
  
  # Check the result
  print("Supervoxels result:")
  print(names(cres))
  print(paste("Class:", class(cres)))
  print(paste("Has centers:", !is.null(cres$centers)))
  print(paste("Centers dim:", paste(dim(cres$centers), collapse="x")))
  
}, error = function(e) {
  print("Error:")
  print(e)
  traceback()
})