library(neurocluster)
library(neuroim2)

# Use test data
test_data_dir <- file.path(getwd(), "testdata")
if (!dir.exists(test_data_dir)) {
  test_data_dir <- file.path(dirname(dirname(getwd())), "testdata")
}

bvec <- read_vec(file.path(test_data_dir, "rscan01.nii.gz"))
mask <- read_vol(file.path(test_data_dir, "mask.nii"))

# Extract the relevant parts from commute_cluster to debug
mask.idx <- which(mask > 0)
grid <- index_to_coord(mask, mask.idx)
feature_mat <- neuroim2::series(bvec, mask.idx)

# Check dimensions
print(paste("Grid dimensions:", paste(dim(grid), collapse="x")))
print(paste("Feature mat dimensions:", paste(dim(feature_mat), collapse="x")))
print(paste("Transposed feature mat:", paste(dim(t(feature_mat)), collapse="x")))

# Try the call
W <- neighborweights::weighted_spatial_adjacency(grid, t(feature_mat),
                                                 dthresh=5*6,
                                                 nnk=27,
                                                 wsigma=.73,
                                                 sigma=5,
                                                 alpha=.5,
                                                 weight_mode="binary",
                                                 include_diagonal=FALSE,
                                                 stochastic=TRUE)

print(paste("W class:", class(W)))
print(paste("W type:", typeof(W)))
if (is.matrix(W)) {
  print(paste("W dimensions:", paste(dim(W), collapse="x")))
} else {
  print("W is not a matrix")
  print(str(W))
}