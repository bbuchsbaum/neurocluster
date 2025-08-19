library(neighborweights)
library(neuroim2)

# Simple test of weighted_spatial_adjacency without weight_mode
dims <- c(5, 5, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
mask.idx <- which(mask > 0)
coords <- index_to_coord(mask, mask.idx)

# Create simple feature matrix (rows = voxels, cols = features)
feature_mat <- matrix(rnorm(length(mask.idx) * 10), nrow=length(mask.idx))

print("Testing weighted_spatial_adjacency without weight_mode:")
tryCatch({
  W <- weighted_spatial_adjacency(coords, feature_mat,
                                  dthresh=6,
                                  nnk=9,
                                  wsigma=1,
                                  sigma=2,
                                  alpha=0.5,
                                  include_diagonal=FALSE,
                                  stochastic=TRUE)
  print("SUCCESS: weighted_spatial_adjacency works without weight_mode!")
  print(paste("W dim:", paste(dim(W), collapse="x")))
}, error = function(e) {
  print("ERROR: Still have error without weight_mode:")
  print(e)
})