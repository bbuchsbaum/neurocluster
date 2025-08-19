library(neurocluster)
library(neuroim2)
library(neighborweights)

# Test different approaches

# Simple test data
dims <- c(5, 5, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
mask.idx <- which(mask > 0)
coords <- index_to_coord(mask, mask.idx)
feature_mat <- matrix(rnorm(length(mask.idx) * 10), nrow=length(mask.idx))

print("Testing different parameter combinations:")

# Test 1: Default weight_mode
print("\n1. With default weight_mode:")
tryCatch({
  W <- weighted_spatial_adjacency(coords, feature_mat,
                                  dthresh=6,
                                  nnk=9,
                                  wsigma=1,
                                  sigma=2,
                                  alpha=0.5,
                                  include_diagonal=FALSE,
                                  stochastic=TRUE)
  print("SUCCESS with defaults!")
}, error = function(e) {
  print(paste("ERROR:", e$message))
})

# Test 2: Explicit weight_mode="heat"
print("\n2. With weight_mode='heat':")
tryCatch({
  W <- weighted_spatial_adjacency(coords, feature_mat,
                                  dthresh=6,
                                  nnk=9,
                                  wsigma=1,
                                  sigma=2,
                                  alpha=0.5,
                                  weight_mode="heat",
                                  include_diagonal=FALSE,
                                  stochastic=TRUE)
  print("SUCCESS with weight_mode='heat'!")
}, error = function(e) {
  print(paste("ERROR:", e$message))
})

# Test 3: Looking at function formals
print("\n3. Function default parameters:")
print(formals(weighted_spatial_adjacency))