library(neurocluster)
library(neuroim2)
library(neighborweights)

# Create minimal test data
dims <- c(5, 5, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
mask.idx <- which(mask > 0)
coords <- index_to_coord(mask, mask.idx)

# Try to call spatial_adjacency directly
tryCatch({
  G <- neighborweights::spatial_adjacency(coords,
                                          dthresh=9, 
                                          nnk=9,
                                          weight_mode="heat",
                                          sigma=.8, 
                                          stochastic=TRUE)
  print("Success!")
}, error = function(e) {
  print("Error occurred:")
  print(e)
  print("")
  print("Let's check the function signature:")
  print(args(neighborweights::spatial_adjacency))
})