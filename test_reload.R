# Unload and reload the packages
if ('neurocluster' %in% loadedNamespaces()) unloadNamespace('neurocluster')
if ('neighborweights' %in% loadedNamespaces()) unloadNamespace('neighborweights')

# Now reload
library(neighborweights)
library(neurocluster)

# Try our debug test again
library(neuroim2)
dims <- c(5, 5, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
mask.idx <- which(mask > 0)
coords <- index_to_coord(mask, mask.idx)

tryCatch({
  G <- neighborweights::spatial_adjacency(coords,
                                          dthresh=9, 
                                          nnk=9,
                                          weight_mode='heat',
                                          sigma=.8, 
                                          stochastic=TRUE)
  cat('SUCCESS: spatial_adjacency works!\n')
}, error = function(e) {
  cat('ERROR: Still have error:\n')
  print(e)
})