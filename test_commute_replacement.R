library(neighborweights)
library(Matrix)

# Create test data
coords <- matrix(rnorm(30), nrow=10, ncol=3)
W <- spatial_adjacency(coords, dthresh=5, nnk=5)

# Test commute_time_distance 
print("Testing commute_time_distance:")
result1 <- commute_time_distance(W, ncomp=3)
print(paste("Result type:", typeof(result1)))
print(paste("Result class:", class(result1)))
print(paste("Result dim:", paste(dim(result1), collapse="x")))

# Test compute_diffusion_map
print("\nTesting compute_diffusion_map:")
result2 <- compute_diffusion_map(W, t=1, k=3)
print(paste("Result type:", typeof(result2)))
print(paste("Result class:", class(result2)))
if (is.list(result2)) {
  print(paste("Result names:", paste(names(result2), collapse=", ")))
} else {
  print(paste("Result dim:", paste(dim(result2), collapse="x")))
}

# The old commute_time likely returned a list with $cds
# commute_time_distance returns a matrix directly
# So we need to wrap it
print("\nLooks like commute_time_distance returns the coordinates directly")