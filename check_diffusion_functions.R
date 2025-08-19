library(neighborweights)

# Check signatures
print("commute_time_distance signature:")
print(args(commute_time_distance))

print("\ncompute_diffusion_map signature:")
print(args(compute_diffusion_map))

# The old commute_time function likely returned something with a $cds component
# Let's see what these return
coords <- matrix(rnorm(30), nrow=10, ncol=3)
feature_mat <- matrix(rnorm(100), nrow=10, ncol=10)

# Create a simple adjacency matrix
W <- spatial_adjacency(coords, dthresh=5, nnk=5)

print("\nTesting compute_diffusion_map:")
tryCatch({
  result <- compute_diffusion_map(W, ndim=3)
  print(paste("Result class:", class(result)))
  print(paste("Result names:", paste(names(result), collapse=", ")))
  if ("cds" %in% names(result)) {
    print("Has 'cds' component!")
  }
}, error = function(e) {
  print(paste("Error:", e$message))
})