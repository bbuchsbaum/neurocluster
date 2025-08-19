library(neighborweights)

# Check the exact parameter order for weighted_spatial_adjacency
print("Function formals:")
print(formals(weighted_spatial_adjacency))

# Test with exact parameter order
coords <- matrix(rnorm(30), nrow=10, ncol=3)
feature_mat <- matrix(rnorm(100), nrow=10, ncol=10)

print("\nTesting with positional arguments in exact order:")
tryCatch({
  W <- weighted_spatial_adjacency(
    coord_mat = coords,
    feature_mat = feature_mat,
    wsigma = 0.73,
    alpha = 0.5,
    nnk = 27,
    weight_mode = "binary",
    sigma = 1,
    dthresh = 2.5,
    include_diagonal = TRUE,
    normalized = FALSE,
    stochastic = FALSE
  )
  print("SUCCESS!")
}, error = function(e) {
  print(paste("ERROR:", e$message))
})