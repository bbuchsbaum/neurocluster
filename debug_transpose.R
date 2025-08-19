library(neighborweights)

# Test with transposed feature matrix like in commute_cluster
coords <- matrix(rnorm(30), nrow=10, ncol=3)
feature_mat <- matrix(rnorm(100), nrow=10, ncol=10)

print("Testing with transposed feature matrix:")
tryCatch({
  W <- weighted_spatial_adjacency(
    coord_mat = coords,
    feature_mat = t(feature_mat),  # Transposed!
    wsigma = 0.73,
    alpha = 0.5,
    nnk = 27,
    weight_mode = "binary",
    sigma = 1,
    dthresh = 2.5,
    include_diagonal = FALSE,
    stochastic = TRUE
  )
  print("SUCCESS with transposed matrix!")
  print(paste("W dimensions:", paste(dim(W), collapse="x")))
}, error = function(e) {
  print(paste("ERROR:", e$message))
})

# Now test the exact call pattern from commute_cluster
print("\nTesting exact commute_cluster pattern (positional args):")
tryCatch({
  W <- neighborweights::weighted_spatial_adjacency(coords, t(feature_mat),
                                                   dthresh=2.5,
                                                   nnk=27,
                                                   wsigma=0.73,
                                                   sigma=1,
                                                   alpha=0.5,
                                                   weight_mode="binary",
                                                   include_diagonal=FALSE,
                                                   stochastic=TRUE)
  print("SUCCESS with commute_cluster pattern!")
}, error = function(e) {
  print(paste("ERROR:", e$message))
})