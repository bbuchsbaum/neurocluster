library(testthat)
library(neurocluster)

test_that("best_candidate functions are available", {
  # Check that the functions exist
  expect_true(exists("best_candidate_sequential", asNamespace("neurocluster")))
  expect_true(exists("best_candidate_parallel", asNamespace("neurocluster")))
  
  # Create minimal test data
  set.seed(123)
  n <- 100
  k <- 5
  
  # Create candidate list - each voxel can switch to neighbors
  candidates <- lapply(1:n, function(i) sample(1:k, size=min(3, k), replace=FALSE))
  
  # Current cluster assignments
  curclus <- sample(1:k, n, replace=TRUE)
  
  # Coordinates (3D)
  coords <- matrix(rnorm(3*n), nrow=3, ncol=n)
  
  # Data (10 features)
  data <- matrix(rnorm(10*n), nrow=10, ncol=n)
  
  # Cluster centroids
  data_centroids <- matrix(rnorm(10*k), nrow=10, ncol=k)
  coord_centroids <- matrix(rnorm(3*k), nrow=3, ncol=k)
  
  # Test sequential version
  res_seq <- neurocluster:::best_candidate_sequential(
    candidates, curclus, coords, data_centroids, 
    coord_centroids, data, 
    sigma1=1, sigma2=1, alpha=0.5
  )
  
  expect_equal(length(res_seq), n)
  expect_true(all(res_seq %in% 1:k))
  
  # Test parallel version - should work but might crash
  # For now, just test that it doesn't error on call
  expect_no_error({
    res_par <- neurocluster:::best_candidate_parallel(
      candidates, curclus, coords, data_centroids, 
      coord_centroids, data, 
      sigma1=1, sigma2=1, alpha=0.5, grain_size=20
    )
  })
})