# Load required packages
library(neurocluster)
library(neuroim2)

test_that("cluster4d validates inputs correctly", {
  # Create small test data
  mask <- NeuroVol(array(1, c(10, 10, 10)), NeuroSpace(c(10, 10, 10)))
  vec <- replicate(20, NeuroVol(array(runif(10*10*10), c(10, 10, 10)),
                                NeuroSpace(c(10, 10, 10))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test invalid vec input
  expect_error(
    cluster4d(list(), mask, n_clusters = 10),
    "vec must be a NeuroVec"
  )
  
  # Test invalid mask input
  expect_error(
    cluster4d(vec, array(1, c(10, 10, 10)), n_clusters = 10),
    "mask must be a NeuroVol"
  )
  
  # Test n_clusters validation
  expect_error(
    cluster4d(vec, mask, n_clusters = -1),
    "n_clusters must be positive"
  )
  
  expect_error(
    cluster4d(vec, mask, n_clusters = 0),
    "n_clusters must be positive"
  )
  
  # Test too many clusters
  expect_error(
    cluster4d(vec, mask, n_clusters = 2000),
    "Cannot create 2000 clusters from 1000 masked voxels"
  )
  
  # Test spatial_weight validation
  expect_error(
    cluster4d(vec, mask, n_clusters = 10, spatial_weight = -0.1),
    "spatial_weight must be between 0 and 1"
  )
  
  expect_error(
    cluster4d(vec, mask, n_clusters = 10, spatial_weight = 1.1),
    "spatial_weight must be between 0 and 1"
  )
  
  # Test connectivity validation
  expect_error(
    cluster4d(vec, mask, n_clusters = 10, connectivity = 5),
    "connectivity must be 6, 18, 26, or 27"
  )
})

test_that("cluster4d works with all methods", {
  # Create small test data
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  methods <- c("supervoxels", "snic", "slic", "slice_msf", "flash3d")
  
  for (method in methods) {
    # Skip methods that might not be available
    if (method == "slic" && !exists("slic4d_core")) {
      skip(paste("Skipping", method, "- C++ implementation not available"))
    }
    if (method == "flash3d" && !exists("flash3d_supervoxels_cpp")) {
      skip(paste("Skipping", method, "- C++ implementation not available"))
    }
    
    result <- cluster4d(vec, mask, n_clusters = 10, method = method,
                       max_iterations = 2, verbose = FALSE)
    
    # Check result structure
    expect_s3_class(result, "cluster4d_result")
    expect_s3_class(result, "cluster_result")
    
    # Check required components
    expect_true("clusvol" %in% names(result))
    expect_true("cluster" %in% names(result))
    expect_true("centers" %in% names(result))
    expect_true("coord_centers" %in% names(result))
    expect_true("n_clusters" %in% names(result))
    expect_true("method" %in% names(result))
    expect_true("parameters" %in% names(result))
    
    # Check dimensions
    expect_equal(length(result$cluster), sum(mask > 0))
    expect_equal(nrow(result$centers), result$n_clusters)
    expect_equal(nrow(result$coord_centers), result$n_clusters)
    expect_equal(ncol(result$coord_centers), 3)
    
    # Check method is recorded correctly
    expect_equal(result$method, method)
  }
})

test_that("parameter mapping works correctly", {
  params_supervoxels <- neurocluster:::map_cluster4d_params("supervoxels", 
                                            K = 100, alpha = 0.3, iterations = 25)
  expect_equal(params_supervoxels$n_clusters, 100)
  expect_equal(params_supervoxels$spatial_weight, 0.7)  # 1 - alpha
  expect_equal(params_supervoxels$max_iterations, 25)
  expect_null(params_supervoxels$K)
  
  params_snic <- neurocluster:::map_cluster4d_params("snic", 
                                      K = 50, compactness = 5, max_iter = 10)
  expect_equal(params_snic$n_clusters, 50)
  expect_equal(params_snic$spatial_weight, 0.5)  # compactness/10
  expect_equal(params_snic$max_iterations, 10)
  
  params_flash <- neurocluster:::map_cluster4d_params("flash3d",
                                       K = 75, lambda_s = 0.8, rounds = 3)
  expect_equal(params_flash$n_clusters, 75)
  expect_equal(params_flash$spatial_weight, 0.8)
  expect_equal(params_flash$max_iterations, 3)
})

test_that("suggest_cluster4d_params provides reasonable suggestions", {
  # Small data
  params_small <- suggest_cluster4d_params(1000, 100, priority = "balanced")
  expect_true(params_small$n_clusters > 0)
  expect_true(params_small$n_clusters < 1000)
  expect_true("recommended_method" %in% names(params_small))
  
  # Large data with speed priority
  params_fast <- suggest_cluster4d_params(100000, 200, priority = "speed")
  expect_true(params_fast$recommended_method %in% c("slice_msf", "flash3d"))
  
  # Quality priority
  params_quality <- suggest_cluster4d_params(5000, 150, priority = "quality")
  expect_true(params_quality$recommended_method %in% c("supervoxels", "slic"))
  
  # Memory priority
  params_memory <- suggest_cluster4d_params(50000, 300, priority = "memory")
  expect_equal(params_memory$recommended_method, "snic")
})

test_that("cluster4d S3 methods work", {
  # Create small test result
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  result <- cluster4d(vec, mask, n_clusters = 5, method = "snic",
                     max_iterations = 1, verbose = FALSE)
  
  # Test print method
  output <- capture.output(print(result))
  expect_true(any(grepl("Cluster4D Result", output)))
  expect_true(any(grepl("Method: snic", output)))
  expect_true(any(grepl("Number of clusters:", output)))
  
  # Test summary method
  summary_output <- capture.output(summary(result))
  expect_true(any(grepl("Cluster4D Analysis Summary", summary_output)))
  expect_true(any(grepl("Parameters:", summary_output)))
  expect_true(any(grepl("Cluster Size Statistics:", summary_output)))
  
  # Test validate method
  validation <- validate_cluster4d(result, vec, mask)
  expect_true(validation$valid)
  expect_equal(validation$summary$n_voxels, sum(mask > 0))
})

test_that("compare_cluster4d works with multiple results", {
  # Create small test data
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Create two results with different methods
  result1 <- cluster4d(vec, mask, n_clusters = 5, method = "snic",
                      max_iterations = 1, verbose = FALSE)
  
  # Use a different method or parameters for result2
  result2 <- cluster4d(vec, mask, n_clusters = 8, method = "snic",
                      spatial_weight = 0.8, max_iterations = 1, verbose = FALSE)
  
  # Compare results
  comparison <- compare_cluster4d(result1, result2, 
                                 metrics = c("summary", "temporal_coherence"))
  
  expect_s3_class(comparison, "data.frame")
  expect_equal(nrow(comparison), 2)
  expect_true("N_Clusters" %in% names(comparison))
  expect_true("Mean_Size" %in% names(comparison))
})

test_that("initialization methods work correctly", {
  # Create test coordinates and features
  coords <- matrix(runif(300), ncol = 3)  # 100 points in 3D
  features <- matrix(rnorm(1000), ncol = 10)  # 100 points, 10 features
  
  # Test kmeans initialization
  init_kmeans <- neurocluster:::initialize_clusters(coords, features, n_clusters = 10,
                                    method = "kmeans")
  expect_equal(length(init_kmeans$seeds), 10)
  expect_equal(length(init_kmeans$initial_labels), 100)
  
  # Test random initialization
  init_random <- neurocluster:::initialize_clusters(coords, features, n_clusters = 10,
                                    method = "random")
  expect_equal(length(init_random$seeds), 10)
  
  # Test grid initialization
  init_grid <- neurocluster:::initialize_clusters(coords, features, n_clusters = 10,
                                  method = "grid")
  expect_true(length(init_grid$seeds) >= 10)
  
  # Test Poisson disk initialization
  init_poisson <- neurocluster:::initialize_clusters(coords, features, n_clusters = 10,
                                     method = "poisson")
  expect_true(length(init_poisson$seeds) >= 1)
})

test_that("compute_cluster_centers works correctly", {
  # Create test data
  labels <- c(1, 1, 2, 2, 2, 3, 3)
  features <- matrix(1:21, nrow = 7, ncol = 3)
  coords <- matrix(runif(21), nrow = 7, ncol = 3)
  
  # Test mean centers
  centers_mean <- neurocluster:::compute_cluster_centers(labels, features, coords, method = "mean")
  expect_equal(centers_mean$n_clusters, 3)
  expect_equal(nrow(centers_mean$centers), 3)
  expect_equal(nrow(centers_mean$coord_centers), 3)
  
  # Test with NA values
  labels_na <- c(1, 1, NA, 2, 2, 3, NA)
  centers_na <- neurocluster:::compute_cluster_centers(labels_na, features, coords, method = "mean")
  expect_equal(centers_na$n_clusters, 3)
  
  # Test medoid centers
  centers_medoid <- neurocluster:::compute_cluster_centers(labels, features, coords, method = "medoid")
  expect_equal(centers_medoid$n_clusters, 3)
})

test_that("backward compatibility is maintained", {
  # Create small test data
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test that old function names still work (skip if not available)
  if (exists("snic")) {
    result_old <- snic(vec, mask, K = 5, compactness = 5, max_iter = 1)
  } else {
    skip("Original snic function not available")
  }
  expect_s3_class(result_old, "cluster_result")
  expect_s3_class(result_old, "snic_cluster_result")
  
  # Test that results have expected structure
  expect_true("clusvol" %in% names(result_old))
  expect_true("cluster" %in% names(result_old))
  
  # New unified interface should produce compatible results
  result_new <- cluster4d(vec, mask, n_clusters = 5, method = "snic",
                         spatial_weight = 0.5, max_iterations = 1)
  
  # Both should have the same essential components
  expect_equal(length(result_old$cluster), length(result_new$cluster))
  expect_equal(dim(result_old$clusvol), dim(result_new$clusvol))
})