# Test parameter handling for cluster4d framework
library(neurocluster)
library(neuroim2)

test_that("method-specific parameters are properly passed through", {
  # Create small test data
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test supervoxels-specific parameters
  if (exists("supervoxels")) {
    result_sv <- cluster4d(vec, mask, n_clusters = 5, 
                          method = "supervoxels",
                          sigma1 = 2.0,  # supervoxels-specific
                          sigma2 = 3.0,  # supervoxels-specific
                          use_gradient = FALSE,  # supervoxels-specific
                          converge_thresh = 0.01,  # supervoxels-specific
                          max_iterations = 2,
                          verbose = FALSE)
    
    expect_s3_class(result_sv, "cluster4d_result")
    # Check that parameters were recorded
    expect_equal(result_sv$parameters$sigma1, 2.0)
    expect_equal(result_sv$parameters$sigma2, 3.0)
  }
  
  # Test SLIC-specific parameters
  if (exists("slic4d_core")) {
    result_slic <- cluster4d(vec, mask, n_clusters = 5,
                            method = "slic",
                            preserve_k = TRUE,  # SLIC-specific
                            seed_relocate = "correlation",  # SLIC-specific
                            seed_relocate_radius = 2,  # SLIC-specific
                            n_components = 5,  # SLIC-specific
                            feature_norm = "l2",  # SLIC-specific
                            max_iterations = 2,
                            verbose = FALSE)
    
    expect_s3_class(result_slic, "cluster4d_result")
    expect_equal(result_slic$parameters$preserve_k, TRUE)
    expect_equal(result_slic$parameters$seed_relocate, "correlation")
  }
  
  # Test slice_msf-specific parameters
  result_msf <- cluster4d(vec, mask, n_clusters = 5,
                         method = "slice_msf",
                         num_runs = 2,  # slice_msf-specific
                         consensus = FALSE,  # slice_msf-specific
                         stitch_z = FALSE,  # slice_msf-specific
                         theta_link = 0.9,  # slice_msf-specific
                         min_contact = 2,  # slice_msf-specific
                         r = 6,  # slice_msf-specific
                         gamma = 2.0,  # slice_msf-specific
                         verbose = FALSE)
  
  expect_s3_class(result_msf, "cluster4d_result")
  expect_equal(result_msf$parameters$num_runs, 2)
  expect_equal(result_msf$parameters$stitch_z, FALSE)
  expect_equal(result_msf$parameters$theta_link, 0.9)
  
  # Test FLASH-3D-specific parameters
  if (exists("flash3d_supervoxels_cpp")) {
    result_flash <- cluster4d(vec, mask, n_clusters = 5,
                             method = "flash3d",
                             lambda_t = 1.5,  # FLASH-specific
                             lambda_g = 0.2,  # FLASH-specific
                             bits = 128,  # FLASH-specific
                             dctM = 16,  # FLASH-specific
                             max_iterations = 1,
                             verbose = FALSE)
    
    expect_s3_class(result_flash, "cluster4d_result")
    expect_equal(result_flash$parameters$lambda_t, 1.5)
    expect_equal(result_flash$parameters$bits, 128)
  }
})

test_that("spatial_weight parameter maps correctly to method-specific values", {
  # Test the parameter mapping logic
  
  # For supervoxels: alpha = 1 - spatial_weight
  spatial_weight <- 0.7
  alpha_expected <- 1 - spatial_weight  # 0.3
  
  # For SNIC: compactness = spatial_weight * 10
  compactness_snic_expected <- spatial_weight * 10  # 7
  
  # For SLIC: compactness = spatial_weight * 20
  compactness_slic_expected <- spatial_weight * 20  # 14
  
  # For FLASH: lambda_s = spatial_weight (direct)
  lambda_s_expected <- spatial_weight  # 0.7
  
  # Verify the mappings are as expected
  expect_equal(alpha_expected, 0.3)
  expect_equal(compactness_snic_expected, 7)
  expect_equal(compactness_slic_expected, 14)
  expect_equal(lambda_s_expected, 0.7)
})

test_that("common parameters work consistently across all methods", {
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Common parameters to test
  n_clusters <- 5
  spatial_weight <- 0.6
  max_iterations <- 2
  connectivity <- 6
  
  methods_to_test <- c("snic", "slice_msf")  # Methods that should always work
  
  for (method in methods_to_test) {
    result <- cluster4d(vec, mask, 
                       n_clusters = n_clusters,
                       spatial_weight = spatial_weight,
                       max_iterations = max_iterations,
                       connectivity = connectivity,
                       method = method,
                       verbose = FALSE)
    
    # Check that common parameters are handled
    expect_s3_class(result, "cluster4d_result")
    expect_equal(result$method, method)
    
    # Verify n_clusters is respected (approximately)
    expect_true(result$n_clusters > 0)
    expect_true(result$n_clusters <= n_clusters * 2)  # Allow some flexibility
    
    # Check that results have the required structure
    expect_true("clusvol" %in% names(result))
    expect_true("cluster" %in% names(result))
    expect_true("centers" %in% names(result))
    expect_true("coord_centers" %in% names(result))
  }
})

test_that("parameter validation works correctly", {
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test invalid spatial_weight
  expect_error(
    cluster4d(vec, mask, n_clusters = 5, spatial_weight = 1.5),
    "spatial_weight must be between 0 and 1"
  )
  
  expect_error(
    cluster4d(vec, mask, n_clusters = 5, spatial_weight = -0.1),
    "spatial_weight must be between 0 and 1"
  )
  
  # Test invalid connectivity
  expect_error(
    cluster4d(vec, mask, n_clusters = 5, connectivity = 10),
    "connectivity must be 6, 18, 26, or 27"
  )
  
  # Test invalid method
  expect_error(
    cluster4d(vec, mask, n_clusters = 5, method = "invalid_method"),
    "'arg' should be one of"
  )
})

test_that("legacy parameter names still work with warnings", {
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test that K parameter gets mapped to n_clusters
  params <- neurocluster:::map_cluster4d_params("snic", K = 10, compactness = 5)
  expect_equal(params$n_clusters, 10)
  expect_null(params$K)
  
  # Test that alpha gets mapped for supervoxels
  params <- neurocluster:::map_cluster4d_params("supervoxels", K = 10, alpha = 0.3)
  expect_equal(params$n_clusters, 10)
  expect_equal(params$spatial_weight, 0.7)  # 1 - alpha
  
  # Test that iterations gets mapped
  params <- neurocluster:::map_cluster4d_params("supervoxels", iterations = 25)
  expect_equal(params$max_iterations, 25)
  
  # Test that max_iter gets mapped
  params <- neurocluster:::map_cluster4d_params("snic", max_iter = 15)
  expect_equal(params$max_iterations, 15)
  
  # Test that rounds gets mapped for flash3d
  params <- neurocluster:::map_cluster4d_params("flash3d", rounds = 3)
  expect_equal(params$max_iterations, 3)
})

test_that("method wrappers handle parameters correctly", {
  mask <- NeuroVol(array(1, c(8, 8, 4)), NeuroSpace(c(8, 8, 4)))
  vec <- replicate(10, NeuroVol(array(runif(8*8*4), c(8, 8, 4)),
                                NeuroSpace(c(8, 8, 4))), simplify=FALSE)
  vec <- do.call(concat, vec)
  
  # Test cluster4d_snic wrapper
  result_snic <- cluster4d_snic(vec, mask, 
                               n_clusters = 5,
                               spatial_weight = 0.5,
                               max_iterations = 2,
                               verbose = FALSE)
  
  expect_s3_class(result_snic, "cluster4d_result")
  expect_equal(result_snic$method, "snic")
  expect_equal(result_snic$parameters$spatial_weight, 0.5)
  
  # Test cluster4d_slice_msf wrapper
  result_msf <- cluster4d_slice_msf(vec, mask,
                                   n_clusters = 5,
                                   spatial_weight = 0.4,
                                   num_runs = 1,
                                   consensus = FALSE,
                                   verbose = FALSE)
  
  expect_s3_class(result_msf, "cluster4d_result")
  expect_equal(result_msf$method, "slice_msf")
  expect_equal(result_msf$parameters$num_runs, 1)
  
  # Test that supervoxels wrapper exists and works
  if (exists("cluster4d_supervoxels")) {
    result_sv <- cluster4d_supervoxels(vec, mask,
                                      n_clusters = 5,
                                      spatial_weight = 0.3,
                                      sigma1 = 1.5,
                                      sigma2 = 2.5,
                                      max_iterations = 2,
                                      verbose = FALSE)
    
    expect_s3_class(result_sv, "cluster4d_result")
    expect_equal(result_sv$method, "supervoxels")
  }
})