library(testthat)
library(neurocluster)
library(neuroim2)

test_that("Gaussian blobs are recovered across major clustering methods", {
  data <- generate_synthetic_volume(
    scenario = "gaussian_blobs",
    dims = c(12, 12, 6),
    n_clusters = 4,
    n_time = 25,
    spread = c(3, 3, 2),
    noise_sd = 0.03,
    seed = 2025
  )

  # Slice-MSF with z-smoothing
  msf <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    target_k_global = data$n_clusters,
    num_runs = 1,
    r = 6,
    compactness = 4,
    min_size = 8,
    stitch_z = TRUE,
    z_mult = 0.3,
    theta_link = 0.82
  )
  expect_gt(cluster_purity(msf$cluster, data$truth), 0.78)

  # Supervoxels via cluster4d wrapper
  sv <- cluster4d_supervoxels(
    vec = data$vec,
    mask = data$mask,
    n_clusters = data$n_clusters,
    spatial_weight = 0.45,
    max_iterations = 8,
    parallel = FALSE,
    sigma1 = 0.8,
    sigma2 = 2.0,
    converge_thresh = 0.01
  )
  expect_gt(cluster_purity(sv$cluster, data$truth), 0.6)

  # SNIC
  sn <- cluster4d_snic(
    vec = data$vec,
    mask = data$mask,
    n_clusters = data$n_clusters,
    spatial_weight = 0.5,
    max_iterations = 60
  )
  expect_gt(cluster_purity(sn$cluster, data$truth), 0.65)

  # FLASH-3D (hash-based supervoxels)
  f3 <- cluster4d_flash3d(
    vec = data$vec,
    mask = data$mask,
    n_clusters = data$n_clusters,
    spatial_weight = 0.55,
    max_iterations = 2,
    lambda_t = 0.9,
    bits = 64,
    dctM = 8,
    verbose = FALSE
  )
  expect_gt(cluster_purity(f3$cluster, data$truth), 0.6)
})

test_that("Layered scenario remains consistent across methods", {
  data <- generate_synthetic_volume(
    scenario = "z_layers",
    dims = c(10, 10, 9),
    n_clusters = 3,
    n_time = 20,
    noise_sd = 0.04,
    seed = 4242
  )

  msf <- slice_msf(
    vec = data$vec,
    mask = data$mask,
    target_k_global = data$n_clusters,
    num_runs = 1,
    r = 5,
    compactness = 3,
    min_size = 5,
    stitch_z = TRUE,
    theta_link = 0.9,
    z_mult = 0.2
  )
  expect_gt(cluster_purity(msf$cluster, data$truth), 0.75)

  sv <- cluster4d_supervoxels(
    vec = data$vec,
    mask = data$mask,
    n_clusters = data$n_clusters,
    spatial_weight = 0.3,
    max_iterations = 12,
    parallel = FALSE,
    sigma1 = 0.6,
    sigma2 = 1.2,
    converge_thresh = 0.005
  )
  expect_gt(cluster_purity(sv$cluster, data$truth), 0.6)
})
