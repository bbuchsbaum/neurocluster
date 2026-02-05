library(testthat)
library(neurocluster)

context("ACSC C++ boundary refinement")

test_that("refine_boundaries_cpp respects normalized dot-product correlation", {
  skip_if_not(exists("refine_boundaries_cpp"), "C++ boundary refinement not available")
  skip_if_not(exists("find_boundary_voxels_cpp"), "C++ boundary voxel finder not available")

  # 4 voxels, 2 timepoints (voxels x time)
  feature_mat <- rbind(
    c(10, 0),  # voxel 1: cluster 1
    c(0, 10),  # voxel 2: should belong to cluster 2 but initially mislabeled
    c(0, 10),  # voxel 3: cluster 2
    c(0, 10)   # voxel 4: cluster 2
  )

  # Normalize rows (as refine_voxel_boundaries() does)
  feature_mat_norm <- feature_mat / sqrt(rowSums(feature_mat^2))

  labels <- c(1L, 1L, 2L, 2L)

  # Neighbor indices (N x 6), 1-based, with 0 for missing
  nbr <- matrix(0L, nrow = 4, ncol = 6)
  nbr[1, 1] <- 2L
  nbr[2, 1] <- 1L
  nbr[2, 2] <- 3L
  nbr[3, 1] <- 2L
  nbr[3, 2] <- 4L
  nbr[4, 1] <- 3L

  boundary <- find_boundary_voxels_cpp(labels, nbr)
  expect_true(2L %in% boundary)

  out <- refine_boundaries_cpp(
    voxel_labels = labels,
    feature_mat_normalized = feature_mat_norm,
    neighbor_indices = nbr,
    boundary_voxels = boundary,
    max_iter = 5L
  )

  expect_length(out$labels, length(labels))
  # Voxel 2 should switch to label 2 given perfect match to cluster 2 centroid
  expect_equal(out$labels[2], 2L)
})

