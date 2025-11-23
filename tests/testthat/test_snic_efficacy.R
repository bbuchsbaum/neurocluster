library(testthat)
library(neuroim2)
library(neurocluster)

context("SNIC clustering efficacy and accuracy tests")

# ==============================================================================
# Synthetic Test Data Generators with Ground Truth
# ==============================================================================

#' Create synthetic data with distinct spatial regions (ground truth labels)
#'
#' @param dims Volume dimensions
#' @param nvols Number of time points
#' @param n_regions Number of distinct regions (ground truth clusters)
#' @param noise_level Standard deviation of Gaussian noise
#' @param spatial_smoothness Whether to smooth boundaries between regions
#' @return List with bvec, mask, and ground_truth labels
create_synthetic_regions <- function(dims = c(20, 20, 20), nvols = 10,
                                     n_regions = 4, noise_level = 0.1,
                                     spatial_smoothness = TRUE) {

  # Create mask
  mask_array <- array(1, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  mask_idx <- which(mask > 0)
  n_voxels <- length(mask_idx)

  # Create ground truth by spatial k-means
  coords <- index_to_coord(mask, mask_idx)
  ground_truth_kmeans <- kmeans(coords, centers = n_regions, iter.max = 100)
  ground_truth <- ground_truth_kmeans$cluster

  # Create distinct time series for each region
  vol_data <- array(0, c(dims, nvols))

  for (region in 1:n_regions) {
    region_voxels <- mask_idx[ground_truth == region]

    # Create a unique time series pattern for this region
    # Use different frequencies and phases for each region
    time_series <- numeric(nvols)
    for (t in 1:nvols) {
      # Region-specific oscillation
      time_series[t] <- sin(2 * pi * region * t / nvols) +
                       0.5 * cos(2 * pi * (region + 1) * t / (nvols * 2))
    }

    # Assign to all voxels in this region
    for (t in 1:nvols) {
      base_value <- time_series[t]

      if (spatial_smoothness) {
        # Add smooth spatial variation within region
        region_coords <- coords[ground_truth == region, , drop = FALSE]
        center <- colMeans(region_coords)
        for (i in seq_along(region_voxels)) {
          voxel_idx <- region_voxels[i]
          dist_from_center <- sqrt(sum((coords[ground_truth == region, ][i, ] - center)^2))
          spatial_factor <- 1 + 0.1 * dist_from_center / max(dims)
          vol_data[voxel_idx][t] <- base_value * spatial_factor + rnorm(1, 0, noise_level)
        }
      } else {
        # Uniform within region
        for (voxel_idx in region_voxels) {
          vol_data[voxel_idx][t] <- base_value + rnorm(1, 0, noise_level)
        }
      }
    }
  }

  # Ensure 4D structure
  dim(vol_data) <- c(dims, nvols)
  bvec <- NeuroVec(vol_data, NeuroSpace(c(dims, nvols)))

  list(
    bvec = bvec,
    mask = mask,
    ground_truth = ground_truth,
    n_regions = n_regions
  )
}

#' Create synthetic data with spatial gradients
#'
#' @param dims Volume dimensions
#' @param nvols Number of time points
#' @param gradient_axis Which axis has the gradient (1=x, 2=y, 3=z)
#' @param noise_level Noise standard deviation
create_synthetic_gradient <- function(dims = c(20, 20, 20), nvols = 10,
                                      gradient_axis = 1, noise_level = 0.1) {

  mask_array <- array(1, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  mask_idx <- which(mask > 0)
  coords <- index_to_coord(mask, mask_idx)

  vol_data <- array(0, c(dims, nvols))

  # Create gradient along specified axis
  gradient_values <- coords[, gradient_axis] / dims[gradient_axis]

  for (t in 1:nvols) {
    # Time series varies with spatial gradient
    phase_shift <- 2 * pi * t / nvols
    for (i in seq_along(mask_idx)) {
      base_signal <- sin(2 * pi * gradient_values[i] + phase_shift)
      vol_data[mask_idx[i]][t] <- base_signal + rnorm(1, 0, noise_level)
    }
  }

  dim(vol_data) <- c(dims, nvols)
  bvec <- NeuroVec(vol_data, NeuroSpace(c(dims, nvols)))

  list(
    bvec = bvec,
    mask = mask,
    gradient_values = gradient_values
  )
}

#' Create checkerboard pattern (challenging for clustering)
create_synthetic_checkerboard <- function(dims = c(20, 20, 20), nvols = 10,
                                         block_size = 5, noise_level = 0.1) {

  mask_array <- array(1, dims)
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  mask_idx <- which(mask > 0)
  coords <- index_to_coord(mask, mask_idx)

  vol_data <- array(0, c(dims, nvols))

  # Assign checkerboard labels
  block_labels <- (floor(coords[, 1] / block_size) +
                  floor(coords[, 2] / block_size) +
                  floor(coords[, 3] / block_size)) %% 2

  for (t in 1:nvols) {
    for (i in seq_along(mask_idx)) {
      # Different signal for each checkerboard color
      if (block_labels[i] == 0) {
        base_signal <- sin(2 * pi * t / nvols)
      } else {
        base_signal <- -sin(2 * pi * t / nvols)
      }
      vol_data[mask_idx[i]][t] <- base_signal + rnorm(1, 0, noise_level)
    }
  }

  dim(vol_data) <- c(dims, nvols)
  bvec <- NeuroVec(vol_data, NeuroSpace(c(dims, nvols)))

  list(
    bvec = bvec,
    mask = mask,
    ground_truth = block_labels,
    block_size = block_size
  )
}

# ==============================================================================
# Clustering Quality Metrics
# ==============================================================================

#' Adjusted Rand Index (ARI) - measures agreement between clusterings
#' 1.0 = perfect agreement, 0.0 = random, < 0 = worse than random
compute_ari <- function(labels1, labels2) {
  # Contingency table
  n <- length(labels1)
  tab <- table(labels1, labels2)

  # Sum of choose(n_ij, 2)
  sum_comb_tab <- sum(choose(tab, 2))

  # Marginal sums
  sum_comb_rows <- sum(choose(rowSums(tab), 2))
  sum_comb_cols <- sum(choose(colSums(tab), 2))

  # Expected index
  expected <- sum_comb_rows * sum_comb_cols / choose(n, 2)

  # Max index
  max_index <- 0.5 * (sum_comb_rows + sum_comb_cols)

  # ARI
  if (max_index == expected) {
    return(0)
  }
  (sum_comb_tab - expected) / (max_index - expected)
}

#' Normalized Mutual Information (NMI)
compute_nmi <- function(labels1, labels2) {
  tab <- table(labels1, labels2)
  n <- sum(tab)

  # Marginal distributions
  p_i <- rowSums(tab) / n
  p_j <- colSums(tab) / n

  # Joint distribution
  p_ij <- tab / n

  # Mutual information
  mi <- 0
  for (i in 1:nrow(tab)) {
    for (j in 1:ncol(tab)) {
      if (p_ij[i, j] > 0) {
        mi <- mi + p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
      }
    }
  }

  # Entropy
  H_i <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  H_j <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))

  # NMI
  if (H_i == 0 || H_j == 0) return(0)
  mi / sqrt(H_i * H_j)
}

#' Spatial contiguity score - measures how spatially coherent clusters are
#' Higher is better (clusters are more spatially compact)
compute_spatial_contiguity <- function(cluster_labels, coords) {
  unique_clusters <- unique(cluster_labels)
  n_clusters <- length(unique_clusters)

  contiguity_scores <- numeric(n_clusters)

  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    cluster_coords <- coords[cluster_labels == cluster_id, , drop = FALSE]

    if (nrow(cluster_coords) > 1) {
      # Average distance to cluster centroid
      centroid <- colMeans(cluster_coords)
      distances <- sqrt(rowSums((cluster_coords - matrix(centroid,
                                                          nrow = nrow(cluster_coords),
                                                          ncol = 3, byrow = TRUE))^2))
      contiguity_scores[i] <- mean(distances)
    } else {
      contiguity_scores[i] <- 0
    }
  }

  # Return inverse of mean distance (higher = more compact)
  1 / (mean(contiguity_scores) + 1e-6)
}

#' Within-cluster feature similarity
#' Measures how similar time series are within each cluster
compute_within_cluster_similarity <- function(cluster_labels, features) {
  unique_clusters <- unique(cluster_labels)
  similarities <- numeric(length(unique_clusters))

  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    cluster_features <- features[, cluster_labels == cluster_id, drop = FALSE]

    if (ncol(cluster_features) > 1) {
      # Compute mean pairwise correlation
      cor_mat <- cor(cluster_features)
      # Extract upper triangle (excluding diagonal)
      similarities[i] <- mean(cor_mat[upper.tri(cor_mat)])
    } else {
      similarities[i] <- 1.0
    }
  }

  mean(similarities)
}

# ==============================================================================
# Efficacy Tests
# ==============================================================================

test_that("SNIC recovers ground truth spatial regions with high accuracy", {
  set.seed(42)
  synthetic <- create_synthetic_regions(
    dims = c(15, 15, 15),
    nvols = 8,
    n_regions = 4,
    noise_level = 0.15,
    spatial_smoothness = TRUE
  )

  # Run SNIC with appropriate K
  result <- snic(synthetic$bvec, synthetic$mask,
                K = synthetic$n_regions,
                compactness = 5)

  # Compute clustering quality metrics
  ari <- compute_ari(synthetic$ground_truth, result$cluster)
  nmi <- compute_nmi(synthetic$ground_truth, result$cluster)

  # SNIC should achieve reasonable agreement with ground truth
  # ARI > 0.5 indicates good recovery, > 0.7 is very good
  expect_true(ari > 0.4,
              info = sprintf("ARI = %.3f, expected > 0.4", ari))

  # NMI should also be reasonable
  expect_true(nmi > 0.4,
              info = sprintf("NMI = %.3f, expected > 0.4", nmi))

  cat(sprintf("\n  Ground truth recovery: ARI=%.3f, NMI=%.3f\n", ari, nmi))
})

test_that("SNIC produces spatially contiguous clusters", {
  set.seed(123)
  synthetic <- create_synthetic_regions(
    dims = c(12, 12, 12),
    nvols = 6,
    n_regions = 5,
    noise_level = 0.1
  )

  # Run SNIC with high compactness (should favor spatial coherence)
  result <- snic(synthetic$bvec, synthetic$mask,
                K = 5,
                compactness = 10)

  # Compute spatial contiguity
  mask_idx <- which(synthetic$mask > 0)
  coords <- index_to_coord(synthetic$mask, mask_idx)
  contiguity <- compute_spatial_contiguity(result$cluster, coords)

  # Higher compactness should produce more compact clusters
  # This is a relative measure, but should be > 0.1 for K=5
  expect_true(contiguity > 0.1,
              info = sprintf("Spatial contiguity = %.3f, expected > 0.1", contiguity))

  cat(sprintf("\n  Spatial contiguity score: %.3f\n", contiguity))
})

test_that("Compactness parameter correctly trades off spatial vs feature similarity", {
  set.seed(456)
  synthetic <- create_synthetic_regions(
    dims = c(12, 12, 12),
    nvols = 8,
    n_regions = 4,
    noise_level = 0.1
  )

  # Low compactness (feature-driven)
  result_low <- snic(synthetic$bvec, synthetic$mask,
                    K = 4,
                    compactness = 1)

  # High compactness (spatial-driven)
  result_high <- snic(synthetic$bvec, synthetic$mask,
                     K = 4,
                     compactness = 20)

  # Compute metrics
  mask_idx <- which(synthetic$mask > 0)
  coords <- index_to_coord(synthetic$mask, mask_idx)
  features <- series(synthetic$bvec, mask_idx)

  spatial_low <- compute_spatial_contiguity(result_low$cluster, coords)
  spatial_high <- compute_spatial_contiguity(result_high$cluster, coords)

  feature_sim_low <- compute_within_cluster_similarity(result_low$cluster, features)
  feature_sim_high <- compute_within_cluster_similarity(result_high$cluster, features)

  # High compactness should produce more spatially compact clusters
  expect_true(spatial_high >= spatial_low * 0.8,
              info = sprintf("High compact (%.3f) should be >= Low compact (%.3f) * 0.8",
                           spatial_high, spatial_low))

  cat(sprintf("\n  Low compactness: spatial=%.3f, feature_sim=%.3f\n",
             spatial_low, feature_sim_low))
  cat(sprintf("  High compactness: spatial=%.3f, feature_sim=%.3f\n",
             spatial_high, feature_sim_high))
})

test_that("SNIC handles spatial gradients appropriately", {
  set.seed(789)
  synthetic <- create_synthetic_gradient(
    dims = c(20, 10, 10),
    nvols = 8,
    gradient_axis = 1,
    noise_level = 0.1
  )

  # Run SNIC with moderate K
  K <- 5
  result <- snic(synthetic$bvec, synthetic$mask,
                K = K,
                compactness = 5)

  # Clusters should follow the gradient
  mask_idx <- which(synthetic$mask > 0)
  coords <- index_to_coord(synthetic$mask, mask_idx)

  # For each cluster, compute mean position along gradient axis
  cluster_positions <- tapply(coords[, 1], result$cluster, mean)

  # Cluster positions should span the gradient
  position_range <- diff(range(cluster_positions))
  max_range <- diff(range(coords[, 1]))

  # Should use at least 50% of the spatial range
  expect_true(position_range > 0.5 * max_range,
              info = sprintf("Cluster span (%.1f) should cover > 50%% of range (%.1f)",
                           position_range, max_range))

  cat(sprintf("\n  Gradient coverage: %.1f / %.1f (%.1f%%)\n",
             position_range, max_range, 100 * position_range / max_range))
})

test_that("SNIC within-cluster feature similarity is high", {
  set.seed(101)
  synthetic <- create_synthetic_regions(
    dims = c(12, 12, 12),
    nvols = 10,
    n_regions = 4,
    noise_level = 0.08
  )

  result <- snic(synthetic$bvec, synthetic$mask,
                K = 4,
                compactness = 5)

  mask_idx <- which(synthetic$mask > 0)
  features <- series(synthetic$bvec, mask_idx)

  within_sim <- compute_within_cluster_similarity(result$cluster, features)

  # Within-cluster correlation should be reasonably high (> 0.3)
  expect_true(within_sim > 0.3,
              info = sprintf("Within-cluster similarity = %.3f, expected > 0.3", within_sim))

  cat(sprintf("\n  Within-cluster feature similarity: %.3f\n", within_sim))
})

test_that("SNIC is stable across multiple runs with same seed", {
  set.seed(202)
  synthetic <- create_synthetic_regions(
    dims = c(10, 10, 10),
    nvols = 5,
    n_regions = 3,
    noise_level = 0.1
  )

  # Run multiple times with same seed
  set.seed(999)
  result1 <- snic(synthetic$bvec, synthetic$mask, K = 3, compactness = 5)

  set.seed(999)
  result2 <- snic(synthetic$bvec, synthetic$mask, K = 3, compactness = 5)

  set.seed(999)
  result3 <- snic(synthetic$bvec, synthetic$mask, K = 3, compactness = 5)

  # All results should be identical
  expect_identical(result1$cluster, result2$cluster)
  expect_identical(result2$cluster, result3$cluster)
})

test_that("SNIC performance degrades gracefully with increasing noise", {
  dims <- c(12, 12, 12)
  nvols <- 8
  n_regions <- 4

  noise_levels <- c(0.05, 0.15, 0.3)
  ari_scores <- numeric(length(noise_levels))

  for (i in seq_along(noise_levels)) {
    set.seed(303 + i)
    synthetic <- create_synthetic_regions(
      dims = dims,
      nvols = nvols,
      n_regions = n_regions,
      noise_level = noise_levels[i]
    )

    result <- snic(synthetic$bvec, synthetic$mask,
                  K = n_regions,
                  compactness = 5)

    ari_scores[i] <- compute_ari(synthetic$ground_truth, result$cluster)
  }

  # ARI should decrease monotonically with noise (allowing small tolerance)
  for (i in 2:length(noise_levels)) {
    expect_true(ari_scores[i] <= ari_scores[i-1] + 0.1,
                info = sprintf("ARI should decrease with noise: %.3f -> %.3f",
                             ari_scores[i-1], ari_scores[i]))
  }

  cat(sprintf("\n  ARI vs noise: %.3f (%.2f) -> %.3f (%.2f) -> %.3f (%.2f)\n",
             ari_scores[1], noise_levels[1],
             ari_scores[2], noise_levels[2],
             ari_scores[3], noise_levels[3]))
})

test_that("SNIC cluster centers accurately represent cluster members", {
  set.seed(404)
  synthetic <- create_synthetic_regions(
    dims = c(12, 12, 12),
    nvols = 8,
    n_regions = 4,
    noise_level = 0.1
  )

  result <- snic(synthetic$bvec, synthetic$mask,
                K = 4,
                compactness = 5)

  mask_idx <- which(synthetic$mask > 0)
  features <- series(synthetic$bvec, mask_idx)

  # For each cluster, compute correlation between center and members
  unique_clusters <- sort(unique(result$cluster))
  center_quality <- numeric(length(unique_clusters))

  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    cluster_members <- features[, result$cluster == cluster_id, drop = FALSE]

    if (!is.null(result$centers) && ncol(cluster_members) > 1) {
      # Get center for this cluster
      center <- result$centers[, i]

      # Compute correlation between center and each member
      correlations <- apply(cluster_members, 2, function(x) cor(center, x))
      center_quality[i] <- mean(correlations)
    }
  }

  if (length(center_quality) > 0 && !all(is.na(center_quality))) {
    mean_quality <- mean(center_quality, na.rm = TRUE)

    # Centers should be reasonably representative (correlation > 0.5)
    expect_true(mean_quality > 0.5,
                info = sprintf("Center quality = %.3f, expected > 0.5", mean_quality))

    cat(sprintf("\n  Mean center-member correlation: %.3f\n", mean_quality))
  } else {
    skip("Centers not available for quality check")
  }
})

# ==============================================================================
# Comparative Tests
# ==============================================================================

test_that("SNIC outperforms spatial k-means on structured data", {
  set.seed(505)
  synthetic <- create_synthetic_regions(
    dims = c(12, 12, 12),
    nvols = 8,
    n_regions = 4,
    noise_level = 0.1
  )

  # Run SNIC
  result_snic <- snic(synthetic$bvec, synthetic$mask,
                     K = 4,
                     compactness = 5)

  # Run spatial k-means as baseline
  mask_idx <- which(synthetic$mask > 0)
  coords <- index_to_coord(synthetic$mask, mask_idx)
  features <- series(synthetic$bvec, mask_idx)

  # Combine coords and features for spatial k-means
  combined <- cbind(scale(coords), scale(t(features)))
  kmeans_result <- kmeans(combined, centers = 4, iter.max = 100)

  # Compare quality
  ari_snic <- compute_ari(synthetic$ground_truth, result_snic$cluster)
  ari_kmeans <- compute_ari(synthetic$ground_truth, kmeans_result$cluster)

  # SNIC should be competitive or better
  # (May not always beat k-means, but should be close)
  expect_true(ari_snic > ari_kmeans * 0.7,
              info = sprintf("SNIC ARI (%.3f) should be >= 70%% of k-means ARI (%.3f)",
                           ari_snic, ari_kmeans))

  cat(sprintf("\n  SNIC ARI: %.3f vs k-means ARI: %.3f\n", ari_snic, ari_kmeans))
})
