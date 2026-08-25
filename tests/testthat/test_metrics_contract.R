metric_contract_fixture <- function(labels = c(1L, 1L, 2L, 2L),
                                    features = rbind(
                                      c(1, 2, 3), c(2, 4, 6),
                                      c(1, 2, 1), c(3, 2, 3)
                                    ),
                                    mask_array = array(1, c(2, 2, 1)),
                                    method = "metric_fixture") {
  space <- neuroim2::NeuroSpace(dim(mask_array))
  mask <- neuroim2::NeuroVol(mask_array, space)
  mask_idx <- which(mask_array > 0)
  volumes <- lapply(seq_len(ncol(features)), function(j) {
    values <- array(0, dim(mask_array))
    values[mask_idx] <- features[, j]
    neuroim2::NeuroVol(values, space)
  })
  vec <- do.call(neuroim2::concat, volumes)
  raw <- structure(
    list(labels = labels, cluster = labels, parameters = list()),
    class = c("cluster4d_result", "cluster_result", "list")
  )
  result <- neurocluster:::finalize_cluster4d_result(
    raw, vec, mask, method, list(n_clusters = length(unique(labels)))
  )
  list(result = result, vec = vec, mask = mask, features = features)
}

test_that("partition estimands match hand-computed tiny oracles", {
  observed <- clustering_accuracy(c(1, 1, 2, 2), c(1, 1, 1, 2))

  expect_equal(observed$ari, 0)
  expect_equal(observed$variation_of_information_bits, 1.18872187554087)
  expect_equal(observed$pairwise_dice, 0.4)
  expect_equal(observed$nmi, 0.343711018485451)
})

test_that("matched accuracy uses the exact sparse assignment optimum", {
  predicted <- c(rep(1, 17), rep(2, 8))
  truth <- c(rep(1, 9), rep(2, 8), rep(1, 8))
  observed <- clustering_accuracy(predicted, truth)

  # Contingency cells are [9, 8; 8, 0]. Greedily taking 9 yields 9/25,
  # whereas the exact cross assignment yields (8 + 8)/25.
  expect_equal(observed$matched_accuracy, 16 / 25)
  expect_equal(observed$accuracy, observed$matched_accuracy)
})

test_that("partition estimands are permutation invariant and handle degeneracy", {
  a <- c(1, 1, 2, 2, 3, 3)
  b <- c(9, 9, 4, 4, 7, 7)
  perfect <- clustering_accuracy(a, b)

  expect_equal(perfect$ari, 1)
  expect_equal(perfect$variation_of_information_bits, 0)
  expect_equal(perfect$pairwise_dice, 1)
  expect_true(perfect$perfect)

  singleton <- clustering_accuracy(1, 99)
  expect_equal(singleton$ari, 1)
  expect_equal(singleton$variation_of_information_bits, 0)
  expect_equal(singleton$pairwise_dice, 1)

  opposite <- clustering_accuracy(rep(1, 4), seq_len(4))
  expect_equal(opposite$ari, 0)
  expect_equal(opposite$variation_of_information_bits, 2)
  expect_equal(opposite$pairwise_dice, 0)
})

test_that("partition estimands are symmetric, bounded, and chance adjusted", {
  set.seed(4201)
  a <- sample(rep(seq_len(8), each = 25))
  b <- sample(rep(seq_len(10), each = 20))
  ab <- clustering_accuracy(a, b)
  ba <- clustering_accuracy(b, a)

  expect_equal(ab$ari, ba$ari)
  expect_equal(ab$variation_of_information_bits,
               ba$variation_of_information_bits)
  expect_equal(ab$pairwise_dice, ba$pairwise_dice)
  expect_gte(ab$ari, -1)
  expect_lte(ab$ari, 1)
  expect_gte(ab$variation_of_information_bits, 0)
  expect_gte(ab$pairwise_dice, 0)
  expect_lte(ab$pairwise_dice, 1)

  chance_ari <- replicate(100, {
    clustering_accuracy(sample(a), sample(b))$ari
  })
  expect_lt(abs(mean(chance_ari)), 0.02)
})

test_that("ARI agrees with the trusted mclust implementation", {
  skip_if_not_installed("mclust")
  set.seed(4202)
  for (i in seq_len(20)) {
    a <- sample(seq_len(sample(2:12, 1)), 200, replace = TRUE)
    b <- sample(seq_len(sample(2:12, 1)), 200, replace = TRUE)
    expect_equal(
      clustering_accuracy(a, b)$ari,
      mclust::adjustedRandIndex(a, b),
      tolerance = 1e-12
    )
  }
})

test_that("partition computation stores only observed contingency cells", {
  n <- 50000L
  stats <- neurocluster:::.partition_statistics(seq_len(n), rev(seq_len(n)))

  expect_equal(length(stats$cell_counts), n)
  expect_equal(stats$n_nonzero_cells, n)
  expect_false(any(vapply(stats, is.matrix, logical(1))))
  expect_lt(as.numeric(object.size(stats)), 20 * n)
})

test_that("spatial dispersion and temporal coherence match exact voxel oracles", {
  labels <- c(1L, 1L, 2L, 2L)
  coords <- cbind(c(0, 2, 10, 14), 0, 0)
  features <- rbind(
    c(1, 2, 3), c(2, 4, 6),
    c(1, 2, 1), c(3, 2, 3)
  )

  expect_equal(
    neurocluster:::.spatial_rms_dispersion(labels, coords),
    sqrt(2.5)
  )
  expect_equal(
    neurocluster:::.temporal_pairwise_correlation(labels, features),
    0,
    tolerance = 1e-14
  )
})

test_that("cluster_metrics reports explicit voxel-level estimands", {
  fixture <- metric_contract_fixture()
  coords <- cbind(c(0, 2, 10, 14), 0, 0)
  metrics <- cluster_metrics(
    fixture$result,
    feature_mat = fixture$features,
    coords = coords,
    truth = c(5, 5, 8, 8)
  )

  expect_named(metrics, c(
    "n_clusters", "size_summary", "ari_truth",
    "variation_of_information_truth_bits", "pairwise_dice_truth",
    "temporal_pairwise_correlation", "spatial_rms_distance_mm"
  ))
  expect_equal(metrics$ari_truth, 1)
  expect_equal(metrics$temporal_pairwise_correlation, 0, tolerance = 1e-14)
  expect_equal(metrics$spatial_rms_distance_mm, sqrt(2.5))
})

test_that("temporal coherence comes from voxel features rather than centers", {
  fixture <- metric_contract_fixture()
  second <- fixture$result
  second$centers[,] <- rev(second$centers)

  comparison <- compare_cluster4d(
    first = fixture$result,
    second = second,
    metrics = "temporal_coherence",
    feature_mat = fixture$features
  )

  expect_equal(comparison$Temporal_Pairwise_Correlation, c(0, 0),
               tolerance = 1e-14)
})

test_that("comparison reports named linear-memory partition estimands", {
  first <- metric_contract_fixture(c(1L, 1L, 2L, 2L))
  second <- metric_contract_fixture(c(1L, 1L, 1L, 2L))

  comparison <- compare_cluster4d(
    first = first$result,
    second = second$result,
    metrics = c("summary", "spatial_dispersion", "partition_agreement")
  )

  expect_s3_class(comparison, "data.frame")
  expect_equal(comparison$Adjusted_Rand_Index, c(0, 0))
  expect_equal(comparison$Variation_of_Information_bits,
               rep(1.18872187554087, 2))
  expect_equal(comparison$Pairwise_Dice, c(0.4, 0.4))
  expect_true(all(is.finite(comparison$Spatial_RMS_Distance_mm)))
})

test_that("comparison rejects incompatible mask and coordinate provenance", {
  features <- rbind(
    c(1, 2, 3), c(2, 4, 6),
    c(1, 2, 1), c(3, 2, 3)
  )
  mask_a <- array(0, c(3, 2, 1))
  mask_b <- mask_a
  mask_a[c(1, 2, 3, 4)] <- 1
  mask_b[c(1, 2, 4, 5)] <- 1
  first <- metric_contract_fixture(mask_array = mask_a, features = features)
  second <- metric_contract_fixture(mask_array = mask_b, features = features)

  expect_error(
    compare_cluster4d(first$result, second$result, metrics = "summary"),
    "same included mask voxels"
  )

  bad_provenance <- first$result
  bad_provenance$provenance$coordinate_space$geometry$spacing[1] <- 9
  expect_error(
    compare_cluster4d(first$result, bad_provenance, metrics = "summary"),
    "coordinate provenance"
  )
})

test_that("metric inputs fail closed instead of silently changing estimands", {
  fixture <- metric_contract_fixture()

  expect_error(clustering_accuracy(c(1, NA), c(1, 2)), "missing")
  expect_error(clustering_accuracy(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(
    cluster_metrics(fixture$result, feature_mat = matrix(1, 4, 3)),
    "non-constant"
  )
  expect_error(
    compare_cluster4d(
      fixture$result, fixture$result,
      metrics = "temporal_coherence"
    ),
    "feature_mat"
  )
  expect_error(
    compare_cluster4d(
      fixture$result, fixture$result,
      metrics = "partition_agreement"
    ),
    NA
  )
  expect_error(
    compare_cluster4d(
      fixture$result, fixture$result, fixture$result,
      metrics = "partition_agreement"
    ),
    "exactly two"
  )
  expect_error(
    compare_cluster4d(fixture$result, metrics = "overlap"),
    "no longer supported"
  )
})
