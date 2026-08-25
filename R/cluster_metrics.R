# Shared, explicit clustering estimands. These helpers deliberately operate on
# observed contingency cells only: no N-by-N co-membership matrix and no dense
# K1-by-K2 contingency table is constructed.

.partition_label_codes <- function(labels, name) {
  if (is.factor(labels)) labels <- as.character(labels)
  if (!is.atomic(labels) || is.matrix(labels) || length(labels) < 1L) {
    stop(name, " must be a non-empty atomic vector", call. = FALSE)
  }
  if (anyNA(labels)) {
    stop(name, " must not contain missing labels", call. = FALSE)
  }
  if (is.numeric(labels) && any(!is.finite(labels))) {
    stop(name, " must contain only finite labels", call. = FALSE)
  }
  as.integer(match(labels, unique(labels)))
}

.partition_statistics <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b)) {
    stop("partition label vectors must have the same length", call. = FALSE)
  }
  a <- .partition_label_codes(labels_a, "labels_a")
  b <- .partition_label_codes(labels_b, "labels_b")
  n <- length(a)
  n_a <- max(a)
  n_b <- max(b)

  row_counts <- tabulate(a, nbins = n_a)
  col_counts <- tabulate(b, nbins = n_b)

  # With at most N distinct labels per partition, these double keys remain
  # exact for all practically representable R vectors. Only observed cells are
  # retained, so singleton partitions remain O(N), not O(N^2).
  keys <- (as.double(b) - 1) * n_a + a
  first <- !duplicated(keys)
  unique_keys <- keys[first]
  cell_id <- match(keys, unique_keys)
  cell_counts <- tabulate(cell_id, nbins = length(unique_keys))
  cell_a <- a[first]
  cell_b <- b[first]

  choose2 <- function(x) as.double(x) * (as.double(x) - 1) / 2
  same_a <- sum(choose2(row_counts))
  same_b <- sum(choose2(col_counts))
  same_both <- sum(choose2(cell_counts))

  entropy <- function(counts) {
    probability <- counts / n
    -sum(probability * log2(probability))
  }
  entropy_a <- entropy(row_counts)
  entropy_b <- entropy(col_counts)
  p_cell <- cell_counts / n
  mutual_information <- sum(
    p_cell * log2(
      p_cell /
        ((row_counts[cell_a] / n) * (col_counts[cell_b] / n))
    )
  )

  list(
    n = n,
    n_a = n_a,
    n_b = n_b,
    n_nonzero_cells = length(cell_counts),
    cell_counts = cell_counts,
    same_a = same_a,
    same_b = same_b,
    same_both = same_both,
    entropy_a = entropy_a,
    entropy_b = entropy_b,
    mutual_information = mutual_information
  )
}

.partition_metrics <- function(labels_a, labels_b) {
  stats <- .partition_statistics(labels_a, labels_b)
  total_pairs <- as.double(stats$n) * (stats$n - 1) / 2
  identical_relation <- stats$same_both == stats$same_a &&
    stats$same_both == stats$same_b

  if (total_pairs == 0) {
    ari <- 1
  } else {
    expected <- stats$same_a * stats$same_b / total_pairs
    maximum <- (stats$same_a + stats$same_b) / 2
    denominator <- maximum - expected
    ari <- if (abs(denominator) <= .Machine$double.eps) {
      if (identical_relation) 1 else 0
    } else {
      (stats$same_both - expected) / denominator
    }
  }

  dice_denominator <- stats$same_a + stats$same_b
  pairwise_dice <- if (dice_denominator == 0) {
    if (identical_relation) 1 else 0
  } else {
    2 * stats$same_both / dice_denominator
  }

  variation_information <- stats$entropy_a + stats$entropy_b -
    2 * stats$mutual_information
  if (abs(variation_information) < 64 * .Machine$double.eps) {
    variation_information <- 0
  }
  entropy_sum <- stats$entropy_a + stats$entropy_b
  nmi <- if (entropy_sum == 0) {
    if (identical_relation) 1 else 0
  } else {
    2 * stats$mutual_information / entropy_sum
  }

  list(
    ari = unname(ari),
    variation_of_information_bits = unname(variation_information),
    pairwise_dice = unname(pairwise_dice),
    nmi = unname(nmi),
    n_a = stats$n_a,
    n_b = stats$n_b,
    perfect = isTRUE(identical_relation)
  )
}

.matched_partition_accuracy <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b)) {
    stop("partition label vectors must have the same length", call. = FALSE)
  }
  a <- .partition_label_codes(labels_a, "labels_a")
  b <- .partition_label_codes(labels_b, "labels_b")
  n_a <- max(a)
  n_b <- max(b)
  keys <- (as.double(b) - 1) * n_a + a
  first <- !duplicated(keys)
  cell_counts <- tabulate(match(keys, keys[first]), nbins = sum(first))
  cell_a <- a[first]
  cell_b <- b[first]

  graph <- igraph::make_empty_graph(n = n_a + n_b, directed = FALSE)
  endpoints <- as.vector(rbind(cell_a, n_a + cell_b))
  graph <- igraph::add_edges(graph, endpoints)
  matching <- igraph::max_bipartite_match(
    graph,
    types = c(rep(FALSE, n_a), rep(TRUE, n_b)),
    weights = cell_counts
  )
  unname(matching$matching_weight / length(a))
}

.voxel_metric_matrix <- function(x, n_voxels, name, min_columns = 1L) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix", call. = FALSE)
  }
  if (nrow(x) != n_voxels) {
    if (ncol(x) == n_voxels && nrow(x) != n_voxels) {
      x <- t(x)
    } else {
      stop(name, " must have one row per labeled voxel", call. = FALSE)
    }
  }
  if (ncol(x) < min_columns) {
    stop(name, " must have at least ", min_columns, " columns", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values", call. = FALSE)
  }
  unname(x)
}

.spatial_rms_dispersion <- function(labels, coords) {
  codes <- .partition_label_codes(labels, "labels")
  coords <- .voxel_metric_matrix(coords, length(codes), "coords")
  counts <- tabulate(codes, nbins = max(codes))
  centroids <- rowsum(coords, codes, reorder = TRUE) / counts
  residuals <- coords - centroids[codes, , drop = FALSE]
  sqrt(mean(rowSums(residuals^2)))
}

.temporal_pairwise_correlation <- function(labels, feature_mat) {
  codes <- .partition_label_codes(labels, "labels")
  feature_mat <- .voxel_metric_matrix(
    feature_mat, length(codes), "feature_mat", min_columns = 2L
  )
  centered <- feature_mat - rowMeans(feature_mat)
  squared_norm <- rowSums(centered^2)
  tolerance <- 64 * .Machine$double.eps *
    pmax(1, rowSums(feature_mat^2))
  if (any(squared_norm <= tolerance)) {
    stop(
      "feature_mat rows must be non-constant for Pearson temporal coherence",
      call. = FALSE
    )
  }
  standardized <- centered / sqrt(squared_norm)
  counts <- tabulate(codes, nbins = max(codes))
  pair_counts <- as.double(counts) * (counts - 1) / 2
  total_pairs <- sum(pair_counts)
  if (total_pairs == 0) return(NA_real_)

  sums <- rowsum(standardized, codes, reorder = TRUE)
  pair_correlation_sums <- (rowSums(sums^2) - counts) / 2
  value <- sum(pair_correlation_sums[pair_counts > 0]) / total_pairs
  max(-1, min(1, unname(value)))
}

.cluster4d_geometry_equal <- function(x, y) {
  is.list(x) && is.list(y) &&
    identical(as.integer(x$dimensions), as.integer(y$dimensions)) &&
    .cluster4d_same_numeric(as.numeric(x$spacing), as.numeric(y$spacing)) &&
    .cluster4d_same_numeric(as.matrix(x$affine), as.matrix(y$affine))
}

.cluster4d_result_support <- function(result, name = "result") {
  if (!inherits(result, "cluster_result") &&
      !inherits(result, "cluster4d_result")) {
    stop(name, " must be a cluster4d_result", call. = FALSE)
  }
  if (!is.list(result$provenance) ||
      !is.list(result$provenance$coordinate_space) ||
      !is.list(result$provenance$coordinate_space$geometry) ||
      !is.list(result$provenance$mask)) {
    stop(name, " lacks coordinate and mask provenance", call. = FALSE)
  }
  if (!inherits(result$clusvol, "ClusteredNeuroVol")) {
    stop(name, "$clusvol must be a ClusteredNeuroVol", call. = FALSE)
  }
  mask <- methods::slot(result$clusvol, "mask")
  geometry <- .cluster4d_geometry(mask)
  if (!.cluster4d_geometry_equal(
    geometry, result$provenance$coordinate_space$geometry
  )) {
    stop(name, " coordinate provenance does not match clusvol geometry",
         call. = FALSE)
  }
  mask_values <- as.array(mask)
  if (any(!is.finite(mask_values))) {
    stop(name, " clusvol mask must be finite", call. = FALSE)
  }
  mask_idx <- which(mask_values > 0)
  labels <- .partition_label_codes(result$cluster, paste0(name, "$cluster"))
  if (length(labels) != length(mask_idx)) {
    stop(name, " labels do not match its included mask voxels", call. = FALSE)
  }
  if (!identical(
    as.integer(result$provenance$mask$n_voxels), as.integer(length(mask_idx))
  )) {
    stop(name, " mask provenance has the wrong voxel count", call. = FALSE)
  }
  list(
    labels = labels,
    mask_idx = mask_idx,
    geometry = geometry,
    coords = .cluster4d_index_to_coord(mask, mask_idx)
  )
}
