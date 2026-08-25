#' Adaptive Correlation Superclustering (ACSC)
#'
#' Clusters fMRI voxels into spatially-coherent groups based on temporal correlation
#' and spatial proximity. Includes optional refinement for boundary corrections.
#' The algorithm supports parallel processing via the future package for improved performance.
#'
#' @param bvec A NeuroVec-like object containing 4D fMRI data.
#' @param mask A NeuroVol-like object (logical or numeric mask).
#' @param block_size Approximate side length of blocks (e.g., 2 or 3). Must be > 0.
#' @param ann_k Number of approximate (or exact) nearest neighbors per block. Must be >= 1.
#' @param alpha Feature-similarity weight in `[0, 1]`. Zero is spatial-only and
#'   one is feature-only; both endpoints are supported.
#' @param correlation_metric Similarity definition: Pearson correlation,
#'   Spearman correlation of average ranks, or robust biweight midcorrelation
#'   with tuning constant 9.
#' @param spatial_weighting Spatial adjacency weighting ("gaussian", "binary").
#' @param refine Logical; whether to refine boundaries.
#' @param max_refine_iter Maximum iterations for boundary refinement. Must be >= 0.
#' @param K (Optional) Desired number of clusters.
#' @param knn_proj_dim Optional integer; if > 0, project block summaries to this
#'   dimension before kNN graph construction (speeds up \code{FNN::get.knn()} for
#'   long time-series).
#' @param knn_proj_method Projection method for \code{knn_proj_dim}: one of
#'   \code{"none"}, \code{"dct"}, or \code{"rp"} (random projection).
#' @param knn_proj_seed Integer seed used when \code{knn_proj_method="rp"} (RNG
#'   state is restored after graph construction).
#'
#' @return A standardized `cluster4d_result` (also classed `acsc_result`) with
#'   final contiguous labels, original-space centers, physical coordinate
#'   centers, provenance, and these ACSC-specific elements:
#'   \describe{
#'     \item{cluster_map}{3D array with cluster labels per voxel.}
#'     \item{graph}{An \code{igraph} object used for clustering.}
#'     \item{init_block_label}{Initial coarse partition (3D array) matching \code{mask} dimensions.}
#'   }
#'
#' @details
#' ACSC builds feature-neighbor edges from the declared correlation metric and
#' unions them with six-neighbor block-grid edges. Edge correlations are mapped
#' to non-negative Louvain similarities as `(correlation + 1) / 2`. At
#' `alpha = 0`, feature-neighbor selection is bypassed entirely, so the graph
#' topology and weights depend only on spatial adjacency. Optional DCT or random
#' projection changes feature-neighbor candidate selection, but edge weights are
#' always recomputed using the declared metric.
#'
#' Boundary refinement uses the same metric representation as graph building.
#' Pearson uses centered unit vectors, Spearman uses centered average ranks, and
#' robust mode uses median-centered Tukey-bisquare residuals. The C++ path
#' accelerates centroid dot products and recomputes the active boundary after
#' every synchronous pass; failures fall back to the equivalent R path.
#'
#' After refinement, labels are flood-filled into connected components. If `K`
#' is supplied, the shared adjacency-preserving exact-K engine merges or splits
#' those components and fails closed when the target is topologically infeasible.
#'
#' @importFrom FNN get.knn
#' @importFrom igraph make_empty_graph graph_from_data_frame cluster_louvain communities membership
#' @importFrom future.apply future_lapply future_apply
#' @export
acsc <- function(bvec, mask,
                 block_size = 2,
                 ann_k = 10,
                 alpha = 0.5,
                 correlation_metric = c("pearson", "spearman", "robust"),
                 spatial_weighting = c("gaussian", "binary"),
                 refine = TRUE,
                 max_refine_iter = 5,
                 K = NULL,
                 knn_proj_dim = 0L,
                 knn_proj_method = c("none", "dct", "rp"),
                 knn_proj_seed = 1L)
{
  ## ------------------------------------------------------------------------
  ## 0. Basic Validation
  ## ------------------------------------------------------------------------
  correlation_metric <- match.arg(correlation_metric)
  spatial_weighting  <- match.arg(spatial_weighting)
  refine <- .cluster4d_scalar_logical(refine, "refine", "acsc")
  block_size <- .cluster4d_scalar_number(
    block_size, "block_size", "acsc", lower = .Machine$double.eps
  )
  ann_k <- .cluster4d_scalar_number(
    ann_k, "ann_k", "acsc", lower = 1, integer = TRUE
  )
  alpha <- .cluster4d_scalar_number(
    alpha, "alpha", "acsc", lower = 0, upper = 1
  )
  max_refine_iter <- .cluster4d_scalar_number(
    max_refine_iter, "max_refine_iter", "acsc", lower = 0, integer = TRUE
  )

  knn_proj_method <- match.arg(knn_proj_method)
  knn_proj_dim <- .cluster4d_scalar_number(
    knn_proj_dim, "knn_proj_dim", "acsc", lower = 0, integer = TRUE
  )
  knn_proj_seed <- .cluster4d_scalar_number(
    knn_proj_seed, "knn_proj_seed", "acsc",
    lower = -.Machine$integer.max, upper = .Machine$integer.max, integer = TRUE
  )
  
  # Validate K parameter
  if (!is.null(K)) {
    K <- .cluster4d_scalar_number(K, "K", "acsc", lower = 1, integer = TRUE)
  }

  if (!inherits(mask, "NeuroVol")) {
    stop("mask must be a NeuroVol-like object.")
  }
  if (!inherits(bvec, "NeuroVec")) {
    stop("bvec must be a NeuroVec-like object.")
  }
  if (!all(dim(mask) >= 1)) {
    stop("mask has invalid dimensions.")
  }

  mask_values <- as.array(mask)
  if (any(!is.finite(mask_values))) {
    stop("acsc: mask values must be finite; included voxels are exactly values > 0")
  }

  mask.idx <- which(mask_values > 0)
  if (length(mask.idx) == 0) {
    stop("No nonzero voxels in mask. Nothing to cluster.")
  }
  if (!is.null(K) && K > length(mask.idx)) {
    stop("acsc: K must not exceed the number of included mask voxels")
  }

  ## ------------------------------------------------------------------------
  ## 1. Preprocess data
  ## ------------------------------------------------------------------------
  feature_mat <- preprocess_time_series(bvec, mask, correlation_metric)

  ## coords is Nx3 integer matrix of voxel (i, j, k)
  coords <- index_to_grid(mask, mask.idx)

  ## ------------------------------------------------------------------------
  ## 2. Coarse partition into blocks
  ## ------------------------------------------------------------------------
  block_id <- block_partition(coords, block_size)
  block_summary <- summarize_blocks(
    feature_mat, coords, block_id, block_size = block_size
  )
  
  # Validate ann_k parameter
  nb <- length(unique(block_id))
  # ann_k must be strictly less than nb for FNN::get.knn
  if (nb > 1L && ann_k >= nb) {
    ann_k <- max(1, nb - 1)
  }

  ## ------------------------------------------------------------------------
  ## 3. Construct adjacency graph
  ## ------------------------------------------------------------------------
  graph <- if (nb == 1L) {
    igraph::make_empty_graph(n = 1L, directed = FALSE)
  } else {
    build_acsc_graph(
      block_summary,
      ann_k,
      alpha,
      spatial_weighting,
      block_size,
      correlation_metric = correlation_metric,
      knn_proj_dim = knn_proj_dim,
      knn_proj_method = knn_proj_method,
      knn_proj_seed = knn_proj_seed
    )
  }

  ## ------------------------------------------------------------------------
  ## 4. Run clustering (Louvain)
  ## ------------------------------------------------------------------------
  if (nb == 1L) {
    cluster_result <- 1L
  } else if (igraph::ecount(graph) == 0L || (!is.null(K) && K >= nb)) {
    cluster_result <- seq_len(nb)
  } else if (!is.null(K)) {
    res_est <- estimate_resolution(min(K, nb), graph)
    cluster_result <- run_louvain_clustering(graph, resolution = res_est)
  } else {
    cluster_result <- run_louvain_clustering(graph)
  }

  ## ------------------------------------------------------------------------
  ## 5. Map block labels back to voxels
  ## ------------------------------------------------------------------------
  voxel_labels <- expand_block_labels(cluster_result, block_id, mask.idx)

  ## ------------------------------------------------------------------------
  ## 6. (Optional) Refine boundaries
  ## ------------------------------------------------------------------------
  if (refine && max_refine_iter > 0) {
    voxel_labels <- refine_voxel_boundaries(
      voxel_labels,
      feature_mat,
      coords,
      max_refine_iter,
      dims = dim(mask),
      correlation_metric = correlation_metric
    )
  }

  # Refinement can split a label into disconnected islands. Canonicalize every
  # connected component before any target-K merge/split operation.
  graph_info <- .exact_k_graph(mask, 6L)
  voxel_labels <- .exact_k_connected_labels(
    voxel_labels, graph_info$graph, graph_info$edges
  )

  ## ------------------------------------------------------------------------
  ## 6b. Force exact K if requested (merge/split as needed for downstream ARI)
  ## ------------------------------------------------------------------------
  if (!is.null(K)) {
    voxel_labels <- force_exact_k(
      voxel_labels, feature_mat, K,
      mask = mask, connectivity = 6L
    )
  }

  ## ------------------------------------------------------------------------
  ## 7. Construct output as cluster_result
  ## ------------------------------------------------------------------------
  cluster_map <- array(0L, dim(mask))
  cluster_map[mask.idx] <- voxel_labels

  # Safety: ensure we produced labels for all voxels
  if (length(voxel_labels) != length(mask.idx)) {
    stop("ACSC: voxel_labels length mismatch; expected ", length(mask.idx),
         " got ", length(voxel_labels))
  }

  # Create cluster_result compatible structure
  result <- structure(
    list(
      cluster = voxel_labels,
      clusvol = suppressWarnings(
        neuroim2::ClusteredNeuroVol(mask, clusters = voxel_labels)
      ),
      n_clusters = length(unique(voxel_labels[voxel_labels > 0])),
      cluster_map = cluster_map,
      method = "acsc",
      parameters = list(
        alpha = alpha,
        block_size = block_size,
        ann_k = ann_k,
        knn_proj_dim = knn_proj_dim,
        knn_proj_method = knn_proj_method,
        correlation_metric = correlation_metric[1],
        spatial_weighting = spatial_weighting[1],
        refine = refine,
        max_refine_iter = max_refine_iter,
        K = K,
        knn_proj_seed = knn_proj_seed
      ),
      metadata = list(
        graph = graph,
        init_block_label = construct_block_label_array(block_id, mask),
        correlation = list(
          metric = correlation_metric,
          similarity = switch(
            correlation_metric,
            pearson = "Pearson correlation",
            spearman = "Pearson correlation of average ranks",
            robust = "biweight midcorrelation with tuning constant 9"
          ),
          graph_weight_mapping = "(correlation + 1) / 2",
          alpha_semantics = "feature similarity weight; 0 is spatial-only and 1 is feature-only"
        ),
        topology = list(
          block_connectivity = 6L,
          final_voxel_connectivity = 6L,
          exact_k = !is.null(K)
        )
      )
    ),
    class = c("acsc_result", "cluster4d_result", "cluster_result", "list")
  )

  finalize_cluster4d_result(result, bvec, mask, "acsc", result$parameters)
}

#' Preprocess fMRI time-series data
#'
#' @keywords internal
preprocess_time_series <- function(bvec, mask, correlation_metric) {
  mask.idx <- which(mask > 0)
  feature_mat <- neuroim2::series(bvec, mask.idx)  # returns T x N
  # Ensure rows = voxels, cols = timepoints
  feature_mat <- t(as.matrix(feature_mat))  # now N x T
  feature_mat0 <- feature_mat
  
  # Handle single voxel case
  if (nrow(feature_mat) == 1) {
    feature_mat <- matrix(feature_mat, nrow = 1)
    feature_mat0 <- feature_mat
  }
  
  # Handle zero-length case
  if (length(feature_mat) == 0) {
    stop("No time series data extracted from the input.")
  }

  # Fast path: use C++ (volume mean-centering + per-voxel linear detrend).
  # This matches the prior behavior of `scale(..., center=TRUE, scale=FALSE)` (column centering)
  # followed by per-row linear detrending.
  if (ncol(feature_mat) > 1 && exists("normalize_detrend_cpp", mode = "function")) {
    feature_mat <- tryCatch(
      normalize_detrend_cpp(feature_mat0),
      error = function(e) {
        # Fallback to pure-R below
        NULL
      }
    )
  } else {
    feature_mat <- NULL
  }

  if (is.null(feature_mat)) {
    ## 1) Mean-center volumes (center each timepoint across voxels)
    feature_mat <- base::scale(feature_mat0, center = TRUE, scale = FALSE)

    ## 2) Detrend each voxel across time
    detrend_one <- function(x) stats::lm(x ~ seq_along(x))$residuals
    n_voxels <- nrow(feature_mat)

    if (requireNamespace("future", quietly = TRUE) && future::nbrOfWorkers() > 1 && n_voxels >= 5000) {
      feature_mat <- suppressWarnings(future.apply::future_apply(feature_mat, 1, detrend_one))
    } else {
      feature_mat <- apply(feature_mat, 1, detrend_one)
    }
    feature_mat <- t(feature_mat)  # (voxels x time)
  }

  feature_mat <- unname(as.matrix(feature_mat))
  if (!identical(dim(feature_mat), dim(feature_mat0))) {
    stop(
      "ACSC preprocessing changed the voxel-by-time shape from ",
      paste(dim(feature_mat0), collapse = "x"), " to ",
      paste(dim(feature_mat), collapse = "x"),
      call. = FALSE
    )
  }
  if (any(!is.finite(feature_mat))) {
    stop("ACSC preprocessing produced non-finite values", call. = FALSE)
  }
  feature_mat
}

# Convert each time series into a unit-vector representation whose dot product
# is the declared similarity. Robust mode is the biweight midcorrelation
# representation: median-centered residuals receive Tukey bisquare weights
# with tuning constant 9 before unit normalization.
acsc_metric_representation <- function(feature_mat,
                                       correlation_metric = c("pearson", "spearman", "robust")) {
  correlation_metric <- match.arg(correlation_metric)
  feature_mat <- unname(as.matrix(feature_mat))
  if (!is.numeric(feature_mat) || nrow(feature_mat) < 1L || ncol(feature_mat) < 2L) {
    stop("ACSC similarity input must be a nonempty numeric N by T matrix with T >= 2")
  }
  if (any(!is.finite(feature_mat))) {
    stop("ACSC similarity input must contain only finite values")
  }

  transformed <- switch(
    correlation_metric,
    pearson = sweep(feature_mat, 1L, rowMeans(feature_mat), "-"),
    spearman = {
      ranked <- t(apply(feature_mat, 1L, rank, ties.method = "average"))
      sweep(ranked, 1L, rowMeans(ranked), "-")
    },
    robust = {
      medians <- apply(feature_mat, 1L, stats::median)
      residuals <- sweep(feature_mat, 1L, medians, "-")
      scales <- apply(abs(residuals), 1L, stats::median)
      safe_scales <- scales
      safe_scales[safe_scales <= sqrt(.Machine$double.eps)] <- 1
      u <- residuals / (9 * safe_scales)
      weights <- (1 - pmin(u^2, 1))^2
      weights[abs(u) >= 1] <- 0
      transformed <- residuals * weights
      transformed[scales <= sqrt(.Machine$double.eps), ] <- 0
      transformed
    }
  )

  norms <- sqrt(rowSums(transformed^2))
  live <- norms > sqrt(.Machine$double.eps)
  transformed[live, ] <- transformed[live, , drop = FALSE] / norms[live]
  transformed[!live, ] <- 0
  unname(transformed)
}

acsc_similarity_matrix <- function(feature_mat,
                                   correlation_metric = c("pearson", "spearman", "robust")) {
  representation <- acsc_metric_representation(feature_mat, correlation_metric)
  similarity <- tcrossprod(representation)
  pmax(pmin(similarity, 1), -1)
}

#' Convert linear indices to 3D grid coordinates
#' @keywords internal
index_to_grid <- function(mask, indices) {
  dims <- dim(mask)
  k <- (indices - 1) %/% (dims[1] * dims[2]) + 1
  j <- ((indices - 1) %% (dims[1] * dims[2])) %/% dims[1] + 1
  i <- (indices - 1) %% dims[1] + 1
  cbind(i, j, k)
}

#' Partition voxel coordinates into coarse blocks
#' @keywords internal
block_partition <- function(coords, block_size) {
  # Shift coords so min is 0
  min_x <- min(coords[, 1])
  min_y <- min(coords[, 2])
  min_z <- min(coords[, 3])

  shifted_x <- coords[, 1] - min_x
  shifted_y <- coords[, 2] - min_y
  shifted_z <- coords[, 3] - min_z

  bx <- shifted_x %/% block_size
  by <- shifted_y %/% block_size
  bz <- shifted_z %/% block_size

  max_bx <- max(bx) + 1
  max_by <- max(by) + 1
  block_id <- bx + by * max_bx + bz * max_bx * max_by
  block_id
}

#' Summarize voxel blocks
#' @keywords internal
summarize_blocks <- function(feature_mat, coords, block_id, block_size = 1) {
  unique_blocks <- sort(unique(block_id))
  nb <- length(unique_blocks)
  # Compute block means in O(N) using rowsum() rather than per-block `which(...)`.
  grp <- factor(block_id, levels = unique_blocks)
  counts <- as.numeric(tabulate(as.integer(grp), nbins = nb))
  counts[counts == 0] <- NA_real_

  block_summaries <- rowsum(feature_mat, grp, reorder = FALSE) / counts
  block_spatial <- rowsum(coords, grp, reorder = FALSE) / counts
  voxel_block_grid <- floor(
    sweep(coords, 2L, apply(coords, 2L, min), "-") / block_size
  )
  block_grid <- rowsum(voxel_block_grid, grp, reorder = FALSE) / counts
  dimnames(block_summaries) <- NULL
  dimnames(block_spatial) <- NULL
  dimnames(block_grid) <- NULL

  list(
    summaries = block_summaries,
    spatial   = block_spatial,
    grid = block_grid,
    unique_blocks = unique_blocks
  )
}

#' Build ACSC adjacency graph
#' @keywords internal
build_acsc_graph <- function(block_summary,
                             ann_k,
                             alpha,
                             spatial_weighting,
                             block_size,
                             correlation_metric = c("pearson", "spearman", "robust"),
                             knn_proj_dim = 0L,
                             knn_proj_method = c("none", "dct", "rp"),
                             knn_proj_seed = 1L) {
  correlation_metric <- match.arg(correlation_metric)
  spatial_weighting <- match.arg(spatial_weighting, c("gaussian", "binary"))
  block_summaries <- unname(as.matrix(block_summary$summaries))
  block_spatial <- unname(as.matrix(block_summary$spatial))
  nb <- nrow(block_summaries)
  if (nb < 1L || nrow(block_spatial) != nb || ncol(block_spatial) != 3L) {
    stop("build_acsc_graph: block summaries and N by 3 spatial centers must align")
  }
  if (any(!is.finite(block_summaries)) || any(!is.finite(block_spatial))) {
    stop("build_acsc_graph: block summaries and spatial centers must be finite")
  }
  ann_k <- min(
    nb - 1L,
    .cluster4d_scalar_number(
      ann_k, "ann_k", "build_acsc_graph", lower = 1, integer = TRUE
    )
  )
  alpha <- .cluster4d_scalar_number(
    alpha, "alpha", "build_acsc_graph", lower = 0, upper = 1
  )
  block_size <- .cluster4d_scalar_number(
    block_size, "block_size", "build_acsc_graph", lower = .Machine$double.eps
  )
  knn_proj_method <- match.arg(knn_proj_method)
  knn_proj_dim <- .cluster4d_scalar_number(
    knn_proj_dim, "knn_proj_dim", "build_acsc_graph", lower = 0, integer = TRUE
  )

  empty_graph <- function() {
    graph <- igraph::make_empty_graph(n = nb, directed = FALSE)
    graph <- igraph::set_graph_attr(graph, "correlation_metric", value = correlation_metric)
    igraph::set_graph_attr(graph, "alpha_feature_weight", value = alpha)
  }
  if (nb == 1L) return(empty_graph())

  # Six-neighbor block adjacency is built independently of feature values. The
  # alpha=0 endpoint uses only these edges, so feature permutations cannot alter
  # either topology or weights.
  block_grid <- block_summary$grid
  spatial_pairs <- matrix(integer(), ncol = 2L)
  if (is.matrix(block_grid) && identical(dim(block_grid), c(nb, 3L))) {
    block_grid <- round(block_grid)
    key <- apply(block_grid, 1L, paste, collapse = ":")
    pairs <- vector("list", 3L)
    for (axis in seq_len(3L)) {
      neighbor_grid <- block_grid
      neighbor_grid[, axis] <- neighbor_grid[, axis] + 1L
      target <- match(apply(neighbor_grid, 1L, paste, collapse = ":"), key)
      keep <- !is.na(target)
      pairs[[axis]] <- cbind(which(keep), target[keep])
    }
    spatial_pairs <- do.call(rbind, pairs)
  } else {
    spatial_distance <- as.matrix(stats::dist(block_spatial))
    spatial_pairs <- which(
      upper.tri(spatial_distance) &
        spatial_distance <= block_size * (1 + 1e-8),
      arr.ind = TRUE
    )
  }
  if (length(spatial_pairs)) {
    spatial_pairs <- cbind(
      pmin(spatial_pairs[, 1L], spatial_pairs[, 2L]),
      pmax(spatial_pairs[, 1L], spatial_pairs[, 2L])
    )
    spatial_pairs <- unique(spatial_pairs)
  } else {
    spatial_pairs <- matrix(integer(), ncol = 2L)
  }

  representation <- NULL
  feature_pairs <- matrix(integer(), ncol = 2L)
  if (alpha > 0) {
    representation <- acsc_metric_representation(
      block_summaries, correlation_metric
    )
    use_projection <- knn_proj_dim > 0L &&
      knn_proj_method != "none" && ncol(representation) > knn_proj_dim
    if (use_projection) {
      search_space <- representation
      if (identical(knn_proj_method, "dct")) {
        basis <- make_dct_basis(
          n_time = ncol(representation), n_basis = knn_proj_dim
        )
        search_space <- representation %*% basis
      } else if (identical(knn_proj_method, "rp")) {
        # Deterministic projection without perturbing the caller's RNG stream.
      seed_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(seed_state)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
          }
        } else {
          assign(".Random.seed", seed_state, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(as.integer(knn_proj_seed))
        projection <- matrix(
          stats::rnorm(ncol(representation) * knn_proj_dim),
          nrow = ncol(representation), ncol = knn_proj_dim
        )
        search_space <- (representation %*% projection) / sqrt(knn_proj_dim)
      }
      nearest <- FNN::get.knn(search_space, k = ann_k)$nn.index
      feature_pairs <- cbind(
        rep(seq_len(nb), each = ann_k),
        as.vector(t(nearest))
      )
    } else {
      similarity <- tcrossprod(representation)
      diag(similarity) <- -Inf
      nearest <- t(vapply(seq_len(nb), function(i) {
        order(-similarity[i, ], seq_len(nb))[seq_len(ann_k)]
      }, integer(ann_k)))
      feature_pairs <- cbind(
        rep(seq_len(nb), each = ann_k),
        as.vector(t(nearest))
      )
    }
    feature_pairs <- cbind(
      pmin(feature_pairs[, 1L], feature_pairs[, 2L]),
      pmax(feature_pairs[, 1L], feature_pairs[, 2L])
    )
    feature_pairs <- unique(feature_pairs[feature_pairs[, 1L] != feature_pairs[, 2L], , drop = FALSE])
  }

  pairs <- if (alpha == 0) spatial_pairs else unique(rbind(spatial_pairs, feature_pairs))
  if (!nrow(pairs)) return(empty_graph())
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  pair_key <- paste(pairs[, 1L], pairs[, 2L], sep = ":")
  spatial_key <- if (nrow(spatial_pairs)) {
    paste(spatial_pairs[, 1L], spatial_pairs[, 2L], sep = ":")
  } else character()
  feature_key <- if (nrow(feature_pairs)) {
    paste(feature_pairs[, 1L], feature_pairs[, 2L], sep = ":")
  } else character()
  is_spatial <- pair_key %in% spatial_key
  is_feature <- pair_key %in% feature_key

  delta <- block_spatial[pairs[, 1L], , drop = FALSE] -
    block_spatial[pairs[, 2L], , drop = FALSE]
  spatial_distance2 <- rowSums(delta^2)
  spatial_similarity <- numeric(nrow(pairs))
  if (spatial_weighting == "gaussian") {
    spatial_similarity[is_spatial] <- exp(
      -spatial_distance2[is_spatial] / (2 * block_size^2)
    )
  } else {
    spatial_similarity[is_spatial] <- 1
  }

  correlation <- numeric(nrow(pairs))
  correlation_similarity <- numeric(nrow(pairs))
  if (alpha > 0) {
    correlation <- rowSums(
      representation[pairs[, 1L], , drop = FALSE] *
        representation[pairs[, 2L], , drop = FALSE]
    )
    correlation <- pmax(-1, pmin(1, correlation))
    correlation_similarity <- (correlation + 1) / 2
  }
  final_weight <- alpha * correlation_similarity +
    (1 - alpha) * spatial_similarity
  keep <- is.finite(final_weight) & final_weight > 0
  if (!any(keep)) return(empty_graph())

  edges_df <- data.frame(
    source = pairs[keep, 1L],
    target = pairs[keep, 2L],
    weight = final_weight[keep],
    correlation = correlation[keep],
    correlation_similarity = correlation_similarity[keep],
    spatial_similarity = spatial_similarity[keep],
    edge_source = ifelse(
      is_spatial[keep] & is_feature[keep], "both",
      ifelse(is_spatial[keep], "spatial", "feature")
    ),
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(
    edges_df, directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(nb)))
  )
  graph <- igraph::set_graph_attr(graph, "correlation_metric", value = correlation_metric)
  igraph::set_graph_attr(graph, "alpha_feature_weight", value = alpha)
}

#' Run Louvain clustering
#' @keywords internal
run_louvain_clustering <- function(graph, resolution = NULL) {
  # Some igraph versions do not support a 'resolution' param in cluster_louvain.
  # If supported, it works. If not, user must rely on default or a different method.
  if (!is.null(resolution)) {
    # Try to pass resolution:
    igraph::cluster_louvain(graph, resolution = resolution)
  } else {
    igraph::cluster_louvain(graph)
  }
}

#' Estimate Louvain resolution parameter
#' @keywords internal
estimate_resolution <- function(K, graph) {
  stopifnot(is.numeric(K), length(K) == 1, K > 0)

  tolerance <- 0.05
  max_attempts <- 12

  # Heuristic assumption (typical for Louvain w/ resolution): higher resolution -> more clusters.
  # First, bracket a range that straddles K as best as we can.
  resolution_lower <- 1e-3
  resolution_upper <- 1.0

  count_at <- function(res) {
    clustering <- igraph::cluster_louvain(graph, resolution = res)
    length(igraph::communities(clustering))
  }

  n_upper <- count_at(resolution_upper)
  attempts <- 0
  while (n_upper < K && attempts < 12 && resolution_upper < 128) {
    resolution_lower <- resolution_upper
    resolution_upper <- resolution_upper * 2
    n_upper <- count_at(resolution_upper)
    attempts <- attempts + 1
  }

  n_lower <- count_at(resolution_lower)
  attempts <- 0
  while (n_lower > K && attempts < 12 && resolution_lower > 1e-6) {
    resolution_upper <- resolution_lower
    resolution_lower <- resolution_lower / 2
    n_lower <- count_at(resolution_lower)
    attempts <- attempts + 1
  }

  # If we still couldn't bracket, fall back to the best we have.
  if (!(n_lower <= K && n_upper >= K)) {
    warning("Could not find a resolution matching K. Returning approximate value.")
    return((resolution_lower + resolution_upper) / 2)
  }

  for (i in seq_len(max_attempts)) {
    mid_res <- (resolution_lower + resolution_upper) / 2
    num_clusters <- count_at(mid_res)

    if (abs(num_clusters - K) <= (tolerance * K)) {
      return(mid_res)
    }

    if (num_clusters < K) {
      # Need more clusters -> increase resolution
      resolution_lower <- mid_res
    } else {
      resolution_upper <- mid_res
    }
  }

  warning("Could not find a resolution matching K. Returning approximate value.")
  (resolution_lower + resolution_upper) / 2
}

#' Expand block-level cluster labels to voxel level
#' @keywords internal
expand_block_labels <- function(cluster_result, block_id, mask.idx) {
  # cluster_result: an igraph clustering object
  # block_id: an integer vector for each voxel
  # We need to map raw block_id -> consecutive index for membership
  block_labels <- if (is.numeric(cluster_result)) {
    as.integer(cluster_result)
  } else {
    igraph::membership(cluster_result)
  }
  nb <- length(block_labels)

  if (nb == 0) {
    stop("ACSC: Louvain returned zero communities; graph may be empty")
  }

  # Must recall how we assigned blocks in summarize_blocks
  # The easiest is to re-get unique_blocks from the cluster object if stored,
  # but we do it the same way:
  unique_blocks <- sort(unique(block_id))
  if (length(unique_blocks) != nb) {
    stop("Mismatch between the number of unique blocks and cluster membership vector.")
  }

  voxel_labels <- integer(length(mask.idx))
  for (i in seq_along(unique_blocks)) {
    b <- unique_blocks[i]
    voxel_labels[block_id == b] <- block_labels[i]
  }
  voxel_labels
}

#' Construct a 3D array of block labels
#' @keywords internal
construct_block_label_array <- function(block_id, mask) {
  block_array <- array(0L, dim(mask))
  mask.idx    <- which(mask > 0)
  block_array[mask.idx] <- match(block_id, sort(unique(block_id)))
  block_array
}

#' Normalize feature matrix to unit-length vectors
#'
#' @keywords internal
normalize_features <- function(feature_mat) {
  row_norms <- sqrt(rowSums(feature_mat^2))
  row_norms[row_norms < 1e-8] <- 1.0  # Avoid division by zero
  feature_mat / row_norms
}

#' Refine voxel boundaries using C++ or R implementation
#'
#' For each boundary voxel, compare correlation with each neighboring cluster's
#' cached centroid. This approach is much faster than comparing against all voxel time-series.
#'
#' @param use_cpp Logical; if TRUE, use C++ accelerated version (default)
#' @keywords internal
refine_voxel_boundaries <- function(voxel_labels, feature_mat, coords, max_iter,
                                    use_cpp = TRUE, dims = NULL,
                                    correlation_metric = c("pearson", "spearman", "robust")) {

  correlation_metric <- match.arg(correlation_metric)
  metric_features <- acsc_metric_representation(
    feature_mat, correlation_metric
  )

  # Try C++ implementation first if requested
  if (use_cpp && requireNamespace("neurocluster", quietly = TRUE)) {
    tryCatch({
      # Build 6-connected neighbors on the voxel grid (much cheaper than kNN)
      neighbor_indices <- if (!is.null(dims)) {
        .grid_neighbors6_from_coords(coords, dims)
      } else {
        FNN::get.knn(coords, k = 6)$nn.index
      }

      # Find boundary voxels (do this once outside the iteration)
      boundary_voxels <- find_boundary_voxels_cpp(voxel_labels, neighbor_indices)

      if (length(boundary_voxels) == 0) {
        return(voxel_labels)  # No boundaries to refine
      }

      # Call C++ refinement
      result <- refine_boundaries_cpp(
        voxel_labels = as.integer(voxel_labels),
        feature_mat_normalized = metric_features,
        neighbor_indices = neighbor_indices,
        boundary_voxels = as.integer(boundary_voxels),
        max_iter = as.integer(max_iter)
      )

      return(result$labels)
    }, error = function(e) {
      warning("C++ refinement failed: ", e$message, ". Falling back to R implementation.")
    })
  }

  # Fallback to R implementation
  refine_voxel_boundaries_r(
    voxel_labels, metric_features, coords, max_iter, dims = dims
  )
}

#' Refine voxel boundaries using R implementation (fallback)
#'
#' Pure R version for backward compatibility and fallback.
#'
#' @keywords internal
refine_voxel_boundaries_r <- function(voxel_labels, feature_mat, coords, max_iter, dims = NULL) {
  # Precompute cluster centroids in time
  cluster_centroids <- compute_cluster_centroids(voxel_labels, feature_mat)

  neighbor_indices <- if (!is.null(dims)) {
    .grid_neighbors6_from_coords(coords, dims)
  } else {
    FNN::get.knn(coords, k = 6)$nn.index
  }

  iter_count <- 0
  repeat {
    changed <- 0
    iter_count <- iter_count + 1

    # We only refine boundary voxels => those that have neighbors with different labels
    boundary_voxels <- find_boundary_voxels(voxel_labels, neighbor_indices)
    if (length(boundary_voxels) == 0) break

    for (i in boundary_voxels) {
      cur_label <- voxel_labels[i]
      nbrs      <- neighbor_indices[i, ]
      nbrs <- nbrs[nbrs > 0]
      nbr_labels <- voxel_labels[nbrs]

      # gather unique neighbor labels
      candidate_lbls <- unique(nbr_labels)
      if (length(candidate_lbls) == 1 && candidate_lbls == cur_label) {
        # all neighbors match current label -> skip
        next
      }

      best_corr <- cor_to_centroid(i, cur_label, feature_mat, cluster_centroids)
      best_lbl  <- cur_label

      # check each distinct label in neighbors
      for (lbl in candidate_lbls) {
        if (lbl != cur_label) {
          alt_corr <- cor_to_centroid(i, lbl, feature_mat, cluster_centroids)
          if (alt_corr > best_corr) {
            best_corr <- alt_corr
            best_lbl  <- lbl
          }
        }
      }
      if (best_lbl != cur_label) {
        voxel_labels[i] <- best_lbl
        changed <- changed + 1
      }
    }

    if (changed > 0) {
      # recompute centroids for updated assignments
      cluster_centroids <- compute_cluster_centroids(voxel_labels, feature_mat)
    }
    if (changed == 0 || iter_count >= max_iter) break
  }
  voxel_labels
}

#' Compute centroids of each cluster
#' @keywords internal
compute_cluster_centroids <- function(voxel_labels, feature_mat) {
  # cluster labels might not be 1..K; compute centroids for positive labels only.
  sel <- voxel_labels > 0 & !is.na(voxel_labels)
  if (!any(sel)) return(list())

  lbls <- as.integer(voxel_labels[sel])
  mat <- feature_mat[sel, , drop = FALSE]  # (#voxels x #time)
  unique_lbls <- sort(unique(lbls))

  grp <- factor(lbls, levels = unique_lbls)
  sums <- rowsum(mat, grp, reorder = FALSE)
  counts <- as.numeric(tabulate(as.integer(grp), nbins = length(unique_lbls)))
  counts[counts == 0] <- NA_real_
  centers <- sums / counts
  norms <- sqrt(rowSums(centers^2))
  live <- norms > sqrt(.Machine$double.eps)
  centers[live, ] <- centers[live, , drop = FALSE] / norms[live]
  centers[!live, ] <- 0
  dimnames(centers) <- NULL

  cluster_centroids <- lapply(seq_along(unique_lbls), function(i) centers[i, ])
  names(cluster_centroids) <- as.character(unique_lbls)
  cluster_centroids
}

#' Correlate a voxel's time-series with a cluster centroid
#' @keywords internal
cor_to_centroid <- function(voxel_idx, lbl, feature_mat, cluster_centroids) {
  lbl_str <- as.character(lbl)
  centroid_vec <- cluster_centroids[[lbl_str]]
  x <- feature_mat[voxel_idx, ]
  if (is.null(centroid_vec) || any(!is.finite(centroid_vec))) return(-Inf)
  correlation <- sum(x * centroid_vec)
  if (!is.finite(correlation)) return(-Inf)
  max(-1, min(1, correlation))
}

#' Identify boundary voxels
#'
#' boundary voxel = has at least one neighbor with a different label
#' @keywords internal
find_boundary_voxels <- function(voxel_labels, nn_index) {
  # nn_index: NxK matrix of neighbor indices
  # voxel_labels: length N
  N <- length(voxel_labels)
  boundary_mask <- logical(N)
  for (i in seq_len(N)) {
    nbrs <- nn_index[i, ]
    if (any(voxel_labels[nbrs] != voxel_labels[i])) {
      boundary_mask[i] <- TRUE
    }
  }
  which(boundary_mask)
}

# Build 6-connected neighbor indices from grid coordinates (no kNN needed).
# Returns an N x 6 integer matrix of 1-based neighbor row indices, with 0 where missing.
.grid_neighbors6_from_coords <- function(coords, dims) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 3) stop("coords must be N x 3")
  if (length(dims) != 3) stop("dims must be length 3")

  nx <- as.integer(dims[1])
  ny <- as.integer(dims[2])
  nz <- as.integer(dims[3])
  N <- nrow(coords)

  # Ensure integer grid coords
  i <- as.integer(coords[, 1])
  j <- as.integer(coords[, 2])
  k <- as.integer(coords[, 3])

  lin <- i + nx * ((j - 1L) + ny * (k - 1L))  # 1..nx*ny*nz
  lookup <- integer(nx * ny * nz)
  lookup[lin] <- seq_len(N)

  nbr <- matrix(0L, nrow = N, ncol = 6)
  stride_y <- nx
  stride_z <- nx * ny

  sel <- i > 1L
  nbr[sel, 1] <- lookup[lin[sel] - 1L]
  sel <- i < nx
  nbr[sel, 2] <- lookup[lin[sel] + 1L]
  sel <- j > 1L
  nbr[sel, 3] <- lookup[lin[sel] - stride_y]
  sel <- j < ny
  nbr[sel, 4] <- lookup[lin[sel] + stride_y]
  sel <- k > 1L
  nbr[sel, 5] <- lookup[lin[sel] - stride_z]
  sel <- k < nz
  nbr[sel, 6] <- lookup[lin[sel] + stride_z]

  nbr
}
