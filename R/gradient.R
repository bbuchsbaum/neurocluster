# ThreeByThreeOffset <- rbind(c(1,0,0),
#                             c(-1,0,0),
#                             c(0,1,0),
#                             c(0,-1,0),
#                             c(0,0,1),
#                             c(0,0,-1))
#
#
#
# TwoByTwoOffset <- rbind(c(1,0,0),
#                         c(-1,0,0),
#                         c(0,1,0),
#                         c(0,-1,0))
#






#' Spatial Gradient Calculation
#'
#' The spatial_gradient function calculates the spatial gradient of a \code{NeuroVol} instance within the specified mask.
#'
#' @param vol A \code{NeuroVol} instance for which the spatial gradient should be calculated.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the spatial gradient calculation.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\link[neuroim2:LogicalNeuroVol-class]{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param sigma A numeric value controlling the spatial weighting function. Default is 0.5.
#'
#' @return A \code{NeuroVol} instance containing the spatial gradient values for the input \code{vol}.
#'
#' @examples
#' \dontrun{
#' mask <- neuroim2::NeuroVol(array(1, c(20,20,20)), neuroim2::NeuroSpace(c(20,20,20)))
#' input_vol <- neuroim2::NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'   neuroim2::NeuroSpace(c(20,20,20)))
#' gradient_vol <- spatial_gradient(input_vol, mask)
#' }
#'
#' @seealso \code{\link[neighborweights:spatial_laplacian]{spatial_laplacian}}, \code{\link[neighborweights:weighted_spatial_adjacency]{weighted_spatial_adjacency}}
#' @importFrom neighborweights spatial_adjacency spatial_laplacian
#' @importFrom neuroim2 NeuroVol
#' @import assertthat
#'
#' @export
spatial_gradient <- function(vol, mask, sigma=.5) {
  mask.idx <- which(mask>0)

  G <- neighborweights::spatial_adjacency(index_to_coord(mask, mask.idx),
                                          dthresh=9, nnk=9,
                                          weight_mode="heat",
                                          sigma=.8, stochastic=TRUE)
  v <- vol[mask.idx]
  vs <- G %*% vol[mask.idx]

  S = neighborweights::spatial_laplacian(index_to_coord(mask, mask.idx),
                                         weight_mode="heat",
                                         nnk=27,
                                         dthresh=6,
                                         sigma=sigma,
                                         normalize=FALSE,
                                         stochastic=FALSE)

  grad <- S %*% vs
  NeuroVol(as.vector(grad), space(mask), indices=mask.idx)


}



#' @keywords internal
#' @noRd
# Function to check if any neighboring 4 voxels have a lower gradient value and update seeds
perturb_seeds <- function(grad, seeds) {
  result <- logical(nrow(seeds))
  new_seeds <- seeds

  for (i in seq_along(result)) {
    seed_i <- seeds[i, 1]
    seed_j <- seeds[i, 2]
    seed_k <- seeds[i, 3]
    seed_gradient <- grad[seed_i, seed_j, seed_k]

    neighbors <- matrix(c(seed_i - 1, seed_j, seed_k,
                          seed_i + 1, seed_j, seed_k,
                          seed_i, seed_j - 1, seed_k,
                          seed_i, seed_j + 1, seed_k), ncol = 3, byrow = TRUE)

    for (n in 1:nrow(neighbors)) {
      ni <- neighbors[n, 1]
      nj <- neighbors[n, 2]
      nk <- neighbors[n, 3]

      if (ni >= 1 && nj >= 1 && nk >= 1 &&
          ni <= dim(grad)[1] && nj <= dim(grad)[2] && nk <= dim(grad)[3]) {

        if (grad[ni, nj, nk] < seed_gradient) {
          result[i] <- TRUE
          new_seeds[i, ] <- c(ni, nj, nk)
          break
        }
      }
    }
  }

  return(list(result = result, new_seeds = new_seeds))
}


#' Find G3S Seeds via Functional Gradient Minima
#'
#' Identifies optimal seed locations for G3S clustering by finding local minima
#' in the functional gradient field. Unlike uniform grid seeding, this approach
#' places seeds in the centers of stable functional regions rather than on boundaries.
#'
#' @param feature_mat Numeric matrix (N x M) of compressed features, where each row
#'   is a voxel and each column is a feature dimension. Should be normalized to
#'   unit length for cosine similarity (as output by \code{\link{compress_features_svd}}).
#' @param coords Numeric matrix (N x 3) of spatial coordinates for each voxel,
#'   typically in mm units from \code{neuroim2::index_to_coord()}.
#' @param K Integer; target number of seeds (clusters). The function may return
#'   fewer seeds if spatial separation constraints cannot be satisfied.
#' @param k_neighbors Integer; number of nearest neighbors to use for gradient
#'   computation. Default: 26 (full 3D connectivity).
#' @param oversample_ratio Numeric; ratio of candidates to K considered before
#'   spatial separation is applied. Default: 3. Candidates are the
#'   lowest-gradient voxels *within each spatial cell* rather than the globally
#'   lowest, so a flat or heavily tied gradient field cannot concentrate the
#'   whole pool in one part of the mask.
#' @param min_separation_factor Numeric; starting minimum separation between
#'   seeds, as a multiple of the expected supervoxel scale. Default: 1. This is
#'   a starting point only: the radius decays until K seeds fit, so the outcome
#'   is the widest separation the mask actually admits and is insensitive to the
#'   exact starting value.
#' @param distance Character; distance metric for gradient computation.
#'   One of `"cosine"` (default) or `"euclidean"`.
#' @param knn Optional pre-computed k-nearest-neighbor object. If NULL (default),
#'   KNN is computed internally.
#' @param spatial_scale Optional positive physical length scale for the seed
#'   inhibition radius, normally the masked-volume scale used for compactness.
#'   When NULL the radius is estimated from the coordinate bounding box, which
#'   can overestimate it on non-convex masks.
#'
#' @return Integer vector of length K containing the row indices of selected
#'   seed voxels in the feature_mat/coords matrices.
#'
#' @details
#' ## Algorithm
#'
#' The seeding process follows these steps:
#'
#' 1. **Functional Gradient Computation**: For each voxel, compute the average
#'    dissimilarity to its k nearest spatial neighbors in feature space:
#'    \deqn{grad(i) = \frac{1}{k} \sum_{j \in N(i)} dist(f_i, f_j)}
#'    where \eqn{dist(f_i, f_j) = 1 - f_i \cdot f_j} (cosine distance).
#'
#' 2. **Candidate Selection**: Rank voxels by gradient value and select the top
#'    K * oversample_ratio voxels with lowest gradient (most stable regions).
#'
#' 3. **Spatial Separation**: Starting with the lowest-gradient candidate, greedily
#'    select seeds that maintain minimum spatial separation from all previously
#'    selected seeds. If that pool cannot provide K seeds at the physical
#'    inhibition radius, deterministic farthest-point sampling over the full
#'    mask fills the shortfall.
#'
#' ## Why Gradient-Based Seeding?
#'
#' Traditional uniform grid seeding often places seeds on functional boundaries
#' (high gradient regions), leading to:
#' - Poor initial cluster centroids
#' - Slower convergence
#' - Less biologically plausible parcellations
#'
#' Gradient-based seeding ensures:
#' - Seeds are in functional cores (low gradient = stable, homogeneous regions)
#' - Better initial centroids = faster, more accurate clustering
#' - Biologically plausible: aligns with cortical organization
#'
#' ## Comparison with Other Methods
#'
#' - **vs. Uniform Grid**: G3S gradient seeding reduces boundary artifacts and
#'   improves convergence speed by 2-3x.
#' - **vs. K-means Initialization**: G3S respects spatial structure, preventing
#'   disconnected clusters.
#' - **vs. Random Seeding**: G3S is deterministic given the data, improving
#'   reproducibility.
#'
#' @examples
#' \dontrun{
#' # Simulate compressed fMRI features
#' n_voxels <- 1000
#' n_components <- 15
#' features <- matrix(rnorm(n_voxels * n_components), n_voxels, n_components)
#' features <- t(apply(features, 1, function(x) x / sqrt(sum(x^2))))  # Normalize
#'
#' # Spatial coordinates (10x10x10 grid)
#' coords <- as.matrix(expand.grid(x = 1:10, y = 1:10, z = 1:10))
#'
#' # Find 20 seeds
#' seeds <- find_gradient_seeds_g3s(features, coords, K = 20)
#' print(length(seeds))  # Should be 20 or close to it
#'
#' # Visualize gradient values
#' grad_vals <- compute_functional_gradient(features, coords)
#' print(summary(grad_vals[seeds]))  # Should be low values
#' }
#'
#' @seealso
#' \code{\link{compute_functional_gradient}} for gradient computation details.
#' \code{\link{cluster4d_g3s}} for the full G3S clustering algorithm.
#' \code{\link{spatial_gradient}} for spatial-only gradient (used in SNIC).
#'
#' @export
#' @importFrom FNN get.knn
find_gradient_seeds_g3s <- function(feature_mat,
                                   coords,
                                   K,
                                   k_neighbors = 26,
                                   oversample_ratio = 3,
                                   min_separation_factor = 1,
                                   distance = c("cosine", "euclidean"),
                                   knn = NULL,
                                   spatial_scale = NULL) {

  distance <- match.arg(distance)

  if (!is.matrix(feature_mat)) feature_mat <- as.matrix(feature_mat)
  if (!is.matrix(coords)) coords <- as.matrix(coords)

  N <- nrow(feature_mat)
  if (nrow(coords) != N) {
    stop("feature_mat and coords must have the same number of rows")
  }

  if (K < 1 || K > N) {
    stop("K must be between 1 and ", N)
  }
  if (K == N) return(seq_len(N))

  if (k_neighbors >= N) {
    k_neighbors <- N - 1
  }

  # k-NN for gradient and spacing checks. Allow callers to pass a precomputed
  # structure so we can reuse one graph across seeding + propagation.
  if (is.null(knn)) {
    knn <- FNN::get.knn(coords, k = k_neighbors)
  } else {
    if (!is.list(knn) || is.null(knn$nn.index)) {
      stop("knn must be NULL or a list with nn.index (and optional nn.dist)")
    }
    if (!is.matrix(knn$nn.index) || nrow(knn$nn.index) != N) {
      stop("knn$nn.index must be a matrix with nrow(feature_mat) rows")
    }
    if (ncol(knn$nn.index) < 1) {
      stop("knn$nn.index must have at least one neighbor column")
    }
  }
  knn_idx <- knn$nn.index
  if (ncol(knn_idx) > k_neighbors) {
    knn_idx <- knn_idx[, seq_len(k_neighbors), drop = FALSE]
  }

  # Compute local gradient; use cosine for multi-feature (correlation-like)
  # and Euclidean when working with single-feature 3D data.
  grad_vals <- if (distance == "cosine") {
    calculate_local_gradient(t(feature_mat), knn_idx)
  } else {
    vapply(seq_len(N), function(i) {
      neighbor_indices <- knn_idx[i, ]
      neighbor_indices <- neighbor_indices[neighbor_indices > 0]
      if (length(neighbor_indices) == 0) return(0)

      neighbor_feats <- feature_mat[neighbor_indices, , drop = FALSE]
      diffs <- neighbor_feats -
        matrix(feature_mat[i, ], nrow = length(neighbor_indices),
               ncol = ncol(feature_mat), byrow = TRUE)
      mean(rowSums(diffs * diffs))
    }, numeric(1))
  }

  # Spatial inhibition radius ~ the expected supervoxel scale. The masked
  # measure is the right basis; the bounding box is only a fallback and can
  # exceed the true masked volume several-fold on a non-convex mask.
  if (is.null(spatial_scale)) {
    bbox <- apply(coords, 2, range)
    ranges <- bbox[2, ] - bbox[1, ]
    tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(coords)))
    active_ranges <- ranges[ranges > tolerance]
    effective_dimension <- max(1L, length(active_ranges))
    measure <- if (length(active_ranges)) prod(active_ranges) else 1
    avg_radius <- (measure / K)^(1 / effective_dimension)
  } else {
    avg_radius <- as.numeric(spatial_scale)[[1L]]
    if (!is.finite(avg_radius) || avg_radius <= 0) {
      stop("spatial_scale must be a positive finite number")
    }
  }

  # Candidate pool, stratified by spatial cell. Ranking every voxel by gradient
  # and keeping the globally lowest K * oversample_ratio concentrates the pool
  # wherever the gradient happens to be flattest; worse, when the field is tied
  # (low noise makes it exactly flat inside a parcel) order() resolves the tie
  # by voxel index, which collapses the pool into one corner of the volume and
  # leaves most parcels with no seed at all. Ranking within cells of the
  # expected supervoxel size keeps the pool the same size while guaranteeing it
  # spans the mask.
  candidate_count <- min(N, max(K, round(K * oversample_ratio)))
  cell <- floor((coords - rep(apply(coords, 2L, min), each = N)) / avg_radius)
  n_x <- max(cell[, 1L]) + 1
  n_y <- max(cell[, 2L]) + 1
  cell_id <- (cell[, 3L] * n_y + cell[, 2L]) * n_x + cell[, 1L]
  ranked <- order(cell_id, grad_vals, seq_len(N))
  rank_in_cell <- sequence(rle(cell_id[ranked])$lengths)
  per_cell <- max(1L, as.integer(ceiling(
    candidate_count / max(1L, length(unique(cell_id)))
  )))
  candidates <- ranked[rank_in_cell <= per_cell]
  candidates <- candidates[order(grad_vals[candidates], candidates)]

  # Greedy spatial inhibition in native code, decaying the radius until K seeds
  # fit. Because the pool now spans the mask, decaying is safe: it converges on
  # the widest separation the mask admits instead of letting one concentrated
  # low-gradient region consume every seed.
  selection <- g3s_select_seeds_cpp(
    coords = coords,
    candidates = as.integer(candidates),
    K = as.integer(K),
    radius = avg_radius * min_separation_factor,
    decay = 0.9,
    max_rounds = 60L
  )
  seeds <- selection$seeds

  if (length(seeds) < K) {
    # Only reachable when the stratified pool itself holds fewer than K
    # voxels. Fill the residual by deterministic farthest-point sampling so
    # the extra seeds are still spread across the mask.
    seeds <- .g3s_maximin_pad(seeds, coords, K)
  }

  sort(seeds)
}

#' Extend a seed set to K by repeated farthest-point selection.
#'
#' @keywords internal
#' @noRd
.g3s_maximin_pad <- function(seeds, coords, K) {
  n <- nrow(coords)
  if (!length(seeds)) seeds <- 1L
  nearest <- rep(Inf, n)
  for (seed in seeds) {
    delta <- coords - matrix(coords[seed, ], n, 3L, byrow = TRUE)
    nearest <- pmin(nearest, .rowSums(delta * delta, n, 3L))
  }
  nearest[seeds] <- -Inf
  while (length(seeds) < K) {
    pick <- which.max(nearest)
    seeds <- c(seeds, pick)
    delta <- coords - matrix(coords[pick, ], n, 3L, byrow = TRUE)
    nearest <- pmin(nearest, .rowSums(delta * delta, n, 3L))
    nearest[seeds] <- -Inf
  }
  seeds
}


#' Compute Functional Gradient from Feature Matrix
#'
#' Calculates the functional gradient for each voxel as the average dissimilarity
#' to its spatial neighbors in feature space. Low gradient values indicate stable,
#' homogeneous regions (good seed locations); high values indicate boundaries.
#'
#' @param feature_mat Numeric matrix (N x M) of features (voxels x dimensions).
#' @param coords Numeric matrix (N x 3) of spatial coordinates.
#' @param k_neighbors Integer; number of nearest spatial neighbors to use.
#'   Default: 26.
#'
#' @return Numeric vector of length N containing gradient values for each voxel.
#'
#' @details
#' For each voxel i, the functional gradient is:
#' \deqn{grad(i) = \frac{1}{k} \sum_{j \in N_k(i)} (1 - f_i \cdot f_j)}
#'
#' where \eqn{N_k(i)} are the k nearest spatial neighbors and \eqn{f_i \cdot f_j}
#' is the dot product (cosine similarity for unit-length vectors).
#'
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(500 * 15), 500, 15)
#' features <- t(apply(features, 1, function(x) x / sqrt(sum(x^2))))
#' coords <- matrix(rnorm(500 * 3), 500, 3)
#'
#' grad <- compute_functional_gradient(features, coords, k_neighbors = 26)
#' hist(grad, breaks = 50, main = "Functional Gradient Distribution")
#' }
#'
#' @export
#' @importFrom FNN get.knn
compute_functional_gradient <- function(feature_mat, coords, k_neighbors = 26) {
  N <- nrow(feature_mat)

  # Find k nearest spatial neighbors
  if (k_neighbors >= N) {
    k_neighbors <- N - 1
  }

  neib <- FNN::get.knn(coords, k = k_neighbors)

  # Compute gradient for each voxel
  grad_vals <- vapply(seq_len(N), function(i) {
    # Get feature vectors of neighbors
    neighbor_indices <- neib$nn.index[i, ]
    neighbor_feats <- feature_mat[neighbor_indices, , drop = FALSE]

    # Compute cosine distances (1 - correlation)
    # For normalized vectors, dot product = correlation
    voxel_feat <- feature_mat[i, , drop = FALSE]

    # Vectorized dot product with all neighbors
    dots <- neighbor_feats %*% t(voxel_feat)  # k x 1 matrix
    dists <- 1 - as.vector(dots)

    # Average distance = gradient
    mean(dists)
  }, numeric(1))

  grad_vals
}


#' Find Gradient Seeds (Wrapper for Compatibility)
#'
#' Legacy wrapper function that maintains compatibility with the old
#' \code{find_initial_points()} interface while using the new G3S gradient
#' seeding implementation.
#'
#' @param coords Numeric matrix (N x 3) of grid coordinates (1-indexed).
#' @param grad_vals Numeric vector of length N containing gradient values.
#' @param K Integer; number of seeds to find.
#' @param min_separation_factor Minimum spatial separation as multiple of
#'   expected grid spacing. Default: 1.5.
#'
#' @return Integer vector of seed indices.
#'
#' @keywords internal
#' @noRd
find_gradient_seeds <- function(coords, grad_vals, K, min_separation_factor = 1.5) {
  N <- nrow(coords)

  # Select candidates based on gradient
  n_candidates <- min(K * 3, N)
  candidate_indices <- order(grad_vals)[1:n_candidates]

  # Compute expected spacing
  volume <- prod(apply(coords, 2, function(x) diff(range(x))))
  expected_spacing <- (volume / K)^(1/3)
  min_dist <- expected_spacing * min_separation_factor

  # Greedy selection with spatial separation
  selected_seeds <- integer(K)
  selected_seeds[1] <- candidate_indices[1]
  n_selected <- 1

  for (cand in candidate_indices[-1]) {
    if (n_selected >= K) break

    dists <- sqrt(rowSums((coords[selected_seeds[1:n_selected], , drop = FALSE] -
                           matrix(coords[cand, ], nrow = n_selected, ncol = 3, byrow = TRUE))^2))

    if (min(dists) > min_dist) {
      n_selected <- n_selected + 1
      selected_seeds[n_selected] <- cand
    }
  }

  if (n_selected < K) {
    warning("Could not find ", K, " spatially separated seeds. Found ", n_selected)
    selected_seeds <- selected_seeds[1:n_selected]
  }

  selected_seeds
}
