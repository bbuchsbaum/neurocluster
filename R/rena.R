#' Recursive Nearest Agglomeration (ReNA) Clustering
#'
#' Performs spatially-constrained hierarchical clustering using recursive 1-nearest
#' neighbor aggregation. ReNA is a fast, linear-time agglomerative algorithm that
#' avoids the "percolation" problem of standard hierarchical methods by using
#' 1-NN graphs to ensure balanced cluster sizes.
#'
#' @note Consider using \code{\link{cluster4d}} with \code{method = "rena"} for a
#' standardized interface across all clustering methods.
#'
#' @param bvec A \code{NeuroVec} instance supplying the 4D data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in clustering.
#'   If numeric, nonzero values define included voxels. If logical, TRUE values
#'   define included voxels.
#' @param K The number of clusters to find (default 100).
#' @param connectivity Neighborhood connectivity (6 or 26). Default 26.
#'   6 = face neighbors only, 26 = face + edge + corner neighbors.
#' @param max_iterations Maximum number of recursion iterations (default 50).
#'   Algorithm stops when K clusters are reached or max_iterations is hit.
#' @param verbose Logical; whether to print progress messages. Default FALSE.
#' @param exact_k Logical; whether to use edge pruning to ensure exactly K clusters.
#'   Default TRUE. If FALSE, may produce slightly more or fewer than K clusters.
#'
#' @return A \code{list} of class \code{rena_cluster_result} (inheriting from
#'   \code{cluster_result}) with the following elements:
#' \describe{
#'   \item{clusvol}{An instance of type \linkS4class{ClusteredNeuroVol}.}
#'   \item{cluster}{A vector of cluster indices equal to the number of voxels in the mask.}
#'   \item{centers}{A matrix of cluster centers with each row representing the feature vector for a cluster.}
#'   \item{coord_centers}{A matrix of spatial coordinates with each row corresponding to a cluster.}
#'   \item{n_clusters}{Actual number of clusters found.}
#'   \item{method}{Clustering method used ("rena").}
#'   \item{parameters}{List of all parameters used.}
#'   \item{metadata}{List containing iteration count and convergence information.}
#' }
#'
#' @details
#' ## Algorithm Description
#'
#' ReNA recursively aggregates features using a 1-Nearest Neighbor (1-NN) graph
#' to ensure clusters remain roughly balanced in size and spatially compact.
#' The algorithm proceeds in iterations:
#'
#' 1. **Connectivity Constraint**: Clustering respects the spatial topology graph G.
#'    Distance is only calculated between spatially adjacent voxels.
#'
#' 2. **Similarity Calculation**: For all connected edges, compute squared Euclidean
#'    distance in feature space.
#'
#' 3. **1-NN Subgraph**: For every voxel, identify its single nearest neighbor.
#'    This forms a directed subgraph Q. The 1-NN graph is less likely to "percolate"
#'    (create giant components) compared to k-NN graphs where k >= 2.
#'
#' 4. **Connected Components**: Extract weakly connected components of Q. These
#'    components become the clusters for the current iteration.
#'
#' 5. **Reduction**: Features within each component are averaged to form a single
#'    "super-feature" for the next iteration. The graph G is contracted so edges
#'    between new clusters exist if any constituent voxels were connected.
#'
#' 6. **Stopping Condition**: Repeat recursively until the number of clusters
#'    drops to the desired K. If exact_k = TRUE, edges are pruned by distance
#'    to ensure exactly K clusters.
#'
#' ## Performance Characteristics
#'
#' - **Complexity**: O(N log N) per iteration where N = number of voxels
#' - **Memory**: O(N) for graph and assignments
#' - **Speed**: Generally faster than iterative methods (supervoxels)
#' - **Iterations**: Typically 5-15 iterations to reach K clusters from N voxels
#'
#' ## Advantages over Other Methods
#'
#' - **No percolation**: Unlike single-linkage, doesn't form one giant cluster
#' - **Balanced sizes**: 1-NN graph encourages roughly equal cluster sizes
#' - **Topology-aware**: Respects spatial structure via connectivity graph
#' - **Deterministic**: No random initialization
#' - **Fast convergence**: Linear-time operations per iteration
#'
#' ## Parallelization Status
#'
#' **Currently NOT parallelized.** ReNA uses a sequential recursive algorithm
#' that processes the graph in a hierarchical manner. Each iteration depends on
#' the previous iteration's clustering result.
#'
#' ### Sequential Operations:
#'
#' 1. **Distance Computation**: Pairwise distances for connected voxels (C++)
#' 2. **1-NN Finding**: Each voxel finds its nearest neighbor (C++)
#' 3. **Component Detection**: Union-Find algorithm for connected components (C++)
#' 4. **Feature Aggregation**: Mean pooling within components (C++)
#' 5. **Graph Contraction**: Building reduced graph for next iteration (C++)
#'
#' ### Why Not Parallelized:
#'
#' - **Sequential dependency**: Each iteration requires previous iteration's result
#' - **Fast per-iteration**: Linear-time operations make parallelization overhead high
#' - **Small iteration count**: Typically 5-15 iterations total
#' - **C++ optimization**: Core operations already optimized in C++
#'
#' ### Performance Tips:
#'
#' - **Reduce connectivity**: Use connectivity=6 instead of 26 for larger graphs
#' - **Use smaller K**: Fewer target clusters means fewer iterations
#' - **Pre-smooth data**: Reduces noise, improves cluster coherence
#' - **Alternative for parallelism**: Consider \code{slice_msf()} or \code{flash3d()}
#'
#' @examples
#' \dontrun{
#' # Small synthetic example
#' library(neuroim2)
#' mask <- NeuroVol(array(1, c(10,10,10)), NeuroSpace(c(10,10,10)))
#' vec <- replicate(20,
#'                  NeuroVol(array(runif(10*10*10), c(10,10,10)),
#'                           NeuroSpace(c(10,10,10))),
#'                  simplify=FALSE)
#' vec <- do.call(concat, vec)
#'
#' # Run ReNA clustering
#' rena_res <- rena(vec, mask, K=50, connectivity=26)
#' print(rena_res$n_clusters)
#'
#' # With verbose output
#' rena_res <- rena(vec, mask, K=50, verbose=TRUE)
#' }
#'
#' @references
#' Hoyos-Idrobo, A., Varoquaux, G., Kahn, J., & Thirion, B. (2019).
#' Recursive Nearest Agglomeration (ReNA): Fast clustering for approximation
#' of structured signals. \emph{Pattern Recognition}, 94, 17-28.
#'
#' @seealso
#' \code{\link{cluster4d}} for unified interface,
#' \code{\link{supervoxels}} for iterative heat kernel clustering,
#' \code{\link{snic}} for non-iterative priority queue clustering
#'
#' @importFrom neuroim2 NeuroVec NeuroVol ClusteredNeuroVol series index_to_coord index_to_grid spacing
#' @importFrom Matrix sparseMatrix
#' @export
rena <- function(bvec, mask,
                 K = 100,
                 connectivity = 26,
                 max_iterations = 50,
                 verbose = FALSE,
                 exact_k = TRUE) {

  # Use common validation
  validate_cluster4d_inputs(bvec, mask, K, "rena")

  # Get mask indices
  mask.idx <- which(mask > 0)
  n_voxels <- length(mask.idx)

  if (verbose) {
    message("ReNA: Starting with ", n_voxels, " voxels, target K = ", K)
  }

  # Get coordinates
  coords <- index_to_grid(mask, mask.idx)

  # Extract and scale features (voxels x time)
  feature_mat <- neuroim2::series(bvec, mask.idx)
  feature_mat <- t(as.matrix(feature_mat))  # Ensure matrix, transpose to voxels x time
  feature_mat <- base::scale(feature_mat, center = TRUE, scale = TRUE)
  feature_mat[is.na(feature_mat)] <- 0

  # Transpose for C++ (features x voxels)
  feature_mat <- t(feature_mat)

  # Initialize connectivity graph using neighborweights
  # This gives us a sparse adjacency matrix based on spatial proximity
  adj_matrix <- rena_build_connectivity(coords, mask, connectivity)

  # Track voxel-to-component mapping
  voxel_to_component <- seq_len(n_voxels) - 1  # 0-based for C++

  # Iterative reduction
  current_features <- feature_mat
  current_coords <- coords
  current_adjacency <- adj_matrix
  current_n_clusters <- n_voxels
  iteration <- 0

  while (current_n_clusters > K && iteration < max_iterations) {
    iteration <- iteration + 1

    if (verbose) {
      message("ReNA: Iteration ", iteration, " - ", current_n_clusters, " clusters")
    }

    # 1. Compute distances for connected pairs
    adj_indices <- Matrix::which(current_adjacency > 0, arr.ind = TRUE)
    adj_i <- adj_indices[, 1]
    adj_j <- adj_indices[, 2]

    # Only keep upper triangle (undirected graph)
    upper_tri <- adj_i < adj_j
    adj_i <- adj_i[upper_tri]
    adj_j <- adj_j[upper_tri]

    if (length(adj_i) == 0) {
      warning("ReNA: No edges in graph, stopping early")
      break
    }

    distances <- compute_masked_distances_cpp(current_features, adj_i, adj_j)

    # 2. Find 1-NN subgraph
    nearest_neighbor <- find_1nn_subgraph_cpp(current_n_clusters, adj_i, adj_j, distances)

    # 3. Pruning for exact K (if clustering would otherwise drop below K)
    if (exact_k && current_n_clusters > K) {
      # Estimate how many components we'd get
      test_components <- find_connected_components_cpp(current_n_clusters, nearest_neighbor)
      test_n_components <- length(unique(test_components))

      if (test_n_components < K) {
        # Need to prune edges to maintain K components
        if (verbose) {
          message("ReNA: Pruning edges to achieve exactly K = ", K, " clusters")
        }

        # Create distance vector for nearest neighbors
        nn_distances <- rep(Inf, current_n_clusters)
        for (e in seq_along(adj_i)) {
          i <- adj_i[e]
          j <- adj_j[e]
          if (nearest_neighbor[i] + 1 == j) {  # Convert 0-based to 1-based
            nn_distances[i] <- distances[e]
          }
          if (nearest_neighbor[j] + 1 == i) {
            nn_distances[j] <- distances[e]
          }
        }

        nearest_neighbor <- prune_edges_for_k_cpp(
          current_n_clusters,
          nearest_neighbor,
          nn_distances,
          K
        )
      }
    }

    # 4. Find connected components
    component_labels <- find_connected_components_cpp(current_n_clusters, nearest_neighbor)
    n_components <- length(unique(component_labels))

    if (n_components == current_n_clusters) {
      warning("ReNA: No progress in iteration ", iteration, ", stopping")
      break
    }

    # Update voxel-to-component mapping
    new_voxel_to_component <- component_labels[voxel_to_component + 1]  # R is 1-based
    voxel_to_component <- new_voxel_to_component

    # 5. Aggregate features and coordinates
    current_features <- aggregate_features_cpp(current_features, component_labels, n_components)
    current_coords <- aggregate_coords_cpp(current_coords, component_labels, n_components)

    # 6. Contract graph
    contracted <- contract_graph_cpp(adj_i, adj_j, component_labels, n_components)

    if (length(contracted$i) == 0) {
      # No edges left, each component is isolated
      if (verbose) {
        message("ReNA: Graph fully contracted, stopping")
      }
      break
    }

    # Build new sparse adjacency matrix
    current_adjacency <- Matrix::sparseMatrix(
      i = contracted$i,
      j = contracted$j,
      x = rep(1, length(contracted$i)),
      dims = c(n_components, n_components),
      symmetric = TRUE
    )

    current_n_clusters <- n_components

    # Check if we've reached target
    if (current_n_clusters <= K) {
      if (verbose) {
        message("ReNA: Reached target K = ", K, " in ", iteration, " iterations")
      }
      break
    }
  }

  # Create final cluster assignments (1-based for R)
  final_labels <- voxel_to_component + 1

  # Prepare data for standardized result
  data_prep <- list(
    features = t(feature_mat),  # Transpose back to voxels x time
    coords = coords,
    mask_idx = mask.idx,
    n_voxels = n_voxels,
    dims = dim(mask),
    spacing = spacing(mask)
  )

  # Create standardized result
  result <- create_cluster4d_result(
    labels = final_labels,
    mask = mask,
    data_prep = data_prep,
    method = "rena",
    parameters = list(
      K = K,
      requested_K = K,
      connectivity = connectivity,
      max_iterations = max_iterations,
      exact_k = exact_k
    ),
    metadata = list(
      iterations = iteration,
      final_n_clusters = current_n_clusters
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  # Add ReNA-specific class
  class(result) <- c("rena_cluster_result", class(result))

  result
}

#' Build connectivity graph for ReNA
#'
#' Creates sparse adjacency matrix based on spatial connectivity.
#'
#' @keywords internal
rena_build_connectivity <- function(coords, mask, connectivity) {

  # Use FNN to find k-nearest neighbors based on connectivity
  k <- min(connectivity, nrow(coords) - 1)

  if (k < 1) {
    # Trivial case: single voxel
    return(Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 1)))
  }

  nn_result <- FNN::get.knn(coords, k = k)

  # Build edge list (undirected). Collect all neighbor pairs, then deduplicate.
  all_i <- rep(seq_len(nrow(coords)), each = k)
  all_j <- as.integer(as.vector(t(nn_result$nn.index)))

  # Drop self-edges if any slipped through and canonicalize ordering
  keep <- all_i != all_j
  if (!any(keep)) {
    n <- nrow(coords)
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n, n)))
  }

  edge_matrix <- cbind(
    pmin(all_i[keep], all_j[keep]),
    pmax(all_i[keep], all_j[keep])
  )

  edge_matrix <- unique(edge_matrix)

  if (nrow(edge_matrix) == 0) {
    n <- nrow(coords)
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n, n)))
  }

  # Create sparse adjacency matrix
  adj <- Matrix::sparseMatrix(
    i = c(edge_matrix[, 1], edge_matrix[, 2]),
    j = c(edge_matrix[, 2], edge_matrix[, 1]),
    x = 1,
    dims = c(nrow(coords), nrow(coords))
  )

  adj
}

#' Cluster4d using ReNA method
#'
#' Wrapper for ReNA algorithm with standardized cluster4d interface.
#'
#' @inheritParams cluster4d
#' @param ... Additional parameters (currently unused for ReNA)
#' @param exact_k Logical; if TRUE, apply edge pruning to target exactly \code{n_clusters}.
#'
#' @return A cluster4d_result object
#' @export
cluster4d_rena <- function(vec, mask, n_clusters = 100,
                           spatial_weight = 0.5,
                           max_iterations = 50,
                           connectivity = 26,
                           verbose = FALSE,
                           exact_k = TRUE,
                           ...) {

  # ReNA doesn't use spatial_weight in the same way as other methods,
  # but we accept it for interface consistency
  # Could potentially use it to weight spatial vs feature distances in future

  # Call original rena
  result <- rena(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    connectivity = connectivity,
    max_iterations = max_iterations,
    verbose = verbose,
    exact_k = exact_k
  )

  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }

  # Ensure method info
  result$method <- "rena"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))

  # Store all parameters
  result$parameters <- modifyList(
    result$parameters,
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      max_iterations = max_iterations,
      connectivity = connectivity,
      exact_k = exact_k
    )
  )

  result
}
