#' Randomize Cluster IDs
#'
#' Randomly permute cluster IDs to break any spatial ordering. This is useful
#' when visualizing clusters with continuous color scales, as spatially-ordered
#' cluster IDs create artificial gradients that can be misleading.
#'
#' @param cluster_result A clustering result object (e.g., from cluster4d,
#'   supervoxels_flash3d, slice_msf, etc.) containing at minimum:
#'   - `cluster`: vector of cluster assignments
#'   - `clusvol`: ClusteredNeuroVol object
#'   - Optionally: `centers`, `coord_centers`, etc.
#' @param seed Optional random seed for reproducibility. If NULL, uses R's
#'   current random state. Default is NULL.
#'
#' @return A modified cluster_result object with randomized cluster IDs. All
#'   components (cluster vector, clusvol, centers, coord_centers, etc.) are
#'   updated consistently.
#'
#' @details
#' Many clustering algorithms assign cluster IDs sequentially as clusters are
#' encountered during processing. When voxels are processed in memory order
#' (which is typically z-major in 3D arrays), this creates systematic ordering
#' where cluster ID correlates with z-position (inferior to superior in brain
#' imaging).
#'
#' This spatial ordering creates two problems:
#'
#' 1. **Visualization artifacts**: When viewed with continuous color scales
#'    (e.g., viridis, rainbow), ordered cluster IDs create artificial gradients
#'    that suggest spatial structure where none exists. Randomized IDs show
#'    each cluster as a distinct color patch.
#'
#' 2. **Interpretation bias**: Ordered IDs may mislead users into thinking
#'    cluster numbering has anatomical significance (e.g., "cluster 1 is always
#'    inferior, cluster 100 is always superior").
#'
#' This function breaks the ordering by randomly permuting cluster IDs while
#' maintaining all cluster memberships and properties.
#'
#' @examples
#' \dontrun{
#' # Standard workflow with ordered IDs
#' result <- cluster4d(vec, mask, K = 100, method = "flash3d")
#'
#' # Randomize for better visualization
#' result_random <- randomize_cluster_ids(result, seed = 42)
#'
#' # Compare visualizations
#' # Ordered IDs show gradient from inferior to superior
#' plot(result$clusvol)
#'
#' # Randomized IDs show distinct color patches
#' plot(result_random$clusvol)
#'
#' # Works with any clustering method
#' result_msf <- cluster4d(vec, mask, K = 100, method = "slice_msf")
#' result_msf_random <- randomize_cluster_ids(result_msf)
#' }
#'
#' @export
randomize_cluster_ids <- function(cluster_result, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }

  # Extract cluster vector
  if (!is.list(cluster_result) || !("cluster" %in% names(cluster_result))) {
    stop("cluster_result must contain a 'cluster' component")
  }

  cluster_ids <- cluster_result$cluster
  unique_ids <- sort(unique(cluster_ids[cluster_ids > 0]))
  n_clusters <- length(unique_ids)

  if (n_clusters <= 1) {
    message("Only one cluster found; no randomization needed")
    return(cluster_result)
  }

  # Create random permutation
  shuffle_map <- sample(seq_len(n_clusters))
  names(shuffle_map) <- unique_ids

  # Update cluster vector
  new_cluster_ids <- cluster_ids
  for (i in seq_along(cluster_ids)) {
    if (cluster_ids[i] > 0) {
      old_id <- as.character(cluster_ids[i])
      new_cluster_ids[i] <- shuffle_map[old_id]
    }
  }

  # Update clusvol if present
  if ("clusvol" %in% names(cluster_result)) {
    # Extract mask and recreate ClusteredNeuroVol with new IDs
    clusvol <- cluster_result$clusvol
    if (inherits(clusvol, "ClusteredNeuroVol")) {
      mask <- neuroim2::space(clusvol)
      cluster_result$clusvol <- suppressWarnings(
        neuroim2::ClusteredNeuroVol(mask, clusters = new_cluster_ids)
      )
    }
  }

  # Update centers if present
  if ("centers" %in% names(cluster_result) && !is.null(cluster_result$centers)) {
    if (is.matrix(cluster_result$centers)) {
      # Centers are typically stored as columns or rows
      # Determine orientation by checking dimensions
      centers <- cluster_result$centers
      if (ncol(centers) == n_clusters) {
        # Columns are clusters
        cluster_result$centers <- centers[, shuffle_map, drop = FALSE]
      } else if (nrow(centers) == n_clusters) {
        # Rows are clusters
        cluster_result$centers <- centers[shuffle_map, , drop = FALSE]
      }
    }
  }

  # Update coord_centers if present
  if ("coord_centers" %in% names(cluster_result) && !is.null(cluster_result$coord_centers)) {
    if (is.matrix(cluster_result$coord_centers) && nrow(cluster_result$coord_centers) == n_clusters) {
      cluster_result$coord_centers <- cluster_result$coord_centers[shuffle_map, , drop = FALSE]
    }
  }

  # Update cluster vector
  cluster_result$cluster <- new_cluster_ids

  # Add metadata about randomization
  if (!"parameters" %in% names(cluster_result)) {
    cluster_result$parameters <- list()
  }
  cluster_result$parameters$cluster_ids_randomized <- TRUE
  cluster_result$parameters$randomization_seed <- seed

  cluster_result
}


#' Check for Spatial Ordering in Cluster IDs
#'
#' Diagnostic function to detect whether cluster IDs are spatially ordered
#' along any axis. Useful for identifying clustering algorithms that may
#' produce ordered IDs.
#'
#' @param cluster_result A clustering result object with `coord_centers` component
#'
#' @return A list with correlation statistics:
#'   - `cor_x`: Spearman correlation between cluster ID and x-coordinate
#'   - `cor_y`: Spearman correlation between cluster ID and y-coordinate
#'   - `cor_z`: Spearman correlation between cluster ID and z-coordinate
#'   - `max_abs_cor`: Maximum absolute correlation across all axes
#'   - `is_ordered`: Logical indicating if any axis shows strong ordering (|r| > 0.3)
#'
#' @examples
#' \dontrun{
#' result <- cluster4d(vec, mask, K = 100, method = "flash3d")
#'
#' # Check for ordering
#' ordering <- check_cluster_ordering(result)
#' print(ordering)
#'
#' # If strongly ordered, randomize
#' if (ordering$is_ordered) {
#'   result <- randomize_cluster_ids(result)
#' }
#' }
#'
#' @export
check_cluster_ordering <- function(cluster_result) {

  if (!("coord_centers" %in% names(cluster_result))) {
    stop("cluster_result must contain 'coord_centers' component")
  }

  coords <- cluster_result$coord_centers
  if (!is.matrix(coords) || ncol(coords) < 3) {
    stop("coord_centers must be a matrix with at least 3 columns")
  }

  n_clusters <- nrow(coords)
  cluster_ids <- seq_len(n_clusters)

  # Compute Spearman correlations
  cor_x <- cor(cluster_ids, coords[, 1], method = "spearman")
  cor_y <- cor(cluster_ids, coords[, 2], method = "spearman")
  cor_z <- cor(cluster_ids, coords[, 3], method = "spearman")

  max_abs_cor <- max(abs(c(cor_x, cor_y, cor_z)))

  list(
    cor_x = cor_x,
    cor_y = cor_y,
    cor_z = cor_z,
    max_abs_cor = max_abs_cor,
    is_ordered = max_abs_cor > 0.3,
    interpretation = if (max_abs_cor > 0.7) {
      "Strong spatial ordering detected - randomization strongly recommended"
    } else if (max_abs_cor > 0.3) {
      "Moderate spatial ordering detected - randomization recommended"
    } else {
      "No significant spatial ordering detected"
    }
  )
}
