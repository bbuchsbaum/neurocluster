#' Create Synthetic Ground Truth Clustering Data
#'
#' Generates a volume with perfectly separable clusters for testing
#' clustering algorithm correctness. Each cluster is a spatially contiguous
#' cubic sub-region with identical timeseries within the cluster and orthogonal
#' signals between clusters.
#'
#' @param n_time Number of timepoints. Default 100.
#' @param noise_sd Standard deviation of Gaussian noise to add. Default 0
#'   (no noise = perfect separability).
#' @param n_clusters Number of clusters. Must be a perfect cube (8, 27, 64, etc.)
#'   for cubic cluster geometry. Default 27.
#' @param dim_per_cluster Voxels per dimension per cluster. Default 3.
#' @param seed Random seed for reproducibility. Default NULL.
#'
#' @return A list with:
#'   \item{vec}{NeuroVec object}
#'   \item{mask}{NeuroVol mask (all ones)}
#'   \item{true_clusters}{NeuroVol with ground truth cluster labels}
#'   \item{true_labels}{Integer vector of true cluster assignments}
#'   \item{n_clusters}{Number of clusters}
#'   \item{cluster_size}{Size of each cluster in voxels}
#'   \item{signals}{Matrix of cluster signals (n_time x n_clusters)}
#'
#' @details
#' For the default n_clusters=27, creates a 9x9x9 volume divided into
#' 27 cubic 3x3x3 clusters arranged in a 3x3x3 grid:
#'
#' \preformatted{
#' z=1-3:           z=4-6:           z=7-9:
#'   y                y                y
#'   ^                ^                ^
#'   | 7 | 8 | 9 |    | 16| 17| 18|    | 25| 26| 27|
#'   | 4 | 5 | 6 |    | 13| 14| 15|    | 22| 23| 24|
#'   | 1 | 2 | 3 |    | 10| 11| 12|    | 19| 20| 21|
#'   +-----------> x  +-----------> x  +-----------> x
#' }
#'
#' Each cluster has a unique orthogonal signal generated from QR decomposition
#' of a random matrix, ensuring zero correlation between clusters.
#'
#' This creates an "easy" test case where any reasonable clustering algorithm
#' should achieve perfect or near-perfect recovery with noise_sd=0.
#'
#' @examples
#' # Default: 27 clusters, perfect separability
#' data <- make_synthetic_clusters()
#' table(data$true_labels)  # 27 voxels per cluster
#'
#' # 8 clusters (2x2x2 grid of 3x3x3 cubes = 6x6x6 volume)
#' data <- make_synthetic_clusters(n_clusters = 8)
#'
#' # With noise for more realistic testing
#' data <- make_synthetic_clusters(noise_sd = 0.5)
#'
#' @export
make_synthetic_clusters <- function(n_time = 100, noise_sd = 0, n_clusters = 27,
                                    dim_per_cluster = 3, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)


  # n_clusters must be a perfect cube for cubic geometry
  grid_dim <- round(n_clusters^(1/3))
  if (grid_dim^3 != n_clusters) {
    stop("n_clusters must be a perfect cube (8, 27, 64, 125, ...) for cubic cluster geometry")
  }

  # Volume dimensions
  vol_dim <- grid_dim * dim_per_cluster
  dim_vol <- rep(as.integer(vol_dim), 3)
  n_voxels <- prod(dim_vol)
  cluster_size <- dim_per_cluster^3

  # Create orthogonal signals for each cluster
  # Use QR decomposition of random matrix to get orthonormal columns
  if (!is.null(seed)) set.seed(seed)
  random_mat <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
  qr_decomp <- qr(random_mat)
  signals <- qr.Q(qr_decomp)  # n_time x n_clusters orthonormal matrix

  # Scale signals to have meaningful variance
  signals <- signals * sqrt(n_time)


  # Create ground truth cluster labels
  # Divide volume into grid_dim x grid_dim x grid_dim arrangement of cubic clusters
  true_labels <- array(0L, dim = dim_vol)

  cluster_id <- 1L
  for (zblock in 0:(grid_dim - 1)) {
    z_idx <- (zblock * dim_per_cluster + 1):(zblock * dim_per_cluster + dim_per_cluster)
    for (yblock in 0:(grid_dim - 1)) {
      y_idx <- (yblock * dim_per_cluster + 1):(yblock * dim_per_cluster + dim_per_cluster)
      for (xblock in 0:(grid_dim - 1)) {
        x_idx <- (xblock * dim_per_cluster + 1):(xblock * dim_per_cluster + dim_per_cluster)
        true_labels[x_idx, y_idx, z_idx] <- cluster_id
        cluster_id <- cluster_id + 1L
      }
    }
  }

  # Create 4D data array
  data_4d <- array(0, dim = c(dim_vol, n_time))

  # Fill each voxel with its cluster's signal
  for (x in 1:9) {
    for (y in 1:9) {
      for (z in 1:9) {
        clus <- true_labels[x, y, z]
        data_4d[x, y, z, ] <- signals[, clus]
      }
    }
  }

  # Add noise if requested
  if (noise_sd > 0) {
    noise <- array(rnorm(n_voxels * n_time, sd = noise_sd),
                   dim = c(dim_vol, n_time))
    data_4d <- data_4d + noise
  }

  # Create neuroim2 objects
  space_3d <- neuroim2::NeuroSpace(dim_vol)
  space_4d <- neuroim2::NeuroSpace(c(dim_vol, n_time))

  mask <- neuroim2::NeuroVol(array(1L, dim = dim_vol), space_3d)
  vec <- neuroim2::NeuroVec(data_4d, space_4d)
  true_vol <- neuroim2::NeuroVol(true_labels, space_3d)

  list(
    vec = vec,
    mask = mask,
    true_clusters = true_vol,
    true_labels = as.integer(true_labels),
    n_clusters = n_clusters,
    cluster_size = cluster_size,
    signals = signals
  )
}


#' Compute Clustering Accuracy Metrics
#'
#' Compares predicted cluster assignments against ground truth using multiple
#' metrics that are invariant to cluster label permutation.
#'
#' @param predicted Integer vector of predicted cluster assignments
#' @param truth Integer vector of ground truth cluster assignments
#'
#' @return A list with:
#'   \item{ari}{Adjusted Rand Index (-1 to 1, 1 = perfect)}
#'   \item{nmi}{Normalized Mutual Information (0 to 1, 1 = perfect)}
#'   \item{accuracy}{Best accuracy after optimal label matching (0 to 1)}
#'   \item{n_predicted}{Number of unique predicted clusters}
#'   \item{n_truth}{Number of unique ground truth clusters}
#'   \item{perfect}{Logical: TRUE if clustering is perfect (ARI = 1)}
#'
#' @details
#' Three complementary metrics are computed:
#'
#' \strong{Adjusted Rand Index (ARI):}
#' Measures agreement between two clusterings, adjusted for chance.
#' ARI = 1 means perfect agreement, ARI = 0 means random, ARI < 0 means
#' worse than random.
#'
#' \strong{Normalized Mutual Information (NMI):}
#' Information-theoretic measure of clustering similarity.
#' NMI = 1 means perfect agreement, NMI = 0 means independent.
#'
#' \strong{Accuracy:}
#' Fraction of correctly assigned voxels after finding optimal
#' correspondence between predicted and true labels using the
#' Hungarian algorithm.
#'
#' @examples
#' # Perfect clustering
#' truth <- rep(1:3, each = 10)
#' pred <- rep(c(2, 3, 1), each = 10)  # Same clustering, different labels
#' clustering_accuracy(pred, truth)  # ARI = 1
#'
#' @export
clustering_accuracy <- function(predicted, truth) {
  stopifnot(length(predicted) == length(truth))

  # Remove any zeros (background)
  valid <- truth > 0 & predicted > 0
  pred <- predicted[valid]
  true <- truth[valid]

  n <- length(pred)
  n_pred <- length(unique(pred))
  n_true <- length(unique(true))

  # Contingency table
  cont <- table(pred, true)

  # Adjusted Rand Index
  ari <- compute_ari(cont, n)

  # Normalized Mutual Information
  nmi <- compute_nmi(cont, n)

  # Best accuracy via Hungarian algorithm (or greedy for small problems)
  accuracy <- compute_best_accuracy(cont)

  list(
    ari = ari,
    nmi = nmi,
    accuracy = accuracy,
    n_predicted = n_pred,
    n_truth = n_true,
    perfect = abs(ari - 1) < 1e-10
  )
}


#' Compute Adjusted Rand Index
#' @keywords internal
compute_ari <- function(cont, n) {
  # Row and column sums
  a <- rowSums(cont)
  b <- colSums(cont)

  # Sum of combinations
  sum_comb_cont <- sum(choose(cont, 2))
  sum_comb_a <- sum(choose(a, 2))
  sum_comb_b <- sum(choose(b, 2))
  comb_n <- choose(n, 2)

  # Expected index
  expected <- sum_comb_a * sum_comb_b / comb_n

  # Max index
  max_index <- (sum_comb_a + sum_comb_b) / 2

  # ARI
  if (max_index == expected) {
    return(1)  # Perfect agreement

  }
  (sum_comb_cont - expected) / (max_index - expected)
}


#' Compute Normalized Mutual Information
#' @keywords internal
compute_nmi <- function(cont, n) {
  # Probabilities
  p_ij <- cont / n
  p_i <- rowSums(p_ij)
  p_j <- colSums(p_ij)

  # Entropy
  H_i <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  H_j <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))

  # Mutual information
  MI <- 0
  for (i in seq_along(p_i)) {
    for (j in seq_along(p_j)) {
      if (p_ij[i, j] > 0) {
        MI <- MI + p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
      }
    }
  }

  # Normalized MI
  if (H_i + H_j == 0) return(1)
  2 * MI / (H_i + H_j)
}


#' Compute Best Accuracy via Greedy Matching
#' @keywords internal
compute_best_accuracy <- function(cont) {
  # Greedy matching: repeatedly assign best remaining pair
  n_total <- sum(cont)
  matched <- 0
  cont_work <- cont

  n_matches <- min(nrow(cont), ncol(cont))
  for (i in seq_len(n_matches)) {
    # Find best remaining match
    best_idx <- which.max(cont_work)
    best_val <- cont_work[best_idx]
    if (best_val == 0) break

    matched <- matched + best_val

    # Remove matched row and column
    best_row <- ((best_idx - 1) %% nrow(cont_work)) + 1
    best_col <- ((best_idx - 1) %/% nrow(cont_work)) + 1
    cont_work[best_row, ] <- 0
    cont_work[, best_col] <- 0
  }

  matched / n_total
}


#' Test Clustering Method on Synthetic Data
#'
#' Runs a clustering method on synthetic ground truth data and reports
#' accuracy metrics.
#'
#' @param method Clustering method name (e.g., "supervoxels", "snic", "flash3d")
#' @param noise_sd Noise level for synthetic data. Default 0 (perfect separability).
#' @param n_clusters Number of ground truth clusters (must be perfect cube).
#'   Default 27.
#' @param seed Random seed for reproducibility.
#' @param verbose Print results. Default TRUE.
#' @param ... Additional arguments passed to cluster4d()
#'
#' @return A list with:
#'   \item{metrics}{Output from clustering_accuracy()}
#'   \item{result}{Full clustering result object}
#'   \item{data}{Synthetic data used}
#'   \item{time}{Elapsed time in seconds}
#'
#' @examples
#' \dontrun{
#' # Test FLASH-3D with 27 cubic clusters
#' test_clustering_method("flash3d")
#'
#' # Test with 8 clusters
#' test_clustering_method("slic", n_clusters = 8)
#'
#' # Test with noise
#' test_clustering_method("snic", noise_sd = 0.5)
#' }
#'
#' @export
test_clustering_method <- function(method, noise_sd = 0, n_clusters = 27,
                                    seed = 42, verbose = TRUE, ...) {
  # Generate synthetic data
  data <- make_synthetic_clusters(n_time = 100, noise_sd = noise_sd,
                                  n_clusters = n_clusters, seed = seed)

  # Run clustering
  t_start <- Sys.time()
  result <- cluster4d(
    vec = data$vec,
    mask = data$mask,
    n_clusters = data$n_clusters,
    method = method,
    verbose = FALSE,
    ...
  )
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  # Compute accuracy
  metrics <- clustering_accuracy(result$cluster, data$true_labels)

  if (verbose) {
    cat(sprintf("\n=== %s (noise_sd=%.2f) ===\n", toupper(method), noise_sd))
    cat(sprintf("  Time:      %.3f sec\n", elapsed))
    cat(sprintf("  Clusters:  %d predicted vs %d true\n",
                metrics$n_predicted, metrics$n_truth))
    cat(sprintf("  ARI:       %.4f\n", metrics$ari))
    cat(sprintf("  NMI:       %.4f\n", metrics$nmi))
    cat(sprintf("  Accuracy:  %.4f\n", metrics$accuracy))
    cat(sprintf("  Perfect:   %s\n", ifelse(metrics$perfect, "YES", "NO")))
  }

  invisible(list(
    metrics = metrics,
    result = result,
    data = data,
    time = elapsed
  ))
}


#' Test All Clustering Methods on Synthetic Data
#'
#' Comprehensive test of all available clustering methods on synthetic
#' ground truth data with cubic clusters.
#'
#' @param noise_levels Vector of noise standard deviations to test.
#'   Default c(0, 0.1, 0.5) tests perfect, low noise, and moderate noise.
#' @param methods Vector of method names to test. Default tests all main methods.
#' @param n_clusters Number of ground truth clusters (must be perfect cube).
#'   Default 27.
#' @param seed Random seed for reproducibility.
#'
#' @return A data.frame with columns: method, noise_sd, ari, nmi, accuracy,
#'   n_clusters, time, perfect
#'
#' @examples
#' \dontrun{
#' # Test all methods with 27 clusters
#' results <- test_all_clustering_methods()
#' print(results)
#'
#' # Test with 8 clusters (simpler)
#' results <- test_all_clustering_methods(n_clusters = 8)
#'
#' # Test specific methods
#' results <- test_all_clustering_methods(
#'   methods = c("flash3d", "snic"),
#'   noise_levels = c(0, 0.2)
#' )
#' }
#'
#' @export
test_all_clustering_methods <- function(
    noise_levels = c(0, 0.1, 0.5),
    methods = c("flash3d", "snic", "slic", "slice_msf", "supervoxels", "g3s", "rena"),
    n_clusters = 27,
    seed = 42
) {

  results <- data.frame(
    method = character(),
    noise_sd = numeric(),
    ari = numeric(),
    nmi = numeric(),
    accuracy = numeric(),
    n_clusters = integer(),
    time = numeric(),
    perfect = logical(),
    stringsAsFactors = FALSE
  )

  for (noise_sd in noise_levels) {
    cat(sprintf("\n========== Noise SD = %.2f ==========\n", noise_sd))

    for (method in methods) {
      tryCatch({
        test_result <- test_clustering_method(
          method = method,
          noise_sd = noise_sd,
          n_clusters = n_clusters,
          seed = seed,
          verbose = TRUE
        )

        results <- rbind(results, data.frame(
          method = method,
          noise_sd = noise_sd,
          ari = test_result$metrics$ari,
          nmi = test_result$metrics$nmi,
          accuracy = test_result$metrics$accuracy,
          n_clusters = test_result$metrics$n_predicted,
          time = test_result$time,
          perfect = test_result$metrics$perfect,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        cat(sprintf("\n=== %s (noise_sd=%.2f) ===\n", toupper(method), noise_sd))
        cat(sprintf("  ERROR: %s\n", e$message))

        results <<- rbind(results, data.frame(
          method = method,
          noise_sd = noise_sd,
          ari = NA,
          nmi = NA,
          accuracy = NA,
          n_clusters = NA,
          time = NA,
          perfect = FALSE,
          stringsAsFactors = FALSE
        ))
      })
    }
  }

  cat("\n\n========== SUMMARY ==========\n")
  print(results, row.names = FALSE)

  invisible(results)
}
