#' Prepare commute-clustering features under an explicit bad-voxel policy
#'
#' @keywords internal
#' @noRd
.commute_prepare_features <- function(features,
                                      bad_voxel_policy = c("error", "zero_constant"),
                                      variance_tol = 1e-9) {
  bad_voxel_policy <- match.arg(bad_voxel_policy)
  variance_tol <- .cluster4d_scalar_number(
    variance_tol, "variance_tol", "commute_cluster", lower = 0
  )
  if (!is.matrix(features) || !is.numeric(features) || nrow(features) < 1L ||
      ncol(features) < 1L) {
    stop("commute_cluster: features must be a non-empty numeric matrix",
         call. = FALSE)
  }

  nonfinite <- which(rowSums(!is.finite(features)) > 0L)
  if (length(nonfinite)) {
    stop(
      "commute_cluster: non-finite data in masked voxel(s): ",
      paste(head(nonfinite, 10L), collapse = ", "),
      if (length(nonfinite) > 10L) " ..." else "",
      "; NA/Inf values are never silently imputed",
      call. = FALSE
    )
  }

  centered <- features - rowMeans(features)
  variances <- if (ncol(features) > 1L) {
    rowSums(centered * centered) / (ncol(features) - 1L)
  } else {
    rep(0, nrow(features))
  }
  constant <- which(variances <= variance_tol^2)
  if (length(constant) && identical(bad_voxel_policy, "error")) {
    stop(
      "commute_cluster: ", length(constant),
      " masked voxel(s) have variance at or below variance_tol; ",
      "use bad_voxel_policy='zero_constant' for deterministic zero features",
      call. = FALSE
    )
  }

  scales <- sqrt(variances)
  good <- scales > variance_tol
  standardized <- matrix(0, nrow(features), ncol(features))
  if (any(good)) {
    standardized[good, ] <- centered[good, , drop = FALSE] / scales[good]
  }
  if (any(!is.finite(standardized))) {
    stop("commute_cluster: standardization produced non-finite features",
         call. = FALSE)
  }

  list(
    features = standardized,
    provenance = list(
      policy = bad_voxel_policy,
      nonfinite_policy = "error",
      variance_tol = variance_tol,
      constant_voxels = as.integer(constant),
      constant_replacement = if (length(constant)) {
        "zero standardized feature vector"
      } else {
        "none"
      },
      random_noise_injected = FALSE
    )
  )
}


.commute_knn_graph <- function(coords, features, connectivity,
                               alpha, spatial_sigma, feature_sigma,
                               weight_mode = c("binary", "heat")) {
  weight_mode <- match.arg(weight_mode)
  n <- nrow(coords)
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "commute_cluster",
    lower = 1, upper = n - 1L, integer = TRUE
  )
  if (!is.matrix(coords) || ncol(coords) != 3L ||
      any(!is.finite(coords))) {
    stop("commute_cluster: coords must be a finite N-by-3 matrix",
         call. = FALSE)
  }
  if (!is.matrix(features) || nrow(features) != n ||
      any(!is.finite(features))) {
    stop("commute_cluster: features must be a finite N-row matrix",
         call. = FALSE)
  }

  neighbors <- FNN::get.knn(coords, k = connectivity)
  directed_i <- rep(seq_len(n), each = connectivity)
  directed_j <- as.integer(t(neighbors$nn.index))
  pairs <- cbind(
    pmin(directed_i, directed_j),
    pmax(directed_i, directed_j)
  )
  pairs <- unique(pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE])
  left <- pairs[, 1L]
  right <- pairs[, 2L]

  if (identical(weight_mode, "binary")) {
    weights <- rep(1, length(left))
  } else {
    spatial_d2 <- rowSums(
      (coords[left, , drop = FALSE] - coords[right, , drop = FALSE])^2
    )
    feature_d2 <- rowSums(
      (features[left, , drop = FALSE] - features[right, , drop = FALSE])^2
    )
    spatial_similarity <- exp(-spatial_d2 / (2 * spatial_sigma^2))
    feature_similarity <- exp(-feature_d2 / (2 * feature_sigma^2))
    weights <- (1 - alpha) * spatial_similarity + alpha * feature_similarity
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("commute_cluster: graph weights must be positive and finite",
         call. = FALSE)
  }

  graph <- Matrix::sparseMatrix(
    i = c(left, right),
    j = c(right, left),
    x = c(weights, weights),
    dims = c(n, n)
  )
  components <- igraph::components(
    igraph::graph_from_adjacency_matrix(
      graph, mode = "undirected", weighted = TRUE, diag = FALSE
    )
  )
  list(
    graph = as(graph, "dgCMatrix"),
    edges = pairs,
    weights = weights,
    components = as.integer(components$membership),
    n_components = as.integer(components$no),
    contract = "symmetric union of directed physical-coordinate k-nearest-neighbor edges"
  )
}


.commute_spectral_embedding <- function(graph, ncomp,
                                        eigen_tol = 1e-10) {
  graph <- as(graph, "dgCMatrix")
  n <- nrow(graph)
  if (n < 2L || n != ncol(graph) || any(!is.finite(graph@x)) ||
      any(graph@x < 0) || !Matrix::isSymmetric(graph, checkDN = FALSE)) {
    stop("commute_cluster: graph must be finite, non-negative, symmetric, and at least 2-by-2",
         call. = FALSE)
  }
  ncomp <- .cluster4d_scalar_number(
    ncomp, "ncomp", "commute_cluster", lower = 1,
    upper = n - 1L, integer = TRUE
  )
  eigen_tol <- .cluster4d_scalar_number(
    eigen_tol, "eigen_tol", "commute_cluster", lower = 0
  )
  if (eigen_tol <= 0) {
    stop("commute_cluster: eigen_tol must be positive", call. = FALSE)
  }

  degree <- Matrix::rowSums(graph)
  if (any(degree <= 0)) {
    stop("commute_cluster: graph contains isolated voxels", call. = FALSE)
  }
  laplacian <- diag(as.numeric(degree)) - as.matrix(graph)
  decomposition <- eigen(laplacian, symmetric = TRUE)
  order_idx <- order(decomposition$values)
  values <- decomposition$values[order_idx]
  vectors <- decomposition$vectors[, order_idx, drop = FALSE]
  threshold <- max(
    eigen_tol,
    max(abs(values), 1) * .Machine$double.eps * n * 10
  )
  positive <- which(values > threshold)
  zero_count <- sum(abs(values) <= threshold)
  if (zero_count != 1L) {
    stop(
      "commute_cluster: k-nearest graph has ", zero_count,
      " connected components; increase connectivity",
      call. = FALSE
    )
  }
  if (length(positive) < ncomp) {
    stop(
      "commute_cluster: ncomp exceeds the positive Laplacian eigenspace (",
      length(positive), ")",
      call. = FALSE
    )
  }
  selected <- positive[seq_len(ncomp)]
  volume <- sum(degree)
  embedding <- sweep(
    vectors[, selected, drop = FALSE],
    2L, sqrt(values[selected]), "/"
  ) * sqrt(volume)

  list(
    embedding = unname(embedding),
    eigenvalues = unname(values[selected]),
    eigenvectors = unname(vectors[, selected, drop = FALSE]),
    laplacian = laplacian,
    graph_volume = volume,
    zero_eigenvalues = as.integer(zero_count)
  )
}


.commute_initial_centers <- function(embedding, K) {
  n <- nrow(embedding)
  if (K == 1L) return(1L)
  selected <- integer(K)
  selected[1L] <- order(-rowSums(embedding * embedding), seq_len(n))[1L]
  min_distance <- rep(Inf, n)
  for (index in 2:K) {
    latest <- selected[index - 1L]
    distance <- rowSums(
      (embedding - matrix(
        embedding[latest, ], nrow = n, ncol = ncol(embedding), byrow = TRUE
      ))^2
    )
    min_distance <- pmin(min_distance, distance)
    min_distance[selected[seq_len(index - 1L)]] <- -Inf
    selected[index] <- order(-min_distance, seq_len(n))[1L]
  }
  if (anyDuplicated(selected)) {
    stop("commute_cluster: deterministic center initialization failed",
         call. = FALSE)
  }
  selected
}


.commute_local_rng <- function(expression, seed = 1L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}


#' Commute-time clustering on a physical k-nearest-neighbor graph
#'
#' Builds an explicitly declared symmetric k-nearest-neighbor graph in physical
#' coordinate space, computes a Laplacian commute-time embedding, and applies
#' deterministically initialized k-means.
#'
#' @param bvec A NeuroVec or SparseNeuroVec supplying voxel features.
#' @param mask A NeuroVol mask. Finite values strictly greater than zero are included.
#' @param K Finite integer cluster count from one through the included voxel count.
#' @param ncomp Finite integer embedding dimension from one through N - 1.
#' @param alpha Feature-similarity weight from zero through one.
#' @param sigma1 Positive physical-distance heat-kernel bandwidth.
#' @param sigma2 Positive feature-distance heat-kernel bandwidth.
#' @param connectivity Number of nearest physical neighbors per directed search.
#'   The final undirected graph is the symmetric union of those directed edges;
#'   this is not masked grid adjacency and may bridge mask holes.
#' @param weight_mode `"binary"` or `"heat"` edge weights.
#' @param bad_voxel_policy `"error"` (default), or `"zero_constant"` to map
#'   finite zero-variance voxel series to deterministic zero standardized vectors.
#'   NA and Inf always fail closed. No noise is injected.
#' @param variance_tol Non-negative threshold defining a constant voxel series.
#' @param eigen_tol Positive tolerance for the Laplacian zero eigenspace.
#' @param verbose Logical; print progress.
#'
#' @return A `commute_time_cluster_result` with contiguous labels, K-by-T
#'   original-space feature means, K-by-3 physical coordinate means, the N-by-ncomp
#'   embedding, graph and spectral provenance, and explicit data-policy provenance.
#' @details The implementation is deterministic. Any RNG initialization performed
#'   inside `stats::kmeans()` is scoped to a fixed local seed, after which the
#'   caller's prior RNG state (including absence of `.Random.seed`) is restored
#'   exactly. It forms a dense Laplacian eigendecomposition,
#'   requiring quadratic memory and cubic time; use only for small regions.
#'
#' @importFrom neuroim2 ClusteredNeuroVol
#' @importFrom stats kmeans
#' @export
commute_cluster <- function(bvec,
                            mask,
                            K = 100,
                            ncomp = ceiling(sqrt(K * 2)),
                            alpha = 0.5,
                            sigma1 = 0.73,
                            sigma2 = 5,
                            connectivity = 6L,
                            weight_mode = c("binary", "heat"),
                            bad_voxel_policy = c("error", "zero_constant"),
                            variance_tol = 1e-9,
                            eigen_tol = 1e-10,
                            verbose = FALSE) {
  caller_had_seed <- exists(
    ".Random.seed", envir = .GlobalEnv, inherits = FALSE
  )
  caller_seed <- if (caller_had_seed) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (caller_had_seed) {
      assign(".Random.seed", caller_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  input <- validate_cluster4d_inputs(bvec, mask, K, "commute_cluster")
  K <- input$n_clusters
  n <- input$n_voxels
  if (n < 2L) {
    stop("commute_cluster: at least two masked voxels are required",
         call. = FALSE)
  }
  ncomp <- .cluster4d_scalar_number(
    ncomp, "ncomp", "commute_cluster", lower = 1,
    upper = n - 1L, integer = TRUE
  )
  alpha <- .cluster4d_scalar_number(
    alpha, "alpha", "commute_cluster", lower = 0, upper = 1
  )
  sigma1 <- .cluster4d_scalar_number(
    sigma1, "sigma1", "commute_cluster", lower = 0
  )
  sigma2 <- .cluster4d_scalar_number(
    sigma2, "sigma2", "commute_cluster", lower = 0
  )
  if (sigma1 <= 0 || sigma2 <= 0) {
    stop("commute_cluster: sigma1 and sigma2 must be positive", call. = FALSE)
  }
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "commute_cluster",
    lower = 1, upper = n - 1L, integer = TRUE
  )
  verbose <- .cluster4d_scalar_logical(verbose, "verbose", "commute_cluster")
  weight_mode <- match.arg(weight_mode)
  bad_voxel_policy <- match.arg(bad_voxel_policy)

  original <- .cluster4d_original_data(bvec, mask, "commute_cluster")
  prepared <- .commute_prepare_features(
    original$features, bad_voxel_policy, variance_tol
  )
  if (verbose) message("commute_cluster: building physical k-nearest graph")
  graph <- .commute_knn_graph(
    original$coords, prepared$features, connectivity,
    alpha, sigma1, sigma2, weight_mode
  )
  if (graph$n_components != 1L) {
    stop(
      "commute_cluster: k-nearest graph has ", graph$n_components,
      " connected components; increase connectivity",
      call. = FALSE
    )
  }

  if (verbose) message("commute_cluster: computing dense spectral embedding")
  spectral <- .commute_spectral_embedding(graph$graph, ncomp, eigen_tol)
  if (K == n) {
    labels <- seq_len(n)
    kmeans_diagnostics <- list(iter = 0L, tot.withinss = 0)
  } else if (K == 1L) {
    labels <- rep.int(1L, n)
    kmeans_diagnostics <- list(
      iter = 0L,
      tot.withinss = sum(
        scale(spectral$embedding, center = TRUE, scale = FALSE)^2
      )
    )
  } else {
    seeds <- .commute_initial_centers(spectral$embedding, K)
    fit <- .commute_local_rng(
      stats::kmeans(
        spectral$embedding,
        centers = spectral$embedding[seeds, , drop = FALSE],
        iter.max = 500,
        algorithm = "Lloyd"
      )
    )
    labels <- as.integer(fit$cluster)
    kmeans_diagnostics <- list(
      iter = as.integer(fit$iter),
      tot.withinss = as.numeric(fit$tot.withinss),
      initial_voxels = as.integer(seeds)
    )
  }

  centers <- compute_cluster_centers(
    labels, original$features, original$coords, method = "mean"
  )
  logical_mask <- neuroim2::NeuroVol(
    cluster4d_mask_array(mask, "commute_cluster"),
    space = neuroim2::space(mask)
  )
  clusvol <- suppressWarnings(
    neuroim2::ClusteredNeuroVol(logical_mask, clusters = labels)
  )
  result <- list(
    labels = labels,
    cluster = labels,
    clusvol = clusvol,
    centers = unname(as.matrix(centers$centers)),
    coord_centers = unname(as.matrix(centers$coord_centers)),
    embedding = spectral$embedding,
    eigenvalues = spectral$eigenvalues,
    actual_k = as.integer(K),
    n_clusters = as.integer(K),
    label_ids = seq_len(K),
    method = "commute_time",
    parameters = list(
      n_clusters_requested = as.integer(K),
      ncomp = ncomp,
      alpha = alpha,
      sigma1 = sigma1,
      sigma2 = sigma2,
      connectivity = connectivity,
      weight_mode = weight_mode,
      bad_voxel_policy = bad_voxel_policy,
      variance_tol = variance_tol,
      eigen_tol = eigen_tol
    ),
    metadata = list(
      data_policy = prepared$provenance,
      graph = list(
        contract = graph$contract,
        directed_k = connectivity,
        undirected_edges = as.integer(nrow(graph$edges)),
        components = graph$n_components
      ),
      spectral = list(
        eigenvalues = spectral$eigenvalues,
        graph_volume = spectral$graph_volume,
        zero_eigenvalues = spectral$zero_eigenvalues
      ),
      kmeans = kmeans_diagnostics,
      rng = list(
        used = K > 1L && K < n,
        local_seed = if (K > 1L && K < n) 1L else NULL,
        purpose = if (K > 1L && K < n) "stats::kmeans internals" else "none",
        global_state_preserved = TRUE
      )
    )
  )
  structure(
    result,
    class = c("commute_time_cluster_result", "cluster_result", "list")
  )
}
