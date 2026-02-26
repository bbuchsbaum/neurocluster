#' Sparse Markov Clustering (MCL) Backend
#'
#' Internal sparse MCL implementation for graph clustering of masked voxels.
#' The implementation keeps iterates sparse via per-column pruning to remain
#' practical on neuroimaging graphs.
#'
#' @keywords internal

.mcl_col_normalize <- function(mat) {
  mat <- as(mat, "dgCMatrix")
  cs <- Matrix::colSums(mat)
  nz <- cs > 0
  if (!all(nz)) {
    empties <- which(!nz)
    if (length(empties) > 0) {
      mat <- mat + Matrix::sparseMatrix(
        i = empties,
        j = empties,
        x = rep(1, length(empties)),
        dims = dim(mat)
      )
      cs <- Matrix::colSums(mat)
    }
  }

  inv <- 1 / pmax(cs, .Machine$double.eps)
  nnz_per_col <- diff(mat@p)
  if (length(mat@x) > 0L) {
    mat@x <- mat@x * rep.int(inv, nnz_per_col)
  }
  mat
}

.mcl_prune_sparse <- function(mat, max_per_col = 128L, min_value = 1e-6) {
  mat <- as(mat, "dgCMatrix")
  if (length(mat@x) == 0L) return(mat)
  max_per_col <- as.integer(max_per_col)
  if (max_per_col < 1L) {
    stop(".mcl_prune_sparse: max_per_col must be >= 1")
  }

  min_keep <- if (is.finite(min_value) && min_value > 0) as.numeric(min_value) else -Inf

  if (!exists("mcl_prune_sparse_cpp", mode = "function")) {
    # Conservative fallback if C++ symbol is unavailable.
    if (is.finite(min_value) && min_value > 0) {
      mat <- Matrix::drop0(mat, tol = min_value)
      if (length(mat@x) == 0L) return(mat)
    }
    nnz_per_col <- diff(mat@p)
    if (all(nnz_per_col <= max_per_col)) {
      return(mat)
    }

    n <- ncol(mat)
    p <- mat@p
    i <- mat@i
    x <- mat@x
    new_i <- vector("list", n)
    new_x <- vector("list", n)

    for (col in seq_len(n)) {
      start <- p[col] + 1L
      end <- p[col + 1L]
      if (start > end) next

      idx <- start:end
      if ((end - start + 1L) > max_per_col) {
        vals <- x[idx]
        top <- order(vals, decreasing = TRUE, method = "radix")[seq_len(max_per_col)]
        idx <- idx[top]
      }
      new_i[[col]] <- i[idx] + 1L
      new_x[[col]] <- x[idx]
    }

    lens <- lengths(new_i)
    if (!any(lens > 0L)) {
      return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = mat@Dim))
    }

    return(Matrix::sparseMatrix(
      i = as.integer(unlist(new_i, use.names = FALSE)),
      j = rep.int(seq_len(n), lens),
      x = as.numeric(unlist(new_x, use.names = FALSE)),
      dims = mat@Dim,
      repr = "C"
    ))
  }

  out <- mcl_prune_sparse_cpp(
    p = mat@p,
    i = mat@i,
    x = mat@x,
    ncol = as.integer(ncol(mat)),
    max_per_col = max_per_col,
    min_value = min_keep
  )

  if (length(out$x) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = mat@Dim))
  }

  methods::new(
    "dgCMatrix",
    i = as.integer(out$i),
    p = as.integer(out$p),
    x = as.numeric(out$x),
    Dim = as.integer(mat@Dim),
    Dimnames = mat@Dimnames
  )
}

.mcl_labels_from_flow <- function(flow) {
  flow <- as(flow, "dgCMatrix")
  n <- ncol(flow)

  p <- flow@p
  rows <- flow@i
  vals <- flow@x

  attractor <- integer(n)

  for (col in seq_len(n)) {
    start <- p[col] + 1L
    end <- p[col + 1L]
    if (start > end) {
      attractor[col] <- col
      next
    }

    idx <- start:end
    best <- which.max(vals[idx])
    attractor[col] <- rows[idx][best] + 1L
  }

  uniq <- sort(unique(attractor))
  labels <- match(attractor, uniq)
  as.integer(labels)
}

.mcl_row_zscore <- function(x) {
  mu <- rowMeans(x)
  xc <- x - mu

  if (ncol(x) > 1L) {
    s <- sqrt(rowSums(xc * xc) / (ncol(x) - 1L))
  } else {
    s <- rep(1, nrow(x))
  }

  s[!is.finite(s) | s < 1e-8] <- 1
  xc / s
}

.mcl_edge_similarity <- function(features, edge_i, edge_j,
                                 feature_metric = c("correlation", "euclidean"),
                                 feature_sigma = NULL,
                                 chunk_size = 50000L) {
  feature_metric <- match.arg(feature_metric)

  n_edges <- length(edge_i)
  if (n_edges == 0L) {
    return(numeric(0))
  }

  if (feature_metric == "correlation") {
    if (ncol(features) < 2L) {
      feature_metric <- "euclidean"
    } else {
      z <- .mcl_row_zscore(features)
      sim <- numeric(n_edges)
      denom <- ncol(z) - 1L

      for (start in seq.int(1L, n_edges, by = chunk_size)) {
        end <- min(n_edges, start + chunk_size - 1L)
        idx <- start:end
        sim[idx] <- rowSums(z[edge_i[idx], , drop = FALSE] * z[edge_j[idx], , drop = FALSE]) / denom
      }

      sim <- (sim + 1) / 2
      sim[!is.finite(sim)] <- 0
      sim <- pmin(1, pmax(0, sim))
      return(sim)
    }
  }

  # Euclidean heat-kernel fallback
  d2 <- compute_masked_distances_cpp(
    t(features),
    as.integer(edge_i),
    as.integer(edge_j)
  )

  if (is.null(feature_sigma) || !is.finite(feature_sigma) || feature_sigma <= 0) {
    d <- sqrt(d2)
    feature_sigma <- stats::median(d[d > 0], na.rm = TRUE)
    if (!is.finite(feature_sigma) || feature_sigma <= 0) {
      feature_sigma <- 1
    }
  }

  sim <- exp(-d2 / (2 * feature_sigma * feature_sigma))
  sim[!is.finite(sim)] <- 0
  sim
}

.mcl_build_weighted_graph <- function(features, coords, mask, mask_idx,
                                      connectivity = 6,
                                      spatial_weight = 0.2,
                                      feature_metric = c("correlation", "euclidean"),
                                      feature_sigma = NULL,
                                      spatial_sigma = NULL,
                                      min_edge_weight = 1e-6,
                                      verbose = FALSE) {
  feature_metric <- match.arg(feature_metric)

  if (!connectivity %in% c(6L, 18L, 26L, 27L)) {
    stop("cluster4d_mcl: connectivity must be one of 6, 18, 26, or 27")
  }

  adj <- build_grid_adjacency(mask, mask_idx, connectivity)
  adj_sum <- Matrix::summary(adj)

  keep <- adj_sum$i < adj_sum$j
  edge_i <- adj_sum$i[keep]
  edge_j <- adj_sum$j[keep]

  if (length(edge_i) == 0L) {
    stop("cluster4d_mcl: no edges found in the masked grid graph")
  }

  if (verbose) {
    message("cluster4d_mcl: computing edge similarities for ", length(edge_i), " edges")
  }

  feat_sim <- .mcl_edge_similarity(
    features = features,
    edge_i = edge_i,
    edge_j = edge_j,
    feature_metric = feature_metric,
    feature_sigma = feature_sigma
  )

  dxyz2 <- rowSums((coords[edge_i, , drop = FALSE] - coords[edge_j, , drop = FALSE])^2)

  if (is.null(spatial_sigma) || !is.finite(spatial_sigma) || spatial_sigma <= 0) {
    d <- sqrt(dxyz2)
    spatial_sigma <- stats::median(d[d > 0], na.rm = TRUE)
    if (!is.finite(spatial_sigma) || spatial_sigma <= 0) {
      spatial_sigma <- 1
    }
  }

  spatial_sim <- exp(-dxyz2 / (2 * spatial_sigma * spatial_sigma))
  spatial_sim[!is.finite(spatial_sim)] <- 0

  w <- (1 - spatial_weight) * feat_sim + spatial_weight * spatial_sim
  w <- pmax(w, min_edge_weight)

  n <- nrow(features)
  Matrix::sparseMatrix(
    i = c(edge_i, edge_j),
    j = c(edge_j, edge_i),
    x = c(w, w),
    dims = c(n, n)
  )
}

.mcl_sparse <- function(graph,
                        inflation = 2.0,
                        expansion = 2L,
                        max_iter = 30L,
                        tol = 1e-4,
                        prune_k = 128L,
                        prune_threshold = 1e-6,
                        loop_weight = 1,
                        verbose = FALSE) {
  if (inflation <= 1) {
    stop("cluster4d_mcl: inflation must be > 1")
  }

  if (expansion < 2L) {
    expansion <- 2L
  }

  n <- nrow(graph)
  flow <- graph + Matrix::Diagonal(n = n, x = rep(loop_weight, n))
  flow <- .mcl_col_normalize(flow)

  converged <- FALSE
  iter <- 0L

  for (it in seq_len(max_iter)) {
    iter <- it
    prev <- flow

    expanded <- flow
    if (expansion > 1L) {
      for (k in seq_len(expansion - 1L)) {
        expanded <- expanded %*% flow
      }
    }

    expanded <- Matrix::drop0(as(expanded, "dgCMatrix"), tol = prune_threshold)

    if (length(expanded@x) > 0L) {
      expanded@x <- expanded@x ^ inflation
    }

    flow <- .mcl_col_normalize(expanded)
    flow <- .mcl_prune_sparse(flow, max_per_col = prune_k, min_value = prune_threshold)
    flow <- .mcl_col_normalize(flow)

    delta <- flow - prev
    max_delta <- if (length(delta@x) == 0L) 0 else max(abs(delta@x))

    if (verbose) {
      message(sprintf("cluster4d_mcl: iter=%d max_delta=%.3e nnz=%d", it, max_delta, length(flow@x)))
    }

    if (max_delta < tol) {
      converged <- TRUE
      break
    }
  }

  labels <- .mcl_labels_from_flow(flow)

  list(
    labels = labels,
    flow = flow,
    converged = converged,
    iterations = iter
  )
}

#' Cluster4d using sparse Markov Clustering (MCL)
#'
#' Builds a sparse voxel graph, runs an MCL flow process, and maps attractors to
#' cluster labels. The implementation is optimized for sparse neuroimaging graphs.
#'
#' @inheritParams cluster4d
#' @param inflation MCL inflation parameter (>1). Higher values produce finer clusters.
#' @param expansion MCL expansion power (integer >= 2). Default 2.
#' @param loop_weight Self-loop weight added before normalization.
#' @param prune_k Maximum nonzero entries kept per column during MCL iterations.
#'   Smaller values increase speed and reduce memory.
#' @param prune_threshold Minimum value retained during pruning.
#' @param tol Convergence tolerance on max absolute matrix delta.
#' @param feature_metric Feature similarity metric: `"correlation"` or `"euclidean"`.
#' @param feature_sigma Optional bandwidth for euclidean heat-kernel feature similarity.
#' @param spatial_sigma Optional bandwidth for spatial heat-kernel similarity.
#' @param exact_k If TRUE, force exactly `n_clusters` using `force_exact_k()` post-processing.
#' @param ... Reserved for future options.
#'
#' @return A `cluster4d_result` object.
#' @export
cluster4d_mcl <- function(vec, mask, n_clusters = 100,
                          spatial_weight = 0.2,
                          max_iterations = 8,
                          connectivity = 6,
                          parallel = TRUE,
                          verbose = FALSE,
                          inflation = 1.6,
                          expansion = 2L,
                          loop_weight = 1,
                          prune_k = NULL,
                          prune_threshold = 1e-6,
                          tol = 1e-4,
                          feature_metric = c("correlation", "euclidean"),
                          feature_sigma = NULL,
                          spatial_sigma = NULL,
                          exact_k = FALSE,
                          ...) {

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_mcl")

  if (!is.numeric(spatial_weight) || length(spatial_weight) != 1 ||
      !is.finite(spatial_weight) || spatial_weight < 0 || spatial_weight > 1) {
    stop("cluster4d_mcl: spatial_weight must be in [0, 1]")
  }

  feature_metric <- match.arg(feature_metric)

  data_prep <- prepare_cluster4d_data(
    vec = vec,
    mask = mask,
    scale_features = FALSE,
    scale_coords = FALSE
  )

  if (is.null(prune_k)) {
    # Keep a conservative cap per column for speed; larger values rarely help
    # in these sparse voxel graphs and can significantly increase runtime.
    prune_k <- as.integer(max(16, min(64, n_clusters)))
  } else {
    prune_k <- as.integer(prune_k)
  }

  if (prune_k < 8L) {
    prune_k <- 8L
  }

  graph <- .mcl_build_weighted_graph(
    features = data_prep$features,
    coords = data_prep$coords,
    mask = mask,
    mask_idx = data_prep$mask_idx,
    connectivity = as.integer(connectivity),
    spatial_weight = spatial_weight,
    feature_metric = feature_metric,
    feature_sigma = feature_sigma,
    spatial_sigma = spatial_sigma,
    min_edge_weight = prune_threshold,
    verbose = verbose
  )

  mcl_fit <- .mcl_sparse(
    graph = graph,
    inflation = inflation,
    expansion = as.integer(expansion),
    max_iter = as.integer(max_iterations),
    tol = tol,
    prune_k = prune_k,
    prune_threshold = prune_threshold,
    loop_weight = loop_weight,
    verbose = verbose
  )

  labels <- mcl_fit$labels

  if (exact_k) {
    labels <- force_exact_k(labels, data_prep$features, n_clusters)
  }

  result <- create_cluster4d_result(
    labels = labels,
    mask = mask,
    data_prep = data_prep,
    method = "mcl",
    parameters = list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      max_iterations = as.integer(max_iterations),
      connectivity = as.integer(connectivity),
      parallel = isTRUE(parallel),
      inflation = inflation,
      expansion = as.integer(expansion),
      loop_weight = loop_weight,
      prune_k = prune_k,
      prune_threshold = prune_threshold,
      tol = tol,
      feature_metric = feature_metric,
      feature_sigma = feature_sigma,
      spatial_sigma = spatial_sigma,
      exact_k = isTRUE(exact_k)
    ),
    metadata = list(
      converged = mcl_fit$converged,
      iterations = mcl_fit$iterations,
      nnz_final = length(mcl_fit$flow@x)
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  result
}
