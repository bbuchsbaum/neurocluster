#' Sparse Markov Clustering (MCL) Backend
#'
#' Internal sparse MCL implementation for graph clustering of masked voxels.
#' The implementation keeps iterates sparse via per-column pruning to remain
#' practical on neuroimaging graphs.
#'
#' @keywords internal

.mcl_col_normalize <- function(mat) {
  mat <- as(mat, "dgCMatrix")
  if (nrow(mat) != ncol(mat) || any(!is.finite(mat@x)) ||
      any(mat@x < 0)) {
    stop(".mcl_col_normalize: mat must be square, finite, and non-negative",
         call. = FALSE)
  }
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
  normalized_sums <- Matrix::colSums(mat)
  if (any(!is.finite(mat@x)) || any(mat@x < 0) ||
      any(abs(normalized_sums - 1) > 1e-12)) {
    stop(".mcl_col_normalize: failed to produce stochastic columns",
         call. = FALSE)
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

  if (!is.numeric(min_value) || length(min_value) != 1L ||
      !is.finite(min_value) || min_value < 0) {
    stop(".mcl_prune_sparse: min_value must be finite and non-negative")
  }
  min_keep <- if (min_value > 0) as.numeric(min_value) else -Inf

  if (!exists("mcl_prune_sparse_cpp", mode = "function")) {
    # Conservative fallback if C++ symbol is unavailable.
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
      above <- x[idx] >= min_keep
      if (any(above)) {
        idx <- idx[above]
      } else {
        # Preserve the strongest entry so pruning never creates an empty flow
        # column that would have to be replaced by an invented self-loop.
        vals <- x[idx]
        idx <- idx[order(-vals, i[idx], method = "radix")[1L]]
      }
      if (length(idx) > max_per_col) {
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

.mcl_attractors_from_flow <- function(flow) {
  flow <- as(flow, "dgCMatrix")
  if (nrow(flow) != ncol(flow) || any(!is.finite(flow@x)) ||
      any(flow@x < 0)) {
    stop(".mcl_attractors_from_flow: flow must be square, finite, and non-negative",
         call. = FALSE)
  }
  n <- ncol(flow)

  p <- flow@p
  rows <- flow@i
  vals <- flow@x

  next_node <- integer(n)

  for (col in seq_len(n)) {
    start <- p[col] + 1L
    end <- p[col + 1L]
    if (start > end) {
      next_node[col] <- col
      next
    }

    idx <- start:end
    best <- which.max(vals[idx])
    next_node[col] <- rows[idx][best] + 1L
  }

  # Resolve the functional graph to canonical cycle representatives.  This
  # makes the returned map idempotent even when raw column argmaxes form chains
  # such as 1 -> 2 -> 3.
  attractor <- integer(n)
  for (start_node in seq_len(n)) {
    if (attractor[start_node] != 0L) next
    path <- integer()
    position <- integer(n)
    node <- start_node
    while (attractor[node] == 0L && position[node] == 0L) {
      position[node] <- length(path) + 1L
      path <- c(path, node)
      node <- next_node[node]
    }

    if (attractor[node] != 0L) {
      root <- attractor[node]
    } else {
      cycle <- path[position[node]:length(path)]
      diagonal <- flow[cbind(cycle, cycle)]
      root <- cycle[order(-diagonal, cycle)[1L]]
      attractor[cycle] <- root
      attractor[root] <- root
    }
    attractor[path] <- root
  }

  if (!identical(attractor[attractor], attractor)) {
    stop(".mcl_attractors_from_flow: canonical map is not idempotent",
         call. = FALSE)
  }
  as.integer(attractor)
}

.mcl_labels_from_flow <- function(flow) {
  attractor <- .mcl_attractors_from_flow(flow)
  as.integer(match(attractor, sort(unique(attractor))))
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

  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("cluster4d_mcl: connectivity must be one of 6, 18, or 26")
  }

  adj <- build_grid_adjacency(mask, mask_idx, connectivity)
  adj_sum <- Matrix::summary(adj)

  keep <- adj_sum$i < adj_sum$j
  edge_i <- adj_sum$i[keep]
  edge_j <- adj_sum$j[keep]

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

.mcl_validate_stochastic <- function(flow, context = ".mcl_sparse") {
  flow <- as(flow, "dgCMatrix")
  if (nrow(flow) != ncol(flow) || any(!is.finite(flow@x)) ||
      any(flow@x < 0)) {
    stop(context, ": flow must be square, finite, and non-negative",
         call. = FALSE)
  }
  sums <- Matrix::colSums(flow)
  error <- max(abs(sums - 1), 0)
  if (!is.finite(error) || error > 1e-10) {
    stop(context, ": flow columns are not normalized", call. = FALSE)
  }
  invisible(error)
}


.mcl_step <- function(flow, inflation, expansion, prune_k,
                      prune_threshold) {
  .mcl_validate_stochastic(flow, ".mcl_step input")
  expanded <- flow
  for (power in seq_len(expansion - 1L)) {
    expanded <- expanded %*% flow
  }
  expanded <- as(expanded, "dgCMatrix")
  if (any(!is.finite(expanded@x)) || any(expanded@x < 0)) {
    stop(".mcl_step: expansion produced invalid flow", call. = FALSE)
  }
  if (length(expanded@x) > 0L) {
    expanded@x <- expanded@x ^ inflation
  }
  next_flow <- .mcl_col_normalize(expanded)
  next_flow <- .mcl_prune_sparse(
    next_flow, max_per_col = prune_k, min_value = prune_threshold
  )
  next_flow <- .mcl_col_normalize(next_flow)
  .mcl_validate_stochastic(next_flow, ".mcl_step output")
  next_flow
}


.mcl_sparse <- function(graph,
                        inflation = 2.0,
                        expansion = 2L,
                        max_iter = 30L,
                        tol = 1e-4,
                        prune_k = 128L,
                        prune_threshold = 1e-6,
                        loop_weight = 1,
                        verbose = FALSE,
                        trace = FALSE) {
  graph <- as(graph, "dgCMatrix")
  if (nrow(graph) < 1L || nrow(graph) != ncol(graph) ||
      any(!is.finite(graph@x)) || any(graph@x < 0)) {
    stop("cluster4d_mcl: graph must be non-empty, square, finite, and non-negative",
         call. = FALSE)
  }
  inflation <- .cluster4d_scalar_number(
    inflation, "inflation", "cluster4d_mcl"
  )
  if (inflation <= 1) {
    stop("cluster4d_mcl: inflation must be greater than 1", call. = FALSE)
  }
  expansion <- .cluster4d_scalar_number(
    expansion, "expansion", "cluster4d_mcl", lower = 2, integer = TRUE
  )
  max_iter <- .cluster4d_scalar_number(
    max_iter, "max_iter", "cluster4d_mcl", lower = 1, integer = TRUE
  )
  tol <- .cluster4d_scalar_number(
    tol, "tol", "cluster4d_mcl", lower = 0
  )
  if (tol <= 0) stop("cluster4d_mcl: tol must be positive", call. = FALSE)
  prune_k <- .cluster4d_scalar_number(
    prune_k, "prune_k", "cluster4d_mcl", lower = 1, integer = TRUE
  )
  prune_threshold <- .cluster4d_scalar_number(
    prune_threshold, "prune_threshold", "cluster4d_mcl", lower = 0,
    upper = 1
  )
  if (prune_threshold >= 1) {
    stop("cluster4d_mcl: prune_threshold must be less than 1", call. = FALSE)
  }
  loop_weight <- .cluster4d_scalar_number(
    loop_weight, "loop_weight", "cluster4d_mcl", lower = 0
  )
  verbose <- .cluster4d_scalar_logical(verbose, "verbose", "cluster4d_mcl")
  trace <- .cluster4d_scalar_logical(trace, "trace", "cluster4d_mcl")

  n <- nrow(graph)
  flow <- graph + Matrix::Diagonal(n = n, x = rep(loop_weight, n))
  flow <- .mcl_col_normalize(flow)
  .mcl_validate_stochastic(flow, "cluster4d_mcl initial flow")

  converged <- FALSE
  iter <- 0L
  max_delta <- Inf
  nnz_history <- as.integer(length(flow@x))
  flow_trace <- if (trace) list(flow) else list()

  for (it in seq_len(max_iter)) {
    iter <- it
    prev <- flow
    flow <- .mcl_step(
      prev, inflation, expansion, prune_k, prune_threshold
    )
    delta <- flow - prev
    max_delta <- if (length(delta@x) == 0L) 0 else max(abs(delta@x))
    nnz_history <- c(nnz_history, as.integer(length(flow@x)))
    if (trace) flow_trace[[length(flow_trace) + 1L]] <- flow

    if (verbose) {
      message(sprintf(
        "cluster4d_mcl: iter=%d max_delta=%.3e nnz=%d",
        it, max_delta, length(flow@x)
      ))
    }
    if (max_delta < tol) {
      converged <- TRUE
      break
    }
  }

  attractors <- .mcl_attractors_from_flow(flow)
  labels <- as.integer(match(attractors, sort(unique(attractors))))

  list(
    labels = labels,
    attractors = attractors,
    flow = flow,
    converged = converged,
    iterations = iter,
    max_delta = max_delta,
    nnz_history = nnz_history,
    trace = flow_trace
  )
}


.mcl_inflation_candidates <- function(inflation) {
  candidates <- inflation * c(0.75, 0.9, 1, 1.25, 1.6)
  candidates <- pmax(candidates, 1 + 1e-6)
  sort(unique(as.numeric(candidates)))
}


.mcl_targeted_sparse <- function(graph, target_k, inflation, ...) {
  candidates <- .mcl_inflation_candidates(inflation)
  fits <- lapply(candidates, function(value) {
    .mcl_sparse(graph, inflation = value, ...)
  })
  achieved <- vapply(fits, function(fit) max(fit$labels), integer(1))
  converged <- vapply(fits, `[[`, logical(1), "converged")
  selected <- order(
    abs(achieved - target_k),
    !converged,
    abs(log(candidates / inflation)),
    candidates
  )[1L]
  fit <- fits[[selected]]
  fit$target_search <- data.frame(
    inflation = candidates,
    achieved_k = achieved,
    converged = converged,
    iterations = vapply(fits, `[[`, integer(1), "iterations"),
    nnz_final = vapply(fits, function(x) length(x$flow@x), integer(1))
  )
  fit$selected_inflation <- candidates[selected]
  fit
}

#' Cluster4d using sparse Markov Clustering (MCL)
#'
#' Builds a sparse voxel graph, runs an MCL flow process, and maps attractors to
#' cluster labels. The implementation is optimized for sparse neuroimaging graphs.
#'
#' @param vec A NeuroVec, SparseNeuroVec, or NeuroVol containing the features.
#' @param mask A NeuroVol mask; finite values strictly greater than zero are included.
#' @param n_clusters Requested cluster count. With `exact_k = FALSE`, this is
#'   the target used to select the closest natural MCL partition across a fixed,
#'   deterministic inflation search. With `exact_k = TRUE`, it is an exact
#'   adjacency-preserving postcondition.
#' @param spatial_weight Convex weight from zero through one for physical edge similarity.
#' @param max_iterations Positive MCL iteration limit for each inflation candidate.
#' @param connectivity Exact masked grid connectivity: 6, 18, or 26.
#' @param verbose Logical; report iteration diagnostics.
#' @param inflation Nominal MCL inflation parameter (>1). In natural target mode,
#'   the deterministic search evaluates this value and fixed multipliers around it.
#' @param expansion MCL expansion power (integer >= 2). Default 2.
#' @param loop_weight Finite non-negative self-loop weight added before normalization.
#' @param prune_k Maximum nonzero entries kept per column during MCL iterations.
#'   Must be a positive integer; the default is 64 and is independent of the target K.
#' @param prune_threshold Finite threshold in [0, 1). The strongest entry in
#'   every column is retained even when all entries fall below the threshold.
#' @param tol Positive convergence tolerance on maximum absolute matrix delta.
#' @param feature_metric Feature similarity metric: `"correlation"` or `"euclidean"`.
#' @param feature_sigma Optional bandwidth for euclidean heat-kernel feature similarity.
#' @param spatial_sigma Optional bandwidth for spatial heat-kernel similarity.
#' @param exact_k If TRUE, use the shared adjacency-preserving exact-K engine.
#'
#' @return A `cluster4d_result`. Metadata records natural and final K,
#'   convergence, selected inflation, candidate search results, flow sparsity,
#'   final normalization error, and exact-K policy.
#' @details Every expansion, inflation, and pruning step is checked to remain
#'   finite, non-negative, and column stochastic. The `parallel` argument was
#'   removed because this sparse Matrix implementation has no parallel backend.
#'   Supplying it is therefore an ordinary unused-argument error rather than a
#'   silently sequential execution.
#' @export
cluster4d_mcl <- function(vec, mask, n_clusters = 100,
                          spatial_weight = 0.2,
                          max_iterations = 8,
                          connectivity = 6,
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
                          exact_k = FALSE) {

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_mcl")

  if (!is.numeric(spatial_weight) || length(spatial_weight) != 1 ||
      !is.finite(spatial_weight) || spatial_weight < 0 || spatial_weight > 1) {
    stop("cluster4d_mcl: spatial_weight must be in [0, 1]")
  }

  feature_metric <- match.arg(feature_metric)
  exact_k <- .cluster4d_scalar_logical(
    exact_k, "exact_k", "cluster4d_mcl"
  )
  inflation <- .cluster4d_scalar_number(
    inflation, "inflation", "cluster4d_mcl"
  )
  if (inflation <= 1) {
    stop("cluster4d_mcl: inflation must be greater than 1", call. = FALSE)
  }
  expansion <- .cluster4d_scalar_number(
    expansion, "expansion", "cluster4d_mcl", lower = 2, integer = TRUE
  )
  loop_weight <- .cluster4d_scalar_number(
    loop_weight, "loop_weight", "cluster4d_mcl", lower = 0
  )
  prune_threshold <- .cluster4d_scalar_number(
    prune_threshold, "prune_threshold", "cluster4d_mcl",
    lower = 0, upper = 1
  )
  if (prune_threshold >= 1) {
    stop("cluster4d_mcl: prune_threshold must be less than 1", call. = FALSE)
  }
  tol <- .cluster4d_scalar_number(
    tol, "tol", "cluster4d_mcl", lower = 0
  )
  if (tol <= 0) stop("cluster4d_mcl: tol must be positive", call. = FALSE)
  max_iterations <- .cluster4d_scalar_number(
    max_iterations, "max_iterations", "cluster4d_mcl",
    lower = 1, integer = TRUE
  )
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "cluster4d_mcl", integer = TRUE
  )
  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("cluster4d_mcl: connectivity must be 6, 18, or 26", call. = FALSE)
  }
  for (entry in c("feature_sigma", "spatial_sigma")) {
    value <- get(entry)
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0)) {
      stop("cluster4d_mcl: ", entry,
           " must be NULL or a positive finite scalar", call. = FALSE)
    }
  }

  data_prep <- prepare_cluster4d_data(
    vec = vec,
    mask = mask,
    scale_features = FALSE,
    scale_coords = FALSE
  )

  if (is.null(prune_k)) {
    prune_k <- 64L
  } else {
    prune_k <- .cluster4d_scalar_number(
      prune_k, "prune_k", "cluster4d_mcl", lower = 1, integer = TRUE
    )
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

  fit_args <- list(
    graph = graph,
    expansion = expansion,
    max_iter = max_iterations,
    tol = tol,
    prune_k = prune_k,
    prune_threshold = prune_threshold,
    loop_weight = loop_weight,
    verbose = verbose
  )
  mcl_fit <- if (exact_k) {
    do.call(.mcl_sparse, c(fit_args, list(inflation = inflation)))
  } else {
    do.call(.mcl_targeted_sparse, c(
      fit_args,
      list(target_k = as.integer(n_clusters), inflation = inflation)
    ))
  }

  labels <- mcl_fit$labels
  natural_k <- as.integer(max(labels))

  if (exact_k) {
    labels <- force_exact_k(
      labels, data_prep$features, n_clusters,
      mask = mask, connectivity = connectivity
    )
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
      inflation = inflation,
      expansion = as.integer(expansion),
      loop_weight = loop_weight,
      prune_k = prune_k,
      prune_threshold = prune_threshold,
      tol = tol,
      feature_metric = feature_metric,
      feature_sigma = feature_sigma,
      spatial_sigma = spatial_sigma,
      exact_k = exact_k
    ),
    metadata = list(
      converged = mcl_fit$converged,
      iterations = mcl_fit$iterations,
      max_delta = mcl_fit$max_delta,
      nnz_history = mcl_fit$nnz_history,
      nnz_final = length(mcl_fit$flow@x),
      flow_column_error = max(
        abs(Matrix::colSums(mcl_fit$flow) - 1), 0
      ),
      natural_k = natural_k,
      achieved_k = as.integer(max(labels)),
      selected_inflation = if (exact_k) inflation else mcl_fit$selected_inflation,
      target_search = if (exact_k) NULL else mcl_fit$target_search,
      target_policy = if (exact_k) {
        "adjacency_preserving_exact_k"
      } else {
        "closest_natural_k_over_deterministic_inflation_grid"
      },
      exact_k_applied = exact_k
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  finalize_cluster4d_result(result, vec, mask, "mcl", result$parameters)
}
