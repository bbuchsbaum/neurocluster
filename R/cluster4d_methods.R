#' S3 Methods for cluster4d_result Objects
#'
#' Methods for printing, summarizing, and plotting cluster4d results.

#' Print cluster4d result
#'
#' @param x A cluster4d_result object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#' @export
print.cluster4d_result <- function(x, ...) {
  cat("Cluster4D Result\n")
  cat("================\n")
  cat("Method:", x$method, "\n")
  cat("Number of clusters:", x$n_clusters, "\n")
  
  if (!is.null(x$parameters)) {
    cat("Requested clusters:", x$parameters$n_clusters_requested, "\n")
    if (!is.null(x$parameters$spatial_weight)) {
      cat("Spatial weight:", x$parameters$spatial_weight, "\n")
    }
  }
  
  # Get dimensions from clusvol
  dims <- dim(x$clusvol)
  cat("Volume dimensions:", paste(dims, collapse=" x "), "\n")
  
  # Count voxels
  n_voxels <- length(x$cluster)
  cat("Masked voxels:", n_voxels, "\n")
  
  # Show cluster sizes
  cluster_sizes <- table(x$cluster)
  cat("\nCluster sizes:\n")
  cat("  Min:", min(cluster_sizes), "voxels\n")
  cat("  Max:", max(cluster_sizes), "voxels\n")
  cat("  Mean:", round(mean(cluster_sizes), 1), "voxels\n")
  cat("  Median:", median(cluster_sizes), "voxels\n")
  
  invisible(x)
}

#' Summarize cluster4d result
#'
#' @param object A cluster4d_result object
#' @param ... Additional arguments (ignored)
#'
#' @return A summary list (invisibly)
#' @method summary cluster4d_result
#' @exportS3Method summary cluster4d_result
summary.cluster4d_result <- function(object, ...) {
  # Basic info
  cat("Cluster4D Analysis Summary\n")
  cat("==========================\n\n")
  
  # Method and parameters
  cat("Method:", object$method, "\n")
  
  if (!is.null(object$parameters)) {
    cat("\nParameters:\n")
    params <- object$parameters
    param_names <- names(params)
    
    # Show key parameters
    key_params <- c("n_clusters_requested", "spatial_weight", "max_iterations", 
                   "connectivity", "parallel")
    for (p in key_params) {
      if (p %in% param_names) {
        cat(sprintf("  %-20s: %s\n", p, 
                   ifelse(is.logical(params[[p]]), 
                         ifelse(params[[p]], "TRUE", "FALSE"),
                         as.character(params[[p]]))))
      }
    }
    
    # Show method-specific parameters
    other_params <- setdiff(param_names, key_params)
    if (length(other_params) > 0) {
      cat("\nMethod-specific parameters:\n")
      for (p in other_params) {
        cat(sprintf("  %-20s: %s\n", p, as.character(params[[p]])))
      }
    }
  }
  
  # Clustering results
  cat("\nClustering Results:\n")
  cat("  Clusters found:", object$n_clusters, "\n")
  
  # Cluster statistics
  cluster_sizes <- table(object$cluster)
  cat("\nCluster Size Statistics:\n")
  cat("  Minimum:", min(cluster_sizes), "voxels\n")
  cat("  1st Quartile:", quantile(cluster_sizes, 0.25), "voxels\n")
  cat("  Median:", median(cluster_sizes), "voxels\n")
  cat("  Mean:", round(mean(cluster_sizes), 1), "voxels\n")
  cat("  3rd Quartile:", quantile(cluster_sizes, 0.75), "voxels\n")
  cat("  Maximum:", max(cluster_sizes), "voxels\n")
  cat("  Std Dev:", round(sd(cluster_sizes), 1), "voxels\n")
  
  # Spatial extent
  if (!is.null(object$coord_centers)) {
    cat("\nSpatial Extent of Clusters:\n")
    ranges <- apply(object$coord_centers, 2, range)
    cat("  X range:", sprintf("%.1f to %.1f mm\n", ranges[1,1], ranges[2,1]))
    cat("  Y range:", sprintf("%.1f to %.1f mm\n", ranges[1,2], ranges[2,2]))
    cat("  Z range:", sprintf("%.1f to %.1f mm\n", ranges[1,3], ranges[2,3]))
  }
  
  # Feature space info
  if (!is.null(object$centers)) {
    cat("\nFeature Space:\n")
    cat("  Time points:", ncol(object$centers), "\n")
    
    # Calculate average within-cluster correlation if possible
    if (ncol(object$centers) > 1) {
      center_cors <- cor(t(object$centers))
      diag(center_cors) <- NA
      mean_cor <- mean(center_cors, na.rm = TRUE)
      cat("  Mean between-cluster correlation:", sprintf("%.3f\n", mean_cor))
    }
  }
  
  # Metadata
  if (!is.null(object$metadata) && length(object$metadata) > 0) {
    cat("\nAdditional Information:\n")
    for (name in names(object$metadata)) {
      if (name != "centers" && name != "coord_centers") {  # Skip redundant info
        val <- object$metadata[[name]]
        if (is.numeric(val) && length(val) == 1) {
          cat(sprintf("  %s: %.3f\n", name, val))
        } else if (is.character(val) && length(val) == 1) {
          cat(sprintf("  %s: %s\n", name, val))
        }
      }
    }
  }
  
  # Return summary invisibly
  invisible(list(
    method = object$method,
    n_clusters = object$n_clusters,
    cluster_sizes = as.numeric(cluster_sizes),
    parameters = object$parameters
  ))
}

#' Plot cluster4d result
#'
#' Creates visualization of clustering results. Shows axial, sagittal, and
#' coronal slices through the clustered volume.
#'
#' @param x A cluster4d_result object
#' @param slice Slice specification. Can be:
#'   \itemize{
#'     \item NULL (default): Shows middle slices
#'     \item Numeric vector of length 3: c(x, y, z) coordinates
#'     \item "montage": Shows multiple slices
#'   }
#' @param view Viewing plane: "axial", "sagittal", "coronal", or "all"
#' @param colors Color palette for clusters. Default uses rainbow colors.
#' @param zlim Numeric length-two label range shared by every displayed plane.
#'   The default, `c(0.5, x$n_clusters + 0.5)`, keeps each positive integer
#'   label mapped to the same categorical color even when a slice omits labels.
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the plotted data
#' @method plot cluster4d_result
#' @export
plot.cluster4d_result <- function(x, slice = NULL, view = "all", 
                                 colors = NULL, zlim = NULL, ...) {
  
  # Get the clustered volume and materialize label array
  clusvol <- x$clusvol
  dims <- dim(clusvol)
  mask <- clusvol@mask
  mask_idx <- which(mask > 0)
  label_arr <- array(NA_integer_, dim = dims)
  label_arr[mask_idx] <- clusvol@clusters
  
  # Determine slices to show
  if (is.null(slice)) {
    # Show middle slices
    slice_x <- dims[1] %/% 2
    slice_y <- dims[2] %/% 2
    slice_z <- dims[3] %/% 2
  } else if (is.numeric(slice) && length(slice) == 3) {
    slice_x <- slice[1]
    slice_y <- slice[2]
    slice_z <- slice[3]
  } else if (slice == "montage") {
    # Will implement montage view
    stop("Montage view not yet implemented")
  } else {
    stop("Invalid slice specification")
  }
  
  # Generate colors if not provided
  if (is.null(colors)) {
    n_colors <- x$n_clusters
    colors <- rainbow(n_colors)
  }
  if (is.null(zlim)) {
    zlim <- c(0.5, x$n_clusters + 0.5)
  }
  
  # Save and restore graphics parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Setup plot layout
  if (view == "all") {
    par(mfrow = c(1, 3))
    views_to_plot <- c("axial", "sagittal", "coronal")
  } else {
    par(mfrow = c(1, 1))
    views_to_plot <- view
  }
  
  # Plot each view
  for (v in views_to_plot) {
    if (v == "axial") {
      # Extract axial slice (x-y plane at fixed z)
      slice_data <- label_arr[, , slice_z, drop = FALSE][, , 1]
      main_title <- paste("Axial slice z =", slice_z)
      xlab <- "X"
      ylab <- "Y"
    } else if (v == "sagittal") {
      # Extract sagittal slice (y-z plane at fixed x)
      slice_data <- label_arr[slice_x, , , drop = FALSE][1, , ]
      main_title <- paste("Sagittal slice x =", slice_x)
      xlab <- "Y"
      ylab <- "Z"
    } else if (v == "coronal") {
      # Extract coronal slice (x-z plane at fixed y)
      slice_data <- label_arr[, slice_y, , drop = FALSE][, 1, ]
      main_title <- paste("Coronal slice y =", slice_y)
      xlab <- "X"
      ylab <- "Z"
    }
    
    # Create color-mapped image
    image(slice_data, col = colors, zlim = zlim, main = main_title,
          xlab = xlab, ylab = ylab, axes = FALSE, ...)
    
    # Add axes
    axis(1)
    axis(2)
    box()
  }
  
  # Reset par
  par(mfrow = c(1, 1))
  
  invisible(list(
    slice_x = slice_x,
    slice_y = slice_y,
    slice_z = slice_z,
    zlim = zlim
  ))
}

#' Compare multiple cluster4d results
#'
#' Compares clustering results over exactly the same included mask voxels and
#' physical coordinate system. Ambiguous center-based scores are not used.
#'
#' @param ... cluster4d_result objects to compare
#' @param metrics Character vector selecting explicit estimands:
#'   \itemize{
#'     \item `"summary"`: cluster-count and cluster-size summaries.
#'     \item `"spatial_dispersion"`: root mean squared physical distance in
#'       millimetres from voxels to their final cluster centroid (lower is more
#'       compact).
#'     \item `"temporal_coherence"`: mean Pearson correlation across all
#'       unordered within-cluster voxel pairs, requiring `feature_mat`.
#'     \item `"partition_agreement"`: adjusted Rand index, variation of
#'       information in bits, and pairwise Dice for exactly two results.
#'   }
#' @param feature_mat A shared finite numeric matrix with voxels in rows and
#'   at least two time/features in columns. Required only for
#'   `"temporal_coherence"`. Rows must be non-constant.
#'
#' @details Every result must carry canonical mask and affine provenance. The
#'   exact included voxel indices and physical geometry are checked before any
#'   cross-result comparison. Partition metrics use only observed contingency
#'   cells and never construct an N-by-N co-membership matrix. Pairwise
#'   agreement columns are repeated on the two result rows so the return value
#'   remains a data frame.
#'
#' @return A data frame with one row per result and explicitly named metrics.
#' @export
compare_cluster4d <- function(..., 
                             metrics = c("summary", "spatial_dispersion"),
                             feature_mat = NULL) {
  results <- list(...)
  if (length(results) == 1L && is.list(results[[1L]]) &&
      !inherits(results[[1L]], "cluster_result") &&
      length(results[[1L]]) > 0L &&
      all(vapply(
        results[[1L]],
        function(x) inherits(x, "cluster_result") ||
          inherits(x, "cluster4d_result"),
        logical(1)
      ))) {
    results <- results[[1L]]
  }
  n_results <- length(results)

  if (n_results < 1L) {
    stop("At least one cluster4d_result is required", call. = FALSE)
  }

  if (!is.character(metrics) || anyNA(metrics) || !length(metrics)) {
    stop("metrics must be a non-empty character vector", call. = FALSE)
  }
  obsolete <- intersect(metrics, c("spatial_coherence", "overlap"))
  if (length(obsolete)) {
    stop(
      paste(obsolete, collapse = ", "),
      " is no longer supported; use spatial_dispersion and/or ",
      "partition_agreement for explicitly defined estimands",
      call. = FALSE
    )
  }
  allowed <- c(
    "summary", "spatial_dispersion", "temporal_coherence",
    "partition_agreement"
  )
  unknown <- setdiff(metrics, allowed)
  if (length(unknown)) {
    stop("Unknown metrics: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  metrics <- unique(metrics)
  if ("partition_agreement" %in% metrics && n_results != 2L) {
    stop("partition_agreement requires exactly two results", call. = FALSE)
  }
  if ("temporal_coherence" %in% metrics && is.null(feature_mat)) {
    stop("feature_mat is required for temporal_coherence", call. = FALSE)
  }

  supports <- lapply(seq_along(results), function(i) {
    .cluster4d_result_support(results[[i]], paste0("result ", i))
  })
  if (n_results > 1L) {
    reference <- supports[[1L]]
    for (i in 2:n_results) {
      if (!.cluster4d_geometry_equal(
        reference$geometry, supports[[i]]$geometry
      )) {
        stop(
          "All results must have compatible physical coordinate provenance",
          call. = FALSE
        )
      }
      if (!identical(reference$mask_idx, supports[[i]]$mask_idx)) {
        stop(
          "All results must use the same included mask voxels and ordering",
          call. = FALSE
        )
      }
    }
  }

  result_names <- names(results)
  if (is.null(result_names)) result_names <- rep("", n_results)
  missing_names <- is.na(result_names) | !nzchar(result_names)
  if (any(missing_names)) {
    method_names <- vapply(seq_along(results), function(i) {
      method <- results[[i]]$method
      if (is.character(method) && length(method) == 1L && nzchar(method)) {
        method
      } else {
        paste0("Result", i)
      }
    }, character(1))
    result_names[missing_names] <- method_names[missing_names]
    if (anyDuplicated(result_names)) {
      result_names <- make.unique(result_names, sep = "_")
    }
  }

  comparison <- data.frame(
    Method = result_names,
    stringsAsFactors = FALSE
  )

  if ("summary" %in% metrics) {
    sizes <- lapply(supports, function(x) {
      tabulate(x$labels, nbins = max(x$labels))
    })
    comparison$N_Clusters <- vapply(sizes, length, integer(1))
    comparison$Min_Size <- vapply(sizes, min, numeric(1))
    comparison$Max_Size <- vapply(sizes, max, numeric(1))
    comparison$Mean_Size <- vapply(sizes, mean, numeric(1))
    comparison$SD_Size <- vapply(sizes, stats::sd, numeric(1))
  }

  if ("spatial_dispersion" %in% metrics) {
    comparison$Spatial_RMS_Distance_mm <- vapply(
      supports,
      function(x) .spatial_rms_dispersion(x$labels, x$coords),
      numeric(1)
    )
  }

  if ("temporal_coherence" %in% metrics) {
    comparison$Temporal_Pairwise_Correlation <- vapply(
      supports,
      function(x) .temporal_pairwise_correlation(x$labels, feature_mat),
      numeric(1)
    )
  }

  if ("partition_agreement" %in% metrics) {
    agreement <- .partition_metrics(
      supports[[1L]]$labels, supports[[2L]]$labels
    )
    comparison$Adjusted_Rand_Index <- rep(agreement$ari, 2L)
    comparison$Variation_of_Information_bits <- rep(
      agreement$variation_of_information_bits, 2L
    )
    comparison$Pairwise_Dice <- rep(agreement$pairwise_dice, 2L)
  }

  comparison
}

#' Validate cluster4d result
#'
#' Checks validity and quality of clustering results.
#'
#' @param result A cluster4d_result object
#' @param vec Original NeuroVec data, required to verify feature centers.
#' @param mask Original mask, required to verify geometry, voxel order, and
#'   physical coordinate centers.
#' @param tolerance Relative numeric tolerance used when independently
#'   recomputing centers and coordinate centers.
#'
#' @return A list with `valid`, `errors`, `warnings`, and a compact `summary`.
#'   Any schema or consistency defect sets `valid = FALSE`.
#' @export
validate_cluster4d <- function(result, vec = NULL, mask = NULL,
                               tolerance = 1e-8) {
  validation <- list(
    valid = TRUE,
    warnings = character(),
    errors = character()
  )

  add_error <- function(message) {
    validation$valid <<- FALSE
    validation$errors <<- unique(c(validation$errors, message))
  }
  add_warning <- function(message) {
    validation$warnings <<- unique(c(validation$warnings, message))
  }
  scalar_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
  }
  matrix_equal <- function(x, y) {
    if (!is.matrix(x) || !is.matrix(y) || !identical(dim(x), dim(y)) ||
        any(!is.finite(x)) || any(!is.finite(y))) return(FALSE)
    scale <- max(1, abs(y))
    max(abs(x - y)) <= tolerance * scale
  }
  geometry_equal <- function(x, y) {
    is.list(x) && is.list(y) &&
      identical(as.integer(x$dimensions), as.integer(y$dimensions)) &&
      .cluster4d_same_numeric(as.numeric(x$spacing), as.numeric(y$spacing)) &&
      .cluster4d_same_numeric(as.matrix(x$affine), as.matrix(y$affine))
  }

  if (!inherits(result, "cluster_result") && !inherits(result, "cluster4d_result")) {
    add_error("Not a valid cluster4d_result object")
  }
  if (!is.list(result)) {
    add_error("Result must be a list-like object")
    validation$summary <- list(
      n_clusters = NA_integer_, n_voxels = NA_integer_,
      cluster_size_range = c(NA_integer_, NA_integer_),
      has_centers = FALSE, has_spatial_centers = FALSE
    )
    return(validation)
  }

  required <- c(
    "labels", "cluster", "clusvol", "centers", "coord_centers",
    "actual_k", "n_clusters", "label_ids", "method", "parameters",
    "provenance"
  )
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    add_error(paste("Missing required fields:", paste(missing, collapse = ", ")))
  }

  labels <- result$labels
  labels_valid <- is.numeric(labels) && length(labels) > 0L &&
    all(is.finite(labels)) && all(labels == floor(labels)) && all(labels > 0)
  if (!labels_valid) {
    add_error("labels must be a non-empty vector of finite positive integers")
    label_ids <- integer()
  } else {
    labels <- as.integer(labels)
    label_ids <- sort(unique(labels))
    if (!identical(label_ids, seq_len(length(label_ids)))) {
      add_error("labels must use contiguous positive IDs starting at 1")
    }
    cluster_sizes <- table(labels)
    tiny_clusters <- sum(cluster_sizes < 5)
    if (tiny_clusters > 0) add_warning(
      paste(tiny_clusters, "clusters have fewer than 5 voxels")
    )
  }

  if (!identical(result$labels, result$cluster)) {
    add_error("cluster must be identical to labels")
  }

  actual_k_valid <- scalar_integer(result$actual_k) && result$actual_k >= 1
  if (!actual_k_valid) {
    add_error("actual_k must be a positive integer scalar")
    actual_k <- NA_integer_
  } else {
    actual_k <- as.integer(result$actual_k)
    if (labels_valid && actual_k != length(label_ids)) {
      add_error("actual_k does not match the final labels")
    }
  }

  if (!scalar_integer(result$n_clusters) ||
      !actual_k_valid || as.integer(result$n_clusters) != actual_k) {
    add_error("n_clusters must be identical to actual_k")
  }
  if (!is.integer(result$label_ids) || !actual_k_valid ||
      !identical(result$label_ids, seq_len(actual_k))) {
    add_error("label_ids must map center rows to labels 1 through actual_k")
  }
  if (!is.character(result$method) || length(result$method) != 1L ||
      is.na(result$method) || !nzchar(result$method)) {
    add_error("method must be a non-empty character scalar")
  }
  if (!is.list(result$parameters)) add_error("parameters must be a list")

  centers_valid <- is.matrix(result$centers) && is.numeric(result$centers) &&
    all(is.finite(result$centers))
  if (!centers_valid) add_error("centers must be a finite numeric matrix")
  coord_centers_valid <- is.matrix(result$coord_centers) &&
    is.numeric(result$coord_centers) && all(is.finite(result$coord_centers))
  if (!coord_centers_valid) {
    add_error("coord_centers must be a finite numeric matrix")
  }
  if (actual_k_valid && centers_valid && nrow(result$centers) != actual_k) {
    add_error("centers must have actual_k rows")
  }
  if (actual_k_valid && coord_centers_valid &&
      !identical(dim(result$coord_centers), c(actual_k, 3L))) {
    add_error("coord_centers must have shape actual_k by 3")
  }

  provenance <- result$provenance
  provenance_valid <- is.list(provenance) &&
    is.list(provenance$label_space) &&
    is.list(provenance$feature_space) &&
    is.list(provenance$coordinate_space) &&
    is.list(provenance$mask)
  if (!provenance_valid) {
    add_error(
      "provenance must define label, feature, coordinate, and mask spaces"
    )
  } else {
    if (!identical(provenance$label_space$row_to_label, result$label_ids)) {
      add_error("provenance row_to_label does not match label_ids")
    }
    if (!identical(provenance$feature_space$representation, "original") ||
        !identical(provenance$feature_space$summary, "mean")) {
      add_error("provenance must type centers as means in original feature space")
    }
    if (!identical(provenance$coordinate_space$units, "mm") ||
        !identical(provenance$coordinate_space$summary, "mean")) {
      add_error("provenance must type coordinate centers as means in mm")
    }
  }

  clusvol_valid <- inherits(result$clusvol, "ClusteredNeuroVol")
  if (!clusvol_valid) {
    add_error("clusvol must be a ClusteredNeuroVol")
  } else {
    volume_labels <- tryCatch(
      as.integer(result$clusvol@clusters),
      error = function(e) NULL
    )
    if (!labels_valid || !identical(volume_labels, labels)) {
      add_error("clusvol labels do not round-trip to labels")
    }
  }

  if (xor(is.null(vec), is.null(mask))) {
    add_error("vec and mask must be supplied together for external validation")
  }
  if (!is.null(vec) && !is.null(mask)) {
    data <- tryCatch(
      .cluster4d_original_data(vec, mask, "validate_cluster4d"),
      error = function(e) {
        add_error(conditionMessage(e))
        NULL
      }
    )
    if (!is.null(data)) {
      if (!labels_valid || length(labels) != nrow(data$features)) {
        add_error("labels length does not match the declared mask voxel set")
      }
      if (centers_valid && ncol(result$centers) != ncol(data$features)) {
        add_error("centers must have one column per original feature")
      }
      if (provenance_valid) {
        if (!identical(
          provenance$feature_space$n_features,
          as.integer(ncol(data$features))
        )) {
          add_error("feature provenance has the wrong feature count")
        }
        if (!geometry_equal(
          provenance$coordinate_space$geometry, data$geometry
        )) {
          add_error("coordinate provenance geometry does not match mask")
        }
      }
      if (clusvol_valid) {
        volume_geometry <- tryCatch(
          .cluster4d_geometry(result$clusvol), error = function(e) NULL
        )
        if (is.null(volume_geometry) ||
            !geometry_equal(volume_geometry, data$geometry)) {
          add_error("clusvol geometry does not match mask")
        }
        volume_mask_idx <- tryCatch(
          which(as.array(result$clusvol@mask) > 0), error = function(e) NULL
        )
        if (!identical(volume_mask_idx, data$mask_idx)) {
          add_error("clusvol mask does not match the declared mask voxel set")
        }
      }
      if (labels_valid && length(labels) == nrow(data$features)) {
        oracle <- compute_cluster_centers(
          labels, data$features, data$coords, method = "mean"
        )
        if (!matrix_equal(result$centers, unname(as.matrix(oracle$centers)))) {
          add_error("centers are stale or do not summarize final labels")
        }
        if (!matrix_equal(
          result$coord_centers, unname(as.matrix(oracle$coord_centers))
        )) {
          add_error(
            "coord_centers are stale or do not summarize final labels in mm"
          )
        }
      }
    }
  } else {
    add_error(
      "Original vec and mask are required to validate centers and geometry"
    )
  }

  cluster_size_range <- if (labels_valid) {
    as.integer(range(table(labels)))
  } else c(NA_integer_, NA_integer_)
  validation$summary <- list(
    n_clusters = if (actual_k_valid) actual_k else NA_integer_,
    n_voxels = if (is.atomic(result$labels)) length(result$labels) else NA_integer_,
    cluster_size_range = cluster_size_range,
    has_centers = centers_valid,
    has_spatial_centers = coord_centers_valid
  )

  validation
}
