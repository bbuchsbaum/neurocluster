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
#' @export
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
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the plotted data
#' @method plot cluster4d_result
#' @export
plot.cluster4d_result <- function(x, slice = NULL, view = "all", 
                                 colors = NULL, ...) {
  
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
    image(slice_data, col = colors, main = main_title,
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
    slice_z = slice_z
  ))
}

#' Compare multiple cluster4d results
#'
#' Compares clustering results from different methods or parameters.
#'
#' @param ... cluster4d_result objects to compare
#' @param metrics Comparison metrics to compute. Options:
#'   \itemize{
#'     \item "summary": Basic statistics
#'     \item "spatial_coherence": Spatial compactness measure
#'     \item "temporal_coherence": Feature similarity within clusters
#'     \item "overlap": Dice coefficient between methods (requires exactly 2 results)
#'   }
#'
#' @return A comparison data frame
#' @export
compare_cluster4d <- function(..., 
                             metrics = c("summary", "spatial_coherence", 
                                       "temporal_coherence")) {
  
  results <- list(...)
  n_results <- length(results)
  
  if (n_results < 1) {
    stop("At least one cluster4d_result required")
  }
  
  # Check all are cluster4d_results
  if (!all(sapply(results, function(x) inherits(x, "cluster_result") || 
                                       inherits(x, "cluster4d_result")))) {
    stop("All arguments must be cluster4d_result or cluster_result objects")
  }
  
  # Extract names or create them
  result_names <- names(results)
  if (is.null(result_names)) {
    result_names <- sapply(results, function(x) x$method)
    if (any(is.null(result_names))) {
      result_names <- paste0("Result", 1:n_results)
    }
  }
  
  # Initialize comparison data frame
  comparison <- data.frame(
    Method = result_names,
    stringsAsFactors = FALSE
  )
  
  # Add basic summary metrics
  if ("summary" %in% metrics) {
    comparison$N_Clusters <- sapply(results, function(x) {
      if (!is.null(x$n_clusters)) x$n_clusters
      else length(unique(x$cluster[!is.na(x$cluster)]))
    })
    
    comparison$Min_Size <- sapply(results, function(x) {
      min(table(x$cluster))
    })
    
    comparison$Max_Size <- sapply(results, function(x) {
      max(table(x$cluster))
    })
    
    comparison$Mean_Size <- sapply(results, function(x) {
      round(mean(table(x$cluster)), 1)
    })
    
    comparison$SD_Size <- sapply(results, function(x) {
      round(sd(table(x$cluster)), 1)
    })
  }
  
  # Spatial coherence
  if ("spatial_coherence" %in% metrics) {
    comparison$Spatial_Coherence <- sapply(results, function(x) {
      if (!is.null(x$coord_centers)) {
        # Calculate average within-cluster spatial variance
        cluster_vars <- numeric()
        for (i in unique(x$cluster)) {
          if (!is.na(i) && i > 0) {
            mask <- x$cluster == i
            if (sum(mask) > 1) {
              # Would need coordinates of all voxels, not just centers
              # For now, return NA
              cluster_vars <- c(cluster_vars, NA)
            }
          }
        }
        if (length(cluster_vars) > 0 && !all(is.na(cluster_vars))) {
          round(mean(cluster_vars, na.rm = TRUE), 2)
        } else {
          NA
        }
      } else {
        NA
      }
    })
  }
  
  # Temporal coherence
  if ("temporal_coherence" %in% metrics) {
    comparison$Temporal_Coherence <- sapply(results, function(x) {
      if (!is.null(x$centers) && ncol(x$centers) > 1) {
        # Calculate mean within-cluster correlation
        # This is approximate using centers
        center_cors <- cor(t(x$centers))
        diag(center_cors) <- NA
        round(mean(center_cors, na.rm = TRUE), 3)
      } else {
        NA
      }
    })
  }
  
  # Overlap metrics (only for exactly 2 results)
  if ("overlap" %in% metrics && n_results == 2) {
    # Calculate Dice coefficient
    clus1 <- results[[1]]$cluster
    clus2 <- results[[2]]$cluster
    
    if (length(clus1) == length(clus2)) {
      # Simple overlap: what fraction of voxel pairs are in same/different clusters
      same1 <- outer(clus1, clus1, "==")
      same2 <- outer(clus2, clus2, "==")
      overlap <- sum(same1 == same2) / length(same1)
      comparison$Overlap <- c(overlap, overlap)
    }
  }
  
  comparison
}

#' Validate cluster4d result
#'
#' Checks validity and quality of clustering results.
#'
#' @param result A cluster4d_result object
#' @param vec Original NeuroVec data
#' @param mask Original mask
#'
#' @return A list with validation results
#' @export
validate_cluster4d <- function(result, vec = NULL, mask = NULL) {
  
  validation <- list(
    valid = TRUE,
    warnings = character(),
    errors = character()
  )
  
  # Check structure
  if (!inherits(result, "cluster_result") && !inherits(result, "cluster4d_result")) {
    validation$valid <- FALSE
    validation$errors <- c(validation$errors, "Not a valid cluster result object")
    return(validation)
  }
  
  # Check required components
  required <- c("clusvol", "cluster", "centers", "coord_centers")
  missing <- setdiff(required, names(result))
  if (length(missing) > 0) {
    validation$warnings <- c(validation$warnings, 
                           paste("Missing components:", paste(missing, collapse = ", ")))
  }
  
  # Check cluster assignments
  if (!is.null(result$cluster)) {
    clusters <- result$cluster[!is.na(result$cluster)]
    
    # Check for empty clusters
    unique_clusters <- sort(unique(clusters))
    expected_clusters <- 1:max(unique_clusters)
    missing_clusters <- setdiff(expected_clusters, unique_clusters)
    
    if (length(missing_clusters) > 0) {
      validation$warnings <- c(validation$warnings,
                             paste("Missing cluster IDs:", 
                                   paste(missing_clusters, collapse = ", ")))
    }
    
    # Check for very small clusters
    cluster_sizes <- table(clusters)
    tiny_clusters <- sum(cluster_sizes < 5)
    if (tiny_clusters > 0) {
      validation$warnings <- c(validation$warnings,
                             paste(tiny_clusters, "clusters have fewer than 5 voxels"))
    }
  }
  
  # Check centers match cluster count
  if (!is.null(result$centers) && !is.null(result$n_clusters)) {
    if (nrow(result$centers) != result$n_clusters) {
      validation$warnings <- c(validation$warnings,
                             "Number of centers doesn't match n_clusters")
    }
  }
  
  # If original data provided, check consistency
  if (!is.null(vec) && !is.null(mask)) {
    mask_voxels <- sum(mask > 0)
    if (length(result$cluster) != mask_voxels) {
      validation$errors <- c(validation$errors,
                           paste("Cluster assignments (", length(result$cluster), 
                                 ") don't match mask voxels (", mask_voxels, ")", sep = ""))
      validation$valid <- FALSE
    }
  }
  
  # Check for NA values
  na_count <- sum(is.na(result$cluster))
  if (na_count > 0) {
    validation$warnings <- c(validation$warnings,
                           paste(na_count, "voxels have NA cluster assignment"))
  }
  
  # Summary
  validation$summary <- list(
    n_clusters = result$n_clusters,
    n_voxels = length(result$cluster),
    cluster_size_range = range(table(result$cluster)),
    has_centers = !is.null(result$centers),
    has_spatial_centers = !is.null(result$coord_centers)
  )
  
  validation
}
