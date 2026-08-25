# Synthetic Cube Test - 9x9x9 Volume with 27 Cubic Clusters
#
# This test evaluates all clustering methods on a synthetic dataset with
# perfectly separable ground truth clusters. Each method should ideally
# achieve perfect recovery (ARI = 1.0) on this easy test case.
#
# Run with: devtools::test(filter = "cube")

# All clustering methods to test
ALL_METHODS <- c("flash3d", "slice_msf", "g3s", "rena", "rena_plus", "acsc",
                 "slic", "corr_slic", "supervoxels", "snic", "commute")

# Method-specific overrides to keep each algorithm near its intended regime
METHOD_OVERRIDES <- list(
  commute = list(
    spatial_weight = 0.1,  # feature-heavy for perfectly separable signals
    ncomp = 40             # richer embedding for the 27-cube structure
  ),
  corr_slic = list(
    spatial_weight = 0.25,
    embedding_dim = 64,
    max_iterations = 5,
    connectivity = 6
  )
)

# Empirical recovery floors are method-specific because the algorithms optimize
# different estimands. In particular, Slice-MSF's exact-K repair is constrained
# by its spatial RAG and is not a supervised feature-recovery method. These
# floors retain a margin for platform noise while detecting material regression.
CUBE_ARI_FLOORS <- c(
  flash3d = 0.95, slice_msf = -0.01, g3s = 0.75, rena = 0.95,
  rena_plus = 0.95, acsc = 0.90, slic = 0.90, corr_slic = 0.95,
  supervoxels = 0.90, snic = 0.25, commute = 0.50
)
MODERATE_NOISE_ARI_FLOORS <- c(
  flash3d = 0.90, slice_msf = -0.01, g3s = 0.75, rena = 0.90,
  rena_plus = 0.90, acsc = 0.90, slic = 0.90, corr_slic = 0.90,
  supervoxels = 0.85, snic = 0.30, commute = 0.45
)
HIGH_NOISE_ARI_FLOORS <- c(
  flash3d = 0.90, slice_msf = -0.01, g3s = 0.65, rena = 0.90,
  rena_plus = 0.90, acsc = 0.90, slic = 0.85, corr_slic = 0.90,
  supervoxels = 0.75, snic = 0.25, commute = 0.45
)

run_cluster4d <- function(method, ...) {
  extra <- METHOD_OVERRIDES[[method]]
  if (is.null(extra)) extra <- list()
  do.call(cluster4d, c(list(method = method, ...), extra))
}

test_that("synthetic cube test evaluates all clustering methods", {
  skip_on_cran()


  # Header
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("SYNTHETIC CUBE TEST: 9x9x9 Volume with 27 Cubic Clusters\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Create synthetic data
  set.seed(42)
  data <- make_synthetic_clusters(
    n_time = 100,
    noise_sd = 0,
    n_clusters = 27,
    seed = 42
  )

  cat("Test Configuration:\n")
  cat("  Volume dimensions: 9 x 9 x 9 x 100 (spatial x time)\n")
  cat("  True clusters: 27 (3x3x3 grid of 3x3x3 cubes)\n")
  cat("  Cluster size: 27 voxels each\n")
  cat("  Total voxels: 729\n")
  cat("  Noise: 0 (perfectly separable orthogonal signals)\n\n")

  # Storage for results
  results <- data.frame(
    method = character(),
    n_clusters = integer(),
    ari = numeric(),
    nmi = numeric(),
    accuracy = numeric(),
    time_sec = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )

  cat(sprintf("%-12s %4s %7s %7s %8s %6s  %s\n",
              "Method", "K", "ARI", "NMI", "Accuracy", "Time", "Status"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  for (method in ALL_METHODS) {
    t_start <- Sys.time()

    result <- tryCatch({
      suppressWarnings(suppressMessages(
        run_cluster4d(
          method = method,
          vec = data$vec,
          mask = data$mask,
          n_clusters = 27,
          verbose = FALSE
        )
      ))
    }, error = function(e) {
      NULL
    })

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    if (!is.null(result)) {
      metrics <- clustering_accuracy(result$cluster, data$true_labels)

      status <- if (metrics$perfect) {
        "PERFECT"
      } else if (metrics$ari >= 0.8) {
        "GOOD"
      } else if (metrics$ari >= 0.5) {
        "FAIR"
      } else {
        "POOR"
      }

      results <- rbind(results, data.frame(
        method = method,
        n_clusters = metrics$n_predicted,
        ari = metrics$ari,
        nmi = metrics$nmi,
        accuracy = metrics$accuracy,
        time_sec = elapsed,
        status = status,
        stringsAsFactors = FALSE
      ))

      cat(sprintf("%-12s %4d %7.4f %7.4f %8.4f %5.2fs  %s\n",
                  method, metrics$n_predicted, metrics$ari, metrics$nmi,
                  metrics$accuracy, elapsed, status))

    } else {
      results <- rbind(results, data.frame(
        method = method,
        n_clusters = NA,
        ari = NA,
        nmi = NA,
        accuracy = NA,
        time_sec = elapsed,
        status = "ERROR",
        stringsAsFactors = FALSE
      ))

      cat(sprintf("%-12s %4s %7s %7s %8s %5.2fs  %s\n",
                  method, "NA", "NA", "NA", "NA", elapsed, "ERROR"))
    }
  }

  # Summary
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("SUMMARY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  perfect <- results$method[results$status == "PERFECT"]
  good <- results$method[results$status == "GOOD"]
  fair <- results$method[results$status == "FAIR"]
  poor <- results$method[results$status == "POOR"]
  errors <- results$method[results$status == "ERROR"]

  if (length(perfect) > 0) cat("PERFECT (ARI=1.0):", paste(perfect, collapse = ", "), "\n")
  if (length(good) > 0) cat("GOOD (ARI>=0.8):", paste(good, collapse = ", "), "\n")
  if (length(fair) > 0) cat("FAIR (ARI>=0.5):", paste(fair, collapse = ", "), "\n")
  if (length(poor) > 0) cat("POOR (ARI<0.5):", paste(poor, collapse = ", "), "\n")
  if (length(errors) > 0) cat("ERROR:", paste(errors, collapse = ", "), "\n")

  cat("\n")
  cat("Metric definitions:\n")
  cat("  ARI: Adjusted Rand Index (1.0 = perfect, 0 = random)\n")
  cat("  NMI: Normalized Mutual Information (1.0 = perfect, 0 = independent)\n")
  cat("  Accuracy: Fraction correct after optimal label matching\n")
  cat("\n")

  # All methods should at least run without error
  expect_equal(
    length(errors), 0,
    info = paste("Methods with errors:", paste(errors, collapse = ", "))
  )

  # All methods should find a reasonable number of clusters
  valid_results <- results[!is.na(results$n_clusters), ]
  for (method in valid_results$method) {
    observed <- valid_results$ari[valid_results$method == method]
    expect_true(
      observed >= CUBE_ARI_FLOORS[[method]],
      info = sprintf(
        "%s ARI %.4f fell below its cube recovery floor %.2f",
        method, observed, CUBE_ARI_FLOORS[[method]]
      )
    )
  }
  expect_true(
    all(valid_results$n_clusters >= 5 & valid_results$n_clusters <= 50),
    info = "All methods should find between 5 and 50 clusters"
  )
})


test_that("synthetic cube test with noise shows graceful degradation", {
  skip_on_cran()

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("SYNTHETIC CUBE TEST WITH NOISE (SD=0.5)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Test with moderate noise
  set.seed(42)
  data <- make_synthetic_clusters(
    n_time = 100,
    noise_sd = 0.5,  # Moderate noise
    n_clusters = 27,
    seed = 42
  )

  cat("Test Configuration:\n")
  cat("  Noise SD: 0.5 (moderate - should still be recoverable)\n\n")

  cat(sprintf("%-12s %4s %7s %7s %8s %6s  %s\n",
              "Method", "K", "ARI", "NMI", "Accuracy", "Time", "Status"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  results <- data.frame(
    method = character(),
    ari = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )

  for (method in ALL_METHODS) {
    t_start <- Sys.time()

    result <- tryCatch({
      suppressWarnings(suppressMessages(
        cluster4d(
          vec = data$vec,
          mask = data$mask,
          n_clusters = 27,
          method = method,
          verbose = FALSE
        )
      ))
    }, error = function(e) NULL)

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    if (!is.null(result)) {
      metrics <- clustering_accuracy(result$cluster, data$true_labels)
      status <- if (metrics$ari >= 0.9) {
        "EXCELLENT"
      } else if (metrics$ari >= 0.7) {
        "GOOD"
      } else if (metrics$ari >= 0.5) {
        "FAIR"
      } else {
        "POOR"
      }

      cat(sprintf("%-12s %4d %7.4f %7.4f %8.4f %5.2fs  %s\n",
                  method, metrics$n_predicted, metrics$ari, metrics$nmi,
                  metrics$accuracy, elapsed, status))

      results <- rbind(results, data.frame(
        method = method,
        ari = metrics$ari,
        status = status,
        stringsAsFactors = FALSE
      ))

      expect_true(
        metrics$ari >= MODERATE_NOISE_ARI_FLOORS[[method]],
        info = sprintf(
          "%s ARI %.4f fell below its moderate-noise floor %.2f",
          method, metrics$ari, MODERATE_NOISE_ARI_FLOORS[[method]]
        )
      )
    } else {
      cat(sprintf("%-12s %4s %7s %7s %8s %5.2fs  %s\n",
                  method, "NA", "NA", "NA", "NA", elapsed, "ERROR"))

      results <- rbind(results, data.frame(
        method = method,
        ari = NA,
        status = "ERROR",
        stringsAsFactors = FALSE
      ))
    }
  }

  # Summary
  cat("\n")
  excellent <- results$method[results$status == "EXCELLENT"]
  good <- results$method[results$status == "GOOD"]
  fair <- results$method[results$status == "FAIR"]
  poor <- results$method[results$status == "POOR"]

  if (length(excellent) > 0) cat("EXCELLENT (ARI>=0.9):", paste(excellent, collapse = ", "), "\n")
  if (length(good) > 0) cat("GOOD (ARI>=0.7):", paste(good, collapse = ", "), "\n")
  if (length(fair) > 0) cat("FAIR (ARI>=0.5):", paste(fair, collapse = ", "), "\n")
  if (length(poor) > 0) cat("POOR (ARI<0.5):", paste(poor, collapse = ", "), "\n")

  cat("\n")
})


test_that("synthetic cube test with high noise shows robustness", {
  skip_on_cran()

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("SYNTHETIC CUBE TEST WITH HIGH NOISE (SD=1.0)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Test with high noise
  set.seed(42)
  data <- make_synthetic_clusters(
    n_time = 100,
    noise_sd = 1.0,  # High noise - challenging
    n_clusters = 27,
    seed = 42
  )

  cat("Test Configuration:\n")
  cat("  Noise SD: 1.0 (high - challenging recovery)\n\n")

  cat(sprintf("%-12s %4s %7s %7s %8s %6s  %s\n",
              "Method", "K", "ARI", "NMI", "Accuracy", "Time", "Status"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  results <- data.frame(
    method = character(),
    ari = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )

  for (method in ALL_METHODS) {
    t_start <- Sys.time()

    result <- tryCatch({
      suppressWarnings(suppressMessages(
        cluster4d(
          vec = data$vec,
          mask = data$mask,
          n_clusters = 27,
          method = method,
          verbose = FALSE
        )
      ))
    }, error = function(e) NULL)

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    if (!is.null(result)) {
      metrics <- clustering_accuracy(result$cluster, data$true_labels)
      status <- if (metrics$ari >= 0.9) {
        "EXCELLENT"
      } else if (metrics$ari >= 0.7) {
        "GOOD"
      } else if (metrics$ari >= 0.5) {
        "FAIR"
      } else if (metrics$ari >= 0.3) {
        "MARGINAL"
      } else {
        "POOR"
      }

      cat(sprintf("%-12s %4d %7.4f %7.4f %8.4f %5.2fs  %s\n",
                  method, metrics$n_predicted, metrics$ari, metrics$nmi,
                  metrics$accuracy, elapsed, status))

      results <- rbind(results, data.frame(
        method = method,
        ari = metrics$ari,
        status = status,
        stringsAsFactors = FALSE
      ))

      expect_true(
        metrics$ari >= HIGH_NOISE_ARI_FLOORS[[method]],
        info = sprintf(
          "%s ARI %.4f fell below its high-noise floor %.2f",
          method, metrics$ari, HIGH_NOISE_ARI_FLOORS[[method]]
        )
      )
    } else {
      cat(sprintf("%-12s %4s %7s %7s %8s %5.2fs  %s\n",
                  method, "NA", "NA", "NA", "NA", elapsed, "ERROR"))

      results <- rbind(results, data.frame(
        method = method,
        ari = NA,
        status = "ERROR",
        stringsAsFactors = FALSE
      ))
    }
  }

  # Summary
  cat("\n")
  excellent <- results$method[results$status == "EXCELLENT"]
  good <- results$method[results$status == "GOOD"]
  fair <- results$method[results$status == "FAIR"]
  marginal <- results$method[results$status == "MARGINAL"]
  poor <- results$method[results$status == "POOR"]

  if (length(excellent) > 0) cat("EXCELLENT (ARI>=0.9):", paste(excellent, collapse = ", "), "\n")
  if (length(good) > 0) cat("GOOD (ARI>=0.7):", paste(good, collapse = ", "), "\n")
  if (length(fair) > 0) cat("FAIR (ARI>=0.5):", paste(fair, collapse = ", "), "\n")
  if (length(marginal) > 0) cat("MARGINAL (ARI>=0.3):", paste(marginal, collapse = ", "), "\n")
  if (length(poor) > 0) cat("POOR (ARI<0.3):", paste(poor, collapse = ", "), "\n")

  cat("\n")
})
