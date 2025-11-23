// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <atomic>

using namespace Rcpp;
using namespace RcppParallel;

// =============================================================================
// Helper Functions
// =============================================================================

/**
 * Fast correlation for pre-normalized (unit-length) vectors
 *
 * For unit-length vectors: cor(x, y) = dot(x, y)
 * This eliminates the need for expensive standard deviation calculations
 *
 * @param x Pointer to first normalized vector
 * @param y Pointer to second normalized vector
 * @param n Length of vectors
 * @return Correlation coefficient
 */
inline double fast_correlation_normalized(const double* x, const double* y, int n) {
  double dot = 0.0;
  for (int i = 0; i < n; i++) {
    dot += x[i] * y[i];
  }
  return dot;
}

/**
 * Normalize a vector to unit length in-place
 *
 * @param vec Pointer to vector to normalize
 * @param n Length of vector
 */
inline void normalize_vector(double* vec, int n) {
  double norm_sq = 0.0;
  for (int i = 0; i < n; i++) {
    norm_sq += vec[i] * vec[i];
  }

  if (norm_sq > 1e-12) {
    double norm = std::sqrt(norm_sq);
    for (int i = 0; i < n; i++) {
      vec[i] /= norm;
    }
  }
}

/**
 * Compute cluster centroids from labeled voxels
 *
 * @param feature_mat Matrix of voxel features (voxels x time)
 * @param labels Vector of cluster labels for each voxel
 * @param n_voxels Number of voxels
 * @param n_time Number of time points
 * @return Map from label to normalized centroid vector
 */
std::map<int, std::vector<double>> compute_centroids(
    const double* feature_mat,
    const int* labels,
    int n_voxels,
    int n_time) {

  std::map<int, std::vector<double>> centroid_sums;
  std::map<int, int> counts;

  // Accumulate sums per cluster
  for (int i = 0; i < n_voxels; i++) {
    int label = labels[i];

    if (centroid_sums.find(label) == centroid_sums.end()) {
      centroid_sums[label] = std::vector<double>(n_time, 0.0);
      counts[label] = 0;
    }

    std::vector<double>& sum = centroid_sums[label];
    for (int t = 0; t < n_time; t++) {
      sum[t] += feature_mat[i * n_time + t];
    }
    counts[label]++;
  }

  // Compute means and normalize
  std::map<int, std::vector<double>> centroids;
  for (auto& kv : centroid_sums) {
    int label = kv.first;
    std::vector<double>& sum = kv.second;
    int count = counts[label];

    std::vector<double> centroid(n_time);
    for (int t = 0; t < n_time; t++) {
      centroid[t] = sum[t] / count;
    }

    // Normalize to unit length
    normalize_vector(centroid.data(), n_time);
    centroids[label] = centroid;
  }

  return centroids;
}

// =============================================================================
// Boundary Refinement Worker
// =============================================================================

/**
 * RcppParallel worker for boundary voxel refinement
 *
 * Processes boundary voxels in parallel, reassigning each voxel to the cluster
 * whose centroid has the highest correlation with the voxel's time series.
 */
struct BoundaryRefinementWorker : public Worker {
  // Input data (read-only)
  const RMatrix<double> feature_mat;      // Voxels x Time (normalized)
  const RVector<int> current_labels;      // Current cluster assignments
  const RMatrix<int> neighbor_indices;    // Voxels x 6 nearest neighbors
  const RVector<int> boundary_voxels;     // Indices of boundary voxels

  // Centroid data
  const std::map<int, std::vector<double>>& centroids;
  int n_time;

  // Output data (write to separate indices, no contention)
  RVector<int> output_labels;
  std::atomic<int>& changes_counter;

  // Constructor
  BoundaryRefinementWorker(
      const NumericMatrix& feat_mat,
      const IntegerVector& curr_labels,
      const IntegerMatrix& nbr_indices,
      const IntegerVector& boundary_vox,
      const std::map<int, std::vector<double>>& cents,
      IntegerVector& out_labels,
      std::atomic<int>& counter)
    : feature_mat(feat_mat),
      current_labels(curr_labels),
      neighbor_indices(nbr_indices),
      boundary_voxels(boundary_vox),
      centroids(cents),
      n_time(feat_mat.ncol()),
      output_labels(out_labels),
      changes_counter(counter) {}

  // Parallel operator - processes a range of boundary voxels
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t idx = begin; idx < end; idx++) {
      // Get actual voxel index (boundary_voxels is 1-based from R)
      int voxel_idx = boundary_voxels[idx] - 1;
      int current_label = current_labels[voxel_idx];

      // Get neighbor labels to identify candidates
      std::vector<int> candidate_labels;
      candidate_labels.push_back(current_label);  // Always consider current

      for (int n = 0; n < neighbor_indices.ncol(); n++) {
        int neighbor_idx = neighbor_indices(voxel_idx, n) - 1;  // Convert to 0-based
        if (neighbor_idx >= 0 && neighbor_idx < current_labels.length()) {
          int neighbor_label = current_labels[neighbor_idx];

          // Add unique neighbor labels
          if (std::find(candidate_labels.begin(), candidate_labels.end(), neighbor_label)
              == candidate_labels.end()) {
            candidate_labels.push_back(neighbor_label);
          }
        }
      }

      // Find best cluster by correlation
      double best_correlation = -2.0;
      int best_label = current_label;

      for (int candidate_label : candidate_labels) {
        auto it = centroids.find(candidate_label);
        if (it == centroids.end()) continue;  // Skip if centroid not found

        const std::vector<double>& centroid = it->second;

        // Compute correlation (both vectors are normalized, so just dot product)
        const double* voxel_features = &feature_mat(voxel_idx, 0);
        double correlation = fast_correlation_normalized(
          voxel_features, centroid.data(), n_time);

        if (correlation > best_correlation) {
          best_correlation = correlation;
          best_label = candidate_label;
        }
      }

      // Update label if changed
      output_labels[voxel_idx] = best_label;
      if (best_label != current_label) {
        changes_counter.fetch_add(1, std::memory_order_relaxed);
      }
    }
  }
};

// =============================================================================
// Exported Functions
// =============================================================================

//' Refine boundary voxels using C++ acceleration
//'
//' @param voxel_labels Integer vector of current cluster assignments (0-based internally)
//' @param feature_mat_normalized Numeric matrix (voxels x time) with NORMALIZED features
//' @param neighbor_indices Integer matrix (voxels x K) of nearest neighbor indices (1-based from R)
//' @param boundary_voxels Integer vector of boundary voxel indices (1-based from R)
//' @param max_iter Maximum number of refinement iterations
//'
//' @return List containing updated labels and number of iterations performed
//'
//' @keywords internal
// [[Rcpp::export]]
List refine_boundaries_cpp(
    IntegerVector voxel_labels,
    NumericMatrix feature_mat_normalized,
    IntegerMatrix neighbor_indices,
    IntegerVector boundary_voxels,
    int max_iter = 5) {

  int n_voxels = feature_mat_normalized.nrow();
  int n_time = feature_mat_normalized.ncol();

  // Validate inputs
  if (voxel_labels.length() != n_voxels) {
    stop("Length of voxel_labels must match number of rows in feature_mat_normalized");
  }
  if (neighbor_indices.nrow() != n_voxels) {
    stop("Number of rows in neighbor_indices must match number of voxels");
  }

  // Working copy of labels
  IntegerVector current_labels = clone(voxel_labels);
  int total_iterations = 0;

  // Iterative refinement
  for (int iter = 0; iter < max_iter; iter++) {
    total_iterations++;

    // Compute centroids for current labeling
    std::map<int, std::vector<double>> centroids = compute_centroids(
      feature_mat_normalized.begin(),
      current_labels.begin(),
      n_voxels,
      n_time
    );

    // Output labels (initialized as copy of current)
    IntegerVector output_labels = clone(current_labels);

    // Atomic counter for changes
    std::atomic<int> changes_counter(0);

    // Create worker and run parallel refinement
    BoundaryRefinementWorker worker(
      feature_mat_normalized,
      current_labels,
      neighbor_indices,
      boundary_voxels,
      centroids,
      output_labels,
      changes_counter
    );

    // Parallelize over boundary voxels with grain size control
    int grain_size = std::max(1, (int)boundary_voxels.length() / 100);
    parallelFor(0, boundary_voxels.length(), worker, grain_size);

    // Check convergence
    int n_changes = changes_counter.load();
    current_labels = output_labels;

    if (n_changes == 0) {
      break;  // Converged
    }
  }

  return List::create(
    Named("labels") = current_labels,
    Named("iterations") = total_iterations
  );
}

//' Find boundary voxels (voxels with neighbors having different labels)
//'
//' @param voxel_labels Integer vector of cluster assignments
//' @param neighbor_indices Integer matrix (voxels x K) of nearest neighbor indices (1-based)
//'
//' @return Integer vector of boundary voxel indices (1-based for R)
//'
//' @keywords internal
// [[Rcpp::export]]
IntegerVector find_boundary_voxels_cpp(
    IntegerVector voxel_labels,
    IntegerMatrix neighbor_indices) {

  int n_voxels = voxel_labels.length();
  std::vector<int> boundary_indices;

  for (int i = 0; i < n_voxels; i++) {
    int current_label = voxel_labels[i];
    bool is_boundary = false;

    // Check if any neighbor has different label
    for (int n = 0; n < neighbor_indices.ncol(); n++) {
      int neighbor_idx = neighbor_indices(i, n) - 1;  // Convert to 0-based

      if (neighbor_idx >= 0 && neighbor_idx < n_voxels) {
        if (voxel_labels[neighbor_idx] != current_label) {
          is_boundary = true;
          break;
        }
      }
    }

    if (is_boundary) {
      boundary_indices.push_back(i + 1);  // Convert back to 1-based for R
    }
  }

  return wrap(boundary_indices);
}
