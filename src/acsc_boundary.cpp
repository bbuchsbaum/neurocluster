// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <atomic>

using namespace Rcpp;
using namespace RcppParallel;

// =============================================================================
// Helper Functions
// =============================================================================

// Dot product between a *row* of an R matrix (column-major) and a contiguous centroid.
// The row is strided by n_rows in memory: x[t] lives at base[row + n_rows * t].
inline double fast_row_dot_strided(const double* base, int row, int n_rows,
                                  const double* centroid, int n_time) {
  double dot = 0.0;
  const double* x0 = base + row;
  for (int t = 0; t < n_time; ++t) {
    dot += x0[n_rows * t] * centroid[t];
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

// Dense centroid table for labels in [0, max_label].
struct CentroidTable {
  int max_label;
  int n_time;
  std::vector<double> data; // (max_label+1) * n_time, contiguous per label
  std::vector<int> counts;  // (max_label+1)
};

static CentroidTable compute_centroid_table(const double* feature_mat,
                                            const int* labels,
                                            int n_voxels,
                                            int n_time) {
  int max_label = 0;
  for (int i = 0; i < n_voxels; ++i) {
    if (labels[i] > max_label) max_label = labels[i];
  }

  CentroidTable tab;
  tab.max_label = max_label;
  tab.n_time = n_time;
  tab.data.assign((size_t)(max_label + 1) * (size_t)n_time, 0.0);
  tab.counts.assign((size_t)(max_label + 1), 0);

  // R matrices are column-major: feature_mat is (voxels x time).
  const int stride = n_voxels;
  for (int i = 0; i < n_voxels; ++i) {
    int label = labels[i];
    if (label < 0) continue;
    if (label > max_label) continue;
    tab.counts[(size_t)label] += 1;
    double* sum = &tab.data[(size_t)label * (size_t)n_time];
    const double* x0 = feature_mat + i;
    for (int t = 0; t < n_time; ++t) {
      sum[t] += x0[stride * t];
    }
  }

  // Convert sums to means and normalize in-place.
  for (int label = 0; label <= max_label; ++label) {
    int c = tab.counts[(size_t)label];
    if (c <= 0) continue;
    double inv = 1.0 / (double)c;
    double* centroid = &tab.data[(size_t)label * (size_t)n_time];
    for (int t = 0; t < n_time; ++t) {
      centroid[t] *= inv;
    }
    normalize_vector(centroid, n_time);
  }

  return tab;
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

  // Centroid table (dense)
  const CentroidTable& centroids;
  const int n_time;
  const int n_voxels;
  const double* feat_base;

  // Output data (write to separate indices, no contention)
  RVector<int> output_labels;
  std::atomic<int>& changes_counter;

  // Constructor
  BoundaryRefinementWorker(
      const NumericMatrix& feat_mat,
      const IntegerVector& curr_labels,
      const IntegerMatrix& nbr_indices,
      const IntegerVector& boundary_vox,
      const CentroidTable& cents,
      IntegerVector& out_labels,
      std::atomic<int>& counter)
    : feature_mat(feat_mat),
      current_labels(curr_labels),
      neighbor_indices(nbr_indices),
      boundary_voxels(boundary_vox),
      centroids(cents),
      n_time(feat_mat.ncol()),
      n_voxels(feat_mat.nrow()),
      feat_base(feat_mat.begin()),
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
        if (candidate_label < 0 || candidate_label > centroids.max_label) continue;
        if (centroids.counts[(size_t)candidate_label] <= 0) continue;
        const double* centroid = &centroids.data[(size_t)candidate_label * (size_t)n_time];

        // Compute correlation (both vectors are normalized, so just dot product).
        // feature_mat is column-major, so rows are strided.
        double correlation = fast_row_dot_strided(feat_base, voxel_idx, n_voxels, centroid, n_time);

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

    // Compute centroids for current labeling (dense + stride-correct)
    CentroidTable centroids = compute_centroid_table(
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
