#include <Rcpp.h>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Helper structs
// -----------------------------------------------------------------------------

namespace {

struct G3S_VoxelNode {
  int voxel_idx;
  int cluster_label;
  double cost;

  bool operator>(const G3S_VoxelNode& other) const {
    if (cost != other.cost) return cost > other.cost;
    if (cluster_label != other.cluster_label) {
      return cluster_label > other.cluster_label;
    }
    return voxel_idx > other.voxel_idx;
  }
};

inline double g3s_edge_feature_distance(
    const NumericMatrix& feature_mat, int left, int right) {
  const int M = feature_mat.nrow();
  if (M > 1) {
    double dot = 0.0;
    for (int m = 0; m < M; ++m) {
      dot += feature_mat(m, left) * feature_mat(m, right);
    }
    dot = std::max(-1.0, std::min(1.0, dot));
    return 1.0 - dot;
  }
  return std::abs(feature_mat(0, left) - feature_mat(0, right));
}

} // namespace

// -----------------------------------------------------------------------------
// 1. Local gradient (C++)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector calculate_local_gradient(const NumericMatrix& feature_mat,
                                       const IntegerMatrix& neighbor_indices) {
  int N = feature_mat.ncol();
  int M = feature_mat.nrow();
  int K = neighbor_indices.ncol();
  bool use_cosine = M > 1;

  NumericVector grad(N);

  for (int i = 0; i < N; ++i) {
    const double* v = &feature_mat(0, i);
    double total = 0.0;
    double dist_sum = 0.0;
    int valid = 0;

    for (int k = 0; k < K; ++k) {
      int n_idx = neighbor_indices(i, k) - 1;
      if (n_idx < 0 || n_idx >= N) continue;
      const double* n = &feature_mat(0, n_idx);
      if (use_cosine) {
        double dot = 0.0;
        for (int m = 0; m < M; ++m) dot += v[m] * n[m];
        total += dot;
      } else {
        double sq = 0.0;
        for (int m = 0; m < M; ++m) {
          double diff = v[m] - n[m];
          sq += diff * diff;
        }
        dist_sum += sq;
      }
      valid++;
    }

    if (valid > 0) {
      grad[i] = 1.0 - (total / valid);
    } else {
      grad[i] = 1.0;
    }

    if (!use_cosine) {
      grad[i] = (valid > 0) ? dist_sum / valid : 0.0;
    }
  }

  return grad;
}

// -----------------------------------------------------------------------------
// 2. Geodesic propagation (neighbor-driven, optimized)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector g3s_propagate_cpp(const NumericMatrix& feature_mat,
                                const IntegerVector& seed_indices,
                                const IntegerMatrix& neighbor_indices,
                                const NumericMatrix& neighbor_dists,
                                double alpha,
                                double compactness) {
  const int N = feature_mat.ncol();
  const int K = seed_indices.size();
  const int K_neighbors = neighbor_indices.ncol();
  if (N < 1 || feature_mat.nrow() < 1) stop("feature_mat must be non-empty");
  if (neighbor_indices.nrow() != N || neighbor_dists.nrow() != N ||
      neighbor_dists.ncol() != K_neighbors) {
    stop("neighbor matrices must have one aligned row per voxel");
  }
  if (K < 1) stop("at least one seed is required");
  if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
    stop("alpha must be in [0, 1]");
  }
  if (!std::isfinite(compactness) || compactness <= 0.0) {
    stop("compactness must be positive and finite");
  }
  for (int i = 0; i < N; ++i) {
    for (int m = 0; m < feature_mat.nrow(); ++m) {
      if (!std::isfinite(feature_mat(m, i))) {
        stop("feature_mat must contain only finite values");
      }
    }
  }

  std::vector<int> labels(N, 0);
  std::vector<int> best_label(N, 0);
  std::vector<double> best_cost(N, std::numeric_limits<double>::infinity());
  std::priority_queue<G3S_VoxelNode,
                      std::vector<G3S_VoxelNode>,
                      std::greater<G3S_VoxelNode>> pq;

  for (int k = 0; k < K; ++k) {
    const int idx = seed_indices[k] - 1;
    if (idx < 0 || idx >= N) stop("seed_indices are out of range");
    const int label = k + 1;
    if (best_cost[idx] == 0.0) stop("seed_indices must be unique");
    best_cost[idx] = 0.0;
    best_label[idx] = label;
    pq.push({idx, label, 0.0});
  }

  int finalized = 0;
  const double tolerance = 1e-12;
  while (!pq.empty()) {
    G3S_VoxelNode top = pq.top();
    pq.pop();
    if (labels[top.voxel_idx] != 0) continue;
    if (std::abs(top.cost - best_cost[top.voxel_idx]) > tolerance ||
        top.cluster_label != best_label[top.voxel_idx]) continue;

    labels[top.voxel_idx] = top.cluster_label;
    ++finalized;
    if (finalized % 2000 == 0) Rcpp::checkUserInterrupt();

    for (int slot = 0; slot < K_neighbors; ++slot) {
      const int neighbor = neighbor_indices(top.voxel_idx, slot) - 1;
      if (neighbor < 0) continue;
      if (neighbor >= N) stop("neighbor index is out of range");
      const double physical_distance = neighbor_dists(top.voxel_idx, slot);
      if (!std::isfinite(physical_distance) || physical_distance <= 0.0) {
        stop("active neighbor distances must be positive and finite");
      }
      if (labels[neighbor] != 0) continue;

      const double feature_distance = g3s_edge_feature_distance(
        feature_mat, top.voxel_idx, neighbor
      );
      const double edge_cost = alpha * feature_distance +
        (1.0 - alpha) * (physical_distance / compactness);
      const double candidate = top.cost + edge_cost;
      if (candidate + tolerance < best_cost[neighbor] ||
          (std::abs(candidate - best_cost[neighbor]) <= tolerance &&
           top.cluster_label < best_label[neighbor])) {
        best_cost[neighbor] = candidate;
        best_label[neighbor] = top.cluster_label;
        pq.push({neighbor, top.cluster_label, candidate});
      }
    }
  }

  IntegerVector out(N);
  for (int i = 0; i < N; ++i) out[i] = labels[i];
  return out;
}

// -----------------------------------------------------------------------------
// 3. Boundary refinement (C++)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector refine_boundaries_g3s_cpp(IntegerVector labels,
                                        const NumericMatrix& feature_mat,
                                        const NumericMatrix& coords,
                                        const IntegerMatrix& neighbor_indices,
                                        double alpha,
                                        double compactness,
                                        int max_iter) {

  int N = labels.size();
  int M = feature_mat.nrow();
  int K_neighbors = neighbor_indices.ncol();
  bool use_cosine = M > 1;
  if (feature_mat.ncol() != N || coords.nrow() != N || coords.ncol() != 3 ||
      neighbor_indices.nrow() != N) {
    stop("refinement inputs must have one aligned row or column per voxel");
  }
  if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
    stop("alpha must be in [0, 1]");
  }
  if (!std::isfinite(compactness) || compactness <= 0.0) {
    stop("compactness must be positive and finite");
  }
  if (max_iter < 0) stop("max_iter must be non-negative");

  std::vector<int> cur(labels.begin(), labels.end());
  std::vector<int> next = cur;

  for (int iter = 0; iter < max_iter; ++iter) {
    int max_label = 0;
    for (int l : cur) if (l > max_label) max_label = l;

    const int label_stride = max_label + 1;
    std::vector<int> counts(max_label + 1, 0);
    std::vector<double> centroids((size_t)label_stride * (size_t)M, 0.0);
    std::vector<double> coord_centroids((size_t)label_stride * 3U, 0.0);

    #ifdef _OPENMP
    {
      const int n_threads = omp_get_max_threads();
      std::vector<int> counts_tls((size_t)n_threads * (size_t)label_stride, 0);
      std::vector<double> sums_tls((size_t)n_threads * (size_t)label_stride * (size_t)M, 0.0);

      #pragma omp parallel num_threads(n_threads)
      {
        const int tid = omp_get_thread_num();
        int* count_local = counts_tls.data() + (size_t)tid * (size_t)label_stride;
        double* sum_local = sums_tls.data() + (size_t)tid * (size_t)label_stride * (size_t)M;

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
          const int lab = cur[i];
          if (lab <= 0) continue;

          count_local[lab]++;
          double* row = sum_local + (size_t)lab * (size_t)M;
          for (int m = 0; m < M; ++m) {
            row[m] += feature_mat(m, i);
          }
        }
      }

      for (int t = 0; t < n_threads; ++t) {
        const int* count_local = counts_tls.data() + (size_t)t * (size_t)label_stride;
        const double* sum_local = sums_tls.data() + (size_t)t * (size_t)label_stride * (size_t)M;
        for (int lab = 1; lab <= max_label; ++lab) {
          counts[lab] += count_local[lab];
          double* dst = centroids.data() + (size_t)lab * (size_t)M;
          const double* src = sum_local + (size_t)lab * (size_t)M;
          for (int m = 0; m < M; ++m) {
            dst[m] += src[m];
          }
        }
      }
    }
    #else
    for (int i = 0; i < N; ++i) {
      const int lab = cur[i];
      if (lab <= 0) continue;
      counts[lab]++;
      double* row = centroids.data() + (size_t)lab * (size_t)M;
      for (int m = 0; m < M; ++m) {
        row[m] += feature_mat(m, i);
      }
    }
    #endif

    for (int i = 0; i < N; ++i) {
      const int lab = cur[i];
      if (lab <= 0) continue;
      double* center = coord_centroids.data() + (size_t)lab * 3U;
      for (int axis = 0; axis < 3; ++axis) center[axis] += coords(i, axis);
    }

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int lab = 1; lab <= max_label; ++lab) {
      if (counts[lab] == 0) continue;
      double inv = 1.0 / counts[lab];
      double sq = 0.0;
      double* row = centroids.data() + (size_t)lab * (size_t)M;
      for (int m = 0; m < M; ++m) {
        row[m] *= inv;
        sq += row[m] * row[m];
      }
      if (use_cosine && sq > 0) {
        double norm = 1.0 / std::sqrt(sq);
        for (int m = 0; m < M; ++m) {
          row[m] *= norm;
        }
      }
      double* coord_center = coord_centroids.data() + (size_t)lab * 3U;
      for (int axis = 0; axis < 3; ++axis) coord_center[axis] *= inv;
    }

    int changed = 0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:changed) schedule(static)
    #endif
    for (int i = 0; i < N; ++i) {
      const int lab = cur[i];
      if (lab <= 0) continue;

      bool boundary = false;
      std::vector<int> candidates;
      candidates.reserve(8);
      candidates.push_back(lab);

      for (int k = 0; k < K_neighbors; ++k) {
        int neigh = neighbor_indices(i, k) - 1;
        if (neigh < 0 || neigh >= N) continue;
        int neigh_lab = cur[neigh];
        if (neigh_lab > 0 && neigh_lab != lab) {
          boundary = true;
          if (std::find(candidates.begin(), candidates.end(), neigh_lab) == candidates.end()) {
            candidates.push_back(neigh_lab);
          }
        }
      }

      if (!boundary) {
        next[i] = lab;
        continue;
      }

      double best_cost = std::numeric_limits<double>::infinity();
      int best_lab = lab;

      for (int cand : candidates) {
        if (counts[cand] == 0) continue;
        const double* c_row = centroids.data() + (size_t)cand * (size_t)M;
        double feature_distance = 0.0;
        if (use_cosine) {
          double dot = 0.0;
          for (int m = 0; m < M; ++m) {
            dot += feature_mat(m, i) * c_row[m];
          }
          dot = std::max(-1.0, std::min(1.0, dot));
          feature_distance = 1.0 - dot;
        } else {
          double dist_sq = 0.0;
          for (int m = 0; m < M; ++m) {
            double diff = feature_mat(m, i) - c_row[m];
            dist_sq += diff * diff;
          }
          feature_distance = std::sqrt(dist_sq);
        }
        const double* coord_center =
          coord_centroids.data() + (size_t)cand * 3U;
        double spatial_sq = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
          const double delta = coords(i, axis) - coord_center[axis];
          spatial_sq += delta * delta;
        }
        const double cost = alpha * feature_distance +
          (1.0 - alpha) * (std::sqrt(spatial_sq) / compactness);
        if (cost < best_cost - 1e-12 ||
            (std::abs(cost - best_cost) <= 1e-12 && cand < best_lab)) {
          best_cost = cost;
          best_lab = cand;
        }
      }

      next[i] = best_lab;
      if (best_lab != lab) changed++;
    }

    cur.swap(next);
    if (changed == 0) break;
  }

  return wrap(cur);
}
