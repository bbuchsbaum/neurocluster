#include <Rcpp.h>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Helper structs
// -----------------------------------------------------------------------------

struct G3S_VoxelNode {
  int voxel_idx;
  int cluster_label;
  double cost;

  bool operator>(const G3S_VoxelNode& other) const {
    return cost > other.cost;
  }
};

struct G3S_Centroid {
  std::vector<double> feature_sum;
  std::vector<double> feature_avg;
  double coord_sum_x, coord_sum_y, coord_sum_z;
  double coord_avg_x, coord_avg_y, coord_avg_z;
  int count;

  G3S_Centroid(int n_features,
               double x, double y, double z,
               const double* feat)
    : coord_sum_x(x), coord_sum_y(y), coord_sum_z(z),
      coord_avg_x(x), coord_avg_y(y), coord_avg_z(z),
      count(1) {
    feature_sum.assign(n_features, 0.0);
    feature_avg.assign(n_features, 0.0);
    for (int i = 0; i < n_features; ++i) {
      feature_sum[i] = feat[i];
      feature_avg[i] = feat[i];
    }
  }

  void add(double x, double y, double z,
           const double* feat,
           int n_features) {
    count++;
    coord_sum_x += x; coord_sum_y += y; coord_sum_z += z;
    coord_avg_x = coord_sum_x / count;
    coord_avg_y = coord_sum_y / count;
    coord_avg_z = coord_sum_z / count;

    double sq_sum = 0.0;
    for (int i = 0; i < n_features; ++i) {
      feature_sum[i] += feat[i];
      feature_avg[i] = feature_sum[i] / count;
      sq_sum += feature_avg[i] * feature_avg[i];
    }

    if (sq_sum > 0) {
      double inv = 1.0 / std::sqrt(sq_sum);
      for (int i = 0; i < n_features; ++i) {
        feature_avg[i] *= inv;
      }
    }
  }
};

inline double g3s_cost(const G3S_Centroid& C,
                       double x, double y, double z,
                       const double* feat,
                       int n_features,
                       double alpha,
                       double compactness) {
  double dx = C.coord_avg_x - x;
  double dy = C.coord_avg_y - y;
  double dz = C.coord_avg_z - z;
  double spatial_dist = std::sqrt(dx*dx + dy*dy + dz*dz);

  double dot = 0.0;
  for (int i = 0; i < n_features; ++i) {
    dot += C.feature_avg[i] * feat[i];
  }
  double feature_dist = 1.0 - dot;

  return alpha * feature_dist + (1.0 - alpha) * (spatial_dist / compactness);
}

// -----------------------------------------------------------------------------
// 1. Local gradient (C++)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector calculate_local_gradient(const NumericMatrix& feature_mat,
                                       const IntegerMatrix& neighbor_indices) {
  int N = feature_mat.ncol();
  int M = feature_mat.nrow();
  int K = neighbor_indices.ncol();

  NumericVector grad(N);

  for (int i = 0; i < N; ++i) {
    const double* v = &feature_mat(0, i);
    double total = 0.0;
    int valid = 0;

    for (int k = 0; k < K; ++k) {
      int n_idx = neighbor_indices(i, k) - 1;
      if (n_idx < 0 || n_idx >= N) continue;
      const double* n = &feature_mat(0, n_idx);
      double dot = 0.0;
      for (int m = 0; m < M; ++m) dot += v[m] * n[m];
      total += dot;
      valid++;
    }

    if (valid > 0) {
      grad[i] = 1.0 - (total / valid);
    } else {
      grad[i] = 1.0;
    }
  }

  return grad;
}

// -----------------------------------------------------------------------------
// 2. Geodesic propagation (neighbor-driven, optimized)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector g3s_propagate_cpp(const NumericMatrix& feature_mat,
                                const NumericMatrix& coords,
                                const IntegerVector& seed_indices,
                                const IntegerMatrix& neighbor_indices,
                                const NumericMatrix& neighbor_dists,
                                double alpha,
                                double compactness) {

  int N = feature_mat.ncol();
  int M = feature_mat.nrow();
  int K = seed_indices.size();
  int K_neighbors = neighbor_indices.ncol();

  std::vector<int> labels(N, 0);
  std::vector<G3S_Centroid> centroids;
  centroids.reserve(K);

  std::priority_queue<G3S_VoxelNode,
                      std::vector<G3S_VoxelNode>,
                      std::greater<G3S_VoxelNode>> pq;

  // Initialize seeds
  for (int k = 0; k < K; ++k) {
    int idx = seed_indices[k] - 1;
    if (idx < 0 || idx >= N) continue;

    centroids.emplace_back(M,
                           coords(idx, 0), coords(idx, 1), coords(idx, 2),
                           &feature_mat(0, idx));
    labels[idx] = k + 1;

    for (int n = 0; n < K_neighbors; ++n) {
      int neigh = neighbor_indices(idx, n) - 1;
      if (neigh < 0 || neigh >= N) continue;
      if (labels[neigh] != 0) continue;

      double cost = g3s_cost(centroids[k],
                             coords(neigh, 0), coords(neigh, 1), coords(neigh, 2),
                             &feature_mat(0, neigh),
                             M, alpha, compactness);
      pq.push({neigh, k + 1, cost});
    }
  }

  int assigned = static_cast<int>(centroids.size());

  while (!pq.empty() && assigned < N) {
    if (assigned % 2000 == 0) Rcpp::checkUserInterrupt();

    G3S_VoxelNode top = pq.top();
    pq.pop();

    if (labels[top.voxel_idx] != 0) continue;

    labels[top.voxel_idx] = top.cluster_label;
    assigned++;

    int c_idx = top.cluster_label - 1;
    centroids[c_idx].add(coords(top.voxel_idx, 0),
                         coords(top.voxel_idx, 1),
                         coords(top.voxel_idx, 2),
                         &feature_mat(0, top.voxel_idx), M);

    for (int n = 0; n < K_neighbors; ++n) {
      int neigh = neighbor_indices(top.voxel_idx, n) - 1;
      if (neigh < 0 || neigh >= N) continue;
      if (labels[neigh] != 0) continue;

      double cost = g3s_cost(centroids[c_idx],
                             coords(neigh, 0), coords(neigh, 1), coords(neigh, 2),
                             &feature_mat(0, neigh),
                             M, alpha, compactness);
      pq.push({neigh, top.cluster_label, cost});
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
                                        const IntegerMatrix& neighbor_indices,
                                        int max_iter) {

  int N = labels.size();
  int M = feature_mat.nrow();
  int K_neighbors = neighbor_indices.ncol();

  std::vector<int> cur(labels.begin(), labels.end());
  std::vector<int> next = cur;

  for (int iter = 0; iter < max_iter; ++iter) {
    int max_label = 0;
    for (int l : cur) if (l > max_label) max_label = l;

    std::vector<std::vector<double>> centroids(max_label + 1,
                                               std::vector<double>(M, 0.0));
    std::vector<int> counts(max_label + 1, 0);

    for (int i = 0; i < N; ++i) {
      int lab = cur[i];
      if (lab <= 0) continue;
      counts[lab]++;
      for (int m = 0; m < M; ++m) centroids[lab][m] += feature_mat(m, i);
    }

    for (int lab = 1; lab <= max_label; ++lab) {
      if (counts[lab] == 0) continue;
      double inv = 1.0 / counts[lab];
      double sq = 0.0;
      for (int m = 0; m < M; ++m) {
        centroids[lab][m] *= inv;
        sq += centroids[lab][m] * centroids[lab][m];
      }
      if (sq > 0) {
        double norm = 1.0 / std::sqrt(sq);
        for (int m = 0; m < M; ++m) centroids[lab][m] *= norm;
      }
    }

    int changed = 0;

    for (int i = 0; i < N; ++i) {
      int lab = cur[i];
      if (lab <= 0) continue;

      bool boundary = false;
      std::vector<int> candidates;
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

      double best_dot = -2.0;
      int best_lab = lab;

      for (int cand : candidates) {
        if (counts[cand] == 0) continue;
        double dot = 0.0;
        for (int m = 0; m < M; ++m) {
          dot += feature_mat(m, i) * centroids[cand][m];
        }
        if (dot > best_dot) {
          best_dot = dot;
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
