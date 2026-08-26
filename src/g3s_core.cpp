#include <Rcpp.h>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstddef>
#include <unordered_map>
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
// 3. Label-connectivity repair
//
// Boundary refinement optimises a voxel-to-centroid objective with no
// connectivity constraint, so it routinely splits a cluster into a large piece
// plus stranded fragments. Left alone those fragments become extra connected
// components that force_exact_k() then has to merge away one at a time, and
// single-voxel fragments can survive every Ward merge and end up as 1-voxel
// clusters. Repairing them here keeps each label to a single component: every
// label retains its largest component and each stranded component is absorbed,
// as a unit, into the adjacent label it is cheapest to attach to under the same
// edge objective the propagation used.
// -----------------------------------------------------------------------------

namespace {

// Components of the subgraph induced by equal labels. Ids ascend with the
// smallest member voxel, which makes every downstream tie-break deterministic.
int g3s_label_components(const std::vector<int>& labels,
                         const std::vector<int>& ptr,
                         const std::vector<int>& idx,
                         std::vector<int>& comp) {
  const int n = static_cast<int>(labels.size());
  comp.assign(n, 0);
  std::vector<int> stack;
  int next_id = 0;
  for (int start = 0; start < n; ++start) {
    if (comp[start] != 0) continue;
    ++next_id;
    comp[start] = next_id;
    stack.push_back(start);
    while (!stack.empty()) {
      const int node = stack.back();
      stack.pop_back();
      for (int p = ptr[node]; p < ptr[node + 1]; ++p) {
        const int nb = idx[p];
        if (comp[nb] == 0 && labels[nb] == labels[node]) {
          comp[nb] = next_id;
          stack.push_back(nb);
        }
      }
    }
  }
  return next_id;
}

// Returns the number of voxels moved. `labels` is updated in place and every
// label that had at least one voxel still has exactly one connected component,
// except in the pathological case where a stranded component touches no
// surviving label at all; those voxels keep their original label and the
// exact-K engine remains the backstop.
int g3s_repair_connectivity(std::vector<int>& labels,
                            const std::vector<int>& ptr,
                            const std::vector<int>& idx,
                            const NumericMatrix& feature_mat,
                            const NumericMatrix& coords,
                            double alpha,
                            double compactness) {
  const int n = static_cast<int>(labels.size());
  std::vector<int> comp;
  const int n_comp = g3s_label_components(labels, ptr, idx, comp);
  if (n_comp <= 0) return 0;

  int max_label = 0;
  for (int v = 0; v < n; ++v) if (labels[v] > max_label) max_label = labels[v];
  if (max_label <= 0) return 0;

  std::vector<int> comp_size(static_cast<std::size_t>(n_comp) + 1, 0);
  std::vector<int> comp_label(static_cast<std::size_t>(n_comp) + 1, 0);
  for (int v = 0; v < n; ++v) {
    ++comp_size[comp[v]];
    comp_label[comp[v]] = labels[v];
  }

  // Largest component per label; ties go to the smaller component id.
  std::vector<int> keeper(static_cast<std::size_t>(max_label) + 1, 0);
  for (int c = 1; c <= n_comp; ++c) {
    const int lab = comp_label[c];
    if (lab <= 0) continue;
    const int cur = keeper[lab];
    if (cur == 0 || comp_size[c] > comp_size[cur]) keeper[lab] = c;
  }

  std::vector<int> original(labels);
  std::vector<char> stranded(static_cast<std::size_t>(n_comp) + 1, 0);
  int n_free = 0;
  for (int c = 1; c <= n_comp; ++c) {
    if (comp_label[c] > 0 && keeper[comp_label[c]] != c) {
      stranded[c] = 1;
      n_free += comp_size[c];
    }
  }
  if (n_free == 0) return 0;

  for (int v = 0; v < n; ++v) if (stranded[comp[v]]) labels[v] = 0;

  std::vector<double> best_cost(static_cast<std::size_t>(n_comp) + 1);
  std::vector<int> best_label(static_cast<std::size_t>(n_comp) + 1);
  const int moved_total = n_free;

  while (n_free > 0) {
    std::fill(best_cost.begin(), best_cost.end(),
              std::numeric_limits<double>::infinity());
    std::fill(best_label.begin(), best_label.end(), 0);

    for (int v = 0; v < n; ++v) {
      if (labels[v] != 0) continue;
      const int c = comp[v];
      for (int p = ptr[v]; p < ptr[v + 1]; ++p) {
        const int nb = idx[p];
        const int nb_label = labels[nb];
        if (nb_label <= 0) continue;
        double spatial_sq = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
          const double delta = coords(v, axis) - coords(nb, axis);
          spatial_sq += delta * delta;
        }
        const double cost =
          alpha * g3s_edge_feature_distance(feature_mat, v, nb) +
          (1.0 - alpha) * (std::sqrt(spatial_sq) / compactness);
        if (cost < best_cost[c] - 1e-12 ||
            (std::abs(cost - best_cost[c]) <= 1e-12 &&
             (best_label[c] == 0 || nb_label < best_label[c]))) {
          best_cost[c] = cost;
          best_label[c] = nb_label;
        }
      }
    }

    int assigned = 0;
    for (int v = 0; v < n; ++v) {
      if (labels[v] != 0) continue;
      const int target = best_label[comp[v]];
      if (target > 0) {
        labels[v] = target;
        ++assigned;
      }
    }
    if (assigned == 0) break;
    n_free -= assigned;
  }

  // Anything still unattached touches no surviving label; put it back.
  for (int v = 0; v < n; ++v) if (labels[v] == 0) labels[v] = original[v];

  return moved_total;
}

// Build a CSR view of the zero-padded neighbour matrix the R side already
// carries, so refine_boundaries_g3s_cpp() keeps its existing signature.
void g3s_csr_from_matrix(const IntegerMatrix& neighbor_indices, int n,
                         std::vector<int>& ptr, std::vector<int>& idx) {
  const int width = neighbor_indices.ncol();
  ptr.assign(static_cast<std::size_t>(n) + 1, 0);
  for (int i = 0; i < n; ++i) {
    int degree = 0;
    for (int k = 0; k < width; ++k) {
      const int nb = neighbor_indices(i, k) - 1;
      if (nb >= 0 && nb < n) ++degree;
    }
    ptr[i + 1] = degree;
  }
  for (int i = 0; i < n; ++i) ptr[i + 1] += ptr[i];
  idx.assign(static_cast<std::size_t>(ptr[n]), 0);
  std::vector<int> fill(ptr.begin(), ptr.begin() + n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < width; ++k) {
      const int nb = neighbor_indices(i, k) - 1;
      if (nb >= 0 && nb < n) idx[fill[i]++] = nb;
    }
  }
}

} // namespace

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector enforce_label_connectivity_cpp(IntegerVector labels,
                                             const NumericMatrix& feature_mat,
                                             const NumericMatrix& coords,
                                             const IntegerMatrix& neighbor_indices,
                                             double alpha,
                                             double compactness) {
  const int n = labels.size();
  if (feature_mat.ncol() != n || coords.nrow() != n || coords.ncol() != 3 ||
      neighbor_indices.nrow() != n) {
    stop("connectivity inputs must have one aligned row or column per voxel");
  }
  if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
    stop("alpha must be in [0, 1]");
  }
  if (!std::isfinite(compactness) || compactness <= 0.0) {
    stop("compactness must be positive and finite");
  }
  std::vector<int> ptr, idx;
  g3s_csr_from_matrix(neighbor_indices, n, ptr, idx);
  std::vector<int> work(labels.begin(), labels.end());
  g3s_repair_connectivity(work, ptr, idx, feature_mat, coords, alpha,
                          compactness);
  return wrap(work);
}

// -----------------------------------------------------------------------------
// 4. Boundary refinement (C++)
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector refine_boundaries_g3s_cpp(IntegerVector labels,
                                        const NumericMatrix& feature_mat,
                                        const NumericMatrix& coords,
                                        const IntegerMatrix& neighbor_indices,
                                        double alpha,
                                        double compactness,
                                        int max_iter,
                                        bool enforce_connectivity = true) {

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

  std::vector<int> ptr, idx;
  if (enforce_connectivity) g3s_csr_from_matrix(neighbor_indices, N, ptr, idx);

  for (int iter = 0; iter < max_iter; ++iter) {
    int max_label = 0;
    for (int l : cur) if (l > max_label) max_label = l;
    if (max_label <= 0) break;

    const int label_stride = max_label + 1;
    std::vector<int> counts(label_stride, 0);
    std::vector<double> centroids((size_t)label_stride * (size_t)M, 0.0);
    std::vector<double> coord_centroids((size_t)label_stride * 3U, 0.0);

    // Serial accumulation. A thread-local reduction would need
    // n_threads * label_stride * M doubles, which is hundreds of megabytes at
    // large K, and this pass is memory-bound rather than compute-bound.
    for (int i = 0; i < N; ++i) {
      const int lab = cur[i];
      if (lab <= 0) continue;
      counts[lab]++;
      double* row = centroids.data() + (size_t)lab * (size_t)M;
      for (int m = 0; m < M; ++m) row[m] += feature_mat(m, i);
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
        for (int m = 0; m < M; ++m) row[m] *= norm;
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
      // Always write next[i]: cur and next are swapped every iteration, so a
      // skipped slot would otherwise resurrect a label from two iterations ago.
      next[i] = lab;
      if (lab <= 0) continue;

      // At most one own label plus 26 distinct neighbour labels.
      int candidates[32];
      int n_cand = 0;
      candidates[n_cand++] = lab;

      for (int k = 0; k < K_neighbors; ++k) {
        int neigh = neighbor_indices(i, k) - 1;
        if (neigh < 0 || neigh >= N) continue;
        int neigh_lab = cur[neigh];
        if (neigh_lab <= 0 || neigh_lab == lab) continue;
        bool seen = false;
        for (int c = 0; c < n_cand; ++c) {
          if (candidates[c] == neigh_lab) { seen = true; break; }
        }
        if (!seen && n_cand < 32) candidates[n_cand++] = neigh_lab;
      }

      if (n_cand == 1) continue;  // interior voxel

      double best_cost = std::numeric_limits<double>::infinity();
      int best_lab = lab;

      for (int c = 0; c < n_cand; ++c) {
        const int cand = candidates[c];
        if (counts[cand] == 0) continue;
        const double* c_row = centroids.data() + (size_t)cand * (size_t)M;
        double feature_distance = 0.0;
        if (use_cosine) {
          double dot = 0.0;
          for (int m = 0; m < M; ++m) dot += feature_mat(m, i) * c_row[m];
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

    // Repair inside the loop, not just at the end: a fragmented cluster has a
    // centroid sitting between its pieces, which drives further fragmentation
    // on the next pass.
    if (enforce_connectivity) {
      g3s_repair_connectivity(cur, ptr, idx, feature_mat, coords, alpha,
                              compactness);
    }
    if (changed == 0) break;
  }

  return wrap(cur);
}

// -----------------------------------------------------------------------------
// 5. Greedy spatial-inhibition seed selection
//
// Walks candidates in priority order and keeps one only if no already-accepted
// seed lies within `radius`. A uniform spatial hash with cell size `radius`
// bounds each test to the 27 surrounding cells, so this is O(candidates)
// rather than the O(candidates * K) all-pairs scan. If fewer than K seeds fit,
// the radius decays and the pass repeats: that finds the largest separation
// that still admits K seeds instead of failing and padding the shortfall.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List g3s_select_seeds_cpp(const NumericMatrix& coords,
                          const IntegerVector& candidates,
                          int K,
                          double radius,
                          double decay,
                          int max_rounds) {
  const int n = coords.nrow();
  const int n_cand = candidates.size();
  if (coords.ncol() != 3) stop("coords must have three columns");
  if (K < 1) stop("K must be positive");
  if (!std::isfinite(radius) || radius < 0.0) stop("radius must be non-negative");
  if (!(decay > 0.0 && decay < 1.0)) stop("decay must be in (0, 1)");
  if (max_rounds < 0) stop("max_rounds must be non-negative");
  for (int c = 0; c < n_cand; ++c) {
    if (candidates[c] < 1 || candidates[c] > n) {
      stop("candidates are out of range");
    }
  }

  std::vector<int> selected;
  selected.reserve(K);
  double current = radius;

  for (int round = 0; ; ++round) {
    selected.clear();
    const double r2 = current * current;

    if (current <= 0.0) {
      const int take = std::min(K, n_cand);
      for (int c = 0; c < take; ++c) selected.push_back(candidates[c]);
      break;
    }

    // Bucket accepted seeds on a lattice of side `current`.
    std::unordered_map<long long, std::vector<int> > buckets;
    buckets.reserve(static_cast<std::size_t>(K) * 2u);
    const double inv = 1.0 / current;

    for (int c = 0; c < n_cand && static_cast<int>(selected.size()) < K; ++c) {
      const int v = candidates[c] - 1;
      const double x = coords(v, 0), y = coords(v, 1), z = coords(v, 2);
      const long long cx = static_cast<long long>(std::floor(x * inv));
      const long long cy = static_cast<long long>(std::floor(y * inv));
      const long long cz = static_cast<long long>(std::floor(z * inv));

      bool blocked = false;
      for (long long dx = -1; dx <= 1 && !blocked; ++dx) {
        for (long long dy = -1; dy <= 1 && !blocked; ++dy) {
          for (long long dz = -1; dz <= 1 && !blocked; ++dz) {
            const long long key = ((cx + dx) * 1000003LL + (cy + dy)) * 1000003LL
              + (cz + dz);
            std::unordered_map<long long, std::vector<int> >::const_iterator it =
              buckets.find(key);
            if (it == buckets.end()) continue;
            for (std::size_t s = 0; s < it->second.size(); ++s) {
              const int u = it->second[s];
              const double ddx = coords(u, 0) - x;
              const double ddy = coords(u, 1) - y;
              const double ddz = coords(u, 2) - z;
              if (ddx * ddx + ddy * ddy + ddz * ddz < r2) { blocked = true; break; }
            }
          }
        }
      }
      if (blocked) continue;

      selected.push_back(v + 1);
      const long long key = (cx * 1000003LL + cy) * 1000003LL + cz;
      buckets[key].push_back(v);
    }

    if (static_cast<int>(selected.size()) >= K || round >= max_rounds) break;
    current *= decay;
    if (current < std::sqrt(std::numeric_limits<double>::epsilon())) current = 0.0;
  }

  std::sort(selected.begin(), selected.end());
  return List::create(_["seeds"] = wrap(selected), _["radius"] = current);
}
