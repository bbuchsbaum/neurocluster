// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ---------- helpers ---------------------------------------------------------

// Compute squared L2 distance between two rows of a matrix without creating
// intermediate rowvec temporaries. This avoids large hidden allocations in
// the hot path of ReNA.
inline double squared_l2_rows(const arma::mat &M, arma::uword i, arma::uword j) {
  const arma::uword n = M.n_cols;
  double acc = 0.0;
  for (arma::uword k = 0; k < n; ++k) {
    double d = M(i, k) - M(j, k);
    acc += d * d;
  }
  return acc;
}

// Build symmetric grid adjacency for masked voxels entirely in C++ to avoid
// repeated R allocations and vector growth.
// [[Rcpp::export]]
arma::sp_mat build_grid_adjacency_cpp(const IntegerVector &mask_idx,
                                      const IntegerVector &dims,
                                      const int connectivity) {
  if (dims.size() != 3)
    stop("build_grid_adjacency_cpp: dims must be length 3");

  const int nx = dims[0];
  const int ny = dims[1];
  const int nz = dims[2];

  const int n_vox = mask_idx.size();
  const int total = nx * ny * nz;

  // map from linear index (0-based) to node id (0..n_vox-1), -1 if not masked
  std::vector<int> loc2node(total, -1);
  for (int v = 0; v < n_vox; ++v) {
    int lin = mask_idx[v] - 1; // R is 1-based
    if (lin < 0 || lin >= total)
      stop("build_grid_adjacency_cpp: mask_idx out of range");
    loc2node[lin] = v;
  }

  // Precompute neighbor offsets in (dx, dy, dz)
  std::vector<int> dx, dy, dz;
  dx.reserve(26);
  dy.reserve(26);
  dz.reserve(26);

  for (int ddx = -1; ddx <= 1; ++ddx) {
    for (int ddy = -1; ddy <= 1; ++ddy) {
      for (int ddz = -1; ddz <= 1; ++ddz) {
        if (ddx == 0 && ddy == 0 && ddz == 0) continue;

        int manhattan = std::abs(ddx) + std::abs(ddy) + std::abs(ddz);
        bool keep = false;

        if (connectivity == 6) {
          keep = (manhattan == 1);
        } else if (connectivity == 18) {
          keep = (manhattan > 0 && manhattan <= 2);
        } else if (connectivity == 26) {
          keep = true; // all non-zero offsets in {-1,0,1}^3
        } else {
          stop("build_grid_adjacency_cpp: unsupported connectivity");
        }

        if (keep) {
          dx.push_back(ddx);
          dy.push_back(ddy);
          dz.push_back(ddz);
        }
      }
    }
  }

  std::vector<uword> ii;
  std::vector<uword> jj;
  ii.reserve(static_cast<size_t>(n_vox) * dx.size());
  jj.reserve(static_cast<size_t>(n_vox) * dx.size());

  for (int v = 0; v < n_vox; ++v) {
    int lin = mask_idx[v] - 1;

    int tmp = lin;
    int x = tmp % nx;      tmp /= nx;
    int y = tmp % ny;
    int z = tmp / ny;

    for (size_t k = 0; k < dx.size(); ++k) {
      int xx = x + dx[k];
      int yy = y + dy[k];
      int zz = z + dz[k];

      if (xx < 0 || xx >= nx ||
          yy < 0 || yy >= ny ||
          zz < 0 || zz >= nz) continue;

      int n_lin = xx + nx * (yy + ny * zz);
      int nb = loc2node[n_lin];
      if (nb >= 0 && nb > v) { // upper triangle only
        ii.push_back(static_cast<uword>(v));
        jj.push_back(static_cast<uword>(nb));
      }
    }
  }

  arma::sp_mat G(n_vox, n_vox);

  if (!ii.empty()) {
    const uword m = ii.size();
    umat locations(2, 2 * m);
    vec  values(2 * m, fill::ones);

    for (uword k = 0; k < m; ++k) {
      locations(0, k)      = ii[k];
      locations(1, k)      = jj[k];
      locations(0, k + m)  = jj[k];
      locations(1, k + m)  = ii[k];
    }

    G = sp_mat(locations, values, n_vox, n_vox);
  }

  return G;
}

struct UnionFind {
  std::vector<int> parent;
  std::vector<int> rank;

  UnionFind() {}
  explicit UnionFind(int n) {
    parent.resize(n);
    rank.assign(n, 0);
    for (int i = 0; i < n; ++i) parent[i] = i;
  }

  int find(int x) {
    if (parent[x] != x) parent[x] = find(parent[x]);
    return parent[x];
  }

  void unite(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y) return;
    if (rank[x] < rank[y]) std::swap(x, y);
    parent[y] = x;
    if (rank[x] == rank[y]) rank[x]++;
  }
};

// ---------- single RNN step (internal) -------------------------------------

IntegerVector rena_rnn_step_internal(const arma::mat &X,
                                     const arma::sp_mat &G,
                                     const vec &grad,
                                     double lambda) {
  uword p = X.n_rows;
  IntegerVector labels(p);

  std::vector<int> nearest_idx(p);

  bool use_grad = (grad.n_elem == p && lambda != 0.0);
  double inf = std::numeric_limits<double>::infinity();

  // Parallel nearest-neighbor search per node
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (uword i = 0; i < p; ++i) {
    double best_d = inf;
    int best_j = static_cast<int>(i);

    const double g_i = use_grad ? grad(i) : 0.0; // avoid repeated loads

    for (sp_mat::const_row_iterator it = G.begin_row(i);
         it != G.end_row(i); ++it) {
      uword j = it.col();
      if (i == j) continue;
      double d = squared_l2_rows(X, i, j);
      if (use_grad) {
        double gdiff = std::abs(g_i - grad(j));
        d *= (1.0 + lambda * gdiff);
      }
      if (d < best_d) {
        best_d = d;
        best_j = static_cast<int>(j);
      }
    }
    nearest_idx[i] = best_j;
  }

  // Reciprocal merges + fallback 1-NN for isolated nodes
  UnionFind uf(p);
  std::vector<char> in_rnn(p, 0);

  for (uword i = 0; i < p; ++i) {
    int j = nearest_idx[i];
    if (j < 0 || j >= static_cast<int>(p)) continue;
    if (nearest_idx[j] == static_cast<int>(i) && static_cast<int>(i) < j) {
      uf.unite(i, j);
      in_rnn[i] = 1;
      in_rnn[j] = 1;
    }
  }

  for (uword i = 0; i < p; ++i) {
    if (!in_rnn[i]) {
      int j = nearest_idx[i];
      if (j >= 0 && j < static_cast<int>(p) && j != static_cast<int>(i)) {
        uf.unite(i, j);
      }
    }
  }

  // Contiguous component labels
  std::unordered_map<int, int> root2id;
  root2id.reserve(p);
  int next_id = 0;
  for (uword i = 0; i < p; ++i) {
    int root = uf.find(i);
    auto it = root2id.find(root);
    if (it == root2id.end()) {
      root2id[root] = next_id;
      labels[i] = next_id;
      ++next_id;
    } else {
      labels[i] = it->second;
    }
  }

  return labels; // 0-based
}

// ---------- Stage 1: edge-aware RNN ReNA coarsening ------------------------

// [[Rcpp::export]]
List rena_rnn_coarse_cpp(const arma::mat &X,
                         const arma::sp_mat &G,
                         const NumericVector &grad_img,
                         int stop_at,
                         double lambda = 0.0,
                         int max_iter = 100) {
  int p_orig = X.n_rows;
  int n = X.n_cols;

  mat X_curr = X;
  sp_mat G_curr = G;
  vec grad_curr;

  bool use_grad = (grad_img.size() == p_orig && lambda != 0.0);
  if (use_grad) {
    grad_curr = as<vec>(grad_img);
  } else {
    grad_curr.set_size(0);
  }

  IntegerVector labels_global(p_orig);
  for (int i = 0; i < p_orig; ++i) labels_global[i] = i;

  vec weights_curr(X_curr.n_rows, fill::ones);

  int iter = 0;

  while ((int)X_curr.n_rows > stop_at && iter < max_iter) {
    ++iter;
    int p_curr = X_curr.n_rows;

    if (use_grad && (int)grad_curr.n_elem != p_curr) {
      stop("grad_curr length mismatch");
    }

    IntegerVector labels_step = rena_rnn_step_internal(X_curr, G_curr, grad_curr, lambda);

    int q = 0;
    for (int i = 0; i < p_curr; ++i) if (labels_step[i] + 1 > q) q = labels_step[i] + 1;
    if (q == p_curr) break; // no merges

    for (int i = 0; i < p_orig; ++i) {
      int old_cluster = labels_global[i];
      if (old_cluster < 0 || old_cluster >= p_curr) stop("labels_global out of range");
      labels_global[i] = labels_step[old_cluster];
    }

    mat X_next(q, n, fill::zeros);
    vec weights_next(q, fill::zeros);
    vec grad_next;
    if (use_grad) grad_next.zeros(q);

    for (int i = 0; i < p_curr; ++i) {
      int c = labels_step[i];
      double w = weights_curr(i);
      X_next.row(c)     += w * X_curr.row(i);
      weights_next(c)   += w;
      if (use_grad) grad_next(c) += w * grad_curr(i);
    }

    for (int c = 0; c < q; ++c) {
      if (weights_next(c) > 0) {
        X_next.row(c) /= weights_next(c);
        if (use_grad) grad_next(c) /= weights_next(c);
      }
    }

    std::vector<uword> ii;
    std::vector<uword> jj;
    ii.reserve(G_curr.n_nonzero);
    jj.reserve(G_curr.n_nonzero);

    // use hash set to deduplicate edges
    std::unordered_set<uint64_t> edge_set;
    edge_set.reserve(G_curr.n_nonzero * 2);

    for (sp_mat::const_iterator it = G_curr.begin(); it != G_curr.end(); ++it) {
      uword i = it.row();
      uword j = it.col();
      if (i == j) continue;
      int ci = labels_step[i];
      int cj = labels_step[j];
      if (ci == cj) continue;
      uint64_t key = (static_cast<uint64_t>(ci) << 32) | static_cast<uint32_t>(cj);
      if (edge_set.insert(key).second) {
        ii.push_back(ci);
        jj.push_back(cj);
      }
    }

    sp_mat G_tmp;
    if (!ii.empty()) {
      uword nnz = ii.size();
      umat locations(2, nnz);
      for (uword k = 0; k < nnz; ++k) {
        locations(0, k) = ii[k];
        locations(1, k) = jj[k];
      }
      vec vals(nnz, fill::ones);
      G_tmp = sp_mat(locations, vals, q, q);
      G_tmp = spones(symmatu(G_tmp));
    } else {
      G_tmp = sp_mat(q, q);
    }

    X_curr = std::move(X_next);
    G_curr = std::move(G_tmp);
    weights_curr = std::move(weights_next);
    if (use_grad) grad_curr = std::move(grad_next);
  }

  int K_prime = X_curr.n_rows;

  IntegerVector labels_out = labels_global; // still 0-based
  IntegerVector sizes(K_prime);
  for (int i = 0; i < p_orig; ++i) {
    int c = labels_out[i];
    if (c < 0 || c >= K_prime) stop("labels_out cluster index out of range");
    sizes[c] += 1;
  }

  for (int i = 0; i < p_orig; ++i) labels_out[i] += 1;

  return List::create(
    Named("X_coarse")      = X_curr,
    Named("G_coarse")      = G_curr,
    Named("labels_coarse") = labels_out,
    Named("sizes_coarse")  = sizes
  );
}

// ---------- Stage 2: adjacency-constrained Ward ----------------------------

struct Edge {
  double cost;
  int a;
  int b;
};

struct EdgeCompare {
  bool operator()(Edge const& e1, Edge const& e2) const {
    return e1.cost > e2.cost; // min-heap
  }
};

// [[Rcpp::export]]
List ward_on_supervoxels_cpp(const arma::mat &X_coarse,
                             const arma::sp_mat &G_coarse,
                             const IntegerVector &sizes,
                             int n_clusters) {
  int N = X_coarse.n_rows;
  int n = X_coarse.n_cols;

  if (sizes.size() != N) stop("sizes length must match rows of X_coarse");
  if (n_clusters <= 0 || n_clusters > N) stop("n_clusters out of range");

  mat means = X_coarse;
  std::vector<double> sz(N);
  for (int i = 0; i < N; ++i) sz[i] = std::max(1, sizes[i]);

  std::vector<bool> active(N, true);
  int active_count = N;

  std::vector<std::vector<int>> neighbors(N);
  for (sp_mat::const_iterator it = G_coarse.begin(); it != G_coarse.end(); ++it) {
    int i = it.row();
    int j = it.col();
    if (i == j) continue;
    neighbors[i].push_back(j);
    neighbors[j].push_back(i);
  }

  UnionFind uf(N);
  std::priority_queue<Edge, std::vector<Edge>, EdgeCompare> pq;

  auto push_edge = [&](int a, int b) {
    if (a == b || !active[a] || !active[b]) return;
    double s_a = sz[a], s_b = sz[b];
    double denom = s_a + s_b;
    if (denom <= 0) return;
    double factor = (s_a * s_b) / denom;
    double d2 = squared_l2_rows(means, a, b);
    double cost = factor * d2;
    pq.push(Edge{cost, a, b});
  };

  for (int i = 0; i < N; ++i) {
    for (int j : neighbors[i]) {
      if (i < j) push_edge(i, j);
    }
  }

  while (active_count > n_clusters && !pq.empty()) {
    Edge e = pq.top();
    pq.pop();
    int a = e.a, b = e.b;
    if (a < 0 || a >= N || b < 0 || b >= N) continue;
    if (!active[a] || !active[b]) continue;

    int ra = uf.find(a), rb = uf.find(b);
    if (ra == rb) continue;

    double s_a = sz[a], s_b = sz[b];
    double s_new = s_a + s_b;
    if (s_new <= 0) continue;

    means.row(a) = (s_a * means.row(a) + s_b * means.row(b)) / s_new;
    sz[a] = s_new;

    active[b] = false;
    sz[b] = 0.0;
    active_count--;
    uf.unite(a, b);

    // Merge neighbor lists with deduplication to avoid growth
    std::vector<int> merged;
    merged.reserve(neighbors[a].size() + neighbors[b].size());
    for (int d : neighbors[a]) if (active[d] && d != a) merged.push_back(d);
    for (int d : neighbors[b]) if (active[d] && d != a) merged.push_back(d);
    std::sort(merged.begin(), merged.end());
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

    neighbors[a].swap(merged);
    neighbors[b].clear();

    // Keep adjacency symmetric and push updated edges
    for (int d : neighbors[a]) {
      neighbors[d].push_back(a);
      push_edge(a, d);
    }
  }

  IntegerVector labels_super(N);
  std::unordered_map<int, int> root2id;
  root2id.reserve(N);
  int next_id = 0;
  for (int i = 0; i < N; ++i) {
    int r = uf.find(i);
    auto it = root2id.find(r);
    if (it == root2id.end()) {
      root2id[r] = next_id;
      labels_super[i] = next_id;
      ++next_id;
    } else {
      labels_super[i] = it->second;
    }
  }

  for (int i = 0; i < N; ++i) labels_super[i] += 1;

  return List::create(
    Named("labels_super") = labels_super,
    Named("n_clusters") = next_id
  );
}
