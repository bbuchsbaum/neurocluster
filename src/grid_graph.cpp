#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Masked-grid graph construction.
//
// build_grid_adjacency_cpp() in er_ml_rena.cpp returns an arma::sp_mat whose
// n_rows * n_cols must fit in an Armadillo uword, which caps the mask at 65535
// voxels unless ARMA_64BIT_WORD is set. Every consumer of that matrix in the
// exact-K engine immediately converts it back to an edge list, so build the
// edge list and the CSR adjacency directly instead. There is no N-by-N object
// and therefore no size ceiling beyond the edge count itself.
// -----------------------------------------------------------------------------

namespace {

// Grid offsets for 6-, 18-, and 26-neighbour connectivity.
void g3s_grid_offsets(int connectivity, std::vector<int>& ox,
                      std::vector<int>& oy, std::vector<int>& oz) {
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dz = -1; dz <= 1; ++dz) {
        const int manhattan = std::abs(dx) + std::abs(dy) + std::abs(dz);
        if (manhattan == 0) continue;
        bool keep;
        if (connectivity == 6) {
          keep = (manhattan == 1);
        } else if (connectivity == 18) {
          keep = (manhattan <= 2);
        } else {
          keep = true;
        }
        if (keep) {
          ox.push_back(dx);
          oy.push_back(dy);
          oz.push_back(dz);
        }
      }
    }
  }
}

// Weak components of a CSR graph, numbered in ascending first-vertex order so
// the labelling is canonical and matches .exact_k_relabel() of an igraph
// membership vector.
std::vector<int> g3s_csr_components(const std::vector<int>& ptr,
                                    const std::vector<int>& idx,
                                    int n) {
  std::vector<int> comp(n, 0);
  std::vector<int> stack;
  stack.reserve(n);
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
        if (comp[nb] == 0) {
          comp[nb] = next_id;
          stack.push_back(nb);
        }
      }
    }
  }
  return comp;
}

} // namespace

// [[Rcpp::export]]
List build_grid_edges_cpp(const IntegerVector& mask_idx,
                          const IntegerVector& dims,
                          const int connectivity) {
  if (dims.size() != 3) stop("build_grid_edges_cpp: dims must be length 3");
  if (connectivity != 6 && connectivity != 18 && connectivity != 26) {
    stop("build_grid_edges_cpp: connectivity must be 6, 18, or 26");
  }

  const int nx = dims[0];
  const int ny = dims[1];
  const int nz = dims[2];
  if (nx < 1 || ny < 1 || nz < 1) {
    stop("build_grid_edges_cpp: dims must be positive");
  }
  const std::size_t total =
    static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
    static_cast<std::size_t>(nz);

  const int n = mask_idx.size();
  if (n < 1) stop("build_grid_edges_cpp: mask_idx must be non-empty");

  // Node ids follow the ascending mask index order produced by which().
  std::vector<int> loc2node(total, -1);
  for (int v = 0; v < n; ++v) {
    if (v > 0 && mask_idx[v] <= mask_idx[v - 1]) {
      stop("build_grid_edges_cpp: mask_idx must be strictly increasing");
    }
    const long long lin = static_cast<long long>(mask_idx[v]) - 1;
    if (lin < 0 || static_cast<std::size_t>(lin) >= total) {
      stop("build_grid_edges_cpp: mask_idx out of range");
    }
    loc2node[static_cast<std::size_t>(lin)] = v;
  }

  std::vector<int> ox, oy, oz;
  ox.reserve(26); oy.reserve(26); oz.reserve(26);
  g3s_grid_offsets(connectivity, ox, oy, oz);
  const int n_off = static_cast<int>(ox.size());

  // Pass 1: degrees. Pass 2: fill. Keeps peak memory at one CSR copy.
  std::vector<int> ptr(static_cast<std::size_t>(n) + 1, 0);
  std::vector<int> coord_x(n), coord_y(n), coord_z(n);
  for (int v = 0; v < n; ++v) {
    long long tmp = static_cast<long long>(mask_idx[v]) - 1;
    const int x = static_cast<int>(tmp % nx); tmp /= nx;
    const int y = static_cast<int>(tmp % ny);
    const int z = static_cast<int>(tmp / ny);
    coord_x[v] = x; coord_y[v] = y; coord_z[v] = z;
    int degree = 0;
    for (int k = 0; k < n_off; ++k) {
      const int xx = x + ox[k], yy = y + oy[k], zz = z + oz[k];
      if (xx < 0 || xx >= nx || yy < 0 || yy >= ny || zz < 0 || zz >= nz) continue;
      const std::size_t lin = static_cast<std::size_t>(xx) +
        static_cast<std::size_t>(nx) * (static_cast<std::size_t>(yy) +
        static_cast<std::size_t>(ny) * static_cast<std::size_t>(zz));
      if (loc2node[lin] >= 0) ++degree;
    }
    ptr[static_cast<std::size_t>(v) + 1] = degree;
  }
  for (int v = 0; v < n; ++v) ptr[v + 1] += ptr[v];

  const int nnz = ptr[n];
  std::vector<int> idx(static_cast<std::size_t>(nnz));
  std::vector<int> fill(ptr.begin(), ptr.begin() + n);
  for (int v = 0; v < n; ++v) {
    const int x = coord_x[v], y = coord_y[v], z = coord_z[v];
    for (int k = 0; k < n_off; ++k) {
      const int xx = x + ox[k], yy = y + oy[k], zz = z + oz[k];
      if (xx < 0 || xx >= nx || yy < 0 || yy >= ny || zz < 0 || zz >= nz) continue;
      const std::size_t lin = static_cast<std::size_t>(xx) +
        static_cast<std::size_t>(nx) * (static_cast<std::size_t>(yy) +
        static_cast<std::size_t>(ny) * static_cast<std::size_t>(zz));
      const int node = loc2node[lin];
      if (node >= 0) idx[fill[v]++] = node;
    }
    // Ascending neighbour order makes the CSR slices and the edge list
    // canonical, which keeps every downstream tie-break deterministic.
    std::sort(idx.begin() + ptr[v], idx.begin() + ptr[v + 1]);
  }

  // Upper-triangular edges, emitted in (i, j) order.
  IntegerMatrix edges(nnz / 2, 2);
  int e = 0;
  for (int v = 0; v < n; ++v) {
    for (int p = ptr[v]; p < ptr[v + 1]; ++p) {
      const int u = idx[p];
      if (u > v) {
        edges(e, 0) = v + 1;
        edges(e, 1) = u + 1;
        ++e;
      }
    }
  }
  if (e != nnz / 2) stop("build_grid_edges_cpp: adjacency is not symmetric");

  const std::vector<int> comp = g3s_csr_components(ptr, idx, n);

  IntegerVector ptr_out(ptr.begin(), ptr.end());
  IntegerVector idx_out(static_cast<std::size_t>(nnz));
  for (int p = 0; p < nnz; ++p) idx_out[p] = idx[p] + 1;
  IntegerVector comp_out(comp.begin(), comp.end());

  return List::create(
    _["edges"] = edges,
    _["neighbor_ptr"] = ptr_out,
    _["neighbor_idx"] = idx_out,
    _["components"] = comp_out,
    _["n_voxels"] = n
  );
}

// -----------------------------------------------------------------------------
// Connected components of the subgraph induced by equal labels, from an edge
// list. Union-find keeps this O(E alpha(N)) and avoids constructing an igraph
// object for what is a single pass over the edges.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector edge_label_components_cpp(const IntegerVector& labels,
                                        const IntegerMatrix& edges) {
  const int n = labels.size();
  if (n < 1) stop("edge_label_components_cpp: labels must be non-empty");
  if (edges.ncol() != 2 && edges.nrow() != 0) {
    stop("edge_label_components_cpp: edges must have two columns");
  }

  std::vector<int> parent(n);
  std::vector<int> rank(n, 0);
  for (int i = 0; i < n; ++i) parent[i] = i;

  // Iterative find with path compression.
  auto find = [&parent](int x) {
    int root = x;
    while (parent[root] != root) root = parent[root];
    while (parent[x] != root) {
      const int next = parent[x];
      parent[x] = root;
      x = next;
    }
    return root;
  };

  const int m = edges.nrow();
  for (int e = 0; e < m; ++e) {
    const int a = edges(e, 0) - 1;
    const int b = edges(e, 1) - 1;
    if (a < 0 || a >= n || b < 0 || b >= n) {
      stop("edge_label_components_cpp: edge endpoint out of range");
    }
    if (labels[a] != labels[b]) continue;
    const int ra = find(a), rb = find(b);
    if (ra == rb) continue;
    if (rank[ra] < rank[rb]) parent[ra] = rb;
    else if (rank[ra] > rank[rb]) parent[rb] = ra;
    else { parent[rb] = ra; ++rank[ra]; }
  }

  // Number components in ascending first-vertex order.
  IntegerVector out(n);
  std::vector<int> seen(n, 0);
  int next_id = 0;
  for (int i = 0; i < n; ++i) {
    const int root = find(i);
    if (seen[root] == 0) seen[root] = ++next_id;
    out[i] = seen[root];
  }
  return out;
}
