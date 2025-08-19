// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------

inline int idx3d(int x, int y, int z, int nx, int ny) {
  return x + nx * (y + ny * z);
}

struct Edge {
  int a;
  int b;
  uint16_t wq; // quantized for counting sort
  float w;     // float distance in [0,2]
};

inline void counting_sort_edges(std::vector<Edge> &edges) {
  const int K = 65536;
  std::vector<int> cnt(K, 0);
  for (const auto &e : edges) cnt[e.wq]++;
  int acc = 0;
  for (int i=0;i<K;++i) { int c = cnt[i]; cnt[i] = acc; acc += c; }
  std::vector<Edge> out(edges.size());
  for (const auto &e : edges) out[cnt[e.wq]++] = e;
  edges.swap(out);
}

inline uint16_t quantize16(double d) {
  if (d < 0.0) d = 0.0;
  if (d > 2.0) d = 2.0;
  return static_cast<uint16_t>(std::lround((d/2.0)*65535.0));
}

// ---------------------------------------------------------------
// Union–Find with internal difference
// ---------------------------------------------------------------

struct UF {
  std::vector<int> parent, rnk, sz;
  std::vector<float> intr;
  UF() {}
  UF(int n) { init(n); }
  void init(int n) {
    parent.resize(n);
    rnk.assign(n,0);
    sz.assign(n,1);
    intr.assign(n,0.0f);
    for (int i=0;i<n;++i) parent[i]=i;
  }
  int find(int x) {
    while (x != parent[x]) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
  }
  int unite(int a, int b, float e) {
    a = find(a); b = find(b);
    if (a == b) return a;
    if (rnk[a] < rnk[b]) std::swap(a,b);
    parent[b] = a;
    sz[a] += sz[b];
    intr[a] = std::max({intr[a], intr[b], e});
    if (rnk[a] == rnk[b]) rnk[a]++;
    return a;
  }
};

// ---------------------------------------------------------------
// Time-series preprocessing
// ---------------------------------------------------------------

inline void detrend_and_zscore(const double* y, int T, std::vector<double> &z) {
  z.assign(T,0.0);
  double mean_t = 0.5 * (T - 1);
  double var_t  = (T*(double)T - 1.0)/12.0;
  double sumy=0.0, sumty=0.0;
  for (int t=0;t<T;++t){ sumy+=y[t]; sumty += t*y[t]; }
  double meany = sumy / T;
  double cov_ty = (sumty / T) - mean_t*meany;
  double beta = cov_ty / var_t;
  double alpha = meany - beta*mean_t;

  double m=0.0, s2=0.0;
  for (int t=0;t<T;++t){ double yt = y[t] - (alpha + beta*t); z[t]=yt; m+=yt; }
  m /= T;
  for (int t=0;t<T;++t){ double c = z[t]-m; s2 += c*c; }
  s2 = (s2<=0.0) ? 1.0 : s2;
  double sd = std::sqrt(s2 / std::max(1,T-1));
  if (sd <= 0.0) sd = 1.0;
  for (int t=0;t<T;++t) z[t] = (z[t]-m)/sd;
}

inline double split_half_corr(const std::vector<double> &z) {
  int T = (int)z.size();
  int n = T/2;
  if (n < 2) return 0.0;
  double sa=0,sb=0,saa=0,sbb=0,sab=0;
  for (int i=0;i<n;++i) {
    double a=z[2*i], b=z[2*i+1];
    sa+=a; sb+=b; saa+=a*a; sbb+=b*b; sab+=a*b;
  }
  double ma=sa/n, mb=sb/n;
  double va=saa - n*ma*ma;
  double vb=sbb - n*mb*mb;
  double cab=sab - n*ma*mb;
  if (va<=0.0 || vb<=0.0) return 0.0;
  double r = cab / std::sqrt(va*vb);
  if (!std::isfinite(r)) r = 0.0;
  if (r>1.0) r=1.0; if (r<-1.0) r=-1.0;
  return r;
}

inline void build_dct_basis(int T, int r, std::vector<double> &phi) {
  phi.assign(T*r, 0.0);
  const double s  = std::sqrt(2.0 / (double)T);
  const double pi = std::acos(-1.0);
  for (int t=0; t<T; ++t) {
    for (int k=0; k<r; ++k) {
      phi[t*r + k] = s * std::cos(pi * ((t + 0.5) * (k + 1) / (double)T));
    }
  }
}

inline double dot_col(const NumericMatrix &U, int i, int j, int r) {
  const double *pi = &U(0,i), *pj = &U(0,j);
  double s=0.0; for (int k=0;k<r;++k) s += pi[k]*pj[k];
  return s;
}

// ---------------------------------------------------------------
// Sketch worker (slice-parallel)
// ---------------------------------------------------------------

struct SliceSketchWorker : public Worker {
  const RMatrix<double> TS;
  const RVector<int> mask;
  const std::vector<double> phi;
  const int nx, ny, nz, T, r;
  const bool rows_are_time;
  const double gamma;
  NumericMatrix &U; // r x N
  NumericVector &W; // N

  SliceSketchWorker(const NumericMatrix &TS_, const IntegerVector &mask_,
                    const std::vector<double> &phi_,
                    int nx_, int ny_, int nz_, int T_, int r_, bool rows_are_time_,
                    double gamma_, NumericMatrix &U_, NumericVector &W_)
    : TS(TS_), mask(mask_), phi(phi_), nx(nx_), ny(ny_), nz(nz_), T(T_), r(r_),
      rows_are_time(rows_are_time_), gamma(gamma_), U(U_), W(W_) {}

  inline void get_ts(int v, std::vector<double> &buf) const {
    buf.assign(T,0.0);
    if (rows_are_time) {
      // TS: T x N -> column v is contiguous (column-major)
      const double* col = &TS(0, v);
      for (int t=0; t<T; ++t) buf[t] = col[t];
    } else {
      // TS: N x T -> row v is NON-contiguous; access element-wise
      for (int t=0; t<T; ++t) buf[t] = TS(v, t);
    }
  }

  void operator()(std::size_t z0, std::size_t z1) {
    std::vector<double> ts(T), zsig(T);
    for (std::size_t z=z0; z<z1; ++z) {
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
        int g = idx3d(x,y,z,nx,ny);
        if (!mask[g]) continue;
        get_ts(g, ts);
        detrend_and_zscore(ts.data(), T, zsig);
        double rrel = split_half_corr(zsig);
        double w = std::pow(std::max(0.0, rrel), gamma);
        W[g] = w;

        // DCT sketch (k=1..r), L2-normalize
        double norm2=0.0;
        for (int k=0;k<r;++k) {
          double acc=0.0;
          const double* ph = &phi[k]; // phi[t*r+k]
          for (int t=0;t<T;++t) acc += zsig[t] * ph[t*r];
          U(k,g)=acc; norm2 += acc*acc;
        }
        norm2 = std::sqrt(std::max(1e-12, norm2));
        for (int k=0;k<r;++k) U(k,g) /= norm2;
      }
    }
  }
};

// ---------------------------------------------------------------
// Build in-slice edges with reliability-aware distance
// d_ij = 1 - (w_i w_j) * (U_i . U_j)  (clamped to [0,2])
// Optional mild spatial penalty for anisotropic voxels.
// ---------------------------------------------------------------

inline void build_slice_edges(
  int z,
  const IntegerVector &mask, const NumericMatrix &U, const NumericVector &W,
  int nx, int ny, int r, int nbhd,
  const NumericVector &voxdim, double spatial_beta,
  std::vector<int> &gids, std::vector<int> &loc_of_plane, std::vector<Edge> &edges
) {
  gids.clear(); edges.clear(); loc_of_plane.assign(nx*ny, -1);

  for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
    int g = idx3d(x,y,z,nx,ny);
    if (mask[g]) { int plane = x + nx*y; loc_of_plane[plane] = (int)gids.size(); gids.push_back(g); }
  }
  if (gids.empty()) return;

  std::vector<std::pair<int,int>> offs;
  offs.emplace_back(1,0); offs.emplace_back(0,1);        // 4-conn
  if (nbhd==8) { offs.emplace_back(1,1); offs.emplace_back(1,-1); }

  const double dx = voxdim.size()>=1 ? (double)voxdim[0] : 1.0;
  const double dy = voxdim.size()>=2 ? (double)voxdim[1] : 1.0;
  const double minlen = std::min(dx,dy);

  for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
    int i = loc_of_plane[x + nx*y];
    if (i < 0) continue;
    int gi = gids[i];
    for (auto &p : offs) {
      int xn=x+p.first, yn=y+p.second;
      if (xn<0||xn>=nx||yn<0||yn>=ny) continue;
      int j = loc_of_plane[xn + nx*yn];
      if (j<0) continue;
      int gj = gids[j];

      double wprod = W[gi]*W[gj];
      double sim   = (wprod>0.0) ? dot_col(U, gi, gj, r) : 0.0;
      double d     = 1.0 - (wprod * sim);               // safer than dividing by wprod
      d = std::max(0.0, std::min(2.0, d));

      // mild spatial penalty for longer diagonals (optional)
      double len = std::sqrt((p.first?dx*dx:0.0) + (p.second?dy*dy:0.0));
      double factor = 1.0 + spatial_beta*(len/minlen - 1.0);
      d *= factor;

      Edge e; e.a=i; e.b=j; e.w=(float)d; e.wq=quantize16(d);
      edges.push_back(e);
    }
  }
}

// ---------------------------------------------------------------
// FH segmentation with canonical small-component cleanup
// ---------------------------------------------------------------

inline void segment_slice_fh(int ns, std::vector<Edge> &edges, double fh_scale,
                             int min_size, std::vector<int> &labels_out) {
  labels_out.assign(ns, -1);
  if (ns==0) return;

  counting_sort_edges(edges);
  UF uf(ns);

  // Kruskal with FH predicate
  for (const auto &e : edges) {
    int a=uf.find(e.a), b=uf.find(e.b);
    if (a==b) continue;
    double Ta = uf.intr[a] + fh_scale / uf.sz[a];
    double Tb = uf.intr[b] + fh_scale / uf.sz[b];
    if (e.w <= Ta && e.w <= Tb) uf.unite(a,b,e.w);
  }

  // Canonical cleanup: keep merging undersized components to their best neighbor until none remain
  {
    bool changed = true;
    int iter = 0, max_iter = 32;
    while (changed && iter++ < max_iter) {
      changed = false;
      std::vector<float> bestW(ns, std::numeric_limits<float>::infinity());
      std::vector<int>   bestN(ns, -1);
      for (const auto &e : edges) {
        int a = uf.find(e.a), b = uf.find(e.b);
        if (a == b) continue;
        if (uf.sz[a] < min_size && e.w < bestW[a]) { bestW[a] = e.w; bestN[a] = b; }
        if (uf.sz[b] < min_size && e.w < bestW[b]) { bestW[b] = e.w; bestN[b] = a; }
      }
      for (int v = 0; v < ns; ++v) {
        int a = uf.find(v);
        if (a != v) continue; // only roots once
        if (uf.sz[a] < min_size && bestN[a] >= 0) {
          uf.unite(a, bestN[a], bestW[a]);
          changed = true;
        }
      }
    }
  }

  // Relabel consecutive
  std::unordered_map<int,int> root2lab;
  int lab=1;
  for (int i=0;i<ns;++i) {
    int r = uf.find(i);
    auto it = root2lab.find(r);
    if (it==root2lab.end()) { root2lab.emplace(r, lab); labels_out[i]=lab; lab++; }
    else labels_out[i]=it->second;
  }
}

// ---------------------------------------------------------------
// Region Adjacency Graph (RAG) agglomeration to exact K
// Keeps contiguity: merges only touching components.
// Similarity = cosine(mean sketch A, mean sketch B).
// ---------------------------------------------------------------

struct PairHash {
  size_t operator()(const std::pair<int,int> &p) const noexcept {
    return ( (size_t)p.first << 32 ) ^ (size_t)p.second;
  }
};

struct PQItem {
  int a, b;
  float d; // distance = 1 - cos
  bool operator<(const PQItem &o) const { return d > o.d; } // min-heap
};

// compute cosine( meanA, meanB )
inline double cos_mean(int a, int b, const NumericMatrix &sumU, const std::vector<int> &sz) {
  if (a==b) return 1.0;
  int r = sumU.nrow();
  double na=0.0, nb=0.0, dp=0.0;
  double invA = 1.0 / std::max(1, sz[a]);
  double invB = 1.0 / std::max(1, sz[b]);
  for (int k=0;k<r;++k) {
    double ma = sumU(k,a)*invA;
    double mb = sumU(k,b)*invB;
    dp += ma*mb;
    na += ma*ma;
    nb += mb*mb;
  }
  if (na<=0.0 || nb<=0.0) return 0.0;
  return dp / (std::sqrt(na)*std::sqrt(nb));
}

// Build adjacency pairs (touching) from voxel->group mapping.
// in_slice: 4/8 neighbors within a slice; if allow_vertical==true, also (x,y,z)-(x,y,z+1).
inline void build_adjacency_pairs(
  const IntegerVector &mask, int nx, int ny, int nz,
  const std::vector<int> &voxel2group, int nbhd, bool allow_vertical,
  std::unordered_set<std::pair<int,int>, PairHash> &pairs_out
) {
  pairs_out.clear();
  auto reg = [&](int a, int b){
    if (a<0 || b<0 || a==b) return;
    int u = std::min(a,b), v = std::max(a,b);
    pairs_out.emplace(u,v);
  };

  // in-slice neighbors
  for (int z=0; z<nz; ++z) {
    for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
      int g = idx3d(x,y,z,nx,ny);
      if (!mask[g]) continue;
      int A = voxel2group[g];
      // E,S
      if (x+1<nx) { int gn=idx3d(x+1,y,z,nx,ny); if(mask[gn]) reg(A, voxel2group[gn]); }
      if (y+1<ny) { int gn=idx3d(x,y+1,z,nx,ny); if(mask[gn]) reg(A, voxel2group[gn]); }
      if (nbhd==8) {
        if (x+1<nx && y+1<ny){ int gn=idx3d(x+1,y+1,z,nx,ny); if(mask[gn]) reg(A, voxel2group[gn]); }
        if (x+1<nx && y-1>=0){ int gn=idx3d(x+1,y-1,z,nx,ny); if(mask[gn]) reg(A, voxel2group[gn]); }
      }
    }
  }
  // vertical neighbors
  if (allow_vertical) {
    for (int z=0; z<nz-1; ++z) {
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
        int g0 = idx3d(x,y,z,  nx,ny);
        int g1 = idx3d(x,y,z+1,nx,ny);
        if (!mask[g0] || !mask[g1]) continue;
        reg(voxel2group[g0], voxel2group[g1]);
      }
    }
  }
}

// Agglomerate to exact target K (global).
// Inputs: initial groups 0..G-1, with size[] and sumU(:,g), and adjacency pairs.
// Returns a mapping group_root -> final_group_root (via uf).
inline void rag_agglomerate_to_K_global(
  int G, int targetK, const std::unordered_set<std::pair<int,int>, PairHash> &pairs,
  NumericMatrix &sumU, std::vector<int> &sz, UF &uf
) {
  if (targetK <= 0) return;
  uf.init(G);

  // Active roots set
  std::vector<char> active(G, 0);
  int nActive = 0;
  for (int g=0; g<G; ++g) { active[g]=1; nActive++; }

  // Neighbor sets (for updates)
  std::vector< std::unordered_set<int> > nbr(G);
  for (const auto &p : pairs) { nbr[p.first].insert(p.second); nbr[p.second].insert(p.first); }

  // Priority queue of edges by distance
  std::priority_queue<PQItem> pq;
  for (const auto &p : pairs) {
    double c = cos_mean(p.first, p.second, sumU, sz);
    float d = (float)(1.0 - std::max(-1.0, std::min(1.0, c)));
    pq.push({p.first, p.second, d});
  }

  while (nActive > targetK && !pq.empty()) {
    PQItem it = pq.top(); pq.pop();
    int a = uf.find(it.a), b = uf.find(it.b);
    if (a==b || !active[a] || !active[b]) continue;

    // Check still neighbors
    if (!nbr[a].count(b)) continue;

    // Recompute distance with current means (lazy-key strategy)
    double ccur = cos_mean(a, b, sumU, sz);
    float dcur = (float)(1.0 - std::max(-1.0, std::min(1.0, ccur)));
    if (dcur > it.d && pq.size()>0) {
      // Push updated key and skip this stale one
      pq.push({a,b,dcur});
      continue;
    }

    // Merge b into a (choose larger root as receiver to reduce map ops)
    if (sz[a] < sz[b]) std::swap(a,b);
    int rnew = uf.unite(a,b,0.0f); // union-find
    (void)rnew; // rnew == a

    // Update sums and size
    int rA = a, rB = b;
    for (int k=0;k<sumU.nrow();++k) sumU(k,rA) += sumU(k,rB);
    sz[rA] += sz[rB];

    // Update neighbor sets: N = (Na ∪ Nb) \ {a,b}
    for (int n : nbr[rB]) { if (n!=rA) { nbr[rA].insert(n); nbr[n].erase(rB); nbr[n].insert(rA); } }
    nbr[rA].erase(rA); nbr[rA].erase(rB);
    nbr[rB].clear();

    // Push new edges
    for (int n : nbr[rA]) {
      double c = cos_mean(rA, n, sumU, sz);
      float d = (float)(1.0 - std::max(-1.0, std::min(1.0, c)));
      pq.push({rA, n, d});
    }

    active[rB]=0; nActive--;
  }
}

// Per-slice variant: merge within each slice to targetK_per_slice

inline int count_active_roots(int G, UF &uf) {
  std::unordered_set<int> roots;
  roots.reserve(G);
  for (int i=0; i<G; ++i) roots.insert(uf.find(i));
  return (int)roots.size();
}

inline void rag_agglomerate_to_K_per_slice(
  const IntegerVector &mask, int nx, int ny, int nz,
  const std::vector<int> &voxel2group, int G, int targetK,
  int nbhd,
  NumericMatrix &sumU, std::vector<int> &sz, UF &uf
) {
  if (targetK <= 0) return;
  for (int z=0; z<nz; ++z) {
    // Identify current roots present in this slice
    std::unordered_set<int> present;
    present.reserve(1024);
    for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
      int v = idx3d(x,y,z, nx,ny);
      if (!mask[v]) continue;
      int g = voxel2group[v];
      if (g < 0) continue;
      present.insert( uf.find(g) );
    }
    int S = (int)present.size();
    if (S <= targetK) continue;

    // Build adjacency among current roots in this slice
    std::unordered_set<std::pair<int,int>, PairHash> pairs;
    auto reg_pair = [&](int ra, int rb){
      if (ra == rb) return;
      int u = std::min(ra, rb), v = std::max(ra, rb);
      pairs.emplace(u, v);
    };
    for (int y=0; y<ny; ++y) for (int x=0; x<nx; ++x) {
      int v0 = idx3d(x,y,z, nx,ny);
      if (!mask[v0]) continue;
      int g0 = voxel2group[v0]; if (g0 < 0) continue;
      int r0 = uf.find(g0);
      // neighbors (forward-only)
      auto try_nb = [&](int xn, int yn) {
        int v1 = idx3d(xn, yn, z, nx,ny);
        if (!mask[v1]) return;
        int g1 = voxel2group[v1]; if (g1 < 0) return;
        int r1 = uf.find(g1);
        reg_pair(r0, r1);
      };
      if (x+1 < nx) try_nb(x+1, y);
      if (y+1 < ny) try_nb(x,   y+1);
      if (nbhd == 8) {
        if (x+1 < nx && y+1 < ny) try_nb(x+1, y+1);
        if (x+1 < nx && y-1 >= 0) try_nb(x+1, y-1);
      }
    }

    // Agglomerate within this slice to achieve targetK clusters
    // We need exactly (S - targetK) merges for this slice
    int merges_needed = S - targetK;
    
    // Build priority queue of merges based on similarity within this slice
    using MergePair = std::pair<double, std::pair<int,int>>;
    std::priority_queue<MergePair> merge_queue;
    
    for (const auto &pair : pairs) {
      int r0 = pair.first, r1 = pair.second;
      // Only consider pairs where both roots are present in this slice
      if (present.count(r0) && present.count(r1)) {
        double similarity = cos_mean(r0, r1, sumU, sz);
        merge_queue.emplace(similarity, std::make_pair(r0, r1));
      }
    }
    
    // Perform exactly merges_needed merges for this slice
    int merges_done = 0;
    while (merges_done < merges_needed && !merge_queue.empty()) {
      auto top = merge_queue.top();
      merge_queue.pop();
      
      int r0 = top.second.first, r1 = top.second.second;
      
      // Check if both roots are still separate and in this slice
      int curr_r0 = uf.find(r0), curr_r1 = uf.find(r1);
      if (curr_r0 == curr_r1) continue; // Already merged
      
      // Verify both roots are still present in this slice
      bool r0_in_slice = false, r1_in_slice = false;
      for (int y=0; y<ny && (!r0_in_slice || !r1_in_slice); ++y) {
        for (int x=0; x<nx && (!r0_in_slice || !r1_in_slice); ++x) {
          int v = idx3d(x,y,z, nx,ny);
          if (!mask[v]) continue;
          int g = voxel2group[v];
          if (g < 0) continue;
          int curr_root = uf.find(g);
          if (curr_root == curr_r0) r0_in_slice = true;
          if (curr_root == curr_r1) r1_in_slice = true;
        }
      }
      
      if (r0_in_slice && r1_in_slice) {
        uf.unite(curr_r0, curr_r1, 0.0f);
        merges_done++;
      }
    }
  }
}

// ---------------------------------------------------------------
// RUN-WISE: SLiCE-MSF + (optional) stitching + (optional) exact-K
// ---------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List slice_msf_runwise(
    Rcpp::NumericMatrix TS,     // T x N or N x T
    Rcpp::IntegerVector mask,   // length = nx*ny*nz (non-zero = in mask)
    Rcpp::IntegerVector vol_dim,// c(nx,ny,nz)
    int r = 12,
    double fh_scale = 0.32,     // FH scale (coarser when larger)
    int min_size = 80,
    int nbhd = 8,               // in-slice 4 or 8
    bool stitch_z = false,      // allow 3-D joining (vertical contact)
    double theta_link = 0.85,   // centroid similarity for stitching
    int min_contact = 1,        // min touching voxels across z to consider stitch
    bool rows_are_time = true,
    double gamma = 1.5,
    Rcpp::NumericVector voxel_dim = R_NilValue, // c(dx,dy,dz)
    double spatial_beta = 0.0,
    int target_k_global = -1,   // exact-K across volume (2-D if !stitch_z, else 3-D)
    int target_k_per_slice = -1 // exact-K per slice (ignored if stitch_z==TRUE)
) {
  if (vol_dim.size()!=3) stop("vol_dim must be c(nx,ny,nz)");
  const int nx=vol_dim[0], ny=vol_dim[1], nz=vol_dim[2];
  const int N = nx*ny*nz;
  if ((int)mask.size()!=N) stop("mask length mismatch");

  int T = rows_are_time ? TS.nrow() : TS.ncol();
  int NN= rows_are_time ? TS.ncol() : TS.nrow();
  if (NN != N) stop("TS voxel dimension does not match nx*ny*nz");
  if (nbhd!=4 && nbhd!=8) stop("nbhd must be 4 or 8");
  if (r<1) stop("r must be >=1");
  if (fh_scale<=0.0) stop("fh_scale must be >0");

  NumericVector voxdim(3);
  if (voxel_dim.size()==3) { voxdim=voxel_dim; } else { voxdim[0]=voxdim[1]=voxdim[2]=1.0; }

  // --- sketches & weights (slice-parallel) ---
  NumericVector W(N);
  NumericMatrix U(r, N);
  std::vector<double> phi; build_dct_basis(T, r, phi);
  SliceSketchWorker skw(TS, mask, phi, nx, ny, nz, T, r, rows_are_time, gamma, U, W);
  parallelFor(0, nz, skw);

  // --- per-slice FH segmentation ---
  std::vector< std::vector<int> > slice_gids(nz);
  std::vector< std::vector<int> > slice_labs(nz);
  std::vector<int> slice_maxlab(nz,0);

  for (int z=0; z<nz; ++z) {
    std::vector<int> gids, loc_of_plane, labs_local;
    std::vector<Edge> edges;
    build_slice_edges(z, mask, U, W, nx, ny, r, nbhd, voxdim, spatial_beta,
                      gids, loc_of_plane, edges);
    int ns = (int)gids.size();
    if (ns>0) {
      segment_slice_fh(ns, edges, fh_scale, min_size, labs_local);
      slice_gids[z] = std::move(gids);
      slice_labs[z] = std::move(labs_local);
      slice_maxlab[z] = *std::max_element(slice_labs[z].begin(), slice_labs[z].end());
    }
  }

  // --- map voxels to initial component ids (unique globally) ---
  std::vector<int> comp_base(nz,0);
  int total_comps=0;
  for (int z=0; z<nz; ++z){ comp_base[z]=total_comps; total_comps+=slice_maxlab[z]; }

  std::vector<int> voxel2comp(N, -1);
  for (int z=0; z<nz; ++z) {
    const auto &gids = slice_gids[z];
    const auto &labs = slice_labs[z];
    for (size_t i=0;i<gids.size();++i) if (labs[i]>0) {
      voxel2comp[gids[i]] = comp_base[z] + (labs[i]-1);
    }
  }

  // --- accumulate per-component sums & sizes ---
  NumericMatrix comp_sum(r, total_comps);
  std::vector<int> comp_sz(total_comps, 0);
  for (int g=0; g<N; ++g) {
    int c = voxel2comp[g];
    if (c<0) continue;
    comp_sz[c] += 1;
    for (int k0=0;k0<r;++k0) comp_sum(k0,c) += U(k0,g);
  }

  // --- optional 3-D stitching across slices (vertical contact) ---
  UF ufC(total_comps);
  if (stitch_z && total_comps>0) {
    std::unordered_map<std::pair<int,int>, int, PairHash> contact;
    contact.reserve((size_t)N/10 + 1);
    for (int z=0; z<nz-1; ++z) {
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
        int g0=idx3d(x,y,z,nx,ny), g1=idx3d(x,y,z+1,nx,ny);
        if (!mask[g0]||!mask[g1]) continue;
        int c0=voxel2comp[g0], c1=voxel2comp[g1];
        if (c0<0 || c1<0 || c0==c1) continue;
        int u=std::min(c0,c1), v=std::max(c0,c1);
        contact[{u,v}]++;
      }
    }
    for (const auto &kv : contact) {
      int c0=kv.first.first, c1=kv.first.second, cnt=kv.second;
      if (cnt < min_contact) continue;
      double c = cos_mean(c0, c1, comp_sum, comp_sz); // use centroid cosine
      if (c >= theta_link) ufC.unite(c0,c1,0.0f);
    }
  }

  // --- collapse components to stitched groups ---
  std::vector<int> comp_root(total_comps);
  std::unordered_map<int,int> root_index; root_index.reserve(total_comps);
  int G=0;
  for (int c=0;c<total_comps;++c) {
    int r0 = stitch_z ? ufC.find(c) : c;
    comp_root[c] = r0;
    if (!root_index.count(r0)) root_index[r0] = G++;
  }

  // aggregate sums/sizes to group level
  NumericMatrix grp_sum(r, G); std::vector<int> grp_sz(G,0);
  for (int c=0;c<total_comps;++c) if (comp_sz[c]>0) {
    int gidx = root_index[ comp_root[c] ];
    grp_sz[gidx] += comp_sz[c];
    for (int k0=0;k0<r;++k0) grp_sum(k0,gidx) += comp_sum(k0,c);
  }

  // voxel -> group index (0..G-1)
  std::vector<int> voxel2group(N, -1);
  for (int v=0; v<N; ++v) {
    int c = voxel2comp[v];
    if (c<0) continue;
    int gidx = root_index[ comp_root[c] ];
    voxel2group[v] = gidx;
  }

  // --- optional exact-K agglomeration (RAG) ---
  UF ufG(G);
  if (target_k_per_slice>0 && !stitch_z) {
    rag_agglomerate_to_K_per_slice(mask, nx, ny, nz, voxel2group, G,
                                   target_k_per_slice, nbhd, grp_sum, grp_sz, ufG);
  } else if (target_k_global>0) {
    std::unordered_set<std::pair<int,int>, PairHash> pairs;
    build_adjacency_pairs(mask, nx, ny, nz, voxel2group, nbhd,
                          /*allow_vertical=*/stitch_z, pairs);
    rag_agglomerate_to_K_global(G, target_k_global, pairs, grp_sum, grp_sz, ufG);
  }

  // --- final labeling ---
  IntegerVector labels(N); // 0 outside mask
  std::unordered_map<int,int> root2lab;
  int lab=1;
  for (int v=0; v<N; ++v) {
    if (!mask[v] || voxel2group[v] < 0) { labels[v]=0; continue; }
    int rG = ufG.find( voxel2group[v] );
    auto it = root2lab.find(rG);
    if (it==root2lab.end()) { root2lab.emplace(rG, lab); labels[v]=lab; lab++; }
    else labels[v] = it->second;
  }

  // zero out features outside mask (for cleanliness)
  for (int v=0; v<N; ++v) if (!mask[v]) { W[v]=0.0; for (int k0=0;k0<r;++k0) U(k0,v)=0.0; }

  return List::create(
    _["labels"]  = labels,
    _["weights"] = W,
    _["sketch"]  = U,
    _["params"]  = List::create(
      _["r"]=r, _["fh_scale"]=fh_scale, _["min_size"]=min_size, _["nbhd"]=nbhd,
      _["stitch_z"]=stitch_z, _["theta_link"]=theta_link, _["min_contact"]=min_contact,
      _["gamma"]=gamma, _["target_k_global"]=target_k_global, _["target_k_per_slice"]=target_k_per_slice
    )
  );
}

// ---------------------------------------------------------------
// CONSENSUS: co-association FH → (optional) exact-K via RAG
// If exact-K is requested, set use_features=TRUE so we can build
// a fused per-voxel sketch to drive merges robustly.
// ---------------------------------------------------------------

struct FuseSliceWorker : public Worker {
  const int nx, ny, nz, r, nbhd;
  const std::vector< RMatrix<double> > sketches; // r x N per run (optional)
  const std::vector< RVector<int> > labels;      // length N per run
  const std::vector< RVector<double> > weights;  // length N per run (optional)
  const RVector<int> mask;
  const double fh_scale;
  const int min_size;
  const bool use_features;
  const double lambda;
  const NumericVector voxdim;
  const double spatial_beta;
  RVector<int> out_labels;

  FuseSliceWorker(
    int nx_, int ny_, int nz_, int r_, int nbhd_,
    const std::vector< RMatrix<double> > &sketches_,
    const std::vector< RVector<int> > &labels_,
    const std::vector< RVector<double> > &weights_,
    const IntegerVector &mask_,
    double fh_scale_, int min_size_,
    bool use_features_, double lambda_,
    const NumericVector &voxdim_, double spatial_beta_,
    IntegerVector &out_labels_
  )
  : nx(nx_), ny(ny_), nz(nz_), r(r_), nbhd(nbhd_),
    sketches(sketches_), labels(labels_), weights(weights_), mask(mask_),
    fh_scale(fh_scale_), min_size(min_size_),
    use_features(use_features_), lambda(lambda_),
    voxdim(voxdim_), spatial_beta(spatial_beta_), out_labels(out_labels_) {}

  void operator()(std::size_t z0, std::size_t z1) {
    const int R = (int)labels.size();
    const double dx = (voxdim.size()>=1)?(double)voxdim[0]:1.0;
    const double dy = (voxdim.size()>=2)?(double)voxdim[1]:1.0;
    const double minlen = std::min(dx, dy);

    for (std::size_t z=z0; z<z1; ++z) {
      std::vector<int> gids; gids.reserve(nx*ny);
      std::vector<int> loc(nx*ny, -1);
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
        int g = idx3d(x,y,z,nx,ny);
        if (mask[g]) { loc[x+nx*y]=(int)gids.size(); gids.push_back(g); }
      }
      int ns=(int)gids.size(); if (ns==0) continue;

      // neighbor offsets
      std::vector<std::pair<int,int>> offs; offs.emplace_back(1,0); offs.emplace_back(0,1);
      if (nbhd==8){ offs.emplace_back(1,1); offs.emplace_back(1,-1); }

      std::vector<Edge> edges; edges.reserve((size_t)ns*(nbhd==8?4:2));
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
        int i = loc[x+nx*y]; if (i<0) continue;
        int gi=gids[i];
        for (auto &p:offs) {
          int xn=x+p.first, yn=y+p.second;
          if (xn<0||xn>=nx||yn<0||yn>=ny) continue;
          int j=loc[xn+nx*yn]; if (j<0) continue;
          int gj=gids[j];

          double Z=0.0, C=0.0, Rh=0.0, Zf=0.0;
          for (int rr=0; rr<R; ++rr) {
            double wi = (weights[rr].size()==0) ? 1.0 : weights[rr][gi];
            double wj = (weights[rr].size()==0) ? 1.0 : weights[rr][gj];
            double a = wi*wj; if (a<=0.0) continue;
            Z += a;
            if (labels[rr][gi]!=0 && labels[rr][gj]!=0 &&
                labels[rr][gi]==labels[rr][gj]) C += a;
            if (use_features && sketches[rr].ncol()>gj && sketches[rr].nrow()==r) {
              const double *pi=&sketches[rr](0,gi), *pj=&sketches[rr](0,gj);
              double dot=0.0; for (int k0=0;k0<r;++k0) dot+=pi[k0]*pj[k0];
              Rh += a*dot; Zf += a;
            }
          }
          double s = (Z>0.0) ? (C/Z) : 0.0;
          double d = 1.0 - s;
          if (use_features) {
            double rho = (Zf>0.0) ? (Rh/Zf) : 0.0;
            d = lambda*(1.0 - s) + (1.0 - lambda)*(1.0 - rho);
          }
          // spatial penalty
          double len = std::sqrt((p.first?dx*dx:0.0)+(p.second?dy*dy:0.0));
          double factor=1.0 + spatial_beta*(len/minlen - 1.0);
          d *= factor;

          Edge e; e.a=i; e.b=j; e.w=(float)std::max(0.0, std::min(2.0, d)); e.wq=quantize16(e.w);
          edges.push_back(e);
        }
      }

      std::vector<int> labs_local;
      segment_slice_fh(ns, edges, fh_scale, min_size, labs_local);
      for (int li=0; li<ns; ++li) out_labels[gids[li]] = labs_local[li];
    }
  }
};

// [[Rcpp::export]]
Rcpp::List slice_fuse_consensus(
  Rcpp::List run_results,           // list of results from slice_msf_runwise
  Rcpp::IntegerVector vol_dim,      // c(nx,ny,nz)
  int nbhd = 8,
  double fh_scale = 0.30,
  int min_size = 80,
  bool use_features = false,
  double lambda = 0.7,
  Rcpp::NumericVector voxel_dim = R_NilValue,
  double spatial_beta = 0.0,
  int target_k_global = -1,         // exact-K across volume (2-D if stitch_z=FALSE)
  int target_k_per_slice = -1,      // exact-K per slice
  bool stitch_z = false             // allow vertical merges if target_k_global>0
) {
  if (vol_dim.size()!=3) stop("vol_dim must be c(nx,ny,nz)");
  const int nx=vol_dim[0], ny=vol_dim[1], nz=vol_dim[2];
  const int N = nx*ny*nz;
  if (nbhd!=4 && nbhd!=8) stop("nbhd must be 4 or 8");

  const int R = run_results.size();
  if (R < 2) stop("Need >=2 runs for consensus");

  // Collect views
  std::vector< RVector<int> > Ls;     Ls.reserve(R);
  std::vector< RVector<double> > Ws;  Ws.reserve(R);
  std::vector< RMatrix<double> > Us;  Us.reserve(R);

  // union mask (voxels that are labeled in all runs OR at least one run?)
  // For stability, require presence in at least one run. We'll mark absent as 0 label and exclude.
  IntegerVector mask(N, 0);
  int r = -1;

  for (int rr=0; rr<R; ++rr) {
    List one = run_results[rr];
    IntegerVector L = one["labels"]; if ((int)L.size()!=N) stop("labels length mismatch");
    Ls.emplace_back(L);

    if (one.containsElementNamed("weights")) {
      NumericVector W = one["weights"]; if ((int)W.size()!=N) stop("weights length mismatch");
      Ws.emplace_back(W);
    } else {
      NumericVector empty(0); Ws.emplace_back(empty);
    }
    if (use_features) {
      if (!one.containsElementNamed("sketch")) stop("use_features=TRUE requires 'sketch' in each run");
      NumericMatrix U = one["sketch"]; if ((int)U.ncol()!=N) stop("sketch ncol mismatch");
      if (r < 0) r = U.nrow(); if (U.nrow()!=r) stop("sketch rank mismatch across runs");
      Us.emplace_back(U);
    } else {
      NumericMatrix empty(0,0); Us.emplace_back(empty);
    }
    for (int i=0;i<N;++i) if (L[i]!=0) mask[i]=1; // union
  }

  NumericVector voxdim(3);
  if (voxel_dim.size()==3) voxdim=voxel_dim; else { voxdim[0]=voxdim[1]=voxdim[2]=1.0; }

  IntegerVector fused(N); for (int i=0;i<N;++i) fused[i]=0;

  // Slice-parallel FH consensus
  FuseSliceWorker worker(
    nx, ny, nz, std::max(0,r), nbhd,
    Us, Ls, Ws, mask, fh_scale, min_size,
    use_features, lambda, voxdim, spatial_beta, fused
  );
  parallelFor(0, nz, worker);

  // Optional exact-K agglomeration (requires features for robust similarity)
  if ((target_k_global>0 || target_k_per_slice>0) && !use_features) {
    stop("Exact-K in consensus requires use_features=TRUE (to drive merges).");
  }
  if (use_features && (target_k_global>0 || target_k_per_slice>0)) {
    // Build a fused per-voxel sketch: weighted average over runs, then L2-normalize
    NumericMatrix Ubar(r, N);
    for (int v=0; v<N; ++v) {
      if (!mask[v]) continue;
      double wsum = 0.0;
      for (int rr=0; rr<R; ++rr) {
        double wv = (Ws[rr].size()==0) ? 1.0 : Ws[rr][v];
        wsum += wv;
        const double *p = &Us[rr](0,v);
        for (int k0=0;k0<r;++k0) Ubar(k0,v) += wv * p[k0];
      }
      if (wsum > 0.0) for (int k0=0;k0<r;++k0) Ubar(k0,v) /= wsum;
      // L2 normalize
      double n2=0.0; for (int k0=0;k0<r;++k0) n2 += Ubar(k0,v)*Ubar(k0,v);
      n2 = std::sqrt(std::max(1e-12, n2));
      for (int k0=0;k0<r;++k0) Ubar(k0,v) /= n2;
    }

    // Map voxels -> components from fused FH labels (per slice)
    // First, make global unique component ids from fused labels per slice.
    std::vector<int> comp_id(N, -1);
    int G=0;
    for (int z=0; z<nz; ++z) {
      // Gather slice IDs
      std::unordered_map<int,int> local2global; local2global.reserve(1024);
      for (int y=0;y<ny;++y) for (int x=0;x<nx;++x){
        int g=idx3d(x,y,z,nx,ny); if(!mask[g]) continue;
        int lab = fused[g]; if (lab<=0) continue;
        auto it = local2global.find(lab);
        if (it==local2global.end()) { local2global.emplace(lab, G); comp_id[g]=G; G++; }
        else comp_id[g] = it->second;
      }
    }
    if (G>0) {
      // Component sums & sizes
      NumericMatrix sumU(r, G); std::vector<int> sz(G,0);
      for (int v=0; v<N; ++v) {
        int c = comp_id[v]; if (c<0) continue;
        sz[c] += 1;
        for (int k0=0;k0<r;++k0) sumU(k0,c) += Ubar(k0,v);
      }
      // Build voxel->group initially == comp_id
      std::vector<int> voxel2group(N, -1);
      for (int v=0; v<N; ++v) if (comp_id[v]>=0) voxel2group[v]=comp_id[v];

      UF ufG(G);
      if (target_k_per_slice>0) {
        rag_agglomerate_to_K_per_slice(mask, nx, ny, nz, voxel2group, G,
                                       target_k_per_slice, nbhd, sumU, sz, ufG);
      } else if (target_k_global>0) {
        std::unordered_set<std::pair<int,int>, PairHash> pairs;
        build_adjacency_pairs(mask, nx, ny, nz, voxel2group, nbhd, stitch_z, pairs);
        rag_agglomerate_to_K_global(G, target_k_global, pairs, sumU, sz, ufG);
      }

      // Write back labels
      std::unordered_map<int,int> root2lab; int lab=1;
      for (int v=0; v<N; ++v) {
        if (!mask[v] || comp_id[v]<0) { fused[v]=0; continue; }
        int rG = ufG.find(comp_id[v]);
        auto it = root2lab.find(rG);
        if (it==root2lab.end()) { root2lab.emplace(rG, lab); fused[v]=lab; lab++; }
        else fused[v] = it->second;
      }
    }
  }

  // Outside union mask -> 0
  for (int i=0;i<N;++i) if (!mask[i]) fused[i]=0;

  return List::create(
    _["labels"] = fused,
    _["params"] = List::create(
      _["nbhd"]=nbhd, _["fh_scale"]=fh_scale, _["min_size"]=min_size,
      _["use_features"]=use_features, _["lambda"]=lambda,
      _["target_k_global"]=target_k_global, _["target_k_per_slice"]=target_k_per_slice,
      _["stitch_z"]=stitch_z
    )
  );
}