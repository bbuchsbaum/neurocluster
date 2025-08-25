
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>

using namespace Rcpp;
using namespace RcppParallel;

struct IntHash { std::size_t operator()(const int &x) const noexcept { return std::hash<int>()(x); } };

inline void lin_to_ijk(int lin, const IntegerVector& dims, int &i, int &j, int &k) {
  int nx = dims[0], ny = dims[1];
  i = lin % nx;
  int tmp = lin / nx;
  j = tmp % ny;
  k = tmp / ny;
}
inline int ijk_to_lin(int i, int j, int k, const IntegerVector& dims) {
  return i + j * dims[0] + k * dims[0] * dims[1];
}

static std::vector<std::array<int,3>> make_offsets6() {
  return {{ {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} }};
}
static std::vector<std::array<int,3>> make_offsets26() {
  std::vector<std::array<int,3>> off;
  for (int dz=-1; dz<=1; ++dz) for (int dy=-1; dy<=1; ++dy) for (int dx=-1; dx<=1; ++dx) {
    if (dx==0 && dy==0 && dz==0) continue;
    off.push_back({dx,dy,dz});
  }
  return off;
}

// Spatial grid index for centers
struct GridIndex {
  std::vector<int> first;
  std::vector<int> next;
  int nx, ny, nz;
  double step;
  NumericMatrix center_coords;
  NumericVector mins;
  GridIndex(int nx_, int ny_, int nz_, double cellStep, NumericMatrix center_coords_, NumericVector mins_)
    : nx(nx_), ny(ny_), nz(nz_), step(cellStep), center_coords(center_coords_), mins(mins_) {}
  inline int cell_of_center(const int k) const {
    double x = center_coords(k,0), y = center_coords(k,1), z = center_coords(k,2);
    int gx = std::max(0, std::min(nx-1, (int)std::floor((x - mins[0]) / step)));
    int gy = std::max(0, std::min(ny-1, (int)std::floor((y - mins[1]) / step)));
    int gz = std::max(0, std::min(nz-1, (int)std::floor((z - mins[2]) / step)));
    return gx + nx * (gy + ny * gz);
  }
  inline int cell_of_point(const double x, const double y, const double z) const {
    int gx = std::max(0, std::min(nx-1, (int)std::floor((x - mins[0]) / step)));
    int gy = std::max(0, std::min(ny-1, (int)std::floor((y - mins[1]) / step)));
    int gz = std::max(0, std::min(nz-1, (int)std::floor((z - mins[2]) / step)));
    return gx + nx * (gy + ny * gz);
  }
};

struct AssignWorker : public Worker {
  RMatrix<double> feats, coords, center_feats, center_coords;
  GridIndex &gindex;
  double compact2, step_mm;
  RVector<int> labels;
  RVector<double> dists;
  AssignWorker(NumericMatrix feats_,
               NumericMatrix coords_,
               NumericMatrix center_feats_,
               NumericMatrix center_coords_,
               GridIndex &gindex_,
               double compactness,
               double step_mm_,
               IntegerVector labels_,
               NumericVector dists_)
    : feats(feats_), coords(coords_), center_feats(center_feats_), center_coords(center_coords_),
      gindex(gindex_), compact2(compactness*compactness), step_mm(step_mm_),
      labels(labels_), dists(dists_) {}
  void operator()(std::size_t begin, std::size_t end) {
    int F = feats.ncol();
    int nx = gindex.nx, ny = gindex.ny, nz = gindex.nz;
    for (std::size_t i = begin; i < end; ++i) {
      double xi = coords(i,0), yi = coords(i,1), zi = coords(i,2);
      int c = gindex.cell_of_point(xi, yi, zi);
      int gz = c / (nx*ny);
      int gy = (c - gz*nx*ny) / nx;
      int gx = c - gz*nx*ny - gy*nx;
      double best = dists[i]; int bestk = labels[i];
      for (int dz = -1; dz <= 1; ++dz) {
        int cz = gz + dz; if (cz < 0 || cz >= nz) continue;
        for (int dy = -1; dy <= 1; ++dy) {
          int cy = gy + dy; if (cy < 0 || cy >= ny) continue;
          for (int dx = -1; dx <= 1; ++dx) {
            int cx = gx + dx; if (cx < 0 || cx >= nx) continue;
            int cellIndex = cx + nx * (cy + ny * cz);
            int k = gindex.first[cellIndex];
            while (k != -1) {
              double df2 = 0.0;
              for (int f = 0; f < F; ++f) {
                double d = feats(i,f) - center_feats(k,f);
                df2 += d*d;
              }
              double dxs = xi - center_coords(k,0);
              double dys = yi - center_coords(k,1);
              double dzs = zi - center_coords(k,2);
              double ds2 = (dxs*dxs + dys*dys + dzs*dzs) / (step_mm*step_mm);
              double D = df2 + compact2 * ds2;
              if (D < best) { best = D; bestk = k; }
              k = gindex.next[k];
            }
          }
        }
      }
      labels[i] = bestk;
      dists[i]  = best;
    }
  }
};

double centers_shift2(const NumericMatrix &a, const NumericMatrix &b) {
  int nrow = a.nrow(), ncol = a.ncol();
  double acc = 0.0;
  for (int i=0;i<nrow;++i) for (int j=0;j<ncol;++j) { double d = a(i,j) - b(i,j); acc += d*d; }
  return acc / (double)std::max(1, nrow);
}

// Strict connectivity: keep largest component per label, merge others to boundary-majority label
void enforce_strict_connectivity(IntegerVector &labels,
                                 const IntegerVector &mask_lin_idx,
                                 const IntegerVector &dims,
                                 const std::vector<std::array<int,3>> &offsets) {
  int N = labels.size();
  std::unordered_map<int,int, IntHash> idxmap; idxmap.reserve(N*1.3+64);
  for (int i=0;i<N;++i) idxmap[ mask_lin_idx[i] ] = i;

  std::vector<char> visited(N, 0);
  std::unordered_map<int, int> largest_comp_size;
  std::unordered_map<int, int> largest_seen_id;

  // First pass: compute component sizes and track largest per label
  int comp_id = 0;
  std::vector<int> comp_size_list;
  std::vector<int> comp_label_list;
  for (int i=0;i<N;++i) if (!visited[i]) {
    int lab = labels[i];
    int size = 0;
    std::vector<int> q; q.reserve(256); q.push_back(i); visited[i]=1;
    for (size_t h=0; h<q.size(); ++h) {
      int u = q[h];
      int lin = mask_lin_idx[u];
      int x,y,z; lin_to_ijk(lin, dims, x,y,z);
      ++size;
      for (auto &d : offsets) {
        int nx = x + d[0], ny = y + d[1], nz = z + d[2];
        if (nx<0 || ny<0 || nz<0 || nx>=dims[0] || ny>=dims[1] || nz>=dims[2]) continue;
        int nlin = ijk_to_lin(nx,ny,nz,dims);
        auto it = idxmap.find(nlin);
        if (it==idxmap.end()) continue;
        int v = it->second;
        if (!visited[v] && labels[v]==lab) { visited[v]=1; q.push_back(v); }
      }
    }
    comp_size_list.push_back(size);
    comp_label_list.push_back(lab);
    auto itL = largest_comp_size.find(lab);
    if (itL == largest_comp_size.end() || size > itL->second) {
      largest_comp_size[lab] = size;
      largest_seen_id[lab] = comp_id;
    }
    comp_id++;
  }

  // Second pass: relabel non-largest components to boundary-majority
  std::fill(visited.begin(), visited.end(), 0);
  comp_id = 0;
  for (int i=0;i<N;++i) if (!visited[i]) {
    int lab = labels[i];
    int target_keep_id = largest_seen_id[lab];
    bool keep = (comp_id == target_keep_id);
    std::vector<int> q; q.reserve(256); q.push_back(i); visited[i]=1;
    std::vector<int> comp; comp.reserve(256);
    std::unordered_map<int,int> boundary_counts;
    while (!q.empty()) {
      int u = q.back(); q.pop_back();
      comp.push_back(u);
      int lin = mask_lin_idx[u];
      int x,y,z; lin_to_ijk(lin, dims, x,y,z);
      for (auto &d : offsets) {
        int nx = x + d[0], ny = y + d[1], nz = z + d[2];
        if (nx<0 || ny<0 || nz<0 || nx>=dims[0] || ny>=dims[1] || nz>=dims[2]) continue;
        int nlin = ijk_to_lin(nx,ny,nz,dims);
        auto it = idxmap.find(nlin);
        if (it==idxmap.end()) continue;
        int v = it->second;
        if (!visited[v] && labels[v]==lab) { visited[v]=1; q.push_back(v); }
        else if (labels[v] != lab) { boundary_counts[ labels[v] ]++; }
      }
    }
    if (!keep) {
      int bestLab = lab, bestCnt = -1;
      for (auto &kv : boundary_counts) if (kv.second > bestCnt) { bestCnt = kv.second; bestLab = kv.first; }
      if (bestCnt >= 0) for (int u : comp) labels[u] = bestLab;
    }
    comp_id++;
  }
}

// Seeding: mask-aware grid + farthest fill
std::vector<int> seed_mask_grid(const NumericMatrix &coords,
                                int K, double step_mm) {
  int N = coords.nrow();
  if (N==0) return {};
  double minx = coords(0,0), maxx = coords(0,0);
  double miny = coords(0,1), maxy = coords(0,1);
  double minz = coords(0,2), maxz = coords(0,2);
  for (int i=1;i<N;++i) {
    minx = std::min(minx, coords(i,0)); maxx = std::max(maxx, coords(i,0));
    miny = std::min(miny, coords(i,1)); maxy = std::max(maxy, coords(i,1));
    minz = std::min(minz, coords(i,2)); maxz = std::max(maxz, coords(i,2));
  }
  int nx = std::max(1, (int)std::round((maxx-minx)/step_mm));
  int ny = std::max(1, (int)std::round((maxy-miny)/step_mm));
  int nz = std::max(1, (int)std::round((maxz-minz)/step_mm));
  std::vector<int> seed_idx; seed_idx.reserve((size_t)nx*ny*nz+1);

  std::vector<int> order(N); for (int i=0;i<N;++i) order[i]=i;
  std::mt19937 rng(123); std::shuffle(order.begin(), order.end(), rng);

  for (int iz=0; iz<nz; ++iz) for (int iy=0; iy<ny; ++iy) for (int ix=0; ix<nx; ++ix) {
    double gx = minx + (ix + 0.5) * ((maxx-minx)/nx);
    double gy = miny + (iy + 0.5) * ((maxy-miny)/ny);
    double gz = minz + (iz + 0.5) * ((maxz-minz)/nz);
    double best = std::numeric_limits<double>::infinity(); int besti=-1;
    int L = std::min(N, 4096);
    for (int t=0; t<L; ++t) {
      int i = order[t];
      double dx = coords(i,0)-gx, dy = coords(i,1)-gy, dz = coords(i,2)-gz;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < best) { best = d2; besti=i; }
    }
    if (besti>=0) seed_idx.push_back(besti);
  }
  std::sort(seed_idx.begin(), seed_idx.end());
  seed_idx.erase(std::unique(seed_idx.begin(), seed_idx.end()), seed_idx.end());

  if ((int)seed_idx.size() > K) seed_idx.resize(K);
  if ((int)seed_idx.size() < K) {
    std::vector<char> chosen(N, 0);
    for (int id: seed_idx) chosen[id]=1;
    while ((int)seed_idx.size() < K) {
      int besti=-1; double bestd=-1.0;
      for (int trial=0; trial<std::min(N, 8192); ++trial) {
        int i = order[trial % N]; if (chosen[i]) continue;
        double mind=std::numeric_limits<double>::infinity();
        for (int id: seed_idx) {
          double dx = coords(i,0)-coords(id,0);
          double dy = coords(i,1)-coords(id,1);
          double dz = coords(i,2)-coords(id,2);
          double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < mind) mind = d2;
        }
        if (mind > bestd) { bestd = mind; besti = i; }
      }
      if (besti>=0) { seed_idx.push_back(besti); chosen[besti]=1; } else break;
    }
  }
  return seed_idx;
}

// Poisson-disk inside mask (approximate, greedy)
std::vector<int> poisson_seeds_in_mask(const IntegerVector &mask_lin_idx,
                                       const IntegerVector &dims,
                                       const NumericVector &voxmm,
                                       double radius_mm, int K_target) {
  int N = mask_lin_idx.size(); if (N==0) return {};
  double vstep = (voxmm[0] + voxmm[1] + voxmm[2]) / 3.0;
  double rvox = std::max(1.0, radius_mm / std::max(1e-6, vstep));

  std::vector<int> order(N); for (int i=0;i<N;++i) order[i]=i;
  std::mt19937 rng(12345); std::shuffle(order.begin(), order.end(), rng);

  auto lin2xyz = [&](int lin){ int x,y,z; lin_to_ijk(lin, dims, x,y,z); return std::array<int,3>{x,y,z}; };
  std::unordered_set<int, IntHash> selected_lin;
  std::vector<std::array<int,3>> selected_xyz; selected_xyz.reserve(K_target*2+64);

  for (int idx : order) {
    if ((int)selected_xyz.size() >= K_target) break;
    int lin = mask_lin_idx[idx];
    auto xyz = lin2xyz(lin);
    bool ok = true;
    for (auto &p : selected_xyz) {
      double dx = xyz[0]-p[0], dy=xyz[1]-p[1], dz=xyz[2]-p[2];
      if (dx*dx+dy*dy+dz*dz < rvox*rvox) { ok=false; break; }
    }
    if (ok) { selected_xyz.push_back(xyz); selected_lin.insert(lin); }
  }

  // map selected lin -> masked index
  std::unordered_map<int,int, IntHash> lin2mask; lin2mask.reserve(N*1.3+64);
  for (int i=0;i<N;++i) lin2mask[ mask_lin_idx[i] ] = i;

  std::vector<int> seeds; seeds.reserve(selected_xyz.size());
  for (auto &xyz : selected_xyz) {
    int lin = ijk_to_lin(xyz[0], xyz[1], xyz[2], dims);
    auto it = lin2mask.find(lin);
    if (it != lin2mask.end()) seeds.push_back(it->second);
  }
  return seeds;
}

// Relocate seeds to minimum gradient in a local cube
void relocate_seeds_by_gradient(std::vector<int> &seed_idx,
                                const IntegerVector &mask_lin_idx,
                                const IntegerVector &dims,
                                const std::unordered_map<int,int, IntHash> &lin2mask,
                                const NumericVector &grad_masked,
                                int radius) {
  if (grad_masked.size()==0 || radius<=0) return;
  int nx=dims[0], ny=dims[1], nz=dims[2];
  for (size_t s=0; s<seed_idx.size(); ++s) {
    int mi = seed_idx[s];
    int lin = mask_lin_idx[mi];
    int x,y,z; lin_to_ijk(lin, dims, x,y,z);
    double bestg = grad_masked[mi]; int best_mi = mi;
    for (int dz=-radius; dz<=radius; ++dz) {
      int zz = z + dz; if (zz<0 || zz>=nz) continue;
      for (int dy=-radius; dy<=radius; ++dy) {
        int yy = y + dy; if (yy<0 || yy>=ny) continue;
        for (int dx=-radius; dx<=radius; ++dx) {
          int xx = x + dx; if (xx<0 || xx>=nx) continue;
          int nlin = ijk_to_lin(xx,yy,zz,dims);
          auto it = lin2mask.find(nlin);
          if (it==lin2mask.end()) continue;
          int nmi = it->second;
          double g = grad_masked[nmi];
          if (g < bestg) { bestg=g; best_mi=nmi; }
        }
      }
    }
    seed_idx[s] = best_mi;
  }
}

// [[Rcpp::export]]
Rcpp::List slic4d_core(NumericMatrix feats,
                       NumericMatrix coords,
                       IntegerVector mask_lin_idx,
                       IntegerVector dims,
                       NumericVector voxmm,
                       int K,
                       double compactness = 10.0,
                       int max_iter = 10,
                       double step_mm = 0.0,
                       int n_threads = 0,
                       std::string seed_method = "mask_poisson",
                       bool enforce_connectivity = true,
                       int min_size = 0,
                       int connectivity = 26,
                       bool strict_connectivity = true,
                       bool preserve_k = false,
                       int topup_iters = 2,
                       NumericVector grad_masked = NumericVector(),
                       int seed_relocate_radius = 1,
                       bool verbose = false) {

  const int N = feats.nrow();
  const int F = feats.ncol();
  if (K < 2 || K > N) stop("K must be in [2, N]");
  if (N < 2) stop("Need at least 2 voxels");
  if (n_threads > 0) RcppParallel::setThreadOptions(n_threads);

  // bbox & step
  double minx = coords(0,0), maxx = coords(0,0);
  double miny = coords(0,1), maxy = coords(0,1);
  double minz = coords(0,2), maxz = coords(0,2);
  for (int i=1;i<N;++i) {
    minx = std::min(minx, coords(i,0)); maxx = std::max(maxx, coords(i,0));
    miny = std::min(miny, coords(i,1)); maxy = std::max(maxy, coords(i,1));
    minz = std::min(minz, coords(i,2)); maxz = std::max(maxz, coords(i,2));
  }
  double dx = std::max(1e-6, maxx-minx), dy = std::max(1e-6, maxy-miny), dz = std::max(1e-6, maxz-minz);
  double vol = dx*dy*dz;
  if (step_mm <= 0.0) step_mm = std::cbrt(vol / (double)K);

  // Seeds
  std::vector<int> seed_idx;
  if (seed_method == "mask_poisson") {
    seed_idx = poisson_seeds_in_mask(mask_lin_idx, dims, voxmm, step_mm, K);
    if ((int)seed_idx.size() < K) {
      auto extra = seed_mask_grid(coords, K - (int)seed_idx.size(), step_mm);
      seed_idx.insert(seed_idx.end(), extra.begin(), extra.end());
      std::sort(seed_idx.begin(), seed_idx.end());
      seed_idx.erase(std::unique(seed_idx.begin(), seed_idx.end()), seed_idx.end());
      if ((int)seed_idx.size() > K) seed_idx.resize(K);
    }
  } else if (seed_method == "mask_grid" || seed_method == "grid" || seed_method == "farthest") {
    seed_idx = seed_mask_grid(coords, K, step_mm);
  } else {
    stop("Unknown seed_method");
  }

  // Build lin->mask map for relocation
  std::unordered_map<int,int, IntHash> lin2mask; lin2mask.reserve(N*1.3+64);
  for (int i=0;i<N;++i) lin2mask[ mask_lin_idx[i] ] = i;

  // Relocate seeds by gradient
  if (grad_masked.size()==N && seed_relocate_radius>0) {
    relocate_seeds_by_gradient(seed_idx, mask_lin_idx, dims, lin2mask, grad_masked, seed_relocate_radius);
  }

  // Initialize centers
  NumericMatrix center_feats(K, F), center_coords(K, 3);
  for (int k=0;k<K;++k) {
    int idx = seed_idx[k % seed_idx.size()];
    for (int f=0; f<F; ++f) center_feats(k,f) = feats(idx,f);
    center_coords(k,0) = coords(idx,0);
    center_coords(k,1) = coords(idx,1);
    center_coords(k,2) = coords(idx,2);
  }

  // Center grid index
  auto build_index = [&](GridIndex &gi){
    std::fill(gi.first.begin(), gi.first.end(), -1);
    for (int k=0;k<K;++k) {
      int cell = gi.cell_of_center(k);
      gi.next[k] = gi.first[cell];
      gi.first[cell] = k;
    }
  };
  double cellStep = 2.0 * step_mm;
  int cnx = std::max(1, (int)std::floor(dx / cellStep));
  int cny = std::max(1, (int)std::floor(dy / cellStep));
  int cnz = std::max(1, (int)std::floor(dz / cellStep));
  GridIndex gindex(cnx, cny, cnz, cellStep, center_coords, NumericVector::create(minx,miny,minz));
  gindex.first.assign(cnx*cny*cnz, -1);
  gindex.next.assign(K, -1);
  build_index(gindex);

  IntegerVector labels(N);
  NumericVector dists(N);
  for (int i=0;i<N;++i) { labels[i]=0; dists[i]=std::numeric_limits<double>::infinity(); }

  // Outer iterations
  NumericMatrix prev_center_coords = clone(center_coords);
  NumericMatrix prev_center_feats  = clone(center_feats);
  for (int it=0; it<max_iter; ++it) {
    AssignWorker worker(feats, coords, center_feats, center_coords, gindex, compactness, step_mm, labels, dists);
    parallelFor(0, (size_t)N, worker, 2000);

    // Update
    std::vector<double> sumf((size_t)K*F, 0.0);
    std::vector<double> sumx(K, 0.0), sumy(K, 0.0), sumz(K, 0.0);
    std::vector<int>    cnt(K, 0);
    for (int i=0;i<N;++i) {
      int k = labels[i];
      cnt[k]++;
      for (int f=0; f<F; ++f) sumf[(size_t)k*F+f] += feats(i,f);
      sumx[k] += coords(i,0); sumy[k] += coords(i,1); sumz[k] += coords(i,2);
    }
    for (int k=0;k<K;++k) {
      double c = std::max(1, cnt[k]);
      for (int f=0; f<F; ++f) center_feats(k,f) = sumf[(size_t)k*F+f] / c;
      center_coords(k,0) = sumx[k] / c;
      center_coords(k,1) = sumy[k] / c;
      center_coords(k,2) = sumz[k] / c;
    }

    build_index(gindex);
    double shift_c = centers_shift2(center_coords, prev_center_coords);
    double shift_f = centers_shift2(center_feats,  prev_center_feats);
    if (verbose) Rcpp::Rcout << "iter " << it+1 << " shift_coords=" << std::sqrt(shift_c)
                             << " shift_feats=" << std::sqrt(shift_f) << std::endl;
    prev_center_coords = clone(center_coords);
    prev_center_feats  = clone(center_feats);
    if (shift_c + shift_f < 1e-6) break;
  }

  if (enforce_connectivity) {
    std::vector<std::array<int,3>> offsets = (connectivity==6) ? make_offsets6() : make_offsets26();
    enforce_strict_connectivity(labels, mask_lin_idx, dims, offsets);
  }

  // Preserve K labels if requested (split largest clusters and refine a bit)
  if (preserve_k) {
    // counts
    std::vector<int> cnt(K,0); for (int i=0;i<N;++i) cnt[ labels[i] ]++;
    auto distance_D = [&](int i, int k)->double{
      double df2=0.0; for (int f=0; f<F; ++f){ double d=feats(i,f)-center_feats(k,f); df2+=d*d; }
      double dxs=coords(i,0)-center_coords(k,0), dys=coords(i,1)-center_coords(k,1), dzs=coords(i,2)-center_coords(k,2);
      double ds2=(dxs*dxs+dys*dys+dzs*dzs)/(step_mm*step_mm);
      return df2 + compactness*compactness*ds2;
    };
    auto nonempty = [&](){ int L=0; for (int k=0;k<K;++k) if (cnt[k]>0) L++; return L; };
    int L = nonempty();
    if (L < K) {
      // split biggest clusters
      for (int need=0; need < K-L; ++need) {
        int big=0; for (int k=1;k<K;++k) if (cnt[k] > cnt[big]) big=k;
        if (cnt[big] <= 1) break;
        int best_i=-1; double bestD=-1.0;
        for (int i=0;i<N;++i) if (labels[i]==big) {
          double Di = distance_D(i,big);
          if (Di>bestD){ bestD=Di; best_i=i; }
        }
        if (best_i<0) break;
        int newLab=-1; for (int k=0;k<K;++k) if (cnt[k]==0){ newLab=k; break; }
        if (newLab<0) break;
        for (int f=0; f<F; ++f) center_feats(newLab,f)=feats(best_i,f);
        center_coords(newLab,0)=coords(best_i,0);
        center_coords(newLab,1)=coords(best_i,1);
        center_coords(newLab,2)=coords(best_i,2);
        labels[best_i]=newLab; cnt[big]--; cnt[newLab]=1;
      }
      // small refinements
      for (int it2=0; it2<topup_iters; ++it2) {
        // rebuild index
        GridIndex gi(cnx, cny, cnz, cellStep, center_coords, NumericVector::create(minx,miny,minz));
        gi.first.assign(cnx*cny*cnz, -1); gi.next.assign(K, -1);
        for (int k=0;k<K;++k){ int cell=gi.cell_of_center(k); gi.next[k]=gi.first[cell]; gi.first[cell]=k; }
        for (int i=0;i<N;++i) dists[i]=std::numeric_limits<double>::infinity();
        AssignWorker w2(feats, coords, center_feats, center_coords, gi, compactness, step_mm, labels, dists);
        parallelFor(0, (size_t)N, w2, 2000);
        std::fill(cnt.begin(), cnt.end(), 0);
        std::vector<double> sumf2((size_t)K*F, 0.0);
        std::vector<double> sx(K,0.0), sy(K,0.0), sz(K,0.0);
        for (int i=0;i<N;++i) {
          int k = labels[i]; cnt[k]++;
          for (int f=0; f<F; ++f) sumf2[(size_t)k*F+f] += feats(i,f);
          sx[k]+=coords(i,0); sy[k]+=coords(i,1); sz[k]+=coords(i,2);
        }
        for (int k=0;k<K;++k){
          double c=std::max(1,cnt[k]);
          for (int f=0; f<F; ++f) center_feats(k,f)=sumf2[(size_t)k*F+f]/c;
          center_coords(k,0)=sx[k]/c; center_coords(k,1)=sy[k]/c; center_coords(k,2)=sz[k]/c;
        }
      }
    }
  }

  // Output labels as 1..K
  IntegerVector out_labels(N);
  for (int i=0;i<N;++i) out_labels[i] = labels[i] + 1;

  return List::create(_["labels"]=out_labels,
                      _["center_feats"]=center_feats,
                      _["center_coords"]=center_coords);
}
