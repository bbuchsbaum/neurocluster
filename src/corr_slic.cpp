// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__aarch64__) && defined(__ARM_FEATURE_DOTPROD)
#include <arm_neon.h>
#endif

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

struct HashProj {
  std::vector<int> h;
  std::vector<float> s;
  std::vector<int> sumS;
};

static inline bool force_portable_dot_i8() {
  // Read once; changing the env var at runtime won't affect behavior.
  static int cached = -1;
  if (cached != -1) return cached == 1;
  const char *val = std::getenv("NEUROCLUSTER_FORCE_PORTABLE_DOT_I8");
  cached = (val && val[0] != '\0' && val[0] != '0') ? 1 : 0;
  return cached == 1;
}

static HashProj make_hash(const int T, const int d, const uint64_t seed) {
  HashProj hp;
  hp.h.resize(T);
  hp.s.resize(T);
  hp.sumS.assign(d, 0);

  for (int t = 0; t < T; ++t) {
    uint64_t x = seed + (uint64_t)t * 0x9e3779b97f4a7c15ULL;
    uint64_t z = splitmix64(x);
    int bin = static_cast<int>(z % static_cast<uint64_t>(d));
    float sign = ((z >> 32) & 1ULL) ? 1.0f : -1.0f;
    hp.h[t] = bin;
    hp.s[t] = sign;
    hp.sumS[bin] += (sign > 0.0f ? 1 : -1);
  }

  return hp;
}

static inline int lin_index(const int x, const int y, const int z,
                            const int X, const int Y) {
  return x + y * X + z * X * Y;
}

template <int D>
static inline float dot_row_fixed(const float *a, const float *b) {
  float s = 0.0f;
#ifdef _OPENMP
#pragma omp simd reduction(+ : s)
#endif
  for (int j = 0; j < D; ++j) {
    s += a[j] * b[j];
  }
  return s;
}

static inline float dot_row(const float *a, const float *b, const int d) {
  if (d == 64) return dot_row_fixed<64>(a, b);
  if (d == 32) return dot_row_fixed<32>(a, b);
  if (d == 16) return dot_row_fixed<16>(a, b);
  float s = 0.0f;
#ifdef _OPENMP
#pragma omp simd reduction(+ : s)
#endif
  for (int j = 0; j < d; ++j) {
    s += a[j] * b[j];
  }
  return s;
}

template <int D>
static inline int32_t dot_row_i8_fixed(const int8_t *a, const int8_t *b) {
  int32_t s = 0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : s)
#endif
  for (int j = 0; j < D; ++j) {
    s += static_cast<int32_t>(a[j]) * static_cast<int32_t>(b[j]);
  }
  return s;
}

#if defined(__aarch64__) && defined(__ARM_FEATURE_DOTPROD)
static inline int32_t dot_row_i8_neon_dotprod(const int8_t *a, const int8_t *b, const int d) {
  int32x4_t acc = vdupq_n_s32(0);
  int j = 0;
  for (; j + 16 <= d; j += 16) {
    const int8x16_t va = vld1q_s8(a + j);
    const int8x16_t vb = vld1q_s8(b + j);
    acc = vdotq_s32(acc, va, vb);
  }
  // Horizontal sum of 4 lanes (portable NEON reduction).
  const int32x2_t sum2 = vadd_s32(vget_low_s32(acc), vget_high_s32(acc));
  int32_t s = vget_lane_s32(sum2, 0) + vget_lane_s32(sum2, 1);
  for (; j < d; ++j) {
    s += static_cast<int32_t>(a[j]) * static_cast<int32_t>(b[j]);
  }
  return s;
}
#endif

static inline int32_t dot_row_i8(const int8_t *a, const int8_t *b, const int d) {
#if defined(__aarch64__) && defined(__ARM_FEATURE_DOTPROD)
  if (!force_portable_dot_i8()) {
    return dot_row_i8_neon_dotprod(a, b, d);
  }
#endif
  if (d == 64) return dot_row_i8_fixed<64>(a, b);
  if (d == 32) return dot_row_i8_fixed<32>(a, b);
  if (d == 16) return dot_row_i8_fixed<16>(a, b);
  int32_t s = 0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : s)
#endif
  for (int j = 0; j < d; ++j) {
    s += static_cast<int32_t>(a[j]) * static_cast<int32_t>(b[j]);
  }
  return s;
}

static inline int8_t quantize_unit_i8(const float x) {
  float y = x;
  if (y > 1.0f) y = 1.0f;
  if (y < -1.0f) y = -1.0f;
  int q = static_cast<int>(std::lround(y * 127.0f));
  if (q > 127) q = 127;
  if (q < -127) q = -127;
  return static_cast<int8_t>(q);
}

struct VoxelGrid {
  int X;
  int Y;
  int Z;
  int V;
  int N;

  std::vector<int32_t> mask_lin;
  std::vector<int16_t> vx;
  std::vector<int16_t> vy;
  std::vector<int16_t> vz;
  std::vector<int32_t> id_map;

  int16_t minx;
  int16_t miny;
  int16_t minz;
  int16_t maxx;
  int16_t maxy;
  int16_t maxz;
};

static VoxelGrid build_voxel_grid(const IntegerVector &mask_lin_idx,
                                  const IntegerVector &dims) {
  if (dims.size() < 3) {
    stop("dims must contain X, Y, Z");
  }

  VoxelGrid G;
  G.X = dims[0];
  G.Y = dims[1];
  G.Z = dims[2];
  G.V = G.X * G.Y * G.Z;
  G.N = mask_lin_idx.size();

  G.mask_lin.resize(G.N);
  G.vx.resize(G.N);
  G.vy.resize(G.N);
  G.vz.resize(G.N);
  G.id_map.assign(G.V, -1);

  G.minx = G.miny = G.minz = std::numeric_limits<int16_t>::max();
  G.maxx = G.maxy = G.maxz = std::numeric_limits<int16_t>::min();

  for (int i = 0; i < G.N; ++i) {
    int lin = mask_lin_idx[i];
    if (lin < 0 || lin >= G.V) {
      stop("mask_lin_idx contains out-of-range linear index");
    }

    int x = lin % G.X;
    int y = (lin / G.X) % G.Y;
    int z = lin / (G.X * G.Y);

    G.mask_lin[i] = static_cast<int32_t>(lin);
    G.vx[i] = static_cast<int16_t>(x);
    G.vy[i] = static_cast<int16_t>(y);
    G.vz[i] = static_cast<int16_t>(z);
    G.id_map[lin] = static_cast<int32_t>(i);

    if (x < G.minx) G.minx = static_cast<int16_t>(x);
    if (y < G.miny) G.miny = static_cast<int16_t>(y);
    if (z < G.minz) G.minz = static_cast<int16_t>(z);
    if (x > G.maxx) G.maxx = static_cast<int16_t>(x);
    if (y > G.maxy) G.maxy = static_cast<int16_t>(y);
    if (z > G.maxz) G.maxz = static_cast<int16_t>(z);
  }

  return G;
}

static int snap_seed_to_mask(const int x, const int y, const int z,
                             const VoxelGrid &G, const int radius) {
  int best = -1;
  int bestd2 = std::numeric_limits<int>::max();

  for (int dz = -radius; dz <= radius; ++dz) {
    int zz = z + dz;
    if (zz < 0 || zz >= G.Z) continue;
    for (int dy = -radius; dy <= radius; ++dy) {
      int yy = y + dy;
      if (yy < 0 || yy >= G.Y) continue;
      for (int dx = -radius; dx <= radius; ++dx) {
        int xx = x + dx;
        if (xx < 0 || xx >= G.X) continue;
        int lin = lin_index(xx, yy, zz, G.X, G.Y);
        int vid = G.id_map[lin];
        if (vid < 0) continue;
        int d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < bestd2) {
          bestd2 = d2;
          best = vid;
        }
      }
    }
  }

  return best;
}

static std::vector<int> axis_points(const int minv, const int maxv, const float S) {
  const int span = std::max(1, maxv - minv + 1);
  const int npts = std::max(1, static_cast<int>(std::round(static_cast<float>(span) / std::max(1.0f, S))));
  std::vector<int> pts;
  pts.reserve(static_cast<size_t>(npts));

  for (int i = 0; i < npts; ++i) {
    float u = (static_cast<float>(i) + 0.5f) / static_cast<float>(npts);
    int p = minv + static_cast<int>(std::floor(u * static_cast<float>(span)));
    if (p < minv) p = minv;
    if (p > maxv) p = maxv;
    if (pts.empty() || pts.back() != p) {
      pts.push_back(p);
    }
  }

  if (pts.empty()) pts.push_back(minv);
  return pts;
}

struct Centers {
  int K;
  int d;
  std::vector<float> cx;
  std::vector<float> cy;
  std::vector<float> cz;
  std::vector<float> cf; // K * d, row-major
};

static Centers init_seeds_grid(const VoxelGrid &G,
                               const std::vector<float> &feat,
                               const int d,
                               const int K,
                               const float S,
                               const int seed) {
  Centers C;
  C.K = K;
  C.d = d;
  C.cx.resize(K);
  C.cy.resize(K);
  C.cz.resize(K);
  C.cf.resize(static_cast<size_t>(K) * static_cast<size_t>(d));

  std::vector<int> seeds;
  seeds.reserve(static_cast<size_t>(K) * 2U);

  const std::vector<int> xpts = axis_points(G.minx, G.maxx, S);
  const std::vector<int> ypts = axis_points(G.miny, G.maxy, S);
  const std::vector<int> zpts = axis_points(G.minz, G.maxz, S);
  const int rmax = std::max(2, static_cast<int>(std::ceil(std::max(1.0f, S))));

  std::vector<uint8_t> seen(static_cast<size_t>(G.N), 0);
  for (int zi : zpts) {
    for (int yi : ypts) {
      for (int xi : xpts) {
        int vid = -1;
        for (int r = 1; r <= rmax; ++r) {
          vid = snap_seed_to_mask(xi, yi, zi, G, r);
          if (vid >= 0) break;
        }
        if (vid >= 0 && !seen[static_cast<size_t>(vid)]) {
          seen[static_cast<size_t>(vid)] = 1;
          seeds.push_back(vid);
        }
      }
    }
  }

  if (seeds.empty()) {
    seeds.push_back(0);
    seen[0] = 1;
  }

  if (static_cast<int>(seeds.size()) > K) {
    // Choose K spatially dispersed seeds (farthest-point sampling) to avoid
    // scan-order bias when reducing dense grid candidates.
    const int M = static_cast<int>(seeds.size());
    std::vector<int> keep;
    keep.reserve(K);
    std::vector<float> min_d2(static_cast<size_t>(M), std::numeric_limits<float>::infinity());

    uint64_t rng = static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 0x243f6a8885a308d3ULL;
    int first = static_cast<int>(splitmix64(rng) % static_cast<uint64_t>(M));
    keep.push_back(seeds[static_cast<size_t>(first)]);

    for (int i = 0; i < M; ++i) {
      const int vid = seeds[static_cast<size_t>(i)];
      const float dx = static_cast<float>(G.vx[vid] - G.vx[keep[0]]);
      const float dy = static_cast<float>(G.vy[vid] - G.vy[keep[0]]);
      const float dz = static_cast<float>(G.vz[vid] - G.vz[keep[0]]);
      min_d2[static_cast<size_t>(i)] = dx * dx + dy * dy + dz * dz;
    }

    while (static_cast<int>(keep.size()) < K) {
      int best_i = -1;
      float best_d2 = -1.0f;
      for (int i = 0; i < M; ++i) {
        const float d2 = min_d2[static_cast<size_t>(i)];
        if (d2 > best_d2) {
          best_d2 = d2;
          best_i = i;
        }
      }

      if (best_i < 0) break;
      const int new_vid = seeds[static_cast<size_t>(best_i)];
      keep.push_back(new_vid);

      for (int i = 0; i < M; ++i) {
        const int vid = seeds[static_cast<size_t>(i)];
        const float dx = static_cast<float>(G.vx[vid] - G.vx[new_vid]);
        const float dy = static_cast<float>(G.vy[vid] - G.vy[new_vid]);
        const float dz = static_cast<float>(G.vz[vid] - G.vz[new_vid]);
        const float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < min_d2[static_cast<size_t>(i)]) {
          min_d2[static_cast<size_t>(i)] = d2;
        }
      }
    }

    seeds.swap(keep);
  } else if (static_cast<int>(seeds.size()) < K) {
    std::vector<int> unseen;
    unseen.reserve(static_cast<size_t>(G.N));
    for (int vid = 0; vid < G.N; ++vid) {
      if (!seen[static_cast<size_t>(vid)]) {
        unseen.push_back(vid);
      }
    }

    uint64_t rng = static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 12345ULL;
    const int need = K - static_cast<int>(seeds.size());
    const int avail = static_cast<int>(unseen.size());
    const int take = std::min(need, avail);
    for (int i = 0; i < take; ++i) {
      const int rem = avail - i;
      const int j = i + static_cast<int>(splitmix64(rng) % static_cast<uint64_t>(rem));
      std::swap(unseen[static_cast<size_t>(i)], unseen[static_cast<size_t>(j)]);
      const int vid = unseen[static_cast<size_t>(i)];
      seeds.push_back(vid);
      seen[static_cast<size_t>(vid)] = 1;
    }

    // Safety fallback if random fill could not satisfy K for any reason.
    for (int vid = 0; static_cast<int>(seeds.size()) < K && vid < G.N; ++vid) {
      if (!seen[static_cast<size_t>(vid)]) {
        seeds.push_back(vid);
        seen[static_cast<size_t>(vid)] = 1;
      }
    }
  }

  for (int k = 0; k < K; ++k) {
    int vid = seeds[k];
    C.cx[k] = static_cast<float>(G.vx[vid]);
    C.cy[k] = static_cast<float>(G.vy[vid]);
    C.cz[k] = static_cast<float>(G.vz[vid]);
    const float *src = feat.data() + static_cast<size_t>(vid) * static_cast<size_t>(d);
    float *dst = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
    std::copy(src, src + d, dst);
  }

  return C;
}

struct Bins {
  int nbx;
  int nby;
  int nbz;
  int nbins;
  std::vector<int32_t> voxel_bin;
  std::vector<uint8_t> nbr_count;
  std::vector<int32_t> nbr_bin; // [nbins, 27]
  std::vector<int32_t> bin_offsets;
  std::vector<int32_t> bin_centers;
};

static inline int bin_id(const int bx, const int by, const int bz,
                         const int nbx, const int nby) {
  return bx + by * nbx + bz * nbx * nby;
}

static Bins init_bins(const VoxelGrid &G, const float S, const int K) {
  Bins B;
  float invS = 1.0f / S;
  B.nbx = static_cast<int>(std::floor((G.X - 1) * invS)) + 1;
  B.nby = static_cast<int>(std::floor((G.Y - 1) * invS)) + 1;
  B.nbz = static_cast<int>(std::floor((G.Z - 1) * invS)) + 1;
  B.nbins = B.nbx * B.nby * B.nbz;

  B.voxel_bin.resize(G.N);
  for (int v = 0; v < G.N; ++v) {
    const int bx = static_cast<int>(std::floor(static_cast<float>(G.vx[v]) * invS));
    const int by = static_cast<int>(std::floor(static_cast<float>(G.vy[v]) * invS));
    const int bz = static_cast<int>(std::floor(static_cast<float>(G.vz[v]) * invS));
    B.voxel_bin[v] = bin_id(bx, by, bz, B.nbx, B.nby);
  }

  B.nbr_count.assign(static_cast<size_t>(B.nbins), 0);
  B.nbr_bin.assign(static_cast<size_t>(B.nbins) * 27U, -1);
  for (int bz = 0; bz < B.nbz; ++bz) {
    for (int by = 0; by < B.nby; ++by) {
      for (int bx = 0; bx < B.nbx; ++bx) {
        const int bid = bin_id(bx, by, bz, B.nbx, B.nby);
        int n = 0;
        for (int dz = -1; dz <= 1; ++dz) {
          const int zz = bz + dz;
          if (zz < 0 || zz >= B.nbz) continue;
          for (int dy = -1; dy <= 1; ++dy) {
            const int yy = by + dy;
            if (yy < 0 || yy >= B.nby) continue;
            for (int dx = -1; dx <= 1; ++dx) {
              const int xx = bx + dx;
              if (xx < 0 || xx >= B.nbx) continue;
              B.nbr_bin[static_cast<size_t>(bid) * 27U + static_cast<size_t>(n)] =
                bin_id(xx, yy, zz, B.nbx, B.nby);
              ++n;
            }
          }
        }
        B.nbr_count[static_cast<size_t>(bid)] = static_cast<uint8_t>(n);
      }
    }
  }

  B.bin_offsets.assign(B.nbins + 1, 0);
  B.bin_centers.clear();
  B.bin_centers.reserve(static_cast<size_t>(K));
  return B;
}

static void rebuild_center_bins(Bins &B, const Centers &C, const float S) {
  float invS = 1.0f / S;
  std::fill(B.bin_offsets.begin(), B.bin_offsets.end(), 0);
  for (int k = 0; k < C.K; ++k) {
    int bx = static_cast<int>(C.cx[k] * invS);
    int by = static_cast<int>(C.cy[k] * invS);
    int bz = static_cast<int>(C.cz[k] * invS);

    if (bx < 0) bx = 0;
    if (by < 0) by = 0;
    if (bz < 0) bz = 0;
    if (bx >= B.nbx) bx = B.nbx - 1;
    if (by >= B.nby) by = B.nby - 1;
    if (bz >= B.nbz) bz = B.nbz - 1;

    const int id = bin_id(bx, by, bz, B.nbx, B.nby);
    ++B.bin_offsets[static_cast<size_t>(id) + 1];
  }
  for (int i = 1; i <= B.nbins; ++i) {
    B.bin_offsets[static_cast<size_t>(i)] += B.bin_offsets[static_cast<size_t>(i - 1)];
  }
  std::vector<int32_t> bin_pos(static_cast<size_t>(B.nbins), 0);
  B.bin_centers.resize(static_cast<size_t>(C.K));
  for (int k = 0; k < C.K; ++k) {
    int bx = static_cast<int>(C.cx[k] * invS);
    int by = static_cast<int>(C.cy[k] * invS);
    int bz = static_cast<int>(C.cz[k] * invS);

    if (bx < 0) bx = 0;
    if (by < 0) by = 0;
    if (bz < 0) bz = 0;
    if (bx >= B.nbx) bx = B.nbx - 1;
    if (by >= B.nby) by = B.nby - 1;
    if (bz >= B.nbz) bz = B.nbz - 1;

    const int id = bin_id(bx, by, bz, B.nbx, B.nby);
    const int offset = B.bin_offsets[static_cast<size_t>(id)];
    const int idx = offset + bin_pos[static_cast<size_t>(id)];
    B.bin_centers[static_cast<size_t>(idx)] = k;
    ++bin_pos[static_cast<size_t>(id)];
  }
}

static void enforce_connectivity(const VoxelGrid &G,
                                 std::vector<int32_t> &labels,
                                 const int K,
                                 const int min_size,
                                 const int connectivity) {
  std::vector<int32_t> comp_id(G.N, -1);
  std::vector<int32_t> comp_label;
  std::vector<int32_t> comp_size;
  std::vector<int32_t> comp_adj;

  comp_label.reserve(K * 2);
  comp_size.reserve(K * 2);
  comp_adj.reserve(K * 2);

  std::vector<int> queue;
  queue.reserve(1024);

  auto add_neighbor = [&](const int xx, const int yy, const int zz, std::vector<int> &out) {
    int lin = lin_index(xx, yy, zz, G.X, G.Y);
    int u = G.id_map[lin];
    if (u >= 0) out.push_back(u);
  };

  auto neighbors = [&](const int v, std::vector<int> &out) {
    out.clear();
    const int x = G.vx[v];
    const int y = G.vy[v];
    const int z = G.vz[v];

    if (connectivity == 6) {
      if (x > 0) add_neighbor(x - 1, y, z, out);
      if (x + 1 < G.X) add_neighbor(x + 1, y, z, out);
      if (y > 0) add_neighbor(x, y - 1, z, out);
      if (y + 1 < G.Y) add_neighbor(x, y + 1, z, out);
      if (z > 0) add_neighbor(x, y, z - 1, out);
      if (z + 1 < G.Z) add_neighbor(x, y, z + 1, out);
    } else {
      for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0 && dz == 0) continue;
            int xx = x + dx;
            int yy = y + dy;
            int zz = z + dz;
            if (xx < 0 || xx >= G.X || yy < 0 || yy >= G.Y || zz < 0 || zz >= G.Z) continue;
            add_neighbor(xx, yy, zz, out);
          }
        }
      }
    }
  };

  std::vector<int> neigh;
  neigh.reserve(32);

  int ncomp = 0;
  for (int v0 = 0; v0 < G.N; ++v0) {
    if (comp_id[v0] != -1) continue;

    int32_t lab = labels[v0];
    queue.clear();
    queue.push_back(v0);
    comp_id[v0] = ncomp;

    int32_t size = 0;
    std::unordered_map<int32_t, int32_t> adj_counts;
    int32_t best_adj = lab;
    int32_t best_cnt = 0;

    for (size_t qi = 0; qi < queue.size(); ++qi) {
      int v = queue[qi];
      ++size;

      neighbors(v, neigh);
      for (int u : neigh) {
        if (labels[u] == lab) {
          if (comp_id[u] == -1) {
            comp_id[u] = ncomp;
            queue.push_back(u);
          }
        } else {
          int32_t nl = labels[u];
          int32_t cnt = ++adj_counts[nl];
          if (cnt > best_cnt) {
            best_cnt = cnt;
            best_adj = nl;
          }
        }
      }
    }

    comp_label.push_back(lab);
    comp_size.push_back(size);
    comp_adj.push_back(best_adj);
    ++ncomp;
  }

  std::vector<int32_t> largest_comp(K, -1);
  std::vector<int32_t> largest_size(K, -1);
  for (int c = 0; c < ncomp; ++c) {
    int lab = comp_label[c];
    if (lab < 0 || lab >= K) continue;
    if (comp_size[c] > largest_size[lab]) {
      largest_size[lab] = comp_size[c];
      largest_comp[lab] = c;
    }
  }

  for (int v = 0; v < G.N; ++v) {
    int c = comp_id[v];
    int lab = labels[v];
    bool keep = (c == largest_comp[lab]);
    bool too_small = (min_size > 0 && comp_size[c] < min_size);
    if (!keep || too_small) {
      labels[v] = comp_adj[c];
    }
  }
}

} // namespace

// [[Rcpp::export]]
List corrslic_core(const NumericMatrix feat,
                   const IntegerVector mask_lin_idx,
                   const IntegerVector dims,
                   const int K,
                   const int d = 64,
                   const int sketch_repeats = 2,
                   const double alpha = 0.5,
                   const int max_iter = 5,
                   const int seed = 1,
                   const int assign_stride = 1,
                   const bool quantize_assign = false,
                   const std::string embed_basis = "hash",
                   const bool whiten_embed = false,
                   const int refine_exact_iters = 0,
                   const bool refine_boundary_only = true,
                   const int refine_stride = 1,
                   const double refine_alpha = -1.0,
                   const int connectivity = 6,
                   const int min_size = 0,
                   int n_threads = 0,
                   const bool verbose = false) {
  const int N = feat.nrow();
  const int T = feat.ncol();

  if (N <= 0) stop("feat has zero rows");
  if (T <= 0) stop("feat has zero columns");
  if (mask_lin_idx.size() != N) stop("mask_lin_idx length must match nrow(feat)");
  if (K < 2) stop("K must be >= 2");
  if (K > N) stop("K must be <= number of masked voxels");
  if (d < 8) stop("embedding dimension d must be >= 8");
  if (sketch_repeats < 1 || sketch_repeats > 8) stop("sketch_repeats must be in [1, 8]");
  if (assign_stride < 1 || assign_stride > 16) stop("assign_stride must be in [1, 16]");
  if (embed_basis != "hash" && embed_basis != "dct") stop("embed_basis must be 'hash' or 'dct'");
  if (refine_exact_iters < 0 || refine_exact_iters > 16) stop("refine_exact_iters must be in [0, 16]");
  if (refine_stride < 1 || refine_stride > 64) stop("refine_stride must be in [1, 64]");
  if (!std::isfinite(refine_alpha) || (refine_alpha != -1.0 && refine_alpha <= 0.0)) {
    stop("refine_alpha must be -1 or a positive finite scalar");
  }
  if (connectivity != 6 && connectivity != 26) stop("connectivity must be 6 or 26");

#ifdef _OPENMP
  if (n_threads <= 0) n_threads = omp_get_max_threads();
#else
  n_threads = 1;
#endif

  if (verbose) {
    Rcout << "corrslic_core: N=" << N << " T=" << T << " K=" << K
          << " d=" << d << " repeats=" << sketch_repeats
          << " stride=" << assign_stride
          << " qassign=" << (quantize_assign ? "true" : "false")
          << " basis=" << embed_basis
          << " whiten=" << (whiten_embed ? "true" : "false")
          << " refine_iters=" << refine_exact_iters
          << " refine_boundary=" << (refine_boundary_only ? "true" : "false")
          << " refine_stride=" << refine_stride
          << " iters=" << max_iter << "\n";
  }

  VoxelGrid G = build_voxel_grid(mask_lin_idx, dims);

  // Embedding accumulation: U[d x N] dimension-major.
  std::vector<float> U(static_cast<size_t>(d) * static_cast<size_t>(N), 0.0f);
  std::vector<float> norm2(static_cast<size_t>(N), 0.0f);
  std::vector<float> invn(static_cast<size_t>(N), 0.0f);
  const float invT = 1.0f / static_cast<float>(T);

  if (embed_basis == "hash") {
    std::vector<HashProj> hps;
    hps.reserve(static_cast<size_t>(sketch_repeats));
    for (int r = 0; r < sketch_repeats; ++r) {
      const uint64_t seed_r = static_cast<uint64_t>(seed) + 0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(r + 1);
      hps.push_back(make_hash(T, d, seed_r));
    }
    const float sketch_scale = 1.0f / std::sqrt(static_cast<float>(sketch_repeats));
    std::vector<float> sumS_total(static_cast<size_t>(d), 0.0f);
    for (int r = 0; r < sketch_repeats; ++r) {
      for (int j = 0; j < d; ++j) {
        sumS_total[static_cast<size_t>(j)] += sketch_scale * static_cast<float>(hps[static_cast<size_t>(r)].sumS[static_cast<size_t>(j)]);
      }
    }

    std::vector<float> sum_x(static_cast<size_t>(N), 0.0f);
    // One threaded region to avoid O(T) thread launch overheads.
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
      int tid = 0;
      int nt = 1;
#ifdef _OPENMP
      tid = omp_get_thread_num();
      nt = omp_get_num_threads();
#endif
      const int start = static_cast<int>((static_cast<int64_t>(tid) * N) / nt);
      const int end = static_cast<int>((static_cast<int64_t>(tid + 1) * N) / nt);

      for (int t = 0; t < T; ++t) {
        const double *col = &feat(0, t);
        for (int v = start; v < end; ++v) {
          sum_x[static_cast<size_t>(v)] += static_cast<float>(col[v]);
        }

        for (int r = 0; r < sketch_repeats; ++r) {
          const HashProj &hp = hps[static_cast<size_t>(r)];
          const int bin = hp.h[static_cast<size_t>(t)];
          const float sign = sketch_scale * hp.s[static_cast<size_t>(t)];
          float *urow = U.data() + static_cast<size_t>(bin) * static_cast<size_t>(N);
          for (int v = start; v < end; ++v) {
            urow[static_cast<size_t>(v)] += sign * static_cast<float>(col[v]);
          }
        }
      }
    }

    for (int j = 0; j < d; ++j) {
      float *urow = U.data() + static_cast<size_t>(j) * static_cast<size_t>(N);
      const float ss = sumS_total[static_cast<size_t>(j)];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        const float mean = sum_x[static_cast<size_t>(v)] * invT;
        urow[static_cast<size_t>(v)] -= mean * ss;
      }
    }
  } else {
    // DCT basis on demeaned time series; emphasizes low-frequency temporal structure.
    std::vector<float> dct_basis(static_cast<size_t>(T) * static_cast<size_t>(d), 0.0f);
    const float dct_scale = std::sqrt(2.0f / static_cast<float>(T));
    const float pi = 3.14159265358979323846f;
    for (int t = 0; t < T; ++t) {
      for (int j = 0; j < d; ++j) {
        const int k = j + 1; // skip DC component
        const float ang = pi * (static_cast<float>(t) + 0.5f) * static_cast<float>(k) / static_cast<float>(T);
        dct_basis[static_cast<size_t>(t) * static_cast<size_t>(d) + static_cast<size_t>(j)] = dct_scale * std::cos(ang);
      }
    }

    const double *feat_ptr = REAL(feat);
    std::vector<float> mean_x(static_cast<size_t>(N), 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      const double *xv = feat_ptr + static_cast<size_t>(v);
      float mean = 0.0f;
      for (int t = 0; t < T; ++t) {
        mean += static_cast<float>(xv[static_cast<size_t>(t) * static_cast<size_t>(N)]);
      }
      mean_x[static_cast<size_t>(v)] = mean * invT;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int j = 0; j < d; ++j) {
      float *urow = U.data() + static_cast<size_t>(j) * static_cast<size_t>(N);
      std::fill(urow, urow + static_cast<size_t>(N), 0.0f);
      for (int t = 0; t < T; ++t) {
        const double *col = feat_ptr + static_cast<size_t>(t) * static_cast<size_t>(N);
        const float bt = dct_basis[static_cast<size_t>(t) * static_cast<size_t>(d) + static_cast<size_t>(j)];
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int v = 0; v < N; ++v) {
          urow[static_cast<size_t>(v)] += (static_cast<float>(col[static_cast<size_t>(v)]) - mean_x[static_cast<size_t>(v)]) * bt;
        }
      }
    }
  }

  if (whiten_embed) {
    // Global whitening across voxels to reduce dominant embedding dimensions.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int j = 0; j < d; ++j) {
      float *urow = U.data() + static_cast<size_t>(j) * static_cast<size_t>(N);
      float mu = 0.0f;
      for (int v = 0; v < N; ++v) mu += urow[static_cast<size_t>(v)];
      mu /= static_cast<float>(N);
      float var = 0.0f;
      for (int v = 0; v < N; ++v) {
        const float dv = urow[static_cast<size_t>(v)] - mu;
        var += dv * dv;
      }
      const float invsd = 1.0f / std::sqrt(var / static_cast<float>(N) + 1e-8f);
      for (int v = 0; v < N; ++v) {
        urow[static_cast<size_t>(v)] = (urow[static_cast<size_t>(v)] - mu) * invsd;
      }
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int v = 0; v < N; ++v) {
    norm2[v] = 0.0f;
  }

  for (int j = 0; j < d; ++j) {
    float *urow = U.data() + static_cast<size_t>(j) * static_cast<size_t>(N);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      float u = urow[v];
      norm2[v] += u * u;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int v = 0; v < N; ++v) {
    invn[v] = 1.0f / std::sqrt(norm2[v] + 1e-8f);
  }

  for (int j = 0; j < d; ++j) {
    float *urow = U.data() + static_cast<size_t>(j) * static_cast<size_t>(N);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      urow[v] *= invn[v];
    }
  }

  // Transpose to voxel-major for fast per-voxel dot products.
  std::vector<float> feat_emb(static_cast<size_t>(N) * static_cast<size_t>(d));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int v = 0; v < N; ++v) {
    float *dst = feat_emb.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
    for (int j = 0; j < d; ++j) {
      dst[j] = U[static_cast<size_t>(j) * static_cast<size_t>(N) + static_cast<size_t>(v)];
    }
  }

  std::vector<int8_t> feat_emb_q;
  if (quantize_assign) {
    feat_emb_q.resize(static_cast<size_t>(N) * static_cast<size_t>(d));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      const float *src = feat_emb.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
      int8_t *dst = feat_emb_q.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
      for (int j = 0; j < d; ++j) {
        dst[j] = quantize_unit_i8(src[j]);
      }
    }
  }

  float S = std::cbrt(static_cast<float>(N) / static_cast<float>(K));
  if (S < 1.0f) S = 1.0f;
  float alpha_over_s2 = static_cast<float>(alpha / (static_cast<double>(S) * static_cast<double>(S)));

  Centers C = init_seeds_grid(G, feat_emb, d, K, S, seed);
  Bins B = init_bins(G, S, K);

  std::vector<int32_t> labels(static_cast<size_t>(N), 0);
  std::vector<int32_t> new_labels(static_cast<size_t>(N), 0);

  uint64_t rng_state = static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 99991ULL;
  const int nt = n_threads;
  std::vector<int32_t> count(static_cast<size_t>(nt) * static_cast<size_t>(K), 0);
  std::vector<float> sum_xyz(static_cast<size_t>(nt) * static_cast<size_t>(K) * 3U, 0.0f);
  std::vector<float> sum_f(static_cast<size_t>(nt) * static_cast<size_t>(K) * static_cast<size_t>(d), 0.0f);
  std::vector<int32_t> count_total(static_cast<size_t>(K), 0);
  const size_t thread_cluster_stride = static_cast<size_t>(K);
  const size_t thread_xyz_stride = static_cast<size_t>(K) * 3U;
  const size_t thread_feat_stride = static_cast<size_t>(K) * static_cast<size_t>(d);
  std::vector<int8_t> center_q;
  const float inv_qscale2 = 1.0f / (127.0f * 127.0f);
  const int stride_eff = std::max(1, assign_stride);
  const bool request_subsample = (stride_eff > 1);
  int subsample_iters = 0;
  if (request_subsample) {
    if (max_iter >= 3) {
      // Intermediate plane-subsampled iterations (z-phase rotation) reduce
      // assignment cost without changing the overall SLIC loop structure.
      // Leave at least one full iteration after the subsampled block.
      const int max_subsample = max_iter - 2;
      const int cap = 4;
      subsample_iters = std::min(max_subsample, std::min(stride_eff, cap));
      subsample_iters = std::max(1, subsample_iters);
    }
  }
  const bool use_subsample = (subsample_iters > 0);

  auto refresh_center_quant = [&]() {
    if (!quantize_assign) return;
    if (center_q.size() != static_cast<size_t>(K) * static_cast<size_t>(d)) {
      center_q.resize(static_cast<size_t>(K) * static_cast<size_t>(d));
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int k = 0; k < K; ++k) {
      const float *src = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
      int8_t *dst = center_q.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
      for (int j = 0; j < d; ++j) {
        dst[j] = quantize_unit_i8(src[j]);
      }
    }
  };
  refresh_center_quant();

  auto assignment_pass = [&](const bool subsample, const int phase,
                             int64_t &changes, int64_t &processed) {
    changes = 0;
    processed = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : changes, processed) schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      if (subsample && (G.vz[v] % stride_eff) != phase) {
        new_labels[v] = labels[v];
        continue;
      }
      processed += 1;

      const float x = static_cast<float>(G.vx[v]);
      const float y = static_cast<float>(G.vy[v]);
      const float z = static_cast<float>(G.vz[v]);
      const int bid0 = B.voxel_bin[v];
      const uint8_t nnb = B.nbr_count[static_cast<size_t>(bid0)];
      const int32_t *nbr = B.nbr_bin.data() + static_cast<size_t>(bid0) * 27U;

      int bestk = labels[v];
      float bestD = std::numeric_limits<float>::infinity();

      if (!quantize_assign) {
        const float *fv = feat_emb.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
        int nchecked = 0;
        for (uint8_t ni = 0; ni < nnb; ++ni) {
          const int bid = nbr[static_cast<size_t>(ni)];
          const int start = B.bin_offsets[static_cast<size_t>(bid)];
          const int end = B.bin_offsets[static_cast<size_t>(bid + 1)];
          for (int idx = start; idx < end; ++idx) {
            const int k = B.bin_centers[static_cast<size_t>(idx)];
            ++nchecked;
            const float dxs = x - C.cx[k];
            const float dys = y - C.cy[k];
            const float dzs = z - C.cz[k];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const float *fk = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
            const float dot = dot_row(fv, fk, d);
            const float D = (1.0f - dot) + alpha_over_s2 * dist2;

            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }

        // Fallback for pathological empty candidate list.
        if (nchecked == 0) {
          for (int k = 0; k < K; ++k) {
            const float dxs = x - C.cx[k];
            const float dys = y - C.cy[k];
            const float dzs = z - C.cz[k];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const float *fk = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
            const float dot = dot_row(fv, fk, d);
            const float D = (1.0f - dot) + alpha_over_s2 * dist2;
            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }
      } else {
        const int8_t *fvq = feat_emb_q.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
        int nchecked = 0;
        for (uint8_t ni = 0; ni < nnb; ++ni) {
          const int bid = nbr[static_cast<size_t>(ni)];
          const int start = B.bin_offsets[static_cast<size_t>(bid)];
          const int end = B.bin_offsets[static_cast<size_t>(bid + 1)];
          for (int idx = start; idx < end; ++idx) {
            const int k = B.bin_centers[static_cast<size_t>(idx)];
            ++nchecked;
            const float dxs = x - C.cx[k];
            const float dys = y - C.cy[k];
            const float dzs = z - C.cz[k];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const int8_t *fkq = center_q.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
            const float dotq = static_cast<float>(dot_row_i8(fvq, fkq, d)) * inv_qscale2;
            const float Dq = (1.0f - dotq) + alpha_over_s2 * dist2;

            if (Dq < bestD) {
              bestD = Dq;
              bestk = k;
            }
          }
        }

        // Fallback for pathological empty candidate list.
        if (nchecked == 0) {
          for (int k = 0; k < K; ++k) {
            const float dxs = x - C.cx[k];
            const float dys = y - C.cy[k];
            const float dzs = z - C.cz[k];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const int8_t *fkq = center_q.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
            const float dotq = static_cast<float>(dot_row_i8(fvq, fkq, d)) * inv_qscale2;
            const float Dq = (1.0f - dotq) + alpha_over_s2 * dist2;
            if (Dq < bestD) {
              bestD = Dq;
              bestk = k;
            }
          }
        }
      }

      new_labels[v] = bestk;
      if (bestk != labels[v]) ++changes;
    }

    labels.swap(new_labels);
  };

  auto update_centers = [&]() {
    std::fill(count_total.begin(), count_total.end(), 0);
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      int32_t *cnt = count.data() + static_cast<size_t>(tid) * thread_cluster_stride;
      float *sxyz = sum_xyz.data() + static_cast<size_t>(tid) * thread_xyz_stride;
      float *sf = sum_f.data() + static_cast<size_t>(tid) * thread_feat_stride;
      std::fill(cnt, cnt + thread_cluster_stride, 0);
      std::fill(sxyz, sxyz + thread_xyz_stride, 0.0f);
      std::fill(sf, sf + thread_feat_stride, 0.0f);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int v = 0; v < N; ++v) {
        int k = labels[v];
        cnt[k] += 1;
        sxyz[static_cast<size_t>(k) * 3U + 0U] += static_cast<float>(G.vx[v]);
        sxyz[static_cast<size_t>(k) * 3U + 1U] += static_cast<float>(G.vy[v]);
        sxyz[static_cast<size_t>(k) * 3U + 2U] += static_cast<float>(G.vz[v]);

        const float *fv = feat_emb.data() + static_cast<size_t>(v) * static_cast<size_t>(d);
        float *acc = sf + static_cast<size_t>(k) * static_cast<size_t>(d);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          acc[j] += fv[j];
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int k = 0; k < K; ++k) {
      const size_t kd = static_cast<size_t>(k) * static_cast<size_t>(d);
      const int32_t *cntk = count.data() + static_cast<size_t>(k);
      const float *sxyzk = sum_xyz.data() + static_cast<size_t>(k) * 3U;
      const float *sfk = sum_f.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
      float *dst = C.cf.data() + kd;
      int32_t total = 0;
      float sx = 0.0f;
      float sy = 0.0f;
      float sz = 0.0f;

      if (nt == 1) {
        total = cntk[0];
        sx = sxyzk[0];
        sy = sxyzk[1];
        sz = sxyzk[2];
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          dst[j] = sfk[j];
        }
      } else if (nt == 2) {
        const int32_t *cnt1 = cntk + thread_cluster_stride;
        const float *xyz1 = sxyzk + thread_xyz_stride;
        const float *f1 = sfk + thread_feat_stride;
        total = cntk[0] + cnt1[0];
        sx = sxyzk[0] + xyz1[0];
        sy = sxyzk[1] + xyz1[1];
        sz = sxyzk[2] + xyz1[2];
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          dst[j] = sfk[j] + f1[j];
        }
      } else if (nt == 4) {
        const int32_t *cnt1 = cntk + thread_cluster_stride;
        const int32_t *cnt2 = cnt1 + thread_cluster_stride;
        const int32_t *cnt3 = cnt2 + thread_cluster_stride;
        const float *xyz1 = sxyzk + thread_xyz_stride;
        const float *xyz2 = xyz1 + thread_xyz_stride;
        const float *xyz3 = xyz2 + thread_xyz_stride;
        const float *f1 = sfk + thread_feat_stride;
        const float *f2 = f1 + thread_feat_stride;
        const float *f3 = f2 + thread_feat_stride;
        total = cntk[0] + cnt1[0] + cnt2[0] + cnt3[0];
        sx = sxyzk[0] + xyz1[0] + xyz2[0] + xyz3[0];
        sy = sxyzk[1] + xyz1[1] + xyz2[1] + xyz3[1];
        sz = sxyzk[2] + xyz1[2] + xyz2[2] + xyz3[2];
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          dst[j] = sfk[j] + f1[j] + f2[j] + f3[j];
        }
      } else {
        total = cntk[0];
        sx = sxyzk[0];
        sy = sxyzk[1];
        sz = sxyzk[2];
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          dst[j] = sfk[j];
        }
        for (int t = 1; t < nt; ++t) {
          cntk += thread_cluster_stride;
          sxyzk += thread_xyz_stride;
          sfk += thread_feat_stride;
          total += cntk[0];
          sx += sxyzk[0];
          sy += sxyzk[1];
          sz += sxyzk[2];
#ifdef _OPENMP
#pragma omp simd
#endif
          for (int j = 0; j < d; ++j) {
            dst[j] += sfk[j];
          }
        }
      }
      count_total[static_cast<size_t>(k)] = total;

      if (total > 0) {
        const float inv = 1.0f / static_cast<float>(total);
        C.cx[k] = sx * inv;
        C.cy[k] = sy * inv;
        C.cz[k] = sz * inv;

        float n2 = 0.0f;
#ifdef _OPENMP
#pragma omp simd reduction(+ : n2)
#endif
        for (int j = 0; j < d; ++j) {
          const float v = dst[j] * inv;
          dst[j] = v;
          n2 += v * v;
        }
        const float invn_center = 1.0f / std::sqrt(n2 + 1e-8f);
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int j = 0; j < d; ++j) {
          dst[j] *= invn_center;
        }
      }
    }

    for (int k = 0; k < K; ++k) {
      const int32_t total = count_total[static_cast<size_t>(k)];
      if (total > 0) continue;
      int vid = static_cast<int>(splitmix64(rng_state) % static_cast<uint64_t>(N));
      C.cx[k] = static_cast<float>(G.vx[vid]);
      C.cy[k] = static_cast<float>(G.vy[vid]);
      C.cz[k] = static_cast<float>(G.vz[vid]);
      const float *src = feat_emb.data() + static_cast<size_t>(vid) * static_cast<size_t>(d);
      float *dst = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
      std::copy(src, src + d, dst);
    }

    refresh_center_quant();
  };

  bool did_full_iter = false;
  for (int it = 0; it < max_iter; ++it) {
    rebuild_center_bins(B, C, S);
    int64_t changes = 0;
    int64_t processed = 0;
    const bool subsample_now = use_subsample && (it >= 1) && (it < (1 + subsample_iters));
    const int phase = (subsample_now ? ((it - 1) % stride_eff) : 0);
    assignment_pass(subsample_now, phase, changes, processed);
    update_centers();

    double change_ratio = static_cast<double>(changes) / static_cast<double>(N);
    if (verbose) {
      Rcout << "corrslic_core iter " << (it + 1) << ": changes=" << changes
            << " (" << change_ratio << ")"
            << " processed=" << processed
            << " mode=" << (subsample_now ? "subsample" : "full") << "\n";
    }

    if (subsample_now) {
      continue;
    }

    did_full_iter = true;
    if (it >= 2 && (changes == 0 || change_ratio < 0.001)) {
      break;
    }
  }

  if (use_subsample && !did_full_iter) {
    rebuild_center_bins(B, C, S);
    int64_t final_changes = 0;
    int64_t final_processed = 0;
    assignment_pass(false, 0, final_changes, final_processed);
    update_centers();
    if (verbose) {
      Rcout << "corrslic_core final full pass: changes=" << final_changes
            << " (" << (static_cast<double>(final_changes) / static_cast<double>(N)) << ")"
            << " processed=" << final_processed << "\n";
    }
  }

  if (refine_exact_iters > 0) {
    const float alpha_refine = static_cast<float>(refine_alpha > 0.0 ? refine_alpha : alpha);
    const float alpha_refine_over_s2 =
      alpha_refine / (static_cast<float>(S) * static_cast<float>(S));

    // Exact correlation refinement works on per-voxel z-normalized original time series.
    const int stride_t = std::max(1, std::min(refine_stride, T));
    const int Td = (T + stride_t - 1) / stride_t;
    const double *feat_ptr_exact = feat.begin();
    std::vector<float> feat_exact(static_cast<size_t>(N) * static_cast<size_t>(Td), 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      float mean = 0.0f;
      for (int tt = 0, t = 0; t < T; t += stride_t, ++tt) {
        mean += static_cast<float>(feat_ptr_exact[static_cast<size_t>(v) + static_cast<size_t>(t) * static_cast<size_t>(N)]);
      }
      mean /= static_cast<float>(Td);
      float n2 = 0.0f;
      float *dst = feat_exact.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);
      for (int tt = 0, t = 0; t < T; t += stride_t, ++tt) {
        const float zt = static_cast<float>(feat_ptr_exact[static_cast<size_t>(v) + static_cast<size_t>(t) * static_cast<size_t>(N)]) - mean;
        dst[static_cast<size_t>(tt)] = zt;
        n2 += zt * zt;
      }
      const float invn_t = 1.0f / std::sqrt(n2 + 1e-8f);
      for (int tt = 0; tt < Td; ++tt) {
        dst[static_cast<size_t>(tt)] *= invn_t;
      }
    }

    // Quantized dot products for exact refinement were slower in practice with
    // the current portable implementation (no dedicated int8 dot-product
    // intrinsics). Keep refinement in float for now.
    const bool quantize_refine = false;
    std::vector<int8_t> feat_exact_q;
    if (quantize_refine) {
      feat_exact_q.resize(static_cast<size_t>(N) * static_cast<size_t>(Td));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        const float *src = feat_exact.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);
        int8_t *dst = feat_exact_q.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);
        for (int tt = 0; tt < Td; ++tt) {
          dst[static_cast<size_t>(tt)] = quantize_unit_i8(src[static_cast<size_t>(tt)]);
        }
      }
    }

    std::vector<float> center_exact(static_cast<size_t>(K) * static_cast<size_t>(Td), 0.0f);
    std::vector<int8_t> center_exact_q;
    if (quantize_refine) {
      center_exact_q.resize(static_cast<size_t>(K) * static_cast<size_t>(Td));
    }
    std::vector<int32_t> count_exact(static_cast<size_t>(K), 0);
    std::vector<float> sx_exact(static_cast<size_t>(K), 0.0f);
    std::vector<float> sy_exact(static_cast<size_t>(K), 0.0f);
    std::vector<float> sz_exact(static_cast<size_t>(K), 0.0f);

    // Thread-local buffers to parallelize exact center updates.
    // Memory: O(threads * K * T) floats, which is typically small for the
    // synthetic benchmarks and moderate for real fMRI use cases (refinement is
    // optional; keep it off for very large K*T if needed).
    std::vector<int32_t> count_exact_tl(static_cast<size_t>(nt) * static_cast<size_t>(K), 0);
    std::vector<float> sum_xyz_exact_tl(static_cast<size_t>(nt) * static_cast<size_t>(K) * 3U, 0.0f);
    std::vector<float> sum_exact_tl(static_cast<size_t>(nt) * static_cast<size_t>(K) * static_cast<size_t>(Td), 0.0f);

    auto update_exact_centers = [&]() {
      std::fill(count_exact_tl.begin(), count_exact_tl.end(), 0);
      std::fill(sum_xyz_exact_tl.begin(), sum_xyz_exact_tl.end(), 0.0f);
      std::fill(sum_exact_tl.begin(), sum_exact_tl.end(), 0.0f);

      // Accumulate per-thread (no atomics) then reduce.
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
      {
        int tid = 0;
        int ntt = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        ntt = omp_get_num_threads();
#endif
        const int start = static_cast<int>((static_cast<int64_t>(tid) * N) / ntt);
        const int end = static_cast<int>((static_cast<int64_t>(tid + 1) * N) / ntt);

        int32_t *cnt_tl = count_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K);
        float *xyz_tl = sum_xyz_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K) * 3U;
        float *sum_tl = sum_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K) * static_cast<size_t>(Td);

        for (int v = start; v < end; ++v) {
          const int k = labels[static_cast<size_t>(v)];
          cnt_tl[static_cast<size_t>(k)] += 1;
          const size_t k3 = static_cast<size_t>(k) * 3U;
          xyz_tl[k3 + 0] += static_cast<float>(G.vx[static_cast<size_t>(v)]);
          xyz_tl[k3 + 1] += static_cast<float>(G.vy[static_cast<size_t>(v)]);
          xyz_tl[k3 + 2] += static_cast<float>(G.vz[static_cast<size_t>(v)]);

          const float *src = feat_exact.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);
          float *dst = sum_tl + static_cast<size_t>(k) * static_cast<size_t>(Td);
#ifdef _OPENMP
#pragma omp simd
#endif
          for (int tt = 0; tt < Td; ++tt) {
            dst[static_cast<size_t>(tt)] += src[static_cast<size_t>(tt)];
          }
        }
      }

      // Reduce thread-local buffers into final centers and normalize.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
      for (int k = 0; k < K; ++k) {
        int32_t cnt = 0;
        float sx = 0.0f, sy = 0.0f, sz = 0.0f;
        float *dst = center_exact.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
        std::fill(dst, dst + Td, 0.0f);

        for (int tid = 0; tid < nt; ++tid) {
          const int32_t *cnt_tl = count_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K);
          const float *xyz_tl = sum_xyz_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K) * 3U;
          const float *sum_tl = sum_exact_tl.data() + static_cast<size_t>(tid) * static_cast<size_t>(K) * static_cast<size_t>(Td);

          cnt += cnt_tl[static_cast<size_t>(k)];
          const size_t k3 = static_cast<size_t>(k) * 3U;
          sx += xyz_tl[k3 + 0];
          sy += xyz_tl[k3 + 1];
          sz += xyz_tl[k3 + 2];

          const float *src = sum_tl + static_cast<size_t>(k) * static_cast<size_t>(Td);
#ifdef _OPENMP
#pragma omp simd
#endif
          for (int tt = 0; tt < Td; ++tt) {
            dst[static_cast<size_t>(tt)] += src[static_cast<size_t>(tt)];
          }
        }

        count_exact[static_cast<size_t>(k)] = cnt;
        sx_exact[static_cast<size_t>(k)] = sx;
        sy_exact[static_cast<size_t>(k)] = sy;
        sz_exact[static_cast<size_t>(k)] = sz;

        if (cnt > 0) {
          const float inv = 1.0f / static_cast<float>(cnt);
          C.cx[static_cast<size_t>(k)] = sx * inv;
          C.cy[static_cast<size_t>(k)] = sy * inv;
          C.cz[static_cast<size_t>(k)] = sz * inv;

          float n2 = 0.0f;
          for (int tt = 0; tt < Td; ++tt) {
            dst[static_cast<size_t>(tt)] *= inv;
            n2 += dst[static_cast<size_t>(tt)] * dst[static_cast<size_t>(tt)];
          }
          const float invn_t = 1.0f / std::sqrt(n2 + 1e-8f);
          for (int tt = 0; tt < Td; ++tt) {
            dst[static_cast<size_t>(tt)] *= invn_t;
          }
        } else {
          // Deterministic reinit (avoid shared RNG in parallel region).
          uint64_t local_state = rng_state ^ (0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(k + 1));
          const int vid = static_cast<int>(splitmix64(local_state) % static_cast<uint64_t>(N));
          C.cx[static_cast<size_t>(k)] = static_cast<float>(G.vx[static_cast<size_t>(vid)]);
          C.cy[static_cast<size_t>(k)] = static_cast<float>(G.vy[static_cast<size_t>(vid)]);
          C.cz[static_cast<size_t>(k)] = static_cast<float>(G.vz[static_cast<size_t>(vid)]);
          const float *src = feat_exact.data() + static_cast<size_t>(vid) * static_cast<size_t>(Td);
          std::copy(src, src + Td, dst);
        }

        if (quantize_refine) {
          int8_t *dstq = center_exact_q.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
          for (int tt = 0; tt < Td; ++tt) {
            dstq[static_cast<size_t>(tt)] = quantize_unit_i8(dst[static_cast<size_t>(tt)]);
          }
        }
      }
    };

    auto mark_boundary_voxels = [&](std::vector<uint8_t> &boundary) {
      boundary.assign(static_cast<size_t>(N), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        const int x = static_cast<int>(G.vx[static_cast<size_t>(v)]);
        const int y = static_cast<int>(G.vy[static_cast<size_t>(v)]);
        const int z = static_cast<int>(G.vz[static_cast<size_t>(v)]);
        const int lv = labels[static_cast<size_t>(v)];
        bool is_boundary = false;

        if (x > 0) {
          const int u = G.id_map[lin_index(x - 1, y, z, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }
        if (!is_boundary && x + 1 < G.X) {
          const int u = G.id_map[lin_index(x + 1, y, z, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }
        if (!is_boundary && y > 0) {
          const int u = G.id_map[lin_index(x, y - 1, z, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }
        if (!is_boundary && y + 1 < G.Y) {
          const int u = G.id_map[lin_index(x, y + 1, z, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }
        if (!is_boundary && z > 0) {
          const int u = G.id_map[lin_index(x, y, z - 1, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }
        if (!is_boundary && z + 1 < G.Z) {
          const int u = G.id_map[lin_index(x, y, z + 1, G.X, G.Y)];
          if (u >= 0 && labels[static_cast<size_t>(u)] != lv) is_boundary = true;
        }

        boundary[static_cast<size_t>(v)] = static_cast<uint8_t>(is_boundary ? 1 : 0);
      }
    };

    auto assignment_exact_pass_f = [&](const std::vector<uint8_t> *boundary_mask,
                                       int64_t &changes, int64_t &processed) {
      changes = 0;
      processed = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : changes, processed) schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        if (boundary_mask && ((*boundary_mask)[static_cast<size_t>(v)] == 0)) {
          new_labels[static_cast<size_t>(v)] = labels[static_cast<size_t>(v)];
          continue;
        }
        processed += 1;

        const float x = static_cast<float>(G.vx[static_cast<size_t>(v)]);
        const float y = static_cast<float>(G.vy[static_cast<size_t>(v)]);
        const float z = static_cast<float>(G.vz[static_cast<size_t>(v)]);
        const int bid0 = B.voxel_bin[static_cast<size_t>(v)];
        const uint8_t nnb = B.nbr_count[static_cast<size_t>(bid0)];
        const int32_t *nbr = B.nbr_bin.data() + static_cast<size_t>(bid0) * 27U;
        const float *fv = feat_exact.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);

        int bestk = labels[static_cast<size_t>(v)];
        float bestD = std::numeric_limits<float>::infinity();
        int nchecked = 0;
        for (uint8_t ni = 0; ni < nnb; ++ni) {
          const int bid = nbr[static_cast<size_t>(ni)];
          const int start = B.bin_offsets[static_cast<size_t>(bid)];
          const int end = B.bin_offsets[static_cast<size_t>(bid + 1)];
          for (int idx = start; idx < end; ++idx) {
            const int k = B.bin_centers[static_cast<size_t>(idx)];
            ++nchecked;
            const float dxs = x - C.cx[static_cast<size_t>(k)];
            const float dys = y - C.cy[static_cast<size_t>(k)];
            const float dzs = z - C.cz[static_cast<size_t>(k)];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const float *fk = center_exact.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
            const float dot = dot_row(fv, fk, Td);
            const float D = (1.0f - dot) + alpha_refine_over_s2 * dist2;
            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }

        if (nchecked == 0) {
          for (int k = 0; k < K; ++k) {
            const float dxs = x - C.cx[static_cast<size_t>(k)];
            const float dys = y - C.cy[static_cast<size_t>(k)];
            const float dzs = z - C.cz[static_cast<size_t>(k)];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const float *fk = center_exact.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
            const float dot = dot_row(fv, fk, Td);
            const float D = (1.0f - dot) + alpha_refine_over_s2 * dist2;
            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }

        new_labels[static_cast<size_t>(v)] = bestk;
        if (bestk != labels[static_cast<size_t>(v)]) ++changes;
      }

      labels.swap(new_labels);
    };

    auto assignment_exact_pass_q = [&](const std::vector<uint8_t> *boundary_mask,
                                       int64_t &changes, int64_t &processed) {
      changes = 0;
      processed = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : changes, processed) schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        if (boundary_mask && ((*boundary_mask)[static_cast<size_t>(v)] == 0)) {
          new_labels[static_cast<size_t>(v)] = labels[static_cast<size_t>(v)];
          continue;
        }
        processed += 1;

        const float x = static_cast<float>(G.vx[static_cast<size_t>(v)]);
        const float y = static_cast<float>(G.vy[static_cast<size_t>(v)]);
        const float z = static_cast<float>(G.vz[static_cast<size_t>(v)]);
        const int bid0 = B.voxel_bin[static_cast<size_t>(v)];
        const uint8_t nnb = B.nbr_count[static_cast<size_t>(bid0)];
        const int32_t *nbr = B.nbr_bin.data() + static_cast<size_t>(bid0) * 27U;
        const int8_t *fvq = feat_exact_q.data() + static_cast<size_t>(v) * static_cast<size_t>(Td);

        int bestk = labels[static_cast<size_t>(v)];
        float bestD = std::numeric_limits<float>::infinity();
        int nchecked = 0;
        for (uint8_t ni = 0; ni < nnb; ++ni) {
          const int bid = nbr[static_cast<size_t>(ni)];
          const int start = B.bin_offsets[static_cast<size_t>(bid)];
          const int end = B.bin_offsets[static_cast<size_t>(bid + 1)];
          for (int idx = start; idx < end; ++idx) {
            const int k = B.bin_centers[static_cast<size_t>(idx)];
            ++nchecked;
            const float dxs = x - C.cx[static_cast<size_t>(k)];
            const float dys = y - C.cy[static_cast<size_t>(k)];
            const float dzs = z - C.cz[static_cast<size_t>(k)];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const int8_t *fkq = center_exact_q.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
            const int32_t dot_i8 = dot_row_i8(fvq, fkq, Td);
            const float dot = inv_qscale2 * static_cast<float>(dot_i8);
            const float D = (1.0f - dot) + alpha_refine_over_s2 * dist2;
            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }

        if (nchecked == 0) {
          for (int k = 0; k < K; ++k) {
            const float dxs = x - C.cx[static_cast<size_t>(k)];
            const float dys = y - C.cy[static_cast<size_t>(k)];
            const float dzs = z - C.cz[static_cast<size_t>(k)];
            const float dist2 = dxs * dxs + dys * dys + dzs * dzs;
            const int8_t *fkq = center_exact_q.data() + static_cast<size_t>(k) * static_cast<size_t>(Td);
            const int32_t dot_i8 = dot_row_i8(fvq, fkq, Td);
            const float dot = inv_qscale2 * static_cast<float>(dot_i8);
            const float D = (1.0f - dot) + alpha_refine_over_s2 * dist2;
            if (D < bestD) {
              bestD = D;
              bestk = k;
            }
          }
        }

        new_labels[static_cast<size_t>(v)] = bestk;
        if (bestk != labels[static_cast<size_t>(v)]) ++changes;
      }

      labels.swap(new_labels);
    };

    update_exact_centers();
    std::vector<uint8_t> boundary_mask;
    for (int it = 0; it < refine_exact_iters; ++it) {
      rebuild_center_bins(B, C, S);
      const std::vector<uint8_t> *bm = nullptr;
      if (refine_boundary_only) {
        mark_boundary_voxels(boundary_mask);
        bm = &boundary_mask;
      }
      int64_t changes = 0;
      int64_t processed = 0;
      if (quantize_refine) {
        assignment_exact_pass_q(bm, changes, processed);
      } else {
        assignment_exact_pass_f(bm, changes, processed);
      }
      if (verbose) {
        const double cr = static_cast<double>(changes) / static_cast<double>(N);
        Rcout << "corrslic_core refine " << (it + 1)
              << ": changes=" << changes
              << " (" << cr << ")"
              << " processed=" << processed
              << " mode=" << (refine_boundary_only ? "boundary" : "full")
              << "\n";
      }
      // Only update prototypes if we will run another refinement pass.
      if (changes == 0 || (it + 1) >= refine_exact_iters) break;
      update_exact_centers();
    }

    // Refresh embedding-space centers to match refined labels for outputs/metadata.
    update_centers();
  }

  int min_size_eff = min_size;
  if (min_size_eff <= 0) {
    int approx = static_cast<int>(std::round(0.25 * S * S * S));
    min_size_eff = std::max(10, approx);
  }

  enforce_connectivity(G, labels, K, min_size_eff, connectivity);

  IntegerVector out_labels(N);
  for (int i = 0; i < N; ++i) {
    out_labels[i] = labels[i] + 1;
  }

  NumericMatrix centers_xyz(K, 3);
  for (int k = 0; k < K; ++k) {
    centers_xyz(k, 0) = C.cx[k] + 1.0;
    centers_xyz(k, 1) = C.cy[k] + 1.0;
    centers_xyz(k, 2) = C.cz[k] + 1.0;
  }

  NumericMatrix center_feats(K, d);
  for (int k = 0; k < K; ++k) {
    const float *src = C.cf.data() + static_cast<size_t>(k) * static_cast<size_t>(d);
    for (int j = 0; j < d; ++j) {
      center_feats(k, j) = src[j];
    }
  }

  return List::create(
    _["labels"] = out_labels,
    _["centers"] = center_feats,
    _["centers_xyz"] = centers_xyz,
    _["params"] = List::create(
      _["K"] = K,
      _["d"] = d,
      _["sketch_repeats"] = sketch_repeats,
      _["alpha"] = alpha,
      _["assign_stride"] = assign_stride,
      _["assign_subsample_iters"] = subsample_iters,
      _["quantize_assign"] = quantize_assign,
      _["embed_basis"] = embed_basis,
      _["whiten_embed"] = whiten_embed,
      _["refine_exact_iters"] = refine_exact_iters,
      _["refine_boundary_only"] = refine_boundary_only,
      _["refine_stride"] = refine_stride,
      _["refine_alpha"] = (refine_alpha > 0.0 ? refine_alpha : alpha),
      _["max_iter"] = max_iter,
      _["seed"] = seed,
      _["connectivity"] = connectivity,
      _["min_size"] = min_size_eff
    )
  );
}

// [[Rcpp::export]]
List brs_slic_core(const NumericMatrix feat,
                   const IntegerVector mask_lin_idx,
                   const IntegerVector dims,
                   const int K,
                   const int d = 32,
                   const int sketch_repeats = 1,
                   const double alpha = 0.05,
                   const int coarse_iter = 3,
                   const int boundary_passes = 2,
                   const int global_passes = 0,
                   const double refine_spatial = 0.02,
                   const double refine_l2 = 0.08,
                   int refine_stride = 0,
                   const int seed = 1,
                   const int connectivity = 6,
                   const int min_size = 0,
                   int n_threads = 0,
                   const bool verbose = false) {
  const int N = feat.nrow();
  const int T = feat.ncol();

  if (N <= 0) stop("feat has zero rows");
  if (T <= 1) stop("feat must have at least 2 timepoints");
  if (mask_lin_idx.size() != N) stop("mask_lin_idx length must match nrow(feat)");
  if (K < 2 || K > N) stop("K must be in [2, N]");
  if (d < 8) stop("embedding dimension d must be >= 8");
  if (coarse_iter < 1) stop("coarse_iter must be >= 1");
  if (boundary_passes < 0) stop("boundary_passes must be >= 0");
  if (global_passes < 0) stop("global_passes must be >= 0");
  if (refine_l2 < 0) stop("refine_l2 must be >= 0");
  if (refine_stride < 0) stop("refine_stride must be >= 0");
  if (connectivity != 6 && connectivity != 26) stop("connectivity must be 6 or 26");

#ifdef _OPENMP
  if (n_threads <= 0) n_threads = omp_get_max_threads();
#else
  n_threads = 1;
#endif

  // Coarse pass: fast sketch-SLIC.
  List coarse = corrslic_core(
    feat, mask_lin_idx, dims, K, d, sketch_repeats, alpha, coarse_iter, seed,
    1, false, "hash", false, 0, true, 1, -1.0, connectivity, 0, n_threads, false
  );

  IntegerVector coarse_labels = coarse["labels"];
  NumericMatrix coarse_centers = coarse["centers"];
  NumericMatrix coarse_centers_xyz = coarse["centers_xyz"];
  std::vector<int32_t> labels(static_cast<size_t>(N), 0);
  for (int i = 0; i < N; ++i) labels[static_cast<size_t>(i)] = static_cast<int32_t>(coarse_labels[i] - 1);

  VoxelGrid G = build_voxel_grid(mask_lin_idx, dims);
  const double *feat_ptr = REAL(feat);
  const int refine_stride_eff = (refine_stride > 0 ? refine_stride : (T >= 40 ? 2 : 1));
  const int Teff = (T + refine_stride_eff - 1) / refine_stride_eff;

  const double S = std::max(1.0, std::cbrt(static_cast<double>(N) / static_cast<double>(K)));
  const double refine_alpha_over_s2 = refine_spatial / (S * S);
  const double refine_l2_weight = refine_l2;
  if (boundary_passes > 0 || global_passes > 0) {
    // Voxel-major decimated time-series buffer for cache-friendly refinement loops.
    std::vector<float> xsub(static_cast<size_t>(N) * static_cast<size_t>(Teff), 0.0f);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      float *xv = xsub.data() + static_cast<size_t>(v) * static_cast<size_t>(Teff);
      for (int tt = 0; tt < Teff; ++tt) {
        const int t = tt * refine_stride_eff;
        xv[tt] = static_cast<float>(
          feat_ptr[static_cast<size_t>(v) + static_cast<size_t>(t) * static_cast<size_t>(N)]
        );
      }
    }

    // Voxel-wise raw time-series stats for exact Pearson.
    std::vector<double> sum_x(static_cast<size_t>(N), 0.0);
    std::vector<double> sum_x2(static_cast<size_t>(N), 0.0);
    std::vector<double> mean_x(static_cast<size_t>(N), 0.0);
    std::vector<double> invsd_x(static_cast<size_t>(N), 0.0);
    const double eps = 1e-10;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int v = 0; v < N; ++v) {
      double sx = 0.0, sx2 = 0.0;
      const float *xv = xsub.data() + static_cast<size_t>(v) * static_cast<size_t>(Teff);
      for (int tt = 0; tt < Teff; ++tt) {
        const double x = static_cast<double>(xv[tt]);
        sx += x;
        sx2 += x * x;
      }
      const double mx = sx / static_cast<double>(Teff);
      const double ss = std::max(0.0, sx2 - static_cast<double>(Teff) * mx * mx);
      sum_x[static_cast<size_t>(v)] = sx;
      sum_x2[static_cast<size_t>(v)] = sx2;
      mean_x[static_cast<size_t>(v)] = mx;
      invsd_x[static_cast<size_t>(v)] = 1.0 / std::sqrt(ss + eps);
    }

    std::vector<double> proto(static_cast<size_t>(K) * static_cast<size_t>(Teff), 0.0);
    std::vector<double> proto_sum(static_cast<size_t>(K), 0.0);
    std::vector<double> proto_mean(static_cast<size_t>(K), 0.0);
    std::vector<double> proto_invsd(static_cast<size_t>(K), 0.0);
    std::vector<double> proto_energy(static_cast<size_t>(K), 0.0);
    std::vector<double> cx(static_cast<size_t>(K), 0.0);
    std::vector<double> cy(static_cast<size_t>(K), 0.0);
    std::vector<double> cz(static_cast<size_t>(K), 0.0);
    std::vector<int32_t> counts(static_cast<size_t>(K), 0);
    std::vector<int32_t> new_labels(static_cast<size_t>(N), 0);

    auto recompute_cluster_stats = [&]() {
      std::fill(proto.begin(), proto.end(), 0.0);
      std::fill(proto_sum.begin(), proto_sum.end(), 0.0);
      std::fill(proto_mean.begin(), proto_mean.end(), 0.0);
      std::fill(proto_invsd.begin(), proto_invsd.end(), 0.0);
      std::fill(proto_energy.begin(), proto_energy.end(), 0.0);
      std::fill(cx.begin(), cx.end(), 0.0);
      std::fill(cy.begin(), cy.end(), 0.0);
      std::fill(cz.begin(), cz.end(), 0.0);
      std::fill(counts.begin(), counts.end(), 0);

      for (int v = 0; v < N; ++v) {
        const int k = labels[static_cast<size_t>(v)];
        if (k < 0 || k >= K) continue;
        counts[static_cast<size_t>(k)] += 1;
        cx[static_cast<size_t>(k)] += static_cast<double>(G.vx[static_cast<size_t>(v)]);
        cy[static_cast<size_t>(k)] += static_cast<double>(G.vy[static_cast<size_t>(v)]);
        cz[static_cast<size_t>(k)] += static_cast<double>(G.vz[static_cast<size_t>(v)]);
      }

      // Traverse feature matrix by columns (R storage order) for cache locality.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
      for (int tt = 0; tt < Teff; ++tt) {
        const int t = tt * refine_stride_eff;
        const double *col = feat_ptr + static_cast<size_t>(t) * static_cast<size_t>(N);
        for (int v = 0; v < N; ++v) {
          const int k = labels[static_cast<size_t>(v)];
          if (k < 0 || k >= K) continue;
          proto[static_cast<size_t>(k) * static_cast<size_t>(Teff) + static_cast<size_t>(tt)] += col[static_cast<size_t>(v)];
        }
      }

      for (int k = 0; k < K; ++k) {
        const int32_t c = counts[static_cast<size_t>(k)];
        if (c <= 0) continue;
        const double invc = 1.0 / static_cast<double>(c);
        cx[static_cast<size_t>(k)] *= invc;
        cy[static_cast<size_t>(k)] *= invc;
        cz[static_cast<size_t>(k)] *= invc;

        double *pk = proto.data() + static_cast<size_t>(k) * static_cast<size_t>(Teff);
        double sp = 0.0, sp2 = 0.0;
        for (int tt = 0; tt < Teff; ++tt) {
          pk[tt] *= invc;
          sp += pk[tt];
          sp2 += pk[tt] * pk[tt];
        }
        const double mp = sp / static_cast<double>(Teff);
        const double ss = std::max(0.0, sp2 - static_cast<double>(Teff) * mp * mp);
        proto_sum[static_cast<size_t>(k)] = sp;
        proto_mean[static_cast<size_t>(k)] = mp;
        proto_invsd[static_cast<size_t>(k)] = 1.0 / std::sqrt(ss + eps);
        proto_energy[static_cast<size_t>(k)] = sp2;
      }
    };

    auto add_candidate = [](int lab, int *cand, int &n, const int cap) {
      if (lab < 0 || n >= cap) return;
      for (int i = 0; i < n; ++i) {
        if (cand[i] == lab) return;
      }
      cand[n++] = lab;
    };

    for (int bp = 0; bp < boundary_passes; ++bp) {
      recompute_cluster_stats();
      int64_t changes = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : changes) schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        const int cur = labels[static_cast<size_t>(v)];
        new_labels[static_cast<size_t>(v)] = cur;

        const int x = G.vx[static_cast<size_t>(v)];
        const int y = G.vy[static_cast<size_t>(v)];
        const int z = G.vz[static_cast<size_t>(v)];

        int cand[7];
        int nc = 0;
        add_candidate(cur, cand, nc, 7);
        bool is_boundary = false;

        if (x > 0) {
          int u = G.id_map[static_cast<size_t>(lin_index(x - 1, y, z, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }
        if (x + 1 < G.X) {
          int u = G.id_map[static_cast<size_t>(lin_index(x + 1, y, z, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }
        if (y > 0) {
          int u = G.id_map[static_cast<size_t>(lin_index(x, y - 1, z, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }
        if (y + 1 < G.Y) {
          int u = G.id_map[static_cast<size_t>(lin_index(x, y + 1, z, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }
        if (z > 0) {
          int u = G.id_map[static_cast<size_t>(lin_index(x, y, z - 1, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }
        if (z + 1 < G.Z) {
          int u = G.id_map[static_cast<size_t>(lin_index(x, y, z + 1, G.X, G.Y))];
          if (u >= 0) {
            const int lu = labels[static_cast<size_t>(u)];
            if (lu != cur) is_boundary = true;
            add_candidate(lu, cand, nc, 7);
          }
        }

        if (!is_boundary || nc <= 1) continue;

        const double sx = sum_x[static_cast<size_t>(v)];
        const double mx = mean_x[static_cast<size_t>(v)];
        const double invsx = invsd_x[static_cast<size_t>(v)];

        int best = cur;
        double best_score = std::numeric_limits<double>::infinity();
        int valid_i[7];
        int nvalid = 0;
        for (int i = 0; i < nc; ++i) {
          const int k = cand[i];
          if (k >= 0 && k < K && counts[static_cast<size_t>(k)] > 0) {
            valid_i[nvalid++] = i;
          }
        }
        if (nvalid == 0) continue;

        double dotxp[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const float *xv = xsub.data() + static_cast<size_t>(v) * static_cast<size_t>(Teff);
        for (int tt = 0; tt < Teff; ++tt) {
          const double xvt = static_cast<double>(xv[tt]);
          for (int vi = 0; vi < nvalid; ++vi) {
            const int i = valid_i[vi];
            const int k = cand[i];
            dotxp[i] += xvt * proto[static_cast<size_t>(k) * static_cast<size_t>(Teff) + static_cast<size_t>(tt)];
          }
        }

        for (int vi = 0; vi < nvalid; ++vi) {
          const int i = valid_i[vi];
          const int k = cand[i];
          const double dotxp_i = dotxp[i];

          const double mp = proto_mean[static_cast<size_t>(k)];
          const double sp = proto_sum[static_cast<size_t>(k)];
          const double cov = dotxp_i - mx * sp - mp * sx + static_cast<double>(Teff) * mx * mp;
          double corr = cov * invsx * proto_invsd[static_cast<size_t>(k)];
          if (corr > 1.0) corr = 1.0;
          if (corr < -1.0) corr = -1.0;
          const double l2 = (sum_x2[static_cast<size_t>(v)] +
                             proto_energy[static_cast<size_t>(k)] -
                             2.0 * dotxp_i) / static_cast<double>(Teff);

          const double dx = static_cast<double>(x) - cx[static_cast<size_t>(k)];
          const double dy = static_cast<double>(y) - cy[static_cast<size_t>(k)];
          const double dz = static_cast<double>(z) - cz[static_cast<size_t>(k)];
          const double dist2 = dx * dx + dy * dy + dz * dz;
          const double score = (1.0 - corr) + refine_l2_weight * l2 + refine_alpha_over_s2 * dist2;

          if (score < best_score) {
            best_score = score;
            best = k;
          }
        }

        if (best != cur) {
          new_labels[static_cast<size_t>(v)] = best;
          changes += 1;
        }
      }

      labels.swap(new_labels);
      if (verbose) {
        const double cr = static_cast<double>(changes) / static_cast<double>(N);
        Rcout << "brs_slic_core refine pass " << (bp + 1)
              << ": changes=" << changes << " (" << cr << ")\n";
      }
      if (changes == 0) break;
    }

    const double invS = 1.0 / S;
    const int nbx = static_cast<int>(std::floor((G.X - 1) * invS)) + 1;
    const int nby = static_cast<int>(std::floor((G.Y - 1) * invS)) + 1;
    const int nbz = static_cast<int>(std::floor((G.Z - 1) * invS)) + 1;
    std::vector<int> head(static_cast<size_t>(nbx) * static_cast<size_t>(nby) * static_cast<size_t>(nbz), -1);
    std::vector<int> next(static_cast<size_t>(K), -1);

    for (int gp = 0; gp < global_passes; ++gp) {
      recompute_cluster_stats();
      std::fill(head.begin(), head.end(), -1);
      std::fill(next.begin(), next.end(), -1);

      for (int k = 0; k < K; ++k) {
        if (counts[static_cast<size_t>(k)] <= 0) continue;
        int bx = static_cast<int>(std::floor(cx[static_cast<size_t>(k)] * invS));
        int by = static_cast<int>(std::floor(cy[static_cast<size_t>(k)] * invS));
        int bz = static_cast<int>(std::floor(cz[static_cast<size_t>(k)] * invS));
        if (bx < 0) bx = 0;
        if (by < 0) by = 0;
        if (bz < 0) bz = 0;
        if (bx >= nbx) bx = nbx - 1;
        if (by >= nby) by = nby - 1;
        if (bz >= nbz) bz = nbz - 1;
        const int bid = bin_id(bx, by, bz, nbx, nby);
        next[static_cast<size_t>(k)] = head[static_cast<size_t>(bid)];
        head[static_cast<size_t>(bid)] = k;
      }

      int64_t changes = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : changes) schedule(static) num_threads(n_threads)
#endif
      for (int v = 0; v < N; ++v) {
        const int cur = labels[static_cast<size_t>(v)];
        new_labels[static_cast<size_t>(v)] = cur;
        const int x = G.vx[static_cast<size_t>(v)];
        const int y = G.vy[static_cast<size_t>(v)];
        const int z = G.vz[static_cast<size_t>(v)];

        int cand[96];
        int nc = 0;
        add_candidate(cur, cand, nc, 96);
        const int bx = static_cast<int>(std::floor(static_cast<double>(x) * invS));
        const int by = static_cast<int>(std::floor(static_cast<double>(y) * invS));
        const int bz = static_cast<int>(std::floor(static_cast<double>(z) * invS));

        for (int dz = -1; dz <= 1; ++dz) {
          const int zz = bz + dz;
          if (zz < 0 || zz >= nbz) continue;
          for (int dy = -1; dy <= 1; ++dy) {
            const int yy = by + dy;
            if (yy < 0 || yy >= nby) continue;
            for (int dx = -1; dx <= 1; ++dx) {
              const int xx = bx + dx;
              if (xx < 0 || xx >= nbx) continue;
              const int bid = bin_id(xx, yy, zz, nbx, nby);
              for (int k = head[static_cast<size_t>(bid)]; k != -1; k = next[static_cast<size_t>(k)]) {
                add_candidate(k, cand, nc, 96);
              }
            }
          }
        }
        if (nc <= 1) continue;

        const double sx = sum_x[static_cast<size_t>(v)];
        const double mx = mean_x[static_cast<size_t>(v)];
        const double invsx = invsd_x[static_cast<size_t>(v)];
        int valid_k[96];
        int nvalid = 0;
        for (int i = 0; i < nc; ++i) {
          const int k = cand[i];
          if (k >= 0 && k < K && counts[static_cast<size_t>(k)] > 0) {
            valid_k[nvalid++] = k;
          }
        }
        if (nvalid <= 1) continue;

        double dotxp[96];
        for (int i = 0; i < nvalid; ++i) dotxp[i] = 0.0;
        const float *xv = xsub.data() + static_cast<size_t>(v) * static_cast<size_t>(Teff);
        for (int tt = 0; tt < Teff; ++tt) {
          const double xvt = static_cast<double>(xv[tt]);
          for (int i = 0; i < nvalid; ++i) {
            const int k = valid_k[i];
            dotxp[i] += xvt * proto[static_cast<size_t>(k) * static_cast<size_t>(Teff) + static_cast<size_t>(tt)];
          }
        }

        int best = cur;
        double best_score = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nvalid; ++i) {
          const int k = valid_k[i];
          const double mp = proto_mean[static_cast<size_t>(k)];
          const double sp = proto_sum[static_cast<size_t>(k)];
          const double cov = dotxp[i] - mx * sp - mp * sx + static_cast<double>(Teff) * mx * mp;
          double corr = cov * invsx * proto_invsd[static_cast<size_t>(k)];
          if (corr > 1.0) corr = 1.0;
          if (corr < -1.0) corr = -1.0;
          const double l2 = (sum_x2[static_cast<size_t>(v)] +
                             proto_energy[static_cast<size_t>(k)] -
                             2.0 * dotxp[i]) / static_cast<double>(Teff);

          const double dx = static_cast<double>(x) - cx[static_cast<size_t>(k)];
          const double dy = static_cast<double>(y) - cy[static_cast<size_t>(k)];
          const double dz = static_cast<double>(z) - cz[static_cast<size_t>(k)];
          const double dist2 = dx * dx + dy * dy + dz * dz;
          const double score = (1.0 - corr) + refine_l2_weight * l2 + refine_alpha_over_s2 * dist2;
          if (score < best_score) {
            best_score = score;
            best = k;
          }
        }

        if (best != cur) {
          new_labels[static_cast<size_t>(v)] = best;
          changes += 1;
        }
      }

      labels.swap(new_labels);
      if (verbose) {
        const double cr = static_cast<double>(changes) / static_cast<double>(N);
        Rcout << "brs_slic_core global pass " << (gp + 1)
              << ": changes=" << changes << " (" << cr << ")\n";
      }
      if (changes == 0) break;
    }
  }

  int min_size_eff = min_size;
  if (min_size_eff <= 0) {
    const int approx = static_cast<int>(std::round(0.25 * S * S * S));
    min_size_eff = std::max(10, approx);
  }

  enforce_connectivity(G, labels, K, min_size_eff, connectivity);

  IntegerVector out_labels(N);
  for (int i = 0; i < N; ++i) out_labels[i] = labels[static_cast<size_t>(i)] + 1;

  std::vector<double> sx(static_cast<size_t>(K), 0.0);
  std::vector<double> sy(static_cast<size_t>(K), 0.0);
  std::vector<double> sz(static_cast<size_t>(K), 0.0);
  std::vector<int32_t> cnt(static_cast<size_t>(K), 0);
  for (int v = 0; v < N; ++v) {
    const int k = labels[static_cast<size_t>(v)];
    if (k < 0 || k >= K) continue;
    cnt[static_cast<size_t>(k)] += 1;
    sx[static_cast<size_t>(k)] += static_cast<double>(G.vx[static_cast<size_t>(v)]);
    sy[static_cast<size_t>(k)] += static_cast<double>(G.vy[static_cast<size_t>(v)]);
    sz[static_cast<size_t>(k)] += static_cast<double>(G.vz[static_cast<size_t>(v)]);
  }

  NumericMatrix centers_xyz(K, 3);
  for (int k = 0; k < K; ++k) {
    const int32_t c = cnt[static_cast<size_t>(k)];
    if (c > 0) {
      const double inv = 1.0 / static_cast<double>(c);
      centers_xyz(k, 0) = sx[static_cast<size_t>(k)] * inv + 1.0;
      centers_xyz(k, 1) = sy[static_cast<size_t>(k)] * inv + 1.0;
      centers_xyz(k, 2) = sz[static_cast<size_t>(k)] * inv + 1.0;
    } else {
      centers_xyz(k, 0) = coarse_centers_xyz(k, 0);
      centers_xyz(k, 1) = coarse_centers_xyz(k, 1);
      centers_xyz(k, 2) = coarse_centers_xyz(k, 2);
    }
  }

  return List::create(
    _["labels"] = out_labels,
    _["centers"] = coarse_centers,
    _["centers_xyz"] = centers_xyz,
    _["params"] = List::create(
      _["K"] = K,
      _["d"] = d,
      _["sketch_repeats"] = sketch_repeats,
      _["alpha"] = alpha,
      _["coarse_iter"] = coarse_iter,
      _["boundary_passes"] = boundary_passes,
      _["global_passes"] = global_passes,
      _["refine_spatial"] = refine_spatial,
      _["refine_l2"] = refine_l2,
      _["refine_stride"] = refine_stride_eff,
      _["seed"] = seed,
      _["connectivity"] = connectivity,
      _["min_size"] = min_size_eff
    )
  );
}
