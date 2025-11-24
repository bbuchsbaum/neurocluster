#include <Rcpp.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// =============================================================================
// OPTIMIZED IMPLEMENTATION: Lightweight structs for performance
// =============================================================================

// Lightweight struct for priority queue - avoids Rcpp::List overhead
struct QueueElement {
    int x, y, z;          // 3D coordinates
    int voxel_idx;        // Linear index in feature matrix
    int k_label;          // Cluster label
    double distance;      // Distance metric for priority

    // Priority queue orders by largest, so reverse for smallest distance
    bool operator>(const QueueElement& other) const {
        return distance > other.distance;
    }
};

// Lightweight struct for centroids with in-place updates
struct Centroid {
    std::vector<double> sum_c;   // Sum of feature vectors
    double sum_x, sum_y, sum_z;  // Sum of coordinates
    int count;

    // Cached averages (updated after each addition)
    std::vector<double> avg_c;
    double avg_x, avg_y, avg_z;

    Centroid(int n_features, double x, double y, double z, const double* c_init)
        : sum_x(x), sum_y(y), sum_z(z), count(1),
          avg_x(x), avg_y(y), avg_z(z) {
        sum_c.resize(n_features);
        avg_c.resize(n_features);
        for(int i = 0; i < n_features; ++i) {
            sum_c[i] = c_init[i];
            avg_c[i] = c_init[i];
        }
    }

    // In-place update - no allocations
    void add_pixel(double x, double y, double z, const double* features, int n_features) {
        count++;
        sum_x += x;
        sum_y += y;
        sum_z += z;

        avg_x = sum_x / count;
        avg_y = sum_y / count;
        avg_z = sum_z / count;

        // Update features and normalize (for fMRI correlation-based distance)
        double sq_norm = 0.0;
        for(int i = 0; i < n_features; ++i) {
            sum_c[i] += features[i];
            avg_c[i] = sum_c[i] / count;
            sq_norm += avg_c[i] * avg_c[i];
        }

        // Re-normalize to unit length (appropriate for fMRI time series)
        if (sq_norm > 0) {
            double norm = std::sqrt(sq_norm);
            for(int i = 0; i < n_features; ++i) {
                avg_c[i] /= norm;
            }
        }
    }
};

// Inline distance computation - no allocations, pure C++
// Combined squared distance (spatial + feature); avoids sqrt for PQ ordering
inline double compute_dist_cpp(double ck_x, double ck_y, double ck_z,
                               const double* ck_c,
                               double ni_x, double ni_y, double ni_z,
                               const double* ni_c,
                               int n_features, double s_inv, double compactness_inv) {

    double dx = ck_x - ni_x;
    double dy = ck_y - ni_y;
    double dz = ck_z - ni_z;
    double ds = dx*dx + dy*dy + dz*dz;

    double dc = 0.0;
    for(int i = 0; i < n_features; ++i) {
        double diff = ck_c[i] - ni_c[i];
        dc += diff * diff;
    }

    // Return squared combined distance; monotonic w.r.t. Euclidean distance
    return (ds * s_inv) + (dc * compactness_inv);
}

// =============================================================================
// OLD IMPLEMENTATION: Kept for backward compatibility
// =============================================================================

class IntegerArray3D {
public:
  IntegerArray3D(IntegerVector data) : data_(data) {
    dims_ = data_.attr("dim");
  }

  IntegerVector toRVector() const {
    IntegerVector result(data_);
    result.attr("dim") = dims_;
    return result;
  }

  int operator()(int i, int j, int k) const {
    int index = i + dims_[0] * (j + dims_[1] * k);
    return data_[index];
  }

  void set(int i, int j, int k, int value) {
    int index = i + dims_[0] * (j + dims_[1] * k);
    data_[index] = value;
  }

  int dim(int i) const {
    return dims_[i];
  }

public:
  IntegerVector data_;
  IntegerVector dims_;
};

IntegerMatrix get_26_connected_neighbors(int i, int j, int k, int max_i, int max_j, int max_k) {
  std::vector<IntegerVector> neighbors;

  for (int i_shift = -1; i_shift <= 1; ++i_shift) {
    for (int j_shift = -1; j_shift <= 1; ++j_shift) {
      for (int k_shift = -1; k_shift <= 1; ++k_shift) {
        // Skip the center point (0, 0, 0) shift
        if (i_shift == 0 && j_shift == 0 && k_shift == 0) {
          continue;
        }

        int i_neighbor = i + i_shift;
        int j_neighbor = j + j_shift;
        int k_neighbor = k + k_shift;

        // Check if the neighbor is within the array bounds
        if (i_neighbor >= 0 && i_neighbor < max_i &&
            j_neighbor >= 0 && j_neighbor < max_j &&
            k_neighbor >= 0 && k_neighbor < max_k) {
          IntegerVector neighbor_coords = IntegerVector::create(i_neighbor, j_neighbor, k_neighbor);
          neighbors.push_back(neighbor_coords);
        }
      }
    }
  }

  IntegerMatrix neighbors_matrix(neighbors.size(), 3);
  for (size_t r = 0; r < neighbors.size(); ++r) {
    neighbors_matrix(r, _) = neighbors[r];
  }

  return neighbors_matrix;
}


// [[Rcpp::export]]
List update_centroid_online(const List& centroid, const NumericVector& x_i, const NumericVector& c_i) {
  // centroid: A list containing the current centroid information with elements x (spatial position), c (color), and n (number of pixels)
  // x_i: Spatial position of the current pixel
  // c_i: Color value of the current pixel

  // Extract the current centroid information
  NumericVector current_x = centroid["x"];
  NumericVector current_c = centroid["c"];
  int current_n = centroid["n"];

  // Calculate the new spatial position
  NumericVector new_x = (current_x * current_n + x_i) / (current_n + 1);

  // Calculate the new color value
  NumericVector new_c = (current_c * current_n + c_i) / (current_n + 1);

  // Normalize the new color value
  new_c = new_c / sqrt(sum(pow(new_c, 2)));

  // Update the number of pixels in the superpixel
  int new_n = current_n + 1;

  // Return the updated centroid as a list
  return List::create(_["x"] = new_x, _["c"] = new_c, _["label"] = centroid["label"], _["n"] = new_n);
}


double compute_distance(NumericVector x_i, NumericVector x_k,
                             NumericVector c_i, NumericVector c_k,
                             double s = 5, double compactness = 5) {
  double dc = sum(pow(c_i - c_k, 2)); // squared Euclidean color distance
  double ds = sum(pow(x_i - x_k, 2)); // squared Euclidean spatial distance

  // Incorporate compactness parameter into the distance calculation
  //double distance = sqrt(pow(dc, 2) + pow(ds / s * compactness, 2));
  double distance = sqrt( (ds/s) + (dc/compactness));
  return distance;
}

List create_element(const NumericVector& x_k, const IntegerVector& voxel,
                    int centroid_idx_k, int k, double distance) {

  List e = List::create(Named("x") = x_k,
                        Named("voxel") = voxel,
                        Named("i_k") = centroid_idx_k,
                        Named("label") = k,
                        Named("distance") = distance);
  return e;
}




// [[Rcpp::export]]
IntegerVector snic_main(IntegerVector L_data, const IntegerVector& mask,
                        const NumericMatrix centroids,
                        const IntegerVector centroid_idx,
                        const IntegerMatrix valid_coords,
                        const NumericMatrix norm_coords,
                        const NumericMatrix& vecmat,
                        int K, double s, double compactness,
                        IntegerVector mask_lookup_data) {

  IntegerArray3D L(L_data);
  IntegerArray3D mask_lookup(mask_lookup_data);

  List C(K);
  int nassigned = 0;
  int nvoxels = vecmat.ncol();

  IntegerVector mask_dims = mask.attr("dim");
  //IntegerVector mask_dims = IntegerVector::create(96,96,26);

  struct CompareDist {
    bool operator()(const List& a, const List& b) {
      return as<double>(a["distance"]) > as<double>(b["distance"]);
    }
  };

  // Create a priority_queue with custom comparison function
  std::priority_queue<List, std::vector<List>, CompareDist> Q;

  // Initialize the priority queue with initial elements
  for (int k = 0; k < K; ++k) {
    // ... Initialize the initial elements (e) ...
    NumericVector x_k = centroids(k, _);
    NumericVector c_k = vecmat(_, centroid_idx[k]);
    IntegerVector voxel = valid_coords(centroid_idx[k], _);


    double distance = 0.0;

    List e = create_element(x_k, voxel, centroid_idx[k], k+1, (k*1.0)/(K+1));

    // Update the centroid information
    C[k] = List::create(_["x"] = x_k, _["c"] = c_k, _["label"]=k+1, _["n"] = 1);

    // Push the element into the priority queue
    Q.push(e);


  }

  List e_i;
  // Main loop
  while (!Q.empty() && nassigned <= nvoxels) {
    // Pop the element with the smallest distance
    e_i = Q.top();
    Q.pop();
    NumericVector x_i = e_i["x"];
    IntegerVector voxel = e_i["voxel"];

    int i_k = e_i["i_k"];
    int k_i = e_i["label"];

    if (L(voxel(0), voxel(1), voxel(2)) == 0) {
      L.set(voxel(0), voxel(1), voxel(2), k_i);
      nassigned++;
      C[k_i-1] = update_centroid_online(C[k_i-1], x_i, vecmat(_, i_k));

      List current_centroid = C.at(k_i - 1);
      //Rcout << "new centroid" << std::endl;
      //Rcpp::print(current_centroid);
      NumericVector cen_k = as<NumericVector>(current_centroid["x"]);
      NumericVector ci = as<NumericVector>(current_centroid["c"]);

      // Iterate through the neighbors of the current pixel
      // Use your implementation of the get_26_connected_neighbors_rcpp function
      IntegerMatrix neighbors = get_26_connected_neighbors(voxel(0), voxel(1), voxel(2),
                                                           mask_dims(0), mask_dims(1), mask_dims(2));


      for (int i = 0; i < neighbors.nrow(); ++i) {
        IntegerVector x_j = neighbors(i, _);
        int neighb_idx = mask_lookup(x_j(0), x_j(1), x_j(2));
        if (neighb_idx < 0) {
          continue;
        }
        if (L(x_j(0), x_j(1), x_j(2)) != 0) {
          continue;
        }


        NumericVector c_j = vecmat(_, neighb_idx);
        NumericVector nx_j = norm_coords(neighb_idx, _);

        // Compute the distance
        // Use your implementation of the compute_distance_rcpp function
        double d_j_ki = compute_distance(cen_k, nx_j, ci, c_j, s, compactness);
        List e_j = create_element(nx_j, x_j, neighb_idx, k_i,d_j_ki);
        Q.push(e_j);

      }
    }

  }

  return L.toRVector();

}

// =============================================================================
// OPTIMIZED SNIC IMPLEMENTATION
// =============================================================================

// [[Rcpp::export]]
IntegerVector snic_main_optimized(IntegerVector L_data,
                                  const IntegerVector& mask,
                                  const NumericMatrix centroids,
                                  const IntegerVector centroid_idx,
                                  const IntegerMatrix valid_coords,
                                  const NumericMatrix norm_coords,
                                  const NumericMatrix& vecmat,
                                  int K, double s, double compactness,
                                  IntegerVector mask_lookup_data) {

    // Setup dimensions - store as plain ints for fast access
    IntegerVector mask_dims = mask.attr("dim");
    int dim_x = mask_dims[0];
    int dim_y = mask_dims[1];
    int dim_z = mask_dims[2];

    // Direct pointer access for fastest array operations
    int* L_ptr = L_data.begin();
    const int* mask_lookup_ptr = mask_lookup_data.begin();

    int nvoxels = vecmat.ncol();
    int n_features = vecmat.nrow();
    int coord_rows = valid_coords.nrow();
    int norm_rows = norm_coords.nrow();
    int mask_lookup_size = mask_lookup_data.size();
    const double* norm_ptr = norm_coords.begin();
    int norm_stride = norm_coords.nrow();

    if (coord_rows != nvoxels || norm_rows != nvoxels) {
        Rcpp::stop("SNIC internal error: coordinate matrices do not match voxel count");
    }
    if (mask_lookup_size != dim_x * dim_y * dim_z) {
        Rcpp::stop("SNIC internal error: mask lookup size mismatch");
    }

    // Precompute inverses to replace division inside hot loops
    const double s_inv = (s > 0.0) ? 1.0 / s : 1.0;
    const double comp_inv = (compactness > 0.0) ? 1.0 / compactness : 1.0;

    // Priority queue with lightweight struct
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<QueueElement>> Q;

    // Initialize centroids
    std::vector<Centroid> C_store;
    C_store.reserve(K);

    // Initialization: create initial centroids and queue elements
    for (int k = 0; k < K; ++k) {
        int idx = centroid_idx[k];  // 0-based index from R
        if (idx < 0 || idx >= nvoxels) {
            Rcpp::stop("SNIC internal error: centroid index out of bounds");
        }

        // Get spatial coordinates
        int vx = valid_coords(idx, 0);
        int vy = valid_coords(idx, 1);
        int vz = valid_coords(idx, 2);

        // Get normalized coordinates
        double nx = norm_ptr[idx + norm_stride * 0];
        double ny = norm_ptr[idx + norm_stride * 1];
        double nz = norm_ptr[idx + norm_stride * 2];

        // Get feature vector - access column in matrix
        // NumericMatrix is column-major, so column idx starts at &vecmat[0] + idx * n_features
        const double* c_init = &vecmat[0] + idx * n_features;

        // Create centroid
        C_store.emplace_back(n_features, nx, ny, nz, c_init);

        // Add to queue with tie-breaking distance
        double dist = (double)k / (K + 1);
        Q.push({vx, vy, vz, idx, k + 1, dist});
    }

    int nassigned = 0;

    // Main loop - process voxels in priority order
    while (!Q.empty()) {
        // Check for user interrupts periodically
        if (nassigned % 1000 == 0) {
            Rcpp::checkUserInterrupt();
        }

        QueueElement e = Q.top();
        Q.pop();

        // Calculate linear index for 3D volume
        int l_index = e.x + dim_x * (e.y + dim_y * e.z);

        // Check if already labeled
        if (L_ptr[l_index] == 0) {
            // Assign label
            L_ptr[l_index] = e.k_label;
            nassigned++;

            // Update centroid in-place
            const double* feat_ptr = &vecmat[0] + e.voxel_idx * n_features;
            double nx = norm_ptr[e.voxel_idx + norm_stride * 0];
            double ny = norm_ptr[e.voxel_idx + norm_stride * 1];
            double nz = norm_ptr[e.voxel_idx + norm_stride * 2];

            Centroid& ck = C_store[e.k_label - 1];
            ck.add_pixel(nx, ny, nz, feat_ptr, n_features);

            // Check 26-connected neighbors - inlined for performance
            for (int dz = -1; dz <= 1; ++dz) {
                int nz_pos = e.z + dz;
                if (nz_pos < 0 || nz_pos >= dim_z) continue;

                for (int dy = -1; dy <= 1; ++dy) {
                    int ny_pos = e.y + dy;
                    if (ny_pos < 0 || ny_pos >= dim_y) continue;

                    for (int dx = -1; dx <= 1; ++dx) {
                        // Skip center voxel
                        if (dx == 0 && dy == 0 && dz == 0) continue;

                        int nx_pos = e.x + dx;
                        if (nx_pos < 0 || nx_pos >= dim_x) continue;

                        // Linear index for neighbor
                        int n_l_index = nx_pos + dim_x * (ny_pos + dim_y * nz_pos);

                        // Check if already labeled
                        if (L_ptr[n_l_index] != 0) continue;

                        // Get neighbor voxel index from lookup
                        int neigh_idx = mask_lookup_ptr[n_l_index];

                        if (neigh_idx < 0) continue;
                        if (neigh_idx >= nvoxels) {
                            Rcpp::stop("SNIC internal error: neighbor index out of bounds");
                        }

                        // Get neighbor features and coordinates
                        const double* n_feat_ptr = &vecmat[0] + neigh_idx * n_features;
                        double nnx = norm_ptr[neigh_idx + norm_stride * 0];
                        double nny = norm_ptr[neigh_idx + norm_stride * 1];
                        double nnz = norm_ptr[neigh_idx + norm_stride * 2];

                        // Compute distance
                        double d = compute_dist_cpp(ck.avg_x, ck.avg_y, ck.avg_z, ck.avg_c.data(),
                                                    nnx, nny, nnz, n_feat_ptr,
                                                    n_features, s_inv, comp_inv);

                        // Add to queue
                        Q.push({nx_pos, ny_pos, nz_pos, neigh_idx, e.k_label, d});
                    }
                }
            }
        }
    }

    return L_data;
}

// [[Rcpp::export]]
IntegerVector compute_boundaryscore_3d_cpp(IntegerVector volume, IntegerVector mask) {
  IntegerArray3D volume_array(volume);
  IntegerArray3D mask_array(mask);

  IntegerVector boundaries_vector(volume.size());
  boundaries_vector.attr("dim") = volume.attr("dim");
  IntegerArray3D boundaries_array(boundaries_vector);

  for (int k = 0; k < volume_array.dim(2); ++k) {
    for (int i = 1; i < volume_array.dim(0) - 1; ++i) {
      for (int j = 1; j < volume_array.dim(1) - 1; ++j) {
        if (mask_array(i, j, k)) {
          int boundary_score = 0;
          for (int a = -1; a <= 1; ++a) {
            for (int b = -1; b <= 1; ++b) {
              if (a != 0 || b != 0) {
                if (volume_array(i, j, k) != volume_array(i + a, j + b, k)) {
                  boundary_score++;
                }
              }
            }
          }
          boundaries_array.set(i, j, k, boundary_score);
        }
      }
    }
  }
  return boundaries_array.toRVector();
}

// [[Rcpp::export]]
IntegerVector detect_boundaries_2d_cpp(IntegerVector volume, IntegerVector mask) {
  IntegerArray3D volume_array(volume);
  IntegerArray3D mask_array(mask);

  IntegerVector boundaries_vector(volume.size());
  boundaries_vector.attr("dim") = volume.attr("dim");
  IntegerArray3D boundaries_array(boundaries_vector);

  for (int k = 0; k < volume_array.dim(2); ++k) {
    for (int i = 1; i < volume_array.dim(0) - 1; ++i) {
      for (int j = 1; j < volume_array.dim(1) - 1; ++j) {
        if (mask_array(i, j, k)) {
          bool is_boundary = false;
          for (int a = -1; a <= 1; ++a) {
            for (int b = -1; b <= 1; ++b) {
              if (a != 0 || b != 0) {
                if (volume_array(i, j, k) != volume_array(i + a, j + b, k)) {
                  is_boundary = true;
                  break;
                }
              }
            }
            if (is_boundary) break;
          }
          boundaries_array.set(i, j, k, is_boundary);
        }
      }
    }
  }
  return boundaries_array.toRVector();
}
