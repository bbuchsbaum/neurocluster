#include <Rcpp.h>
#include <queue>
#include <vector>
#include <algorithm>
using namespace Rcpp;

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




