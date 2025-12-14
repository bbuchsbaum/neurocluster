#include <Rcpp.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>

using namespace Rcpp;

// =============================================================================
// UNION-FIND DATA STRUCTURE FOR CONNECTED COMPONENTS
// =============================================================================

struct UnionFind {
    std::vector<int> parent;
    std::vector<int> rank;
    int n_components;

    UnionFind(int n) : parent(n), rank(n, 0), n_components(n) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    // Find with path compression
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    // Union by rank
    bool unite(int x, int y) {
        int px = find(x);
        int py = find(y);

        if (px == py) return false;

        if (rank[px] < rank[py]) {
            parent[px] = py;
        } else if (rank[px] > rank[py]) {
            parent[py] = px;
        } else {
            parent[py] = px;
            rank[px]++;
        }
        n_components--;
        return true;
    }

    // Get component labels
    std::vector<int> get_labels() {
        std::vector<int> labels(parent.size());
        std::map<int, int> root_to_label;
        int next_label = 0;

        for (size_t i = 0; i < parent.size(); ++i) {
            int root = find(i);
            if (root_to_label.find(root) == root_to_label.end()) {
                root_to_label[root] = next_label++;
            }
            labels[i] = root_to_label[root];
        }

        return labels;
    }
};

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

//' Compute masked distances for ReNA
//'
//' Computes squared Euclidean distances only for connected pairs in sparse adjacency.
//'
//' @param feature_mat Numeric matrix (features x voxels)
//' @param adjacency_i Integer vector of row indices for sparse adjacency
//' @param adjacency_j Integer vector of col indices for sparse adjacency
//' @return NumericVector of distances (same length as adjacency_i)
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_masked_distances_cpp(NumericMatrix feature_mat,
                                           IntegerVector adjacency_i,
                                           IntegerVector adjacency_j) {
    int n_edges = adjacency_i.size();
    int n_features = feature_mat.nrow();
    int n_nodes = feature_mat.ncol();
    NumericVector distances(n_edges);

    for (int e = 0; e < n_edges; ++e) {
        int i = adjacency_i[e] - 1;  // Convert from 1-based to 0-based
        int j = adjacency_j[e] - 1;

        // Guard against malformed adjacency indices; skip invalid edges
        if (i < 0 || j < 0 || i >= n_nodes || j >= n_nodes) {
            distances[e] = R_PosInf;
            continue;
        }

        double dist_sq = 0.0;
        for (int f = 0; f < n_features; ++f) {
            double diff = feature_mat(f, i) - feature_mat(f, j);
            dist_sq += diff * diff;
        }

        distances[e] = dist_sq;
    }

    return distances;
}

//' Find 1-Nearest Neighbor subgraph for ReNA
//'
//' For each node, finds its single nearest neighbor to form directed 1-NN graph.
//'
//' @param n_nodes Number of nodes
//' @param adjacency_i Integer vector of row indices for sparse adjacency
//' @param adjacency_j Integer vector of col indices for sparse adjacency
//' @param distances Numeric vector of distances for each edge
//' @return IntegerVector of nearest neighbor indices (0-based, -1 if no neighbors)
//'
//' @keywords internal
// [[Rcpp::export]]
IntegerVector find_1nn_subgraph_cpp(int n_nodes,
                                    IntegerVector adjacency_i,
                                    IntegerVector adjacency_j,
                                    NumericVector distances) {

    IntegerVector nearest_neighbor(n_nodes, -1);
    NumericVector min_distance(n_nodes, R_PosInf);

    int n_edges = adjacency_i.size();

    // Process edges in deterministic order (by i then j)
    std::vector<int> order(n_edges);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (adjacency_i[a] == adjacency_i[b]) return adjacency_j[a] < adjacency_j[b];
        return adjacency_i[a] < adjacency_i[b];
    });

    // For each edge, update nearest neighbor if this is the closest so far
    for (int idx = 0; idx < n_edges; ++idx) {
        int e = order[idx];
        int i = adjacency_i[e] - 1;  // Convert from 1-based to 0-based
        int j = adjacency_j[e] - 1;
        double dist = distances[e];

        // Update i's nearest neighbor if j is closer
        if (dist < min_distance[i] ||
            (dist == min_distance[i] && (nearest_neighbor[i] < 0 || j < nearest_neighbor[i]))) {
            min_distance[i] = dist;
            nearest_neighbor[i] = j;
        }

        // Update j's nearest neighbor if i is closer (undirected graph)
        if (dist < min_distance[j] ||
            (dist == min_distance[j] && (nearest_neighbor[j] < 0 || i < nearest_neighbor[j]))) {
            min_distance[j] = dist;
            nearest_neighbor[j] = i;
        }
    }

    return nearest_neighbor;
}

//' Find connected components using Union-Find
//'
//' Finds weakly connected components in a directed graph defined by 1-NN edges.
//'
//' @param n_nodes Number of nodes
//' @param nearest_neighbor IntegerVector of nearest neighbor for each node (0-based)
//' @return IntegerVector of component labels (0-based, contiguous)
//'
//' @keywords internal
// [[Rcpp::export]]
IntegerVector find_connected_components_cpp(int n_nodes,
                                            IntegerVector nearest_neighbor) {

    UnionFind uf(n_nodes);

    // Union each node with its nearest neighbor
    for (int i = 0; i < n_nodes; ++i) {
        int nn = nearest_neighbor[i];
        if (nn >= 0 && nn < n_nodes) {  // Valid neighbor
            uf.unite(i, nn);
        }
    }

    // Get component labels
    std::vector<int> labels = uf.get_labels();

    // Convert to R vector
    IntegerVector result(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        result[i] = labels[i];
    }

    return result;
}

//' Aggregate features by component (mean pooling)
//'
//' Computes mean feature vector for each component efficiently.
//'
//' @param feature_mat Numeric matrix (features x voxels)
//' @param component_labels IntegerVector of component labels (0-based)
//' @param n_components Number of unique components
//' @return NumericMatrix of aggregated features (features x n_components)
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix aggregate_features_cpp(NumericMatrix feature_mat,
                                     IntegerVector component_labels,
                                     int n_components) {

    int n_features = feature_mat.nrow();
    int n_voxels = feature_mat.ncol();

    // Initialize result and counts
    NumericMatrix result(n_features, n_components);
    std::vector<int> counts(n_components, 0);

    // Accumulate sums
    for (int v = 0; v < n_voxels; ++v) {
        int comp = component_labels[v];
        if (comp >= 0 && comp < n_components) {
            counts[comp]++;
            for (int f = 0; f < n_features; ++f) {
                result(f, comp) += feature_mat(f, v);
            }
        }
    }

    // Divide by counts to get means
    for (int c = 0; c < n_components; ++c) {
        if (counts[c] > 0) {
            for (int f = 0; f < n_features; ++f) {
                result(f, c) /= counts[c];
            }
        }
    }

    return result;
}

//' Aggregate coordinates by component (mean pooling)
//'
//' Computes mean coordinate for each component efficiently.
//'
//' @param coords Numeric matrix (voxels x 3)
//' @param component_labels IntegerVector of component labels (0-based)
//' @param n_components Number of unique components
//' @return NumericMatrix of aggregated coordinates (n_components x 3)
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix aggregate_coords_cpp(NumericMatrix coords,
                                   IntegerVector component_labels,
                                   int n_components) {

    int n_voxels = coords.nrow();

    // Initialize result and counts
    NumericMatrix result(n_components, 3);
    std::vector<int> counts(n_components, 0);

    // Accumulate sums
    for (int v = 0; v < n_voxels; ++v) {
        int comp = component_labels[v];
        if (comp >= 0 && comp < n_components) {
            counts[comp]++;
            for (int d = 0; d < 3; ++d) {
                result(comp, d) += coords(v, d);
            }
        }
    }

    // Divide by counts to get means
    for (int c = 0; c < n_components; ++c) {
        if (counts[c] > 0) {
            for (int d = 0; d < 3; ++d) {
                result(c, d) /= counts[c];
            }
        }
    }

    return result;
}

//' Contract adjacency graph by merging components
//'
//' Builds new adjacency matrix where nodes are components and edges exist if
//' any constituent nodes were connected.
//'
//' @param adjacency_i Integer vector of row indices for sparse adjacency
//' @param adjacency_j Integer vector of col indices for sparse adjacency
//' @param component_labels IntegerVector of component labels (0-based)
//' @param n_components Number of unique components
//' @return List with contracted_i and contracted_j vectors
//'
//' @keywords internal
// [[Rcpp::export]]
List contract_graph_cpp(IntegerVector adjacency_i,
                        IntegerVector adjacency_j,
                        IntegerVector component_labels,
                        int n_components) {

    int n_edges = adjacency_i.size();

    // Use set to track unique edges between components
    std::set<std::pair<int, int>> edge_set;

    for (int e = 0; e < n_edges; ++e) {
        int i = adjacency_i[e] - 1;  // Convert from 1-based to 0-based
        int j = adjacency_j[e] - 1;

        int comp_i = component_labels[i];
        int comp_j = component_labels[j];

        // Skip self-edges
        if (comp_i == comp_j) continue;

        // Add edge (ensure i < j for undirected graph)
        if (comp_i < comp_j) {
            edge_set.insert(std::make_pair(comp_i, comp_j));
        } else {
            edge_set.insert(std::make_pair(comp_j, comp_i));
        }
    }

    // Convert set to vectors
    int n_new_edges = edge_set.size();
    IntegerVector new_i(n_new_edges);
    IntegerVector new_j(n_new_edges);

    int idx = 0;
    for (const auto& edge : edge_set) {
        new_i[idx] = edge.first + 1;   // Convert back to 1-based
        new_j[idx] = edge.second + 1;
        idx++;
    }

    return List::create(
        Named("i") = new_i,
        Named("j") = new_j
    );
}

//' Prune 1-NN edges to achieve target number of components
//'
//' Sorts 1-NN edges by distance and adds them incrementally until K components remain.
//' Used for exact-K stopping condition.
//'
//' @param n_nodes Number of nodes
//' @param nearest_neighbor IntegerVector of nearest neighbor for each node (0-based)
//' @param distances NumericVector of distances to nearest neighbors
//' @param target_k Target number of components
//' @return IntegerVector of pruned nearest neighbors (-1 for pruned edges)
//'
//' @keywords internal
// [[Rcpp::export]]
IntegerVector prune_edges_for_k_cpp(int n_nodes,
                                     IntegerVector nearest_neighbor,
                                     NumericVector distances,
                                     int target_k) {

    // Create sorted list of edges by distance
    struct Edge {
        int from;
        int to;
        double dist;
        bool operator<(const Edge& other) const {
            return dist < other.dist;
        }
    };

    std::vector<Edge> edges;
    for (int i = 0; i < n_nodes; ++i) {
        int nn = nearest_neighbor[i];
        if (nn >= 0 && nn < n_nodes) {
            edges.push_back({i, nn, distances[i]});
        }
    }

    // Sort by distance
    std::sort(edges.begin(), edges.end());

    // Add edges incrementally until we reach target_k components
    UnionFind uf(n_nodes);
    IntegerVector pruned_neighbors(n_nodes, -1);

    for (const auto& edge : edges) {
        if (uf.unite(edge.from, edge.to)) {
            pruned_neighbors[edge.from] = edge.to;

            // Stop if we've reached target_k components
            if (uf.n_components <= target_k) {
                break;
            }
        }
    }

    return pruned_neighbors;
}
