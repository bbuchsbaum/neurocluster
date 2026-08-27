# Transform New Data Using Existing SVD Compression

Applies a previously computed SVD compression to new data, ensuring
consistency across multiple runs or when adding new voxels.

## Usage

``` r
transform_new_data_svd(new_data, compression_result)
```

## Arguments

- new_data:

  Numeric matrix (N_new x T) with same number of columns as the original
  training data.

- compression_result:

  List returned by
  [`compress_features_svd`](https://bbuchsbaum.github.io/neurocluster/reference/compress_features_svd.md)
  containing the rotation matrix and scaling parameters.

## Value

Compressed and normalized feature matrix (N_new x M) where M is the
number of components in the original compression.

## Examples

``` r
# Train on subset of data
train_data <- matrix(rnorm(500 * 300), 500, 300)
compression <- compress_features_svd(train_data, n_components = 15)
#> Requested n_components (15) explained 13.9% variance; increased to 223 components to retain 95%.

# Apply to new data
test_data <- matrix(rnorm(100 * 300), 100, 300)
compressed_test <- transform_new_data_svd(test_data, compression)
```
