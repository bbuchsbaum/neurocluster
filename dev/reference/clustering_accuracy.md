# Compute Clustering Accuracy Metrics

Compares predicted cluster assignments against ground truth using
multiple metrics that are invariant to cluster label permutation.

## Usage

``` r
clustering_accuracy(predicted, truth)
```

## Arguments

- predicted:

  Integer vector of predicted cluster assignments

- truth:

  Integer vector of ground truth cluster assignments

## Value

A list with adjusted Rand index (`ari`), variation of information in
bits (`variation_of_information_bits`), pairwise Dice (`pairwise_dice`),
arithmetic normalized mutual information (`nmi`), exact maximum
one-to-one label-matched accuracy (`matched_accuracy`), its
compatibility alias `accuracy`, partition counts, and a
permutation-invariant `perfect` flag.

## Details

The metrics compare the two induced equivalence relations; cluster names
do not need to match. All observations are included, including a label
of zero. Missing or non-finite labels are rejected rather than silently
removed. Computation uses only observed contingency cells and is linear
in the number of observations for singleton-heavy partitions.

**Adjusted Rand Index (ARI):** Measures agreement between two
clusterings, adjusted for chance. ARI = 1 means perfect agreement, ARI =
0 means random, ARI \< 0 means worse than random.

**Variation of information (VI):** The sum of the two conditional
entropies, measured with base-2 logarithms. VI is zero exactly for
identical partitions and increases with disagreement.

**Pairwise Dice:** `2 * TP / (2 * TP + FP + FN)` over unordered pairs of
observations, where a positive pair is co-clustered. If both partitions
contain only singleton clusters, their pairwise Dice is defined as one.

**Matched accuracy:** The greatest fraction of matched observations
under a one-to-one mapping between predicted and truth labels. It is
solved exactly on the sparse observed contingency graph. `accuracy` is
retained as an alias for existing callers; it is not the former greedy
approximation.

## Examples

``` r
# Perfect clustering
truth <- rep(1:3, each = 10)
pred <- rep(c(2, 3, 1), each = 10)  # Same clustering, different labels
clustering_accuracy(pred, truth)  # ARI = 1
#> $ari
#> [1] 1
#> 
#> $variation_of_information_bits
#> [1] 0
#> 
#> $pairwise_dice
#> [1] 1
#> 
#> $nmi
#> [1] 1
#> 
#> $matched_accuracy
#> [1] 1
#> 
#> $accuracy
#> [1] 1
#> 
#> $n_predicted
#> [1] 3
#> 
#> $n_truth
#> [1] 3
#> 
#> $perfect
#> [1] TRUE
#> 
```
