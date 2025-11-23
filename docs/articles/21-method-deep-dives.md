# Method deep dives

## Supervoxels

- Iterative, heat-kernel similarities; parallelized updates available.

## SNIC

- Non-iterative; sequential priority queue; connected components by
  design.

## SLIC

- Iterative with local windows; preserves K; good speed/quality balance.

## Slice-MSF

- Slice-wise minimum spanning forest; supports consensus and stitching.

## FLASH-3D

- DCT-based hashing/compression; fast approximation; temporal smoothing
  via `lambda_t`. includes: in_header: \|-
