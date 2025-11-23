# Performance and memory tradeoffs

## Scaling

- Complexity grows with voxels × clusters × iterations; prefer
  slice-wise or hashing methods as N grows.

## Memory

- Rough estimate: `n_vox * n_time * 8 bytes` for raw doubles.

## Rules of thumb

- Small (\<10k vox): supervoxels or SNIC
- Medium (10k–100k): SLIC or slice-MSF
- Large (\>100k): slice-MSF or FLASH-3D includes: in_header: \|-
