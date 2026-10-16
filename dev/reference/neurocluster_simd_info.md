# SIMD / CPU Feature Report (Development Helper)

Reports build-time SIMD feature macros and (when available) a small set
of runtime CPU feature flags, plus the currently-selected int8
dot-product kernel used by neurocluster's correlation-SLIC code paths.

## Usage

``` r
neurocluster_simd_info()
```

## Value

A named list with fields `compiler`, `arch`, `compile`, `runtime`,
`dot_i8_kernel`, `shared_library`, and `notes`.

## Details

This is primarily a development/benchmarking aid. For definitive
evidence that an ISA-specific kernel is in use, disassemble the package
shared library and look for the corresponding instructions.
