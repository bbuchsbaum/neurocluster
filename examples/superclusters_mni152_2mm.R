# Example: Generate 1,000 superclusters on MNI152NLin2009cAsym (2mm) for multiple methods
#
# This script downloads the MNI152 2mm template + brainmask via neuroatlas,
# runs several clustering methods at K=1000, and writes axial slice PNGs for
# quick visual comparison.
#
# Notes:
# - The template is 3D; cluster4d() will wrap it to a single-timepoint NeuroVec.
# - Output PNGs are saved under "examples/output_superclusters".
# - Running all methods is computationally heavy; feel free to subset `methods`.

library(neurocluster)
library(neuroim2)
library(neuroatlas)

set.seed(123)  # reproducible seeds for methods that use RNG

# ---------------------------------------------------------------------------
# Fetch template and mask (2mm resolution)
# ---------------------------------------------------------------------------
t1_path <- neuroatlas::get_template(
  space = "MNI152NLin2009cAsym",
  modality = "T1w",
  desc = "brain",
  resolution = 2,
  path_only = TRUE
)

mask_path <- neuroatlas::get_template(
  space = "MNI152NLin2009cAsym",
  desc = "brain",
  suffix = "mask",
  resolution = 2,
  path_only = TRUE
)

cat("Template:", t1_path, "\n")
cat("Mask    :", mask_path, "\n")

brain_vol <- neuroim2::read_vol(t1_path)
mask_vol  <- neuroim2::read_vol(mask_path)

# ---------------------------------------------------------------------------
# Helper: render axial slices for a ClusteredNeuroVol
# ---------------------------------------------------------------------------
save_axial_slices <- function(clusvol, outfile, n_slices = 10) {
  arr <- as.array(clusvol)
  dims <- dim(arr)
  z_idx <- round(seq(1, dims[3], length.out = n_slices))

  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  png(outfile, width = 1600, height = 800)
  par(mfrow = c(2, ceiling(n_slices / 2)), mar = c(1, 1, 2, 1))

  # Use a moderately sized palette to keep colors distinguishable
  pal <- grDevices::hcl.colors(256, "Zissou 1", alpha = 1)

  for (i in seq_along(z_idx)) {
    z <- z_idx[i]
    sl <- arr[, , z]
    # Simple axial display; flip y for radiological-style view
    image(
      t(apply(sl, 2, rev)),
      axes = FALSE,
      col = pal,
      main = paste0("z = ", z)
    )
  }

  dev.off()
}

# ---------------------------------------------------------------------------
# Method configuration
# ---------------------------------------------------------------------------
methods <- c("supervoxels", "snic", "slic", "slice_msf", "flash3d", "g3s", "rena")

method_params <- list(
  supervoxels = list(
    max_iterations = 20,       # lighter for full-brain 2mm
    connectivity   = 18,
    use_gradient   = FALSE,    # speeds up seeding on 3D structural data
    parallel       = TRUE,
    grain_size     = 200
  ),
  snic        = list(spatial_weight = 0.5),
  slic        = list(spatial_weight = 0.5, max_iterations = 6),
  slice_msf   = list(spatial_weight = 0.6, max_iterations = 4),
  flash3d     = list(spatial_weight = 0.6, max_iterations = 2),
  g3s         = list(spatial_weight = 0.5, max_iterations = 2),
  rena        = list(spatial_weight = 0.5, max_iterations = 20)
)

# ---------------------------------------------------------------------------
# Run clustering and save slices
# ---------------------------------------------------------------------------
output_dir <- file.path("examples", "output_superclusters")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (m in methods) {
  cat("\nRunning method:", m, "\n")
  args <- method_params[[m]]

  res <- do.call(
    cluster4d,
    c(
      list(
        vec = brain_vol,
        mask = mask_vol,
        n_clusters = 1000,
        method = m,
        verbose = FALSE
      ),
      args
    )
  )

  # Save slices
  outfile <- file.path(output_dir, paste0("mni152_2mm_", m, "_k1000.png"))
  save_axial_slices(res$clusvol, outfile)
  cat("  Saved:", outfile, "\n")
}

cat("\nDone. Check PNGs in", output_dir, "\n")
