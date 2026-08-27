#!/usr/bin/env Rscript

# Rebuild the README's synthetic workflow figure from the current source tree.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

syn <- generate_synthetic_volume(
  "gaussian_blobs",
  dims = c(16, 16, 8),
  n_clusters = 4,
  n_time = 40,
  noise_sd = 0.08,
  seed = 42
)
fit <- cluster4d(syn$vec, syn$mask, n_clusters = 4, method = "rena")
audit <- validate_cluster4d(fit, syn$vec, syn$mask)
quality <- cluster_metrics(fit, truth = syn$truth)

stopifnot(
  isTRUE(audit$valid),
  fit$actual_k == 4L,
  is.finite(quality$ari_truth),
  quality$ari_truth > 0.99
)

# Cluster IDs are arbitrary. Align them to the known synthetic labels only for
# side-by-side display; the reported agreement metrics are permutation-invariant.
label_ids <- sort(unique(fit$labels))
display_map <- vapply(label_ids, function(id) {
  counts <- table(syn$truth[fit$labels == id])
  as.integer(names(counts)[which.max(counts)])
}, integer(1))
aligned_labels <- display_map[fit$labels]

truth <- array(syn$truth, dim = syn$dims)
recovered <- array(aligned_labels, dim = syn$dims)
z <- ceiling(syn$dims[[3L]] / 2)
palette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
ink <- "#212121"
muted <- "#6B625B"

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
png(
  "man/figures/neurocluster-workflow.png",
  width = 1800,
  height = 640,
  res = 180,
  bg = "#FBF8F2"
)
old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  dev.off()
}, add = TRUE)

layout(matrix(1:3, nrow = 1), widths = c(1.25, 1, 1))
par(mar = c(3.6, 3.8, 3.2, 1.2), family = "sans", col.axis = muted,
    col.lab = ink, col.main = ink, bty = "n")

matplot(
  t(syn$patterns),
  type = "l",
  lty = 1:4,
  lwd = 2.5,
  col = palette,
  xlab = "Time point",
  ylab = "Standardized signal",
  main = "A  Four temporal signatures"
)
abline(h = 0, col = "#D8D0C4", lwd = 1)
legend(
  "topright",
  legend = paste("Region", 1:4),
  col = palette,
  lty = 1:4,
  lwd = 2.5,
  ncol = 2,
  bty = "n",
  text.col = ink,
  cex = 0.78
)

draw_slice <- function(x, title) {
  image(
    t(x[, ncol(x):1, drop = FALSE]),
    col = palette,
    breaks = seq(0.5, 4.5, by = 1),
    axes = FALSE,
    asp = 1,
    main = title
  )
  box(col = "#D8D0C4")
}

par(mar = c(1.5, 1.5, 3.2, 1.2))
draw_slice(truth[, , z], "B  Known spatial regions")
draw_slice(
  recovered[, , z],
  sprintf("C  ReNA recovery  |  ARI %.2f", quality$ari_truth)
)
