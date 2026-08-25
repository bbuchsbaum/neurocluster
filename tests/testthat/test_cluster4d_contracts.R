contract_fixture <- function(mask_values = NULL, bad_feature = NULL,
                             mask_spacing = c(1, 1, 1),
                             mask_origin = c(0, 0, 0)) {
  dims <- c(3L, 3L, 2L)
  vec_space <- neuroim2::NeuroSpace(
    dims, spacing = c(1, 1, 1), origin = c(0, 0, 0)
  )
  frames <- lapply(seq_len(4L), function(i) {
    voxel <- seq_len(prod(dims))
    values <- array(
      sin(voxel * (i + 1) * 0.37) + cos(voxel * 0.11 + i),
      dims
    )
    if (!is.null(bad_feature) && i == 2L) values[2] <- bad_feature
    neuroim2::NeuroVol(values, vec_space)
  })
  vec <- do.call(neuroim2::concat, frames)

  if (is.null(mask_values)) mask_values <- array(1, dims)
  mask_space <- neuroim2::NeuroSpace(
    dims, spacing = mask_spacing, origin = mask_origin
  )
  mask <- neuroim2::NeuroVol(array(mask_values, dims), mask_space)
  list(vec = vec, mask = mask)
}

cluster4d_methods <- c(
  "supervoxels", "snic", "slic", "corr_slic", "brs_slic", "slice_msf",
  "flash3d", "g3s", "rena", "rena_plus", "mcl", "acsc", "commute"
)

direct_cluster4d_methods <- list(
  supervoxels = cluster4d_supervoxels,
  snic = cluster4d_snic,
  slic = cluster4d_slic,
  corr_slic = cluster4d_corrslic,
  brs_slic = cluster4d_brsslic,
  slice_msf = cluster4d_slice_msf,
  flash3d = cluster4d_flash3d,
  g3s = cluster4d_g3s,
  rena = cluster4d_rena,
  rena_plus = cluster4d_rena_plus,
  mcl = cluster4d_mcl,
  acsc = cluster4d_acsc,
  commute = cluster4d_commute
)

mock_cluster4d_dispatch <- function() {
  sentinel <- function(...) stop("CONTRACT_REACHED_DISPATCH", call. = FALSE)
  testthat::local_mocked_bindings(
    cluster4d_supervoxels = sentinel,
    cluster4d_snic = sentinel,
    cluster4d_slic = sentinel,
    cluster4d_corrslic = sentinel,
    cluster4d_brsslic = sentinel,
    cluster4d_slice_msf = sentinel,
    cluster4d_flash3d = sentinel,
    cluster4d_g3s = sentinel,
    cluster4d_rena = sentinel,
    cluster4d_rena_plus = sentinel,
    cluster4d_mcl = sentinel,
    cluster4d_acsc = sentinel,
    cluster4d_commute = sentinel,
    .package = "neurocluster",
    .env = parent.frame()
  )
}

test_that("every unified method shares the same pre-dispatch input contract", {
  fixture <- contract_fixture()
  mock_cluster4d_dispatch()

  for (method in cluster4d_methods) {
    expect_error(
      cluster4d(fixture$vec, fixture$mask, 2, method = method),
      "CONTRACT_REACHED_DISPATCH",
      info = paste(method, "accepts the common valid input")
    )

    for (bad_k in list(1.5, NaN, Inf, 0, 19)) {
      expect_error(
        cluster4d(fixture$vec, fixture$mask, bad_k, method = method),
        "n_clusters|Cannot create",
        info = paste(method, "rejects invalid K", format(bad_k))
      )
    }
  }
})

test_that("geometry equality includes spacing and affine, not dimensions alone", {
  spacing_bad <- contract_fixture(mask_spacing = c(2, 1, 1))
  origin_bad <- contract_fixture(mask_origin = c(5, 0, 0))

  for (method in cluster4d_methods) {
    expect_error(
      cluster4d(spacing_bad$vec, spacing_bad$mask, 2, method = method),
      "identical spatial spacing",
      info = method
    )
    expect_error(
      cluster4d(origin_bad$vec, origin_bad$mask, 2, method = method),
      "identical spatial affine",
      info = method
    )
  }
})

test_that("mask inclusion is finite and strictly positive", {
  values <- array(0, c(3, 3, 2))
  values[c(1, 4, 7)] <- c(-2, 0.25, 3)
  fixture <- contract_fixture(mask_values = values)

  for (method in cluster4d_methods) {
    contract <- validate_cluster4d_inputs(
      fixture$vec, fixture$mask, 2, paste0("cluster4d:", method)
    )
    expect_identical(contract$mask_idx, c(4L, 7L), info = method)
    expect_identical(which(contract$mask), c(4L, 7L), info = method)
  }

  for (bad_mask in list(NaN, Inf, -Inf)) {
    bad_values <- values
    bad_values[1] <- bad_mask
    bad <- contract_fixture(mask_values = bad_values)
    expect_error(
      validate_cluster4d_inputs(bad$vec, bad$mask, 2),
      "mask values must be finite"
    )
  }
})

test_that("non-finite feature data inside the declared mask fail closed", {
  for (bad_value in list(NaN, Inf, -Inf)) {
    fixture <- contract_fixture(bad_feature = bad_value)
    for (method in cluster4d_methods) {
      expect_error(
        cluster4d(fixture$vec, fixture$mask, 2, method = method),
        "feature data contain 1 non-finite value",
        info = paste(method, format(bad_value))
      )
    }
  }
})

test_that("common parameter capabilities and endpoints are fail-closed", {
  fixture <- contract_fixture()
  mock_cluster4d_dispatch()

  spatial_methods <- cluster4d_methods[!cluster4d_methods %in% c("rena", "rena_plus")]
  for (method in spatial_methods) {
    for (endpoint in c(0, 1)) {
      expect_error(
        cluster4d(
          fixture$vec, fixture$mask, 2, method = method,
          spatial_weight = endpoint
        ),
        "CONTRACT_REACHED_DISPATCH",
        info = paste(method, "accepts spatial_weight", endpoint)
      )
    }
  }
  for (method in c("rena", "rena_plus")) {
    expect_error(
      cluster4d(
        fixture$vec, fixture$mask, 2, method = method, spatial_weight = 0.5
      ),
      "spatial_weight is not supported",
      info = method
    )
  }

  for (method in c("snic", "slice_msf", "commute")) {
    expect_error(
      cluster4d(
        fixture$vec, fixture$mask, 2, method = method, max_iterations = 2
      ),
      "max_iterations is not supported",
      info = method
    )
  }
  for (method in c("snic", "flash3d", "acsc", "commute")) {
    expect_error(
      cluster4d(
        fixture$vec, fixture$mask, 2, method = method, connectivity = 6
      ),
      "connectivity is not supported",
      info = method
    )
  }
  for (method in c(
    "snic", "corr_slic", "brs_slic", "slice_msf", "flash3d", "g3s",
    "rena", "rena_plus", "mcl", "acsc", "commute"
  )) {
    expect_error(
      cluster4d(fixture$vec, fixture$mask, 2, method = method, parallel = TRUE),
      "parallel is not supported",
      info = method
    )
  }

  expect_error(
    cluster4d(fixture$vec, fixture$mask, 2, method = "slic", connectivity = 18),
    "connectivity must be one of 6, 26"
  )
  expect_error(
    cluster4d(fixture$vec, fixture$mask, 2, method = "mcl", connectivity = 6.5),
    "integer-valued"
  )
  expect_error(
    cluster4d(fixture$vec, fixture$mask, 2, method = "g3s", max_iterations = Inf),
    "finite numeric scalar"
  )
  expect_error(
    cluster4d(fixture$vec, fixture$mask, 2, method = "supervoxels", verbose = NA),
    "must be TRUE or FALSE"
  )
})

test_that("exported wrappers reject common parameters they do not implement", {
  fixture <- contract_fixture()
  expect_error(
    cluster4d_snic(fixture$vec, fixture$mask, 2, max_iterations = 2),
    "max_iterations is not supported"
  )
  expect_error(
    cluster4d_slice_msf(fixture$vec, fixture$mask, 2, parallel = TRUE),
    "parallel is not supported"
  )
  expect_error(
    cluster4d_rena(fixture$vec, fixture$mask, 2, spatial_weight = 0.5),
    "spatial_weight is not supported"
  )
  expect_error(
    cluster4d_rena_plus(fixture$vec, fixture$mask, 2, spatial_weight = 0.5),
    "spatial_weight is not supported"
  )
  expect_error(
    cluster4d_mcl(fixture$vec, fixture$mask, 2, parallel = TRUE),
    "unused argument.*parallel"
  )
  expect_error(
    cluster4d_corrslic(fixture$vec, fixture$mask, 2, parallel = TRUE),
    "parallel is not supported"
  )
  expect_error(
    cluster4d_brsslic(fixture$vec, fixture$mask, 2, parallel = TRUE),
    "parallel is not supported"
  )
})

test_that("exported wrappers all apply the shared input contract", {
  fixture <- contract_fixture()
  geometry_bad <- contract_fixture(mask_origin = c(5, 0, 0))
  feature_bad <- contract_fixture(bad_feature = NaN)

  for (method in names(direct_cluster4d_methods)) {
    wrapper <- direct_cluster4d_methods[[method]]
    expect_error(
      wrapper(fixture$vec, fixture$mask, 1.5),
      "integer-valued",
      info = method
    )
    expect_error(
      wrapper(geometry_bad$vec, geometry_bad$mask, 2),
      "identical spatial affine",
      info = method
    )
    expect_error(
      wrapper(feature_bad$vec, feature_bad$mask, 2),
      "feature data contain 1 non-finite value",
      info = method
    )
  }
})

test_that("advertised spatial_weight endpoints reach real implementations", {
  fixture <- contract_fixture()
  endpoint_args <- list(
    supervoxels = list(max_iterations = 1, parallel = FALSE),
    snic = list(),
    slic = list(
      max_iterations = 1, parallel = FALSE, seed_relocate = "none",
      preserve_k = TRUE
    ),
    corr_slic = list(
      max_iterations = 1, parallel = FALSE, embedding_dim = 8,
      refine_exact_iters = 0
    ),
    brs_slic = list(
      max_iterations = 1, parallel = FALSE, embedding_dim = 8,
      boundary_passes = 0, global_passes = 0
    ),
    slice_msf = list(
      num_runs = 1, consensus = FALSE, stitch_z = TRUE, r = 2
    ),
    flash3d = list(max_iterations = 1, dctM = 4),
    g3s = list(
      max_iterations = 1, n_components = 2,
      use_irlba = FALSE, use_rsvd = FALSE
    ),
    mcl = list(max_iterations = 2, connectivity = 6),
    acsc = list(max_iterations = 1, block_size = 2, refine = FALSE),
    commute = list(ncomp = 2)
  )

  for (method in names(endpoint_args)) {
    for (endpoint in c(0, 1)) {
      args <- c(
        list(
          vec = fixture$vec,
          mask = fixture$mask,
          n_clusters = 2,
          method = method,
          spatial_weight = endpoint,
          verbose = FALSE
        ),
        endpoint_args[[method]]
      )
      result <- suppressWarnings(do.call(cluster4d, args))
      expect_s3_class(result, "cluster4d_result")
      expect_identical(
        length(result$cluster), 18L,
        info = paste(method, "preserves the declared mask at", endpoint)
      )
    }
  }
})
