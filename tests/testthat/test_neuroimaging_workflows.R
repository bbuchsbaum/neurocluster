library(neurocluster)
library(neuroim2)

# Real-World Neuroimaging Workflow Tests
# Tests for realistic neuroimaging analysis scenarios

test_that("resting state network parcellation workflow", {
  # Simulate resting-state fMRI data with network structure
  dims <- c(16, 16, 3)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 200  # Typical resting-state length
  
  # Create network-like structure
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  # Simulate resting-state networks with realistic properties
  t_seq <- seq(0, 20*pi, length.out = ntime)  # Multiple cycles
  network_seeds <- list(
    # Default Mode Network (DMN) - low frequency
    dmn = list(centers = rbind(c(4, 4, 2), c(12, 12, 2)), freq = 0.5),
    # Visual Network - medium frequency  
    visual = list(centers = rbind(c(8, 2, 1), c(8, 14, 1)), freq = 1.0),
    # Motor Network - higher frequency
    motor = list(centers = rbind(c(2, 8, 3), c(14, 8, 3)), freq = 1.5)
  )
  
  # Assign voxels to networks based on distance to seeds
  for (i in 1:nvox) {
    voxel_coord <- coords[i, ]
    
    # Find closest network
    min_dist <- Inf
    closest_network <- "dmn"
    
    for (net_name in names(network_seeds)) {
      net_centers <- network_seeds[[net_name]]$centers
      for (j in 1:nrow(net_centers)) {
        dist <- sqrt(sum((voxel_coord - net_centers[j, ])^2))
        if (dist < min_dist) {
          min_dist <- dist
          closest_network <- net_name
        }
      }
    }
    
    # Generate network-specific signal
    freq <- network_seeds[[closest_network]]$freq
    base_signal <- sin(freq * t_seq)
    
    # Add distance-based attenuation and noise
    attenuation <- exp(-min_dist / 3)  # Signal strength decreases with distance
    noise_level <- 0.3 + (1 - attenuation) * 0.5  # More noise for distant voxels
    
    ts_data[i, ] <- attenuation * base_signal + rnorm(ntime, sd = noise_level)
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test resting-state parcellation workflow
  expect_silent(rs_result <- slice_msf(vec, mask, 
                                      r = 12,  # Higher rank for resting-state
                                      min_size = 20,
                                      compactness = 2,  # Allow more irregular shapes
                                      num_runs = 3,
                                      consensus = TRUE))
  
  expect_true(!is.null(rs_result))
  expect_true(length(rs_result$cluster) == nvox)
  
  # Should find reasonable number of networks (3-10)
  n_networks <- length(unique(rs_result$cluster))
  expect_true(n_networks >= 3 && n_networks <= 10,
              info = sprintf("Should find 3-10 networks, found %d", n_networks))
  
  # Test network connectivity analysis
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- rs_result$cluster
  
  # Calculate average within-network correlation vs between-network correlation
  within_corr <- c()
  between_corr <- c()
  
  unique_clusters <- sort(unique(rs_result$cluster))
  for (i in 1:(length(unique_clusters)-1)) {
    for (j in (i+1):length(unique_clusters)) {
      cluster_i_voxels <- which(rs_result$cluster == unique_clusters[i])
      cluster_j_voxels <- which(rs_result$cluster == unique_clusters[j])
      
      if (length(cluster_i_voxels) > 1 && length(cluster_j_voxels) > 1) {
        # Sample voxels to avoid excessive computation
        sample_i <- sample(cluster_i_voxels, min(5, length(cluster_i_voxels)))
        sample_j <- sample(cluster_j_voxels, min(5, length(cluster_j_voxels)))
        
        for (vi in sample_i) {
          for (vj in sample_j) {
            corr_val <- cor(ts_data[vi, ], ts_data[vj, ])
            if (!is.na(corr_val)) {
              between_corr <- c(between_corr, corr_val)
            }
          }
        }
        
        # Within-network correlations for cluster i
        if (length(sample_i) > 1) {
          for (vi1 in 1:(length(sample_i)-1)) {
            for (vi2 in (vi1+1):length(sample_i)) {
              corr_val <- cor(ts_data[sample_i[vi1], ], ts_data[sample_i[vi2], ])
              if (!is.na(corr_val)) {
                within_corr <- c(within_corr, corr_val)
              }
            }
          }
        }
      }
    }
  }
  
  if (length(within_corr) > 0 && length(between_corr) > 0) {
    mean_within_corr <- mean(within_corr)
    mean_between_corr <- mean(between_corr)
    
    # Within-network correlations should be higher than between-network
    expect_true(mean_within_corr > mean_between_corr,
                info = sprintf("Within-network correlation (%.3f) should be > between-network (%.3f)",
                              mean_within_corr, mean_between_corr))
  }
})

test_that("task-related clustering workflow", {
  # Simulate task fMRI with event-related responses
  dims <- c(12, 12, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 120  # Task run length
  
  # Create task design
  TR <- 2  # 2-second repetition time
  task_onsets <- seq(20, 100, by = 20)  # Task events every 20 TRs
  task_duration <- 4  # 4 TR task blocks
  
  # Create canonical HRF
  create_hrf <- function(time_points, peak_delay = 6, undershoot_delay = 16) {
    t <- seq(0, 30, length.out = time_points)
    hrf <- (t / peak_delay)^2 * exp(-t / peak_delay) - 
           0.1 * (t / undershoot_delay)^2 * exp(-t / undershoot_delay)
    hrf / max(hrf)  # Normalize
  }
  
  hrf <- create_hrf(15)  # 15 time point HRF
  
  # Create task design matrix
  task_signal <- numeric(ntime)
  for (onset in task_onsets) {
    for (dur in 1:task_duration) {
      if ((onset + dur - 1) <= ntime) {
        task_signal[onset + dur - 1] <- 1
      }
    }
  }
  
  # Convolve with HRF
  task_convolved <- convolve(task_signal, rev(hrf), type = "open")[1:ntime]
  
  # Create task-related regions
  coords <- arrayInd(1:nvox, dims)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)
  
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    z <- coords[i, 3]
    
    # Define task-responsive regions
    if (x <= 6 && y <= 6) {
      # Strong task response
      activation_strength <- 2.0
      ts_data[i, ] <- activation_strength * task_convolved + 
                     rnorm(ntime, sd = 0.5) + 
                     sin(seq(0, 4*pi, length.out = ntime)) * 0.3  # baseline oscillation
    } else if (x > 6 && y <= 6) {
      # Moderate task response  
      activation_strength <- 1.0
      ts_data[i, ] <- activation_strength * task_convolved + 
                     rnorm(ntime, sd = 0.5) +
                     cos(seq(0, 2*pi, length.out = ntime)) * 0.3
    } else if (x <= 6 && y > 6) {
      # Weak task response
      activation_strength <- 0.3
      ts_data[i, ] <- activation_strength * task_convolved + 
                     rnorm(ntime, sd = 0.5) +
                     sin(seq(0, 6*pi, length.out = ntime)) * 0.2
    } else {
      # No task response (control region)
      ts_data[i, ] <- rnorm(ntime, sd = 0.5) +
                     0.1 * sin(seq(0, 8*pi, length.out = ntime))
    }
  }
  
  # Create NeuroVec
  vec_list <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Test task-related clustering
  expect_silent(task_result <- slice_msf(vec, mask,
                                        r = 8,
                                        min_size = 15,
                                        compactness = 3,
                                        num_runs = 1))
  
  expect_true(!is.null(task_result))
  
  # Should find regions with different activation profiles
  n_clusters <- length(unique(task_result$cluster))
  expect_true(n_clusters >= 2 && n_clusters <= 8,
              info = sprintf("Should find 2-8 task-related clusters, found %d", n_clusters))
  
  # Validate that clustering captures task-related variance
  cluster_vol <- array(0, dims)
  cluster_vol[mask > 0] <- task_result$cluster
  
  # Calculate task-relatedness for each cluster
  cluster_task_correlation <- c()
  unique_clusters <- sort(unique(task_result$cluster))
  
  for (cluster_id in unique_clusters) {
    cluster_voxels <- which(task_result$cluster == cluster_id)
    if (length(cluster_voxels) > 0) {
      # Average signal in this cluster
      cluster_mean_signal <- colMeans(ts_data[cluster_voxels, , drop = FALSE])
      
      # Correlation with task design
      task_corr <- cor(cluster_mean_signal, task_convolved)
      if (!is.na(task_corr)) {
        cluster_task_correlation <- c(cluster_task_correlation, abs(task_corr))
      }
    }
  }
  
  # At least some clusters should show task-relatedness
  if (length(cluster_task_correlation) > 0) {
    max_task_corr <- max(cluster_task_correlation)
    expect_true(max_task_corr > 0.3,
                info = sprintf("At least one cluster should be task-related (max corr=%.3f)", max_task_corr))
  }
})

test_that("multi-subject clustering workflow", {
  # Simulate multi-subject data for group-level analysis
  dims <- c(10, 10, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 60
  n_subjects <- 3
  
  # Create consistent network structure across subjects with individual variations
  coords <- arrayInd(1:nvox, dims)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Define group-level network centers
  group_networks <- list(
    net1 = rbind(c(3, 3, 1), c(7, 7, 1)),
    net2 = rbind(c(3, 7, 2), c(7, 3, 2))
  )
  
  subject_results <- list()
  
  for (subj in 1:n_subjects) {
    # Individual subject variations
    subject_shift <- rnorm(2, mean = 0, sd = 0.5)  # Random spatial shifts
    subject_amplitude <- runif(2, 0.8, 1.2)       # Random amplitude scaling
    
    ts_data_subj <- matrix(0, nrow = nvox, ncol = ntime)
    
    for (i in 1:nvox) {
      voxel_coord <- coords[i, ]
      
      # Find closest network with subject-specific variations
      min_dist <- Inf
      closest_net <- 1
      
      for (net_idx in 1:length(group_networks)) {
        net_centers <- group_networks[[net_idx]]
        # Add subject-specific spatial variations
        net_centers[, 1:2] <- net_centers[, 1:2] + 
          matrix(rep(subject_shift, nrow(net_centers)), nrow = nrow(net_centers), byrow = TRUE)
        
        for (j in 1:nrow(net_centers)) {
          dist <- sqrt(sum((voxel_coord - net_centers[j, ])^2))
          if (dist < min_dist) {
            min_dist <- dist
            closest_net <- net_idx
          }
        }
      }
      
      # Generate network signal with subject variations
      base_freq <- closest_net * 0.8  # Different frequencies for networks
      amplitude <- subject_amplitude[closest_net]
      signal_strength <- exp(-min_dist / 2.5)  # Distance-based strength
      
      ts_data_subj[i, ] <- amplitude * signal_strength * sin(base_freq * t_seq) + 
                          rnorm(ntime, sd = 0.3)
    }
    
    # Create subject NeuroVec
    vec_list_subj <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data_subj[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec_subj <- do.call(concat, vec_list_subj)
    
    # Individual subject clustering
    expect_silent(subject_results[[subj]] <- slice_msf(vec_subj, mask,
                                                      r = 6,
                                                      min_size = 10,
                                                      compactness = 3,
                                                      num_runs = 1))
  }
  
  # All subjects should have valid results
  for (subj in 1:n_subjects) {
    expect_true(!is.null(subject_results[[subj]]))
    expect_true(length(subject_results[[subj]]$cluster) == nvox)
    n_clusters_subj <- length(unique(subject_results[[subj]]$cluster))
    expect_true(n_clusters_subj > 0,
                info = sprintf("Subject %d should find clusters", subj))
  }
  
  # Test group-level consistency
  # Calculate overlap in cluster assignments across subjects
  cluster_overlaps <- c()
  
  for (i in 1:(n_subjects-1)) {
    for (j in (i+1):n_subjects) {
      # Simple overlap measure: fraction of voxel pairs with consistent clustering
      clusters_i <- subject_results[[i]]$cluster
      clusters_j <- subject_results[[j]]$cluster
      
      # Create co-assignment matrices
      coassign_i <- outer(clusters_i, clusters_i, "==")
      coassign_j <- outer(clusters_j, clusters_j, "==")
      
      # Calculate overlap
      agreement <- sum(coassign_i == coassign_j) / length(coassign_i)
      cluster_overlaps <- c(cluster_overlaps, agreement)
    }
  }
  
  mean_overlap <- mean(cluster_overlaps)
  
  # Should have some consistency across subjects (>60%)
  expect_true(mean_overlap > 0.6,
              info = sprintf("Cross-subject consistency should be >0.6, got %.3f", mean_overlap))
  
  # Test group average clustering
  # Simple approach: average time series across subjects
  group_ts_data <- array(0, dim = c(nvox, ntime))
  
  for (subj in 1:n_subjects) {
    subj_data <- series(vec_subj, which(mask > 0))
    group_ts_data <- group_ts_data + t(subj_data)
  }
  group_ts_data <- group_ts_data / n_subjects
  
  # Create group-level NeuroVec
  vec_list_group <- lapply(1:ntime, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- group_ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec_group <- do.call(concat, vec_list_group)
  
  expect_silent(group_result <- slice_msf(vec_group, mask,
                                         r = 6,
                                         min_size = 10,
                                         compactness = 3,
                                         num_runs = 1))
  
  expect_true(!is.null(group_result))
  group_n_clusters <- length(unique(group_result$cluster))
  expect_true(group_n_clusters > 0,
              info = "Group-level clustering should find clusters")
})

test_that("high temporal resolution workflow", {
  syn <- make_block_synthetic(dims = c(10, 10, 1), ntime = 160, noise = 0.12, seed = 123)

  # Higher rank to reflect faster TR; single run to keep runtime low
  expect_silent(htr_result <- slice_msf(
    vec = syn$vec,
    mask = syn$mask,
    r = 16,
    min_size = 6,
    compactness = 2,
    num_runs = 1,
    stitch_z = FALSE
  ))

  n_clusters <- length(unique(htr_result$cluster))
  expect_true(n_clusters >= 2 && n_clusters <= 6,
              info = sprintf("High TR synthetic should find 2-6 clusters, found %d", n_clusters))

  # Temporal structure: each center should vary
  cluster_centers <- htr_result$centers
  cluster_vars <- apply(cluster_centers, 2, var)
  expect_true(all(cluster_vars > 0),
              info = "All cluster centers should show temporal variation")
})

test_that("clinical neuroimaging workflow robustness", {
  # Test robustness to clinical data characteristics (artifacts, motion, etc.)
  dims <- c(12, 12, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 80
  
  # Simulate clinical data artifacts
  coords <- arrayInd(1:nvox, dims)
  ts_data_clean <- matrix(0, nrow = nvox, ncol = ntime)
  t_seq <- seq(0, 4*pi, length.out = ntime)
  
  # Create base neural signals
  for (i in 1:nvox) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    
    if (x <= 6) {
      ts_data_clean[i, ] <- sin(t_seq) + rnorm(ntime, sd = 0.1)
    } else {
      ts_data_clean[i, ] <- cos(t_seq) + rnorm(ntime, sd = 0.1)
    }
  }
  
  # Add clinical artifacts
  clinical_scenarios <- list(
    # Motion artifacts (sudden signal changes)
    motion = {
      ts_data_motion <- ts_data_clean
      motion_timepoints <- c(15, 35, 55)  # Sudden motion events
      for (tp in motion_timepoints) {
        if (tp <= ntime) {
          motion_magnitude <- rnorm(nvox, mean = 0, sd = 1.5)
          ts_data_motion[, tp] <- ts_data_motion[, tp] + motion_magnitude
          # Gradual return to baseline over next 3 timepoints
          if (tp + 3 <= ntime) {
            for (decay in 1:3) {
              ts_data_motion[, tp + decay] <- ts_data_motion[, tp + decay] + 
                motion_magnitude * exp(-decay)
            }
          }
        }
      }
      ts_data_motion
    },
    
    # Scanner drift (slow signal changes)
    drift = {
      ts_data_drift <- ts_data_clean
      drift_pattern <- seq(-0.5, 0.5, length.out = ntime)  # Linear drift
      for (i in 1:nvox) {
        ts_data_drift[i, ] <- ts_data_drift[i, ] + drift_pattern + 
          0.1 * sin(seq(0, pi, length.out = ntime))  # Slow oscillation
      }
      ts_data_drift
    },
    
    # Dropout/signal loss (missing data simulation)
    dropout = {
      ts_data_dropout <- ts_data_clean
      # Simulate signal dropout in a region (e.g., due to susceptibility artifact)
      dropout_region <- coords[, 1] >= 9 & coords[, 1] <= 12 & coords[, 2] >= 6
      dropout_voxels <- which(dropout_region)
      dropout_timepoints <- 25:45  # Period of signal loss
      
      for (tp in dropout_timepoints) {
        if (tp <= ntime) {
          ts_data_dropout[dropout_voxels, tp] <- ts_data_dropout[dropout_voxels, tp] * 0.1 +
            rnorm(length(dropout_voxels), sd = 0.5)  # Severe noise
        }
      }
      ts_data_dropout
    }
  )
  
  for (scenario_name in names(clinical_scenarios)) {
    ts_data_artifact <- clinical_scenarios[[scenario_name]]
    
    vec_list_artifact <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data_artifact[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec_artifact <- do.call(concat, vec_list_artifact)
    
    # Test robustness to artifacts
    expect_silent({
      clinical_result <- slice_msf(vec_artifact, mask,
                                  r = 8,
                                  min_size = 12,
                                  compactness = 4,  # More compact for noisy data
                                  num_runs = 3,     # Multiple runs for stability
                                  consensus = TRUE)
    })
    
    expect_true(!is.null(clinical_result),
                info = sprintf("Should produce result with %s artifacts", scenario_name))
    
    expect_true(length(clinical_result$cluster) == nvox,
                info = sprintf("Should cluster all voxels with %s artifacts", scenario_name))
    
    n_clusters <- length(unique(clinical_result$cluster))
    expect_true(n_clusters > 0 && n_clusters <= 15,
                info = sprintf("Should find reasonable clusters with %s artifacts (found %d)", 
                              scenario_name, n_clusters))
    
    # Test that consensus helps with artifacts
    expect_true(length(clinical_result$runs) == 3,
                info = sprintf("Should store multiple runs for %s scenario", scenario_name))
  }
})
