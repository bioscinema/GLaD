polar_cluster <- function(dist_matrix, k, 
                          max_cluster_size = Inf,
                          first_pole_method = c("median", "max", "random"),
                          subsequent_pole_quantile = 0.75,
                          assignment_method = c("nearest", "balanced"),
                          outlier_filter = TRUE,
                          outlier_threshold = 1.5) {
  #' Polar Ordination Clustering with Robust Pole Selection and Flexible Assignment
  #' 
  #' Performs clustering based on polar ordination with robust pole selection
  #' to avoid outliers and flexible assignment methods.
  #' 
  #' @param dist_matrix A distance matrix (dist object or matrix) with complete dimnames
  #' @param k Number of clusters (poles) to identify
  #' @param max_cluster_size Maximum number of units allowed per cluster (default: Inf)
  #' @param first_pole_method Method for selecting first pole: "median" (median of row sums),
  #'                          "max" (maximum of row sums), or "random" (random selection)
  #' @param subsequent_pole_quantile Quantile to use for subsequent pole selection (default: 0.75).
  #'                                 Use 1.0 for maximum, 0.75 for Q3 (robust), 0.5 for median
  #' @param assignment_method Method for assigning units to poles:
  #'                          "nearest" - assign to nearest pole with size constraints
  #'                          "balanced" - enforce equal cluster sizes (ignores max_cluster_size)
  #' @param outlier_filter If TRUE, exclude outliers when selecting poles (default: TRUE)
  #' @param outlier_threshold IQR multiplier for outlier detection (default: 1.5)
  #' @return A list of class 'polar_cluster' containing clustering results
  #' @examples
  #' dist_mat <- dist(matrix(rnorm(50), ncol = 5))
  #' # Robust method with balanced clusters
  #' result <- polar_cluster(dist_mat, k = 3, 
  #'                         first_pole_method = "median",
  #'                         assignment_method = "balanced",
  #'                         outlier_filter = TRUE)
  
  # Match arguments
  first_pole_method <- match.arg(first_pole_method)
  assignment_method <- match.arg(assignment_method)
  
  # Validate input parameters
  if (!inherits(dist_matrix, "dist") && !is.matrix(dist_matrix)) {
    stop("dist_matrix must be a dist object or matrix")
  }
  
  # Convert to matrix format if necessary
  if (inherits(dist_matrix, "dist")) {
    dist_mat <- as.matrix(dist_matrix)
  } else {
    dist_mat <- dist_matrix
  }
  
  # Check for complete dimnames
  if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
    stop("dist_matrix must have complete dimnames")
  }
  
  n <- nrow(dist_mat)
  unit_names <- rownames(dist_mat)
  
  if (k < 1 || k > n) {
    stop("k must be between 1 and the number of units")
  }
  
  if (max_cluster_size < 1) {
    stop("max_cluster_size must be at least 1")
  }
  
  if (subsequent_pole_quantile < 0 || subsequent_pole_quantile > 1) {
    stop("subsequent_pole_quantile must be between 0 and 1")
  }
  
  if (assignment_method == "nearest" && max_cluster_size * k < n) {
    warning("max_cluster_size * k < n: not all units can be assigned. ",
            "Increasing max_cluster_size to ", ceiling(n / k))
    max_cluster_size <- ceiling(n / k)
  }
  
  # ========== Detect and exclude outliers for pole selection ==========
  outlier_units <- character(0)
  pole_candidates <- unit_names
  
  if (outlier_filter) {
    # Calculate average distance for each unit
    row_sums <- rowSums(dist_mat)
    
    # Use IQR method to detect outliers
    Q1 <- quantile(row_sums, 0.25)
    Q3 <- quantile(row_sums, 0.75)
    IQR_val <- Q3 - Q1
    
    # Units with row sum > Q3 + threshold * IQR are considered outliers
    outlier_mask <- row_sums > (Q3 + outlier_threshold * IQR_val)
    outlier_units <- unit_names[outlier_mask]
    
    # Pole candidates exclude outliers
    pole_candidates <- unit_names[!outlier_mask]
    
    if (length(pole_candidates) < k) {
      warning("Too few non-outlier units (", length(pole_candidates), 
              ") for k = ", k, ". Using all units as candidates.")
      pole_candidates <- unit_names
      outlier_units <- character(0)
    }
  }
  
  # ========== Step 1: Select poles from non-outlier candidates ==========
  poles <- character(k)
  
  # Calculate row sums for candidates only
  candidate_row_sums <- rowSums(dist_mat[pole_candidates, , drop = FALSE])
  
  # Select first pole based on specified method
  if (first_pole_method == "median") {
    # Find candidate closest to median of row sums
    median_val <- median(candidate_row_sums)
    first_pole_idx <- which.min(abs(candidate_row_sums - median_val))
    poles[1] <- pole_candidates[first_pole_idx]
    
  } else if (first_pole_method == "max") {
    # Use maximum row sum among candidates
    poles[1] <- pole_candidates[which.max(candidate_row_sums)]
    
  } else if (first_pole_method == "random") {
    # Random selection from candidates
    set.seed(NULL)
    poles[1] <- sample(pole_candidates, 1)
  }
  
  # Select subsequent poles
  for (i in 2:k) {
    if (i == 2) {
      # Second pole: use quantile of distances to first pole (among candidates)
      dist_to_pole1 <- dist_mat[pole_candidates, poles[1]]
      
      if (subsequent_pole_quantile == 1.0) {
        # Use maximum
        poles[2] <- pole_candidates[which.max(dist_to_pole1)]
      } else {
        # Use quantile-based selection
        threshold <- quantile(dist_to_pole1, subsequent_pole_quantile)
        candidates_above_threshold <- pole_candidates[dist_to_pole1 >= threshold]
        
        if (length(candidates_above_threshold) == 1) {
          poles[2] <- candidates_above_threshold
        } else {
          poles[2] <- candidates_above_threshold[which.max(dist_to_pole1[pole_candidates %in% candidates_above_threshold])]
        }
      }
      
    } else {
      # For 3rd+ poles: compute minimum distance to existing poles
      min_dist_to_poles <- apply(dist_mat[pole_candidates, poles[1:(i-1)], drop = FALSE], 1, min)
      
      # Exclude already selected poles
      available_mask <- !pole_candidates %in% poles[1:(i-1)]
      available_dists <- min_dist_to_poles[available_mask]
      available_candidates <- pole_candidates[available_mask]
      
      if (subsequent_pole_quantile == 1.0) {
        # Use maximum
        poles[i] <- available_candidates[which.max(available_dists)]
      } else {
        # Use quantile-based selection
        threshold <- quantile(available_dists, subsequent_pole_quantile)
        candidates_above_threshold <- available_candidates[available_dists >= threshold]
        
        if (length(candidates_above_threshold) == 1) {
          poles[i] <- candidates_above_threshold
        } else {
          poles[i] <- candidates_above_threshold[which.max(available_dists[available_candidates %in% candidates_above_threshold])]
        }
      }
    }
  }
  
  # ========== Step 2: Assign ALL units (including outliers) to poles ==========
  clustering <- integer(n)
  names(clustering) <- unit_names
  reassignments <- 0
  
  # Calculate distance of each unit to each pole
  dist_to_poles_mat <- dist_mat[, poles, drop = FALSE]
  colnames(dist_to_poles_mat) <- 1:k
  
  if (assignment_method == "nearest") {
    # ===== Method 1: Assign to nearest pole with size constraints =====
    
    # Sort units by their minimum distance to any pole
    min_dist_to_any_pole <- apply(dist_to_poles_mat, 1, min)
    unit_order <- order(min_dist_to_any_pole)
    
    # Track current size of each cluster
    current_sizes <- integer(k)
    
    for (idx in unit_order) {
      unit <- unit_names[idx]
      dists <- dist_to_poles_mat[unit, ]
      pole_preference <- order(dists)
      
      assigned <- FALSE
      for (pole_id in pole_preference) {
        if (current_sizes[pole_id] < max_cluster_size) {
          clustering[unit] <- pole_id
          current_sizes[pole_id] <- current_sizes[pole_id] + 1
          assigned <- TRUE
          
          if (pole_id != pole_preference[1]) {
            reassignments <- reassignments + 1
          }
          break
        }
      }
      
      if (!assigned) {
        # Force assign to cluster with minimum current size
        pole_id <- which.min(current_sizes)
        clustering[unit] <- pole_id
        current_sizes[pole_id] <- current_sizes[pole_id] + 1
        reassignments <- reassignments + 1
        warning("Had to exceed max_cluster_size for cluster ", pole_id)
      }
    }
    
  } else if (assignment_method == "balanced") {
    # ===== Method 2: Force balanced cluster sizes =====
    
    # Target size for each cluster
    target_size <- floor(n / k)
    extra_units <- n %% k  # Some clusters will have one extra unit
    
    # Determine target size for each cluster
    cluster_targets <- rep(target_size, k)
    if (extra_units > 0) {
      cluster_targets[1:extra_units] <- target_size + 1
    }
    
    # Track current size of each cluster
    current_sizes <- integer(k)
    
    # Get all pairwise (unit, pole) combinations with distances
    assignment_options <- data.frame(
      unit = rep(unit_names, each = k),
      pole = rep(1:k, times = n),
      distance = as.vector(t(dist_to_poles_mat))
    )
    
    # Sort by distance (prefer closer assignments)
    assignment_options <- assignment_options[order(assignment_options$distance), ]
    
    # Greedy assignment: go through options in order of distance
    for (i in 1:nrow(assignment_options)) {
      unit <- assignment_options$unit[i]
      pole <- assignment_options$pole[i]
      
      # Skip if unit already assigned
      if (clustering[unit] > 0) next
      
      # Skip if this cluster is full
      if (current_sizes[pole] >= cluster_targets[pole]) next
      
      # Assign
      clustering[unit] <- pole
      current_sizes[pole] <- current_sizes[pole] + 1
      
      # Stop if all units assigned
      if (all(clustering > 0)) break
    }
    
    # If any units still unassigned (shouldn't happen), force assign
    unassigned <- which(clustering == 0)
    if (length(unassigned) > 0) {
      for (unit_idx in unassigned) {
        # Assign to cluster with minimum current size
        pole_id <- which.min(current_sizes)
        clustering[unit_names[unit_idx]] <- pole_id
        current_sizes[pole_id] <- current_sizes[pole_id] + 1
      }
    }
  }
  
  # ========== Step 3: Calculate additional statistics ==========
  
  # Distance of each unit to its assigned pole
  dist_to_pole <- numeric(n)
  names(dist_to_pole) <- unit_names
  for (i in 1:n) {
    unit <- unit_names[i]
    pole_id <- clustering[i]
    dist_to_pole[i] <- dist_mat[unit, poles[pole_id]]
  }
  
  # Size of each cluster
  cluster_sizes <- table(factor(clustering, levels = 1:k))
  
  # Within-cluster cohesion: average distance to pole for each cluster
  within_cluster_dist <- sapply(1:k, function(cluster_id) {
    units_in_cluster <- names(clustering)[clustering == cluster_id]
    if (length(units_in_cluster) == 0) return(NA)
    mean(dist_to_pole[units_in_cluster])
  })
  names(within_cluster_dist) <- paste0("Cluster_", 1:k)
  
  # Calculate silhouette coefficient for clustering quality assessment
  silhouette <- numeric(n)
  names(silhouette) <- unit_names
  
  for (i in 1:n) {
    unit <- unit_names[i]
    cluster_id <- clustering[i]
    
    # a(i): average distance to other units in same cluster
    same_cluster <- names(clustering)[clustering == cluster_id]
    if (length(same_cluster) > 1) {
      a_i <- mean(dist_mat[unit, setdiff(same_cluster, unit)])
    } else {
      a_i <- 0
    }
    
    # b(i): minimum average distance to units in other clusters
    other_clusters <- setdiff(1:k, cluster_id)
    if (length(other_clusters) > 0) {
      b_i <- min(sapply(other_clusters, function(other_id) {
        other_cluster <- names(clustering)[clustering == other_id]
        if (length(other_cluster) == 0) return(Inf)
        mean(dist_mat[unit, other_cluster])
      }))
    } else {
      b_i <- 0
    }
    
    # Silhouette coefficient: (b - a) / max(a, b)
    if (max(a_i, b_i) > 0 && is.finite(b_i)) {
      silhouette[i] <- (b_i - a_i) / max(a_i, b_i)
    } else {
      silhouette[i] <- 0
    }
  }
  
  # Average silhouette width across all units
  avg_silhouette <- mean(silhouette)
  
  # ========== Return results ==========
  result <- list(
    clustering = clustering,           # Named vector of cluster assignments
    poles = poles,                     # Names of pole units
    pole_indices = match(poles, unit_names),  # Numeric indices of poles
    dist_to_pole = dist_to_pole,      # Distance to assigned pole for each unit
    cluster_sizes = as.vector(cluster_sizes), # Size of each cluster
    within_cluster_dist = within_cluster_dist, # Average dist to pole per cluster
    silhouette = silhouette,           # Silhouette coefficient per unit
    avg_silhouette = avg_silhouette,   # Overall clustering quality metric
    reassignments = reassignments,     # Number of units reassigned
    max_cluster_size = max_cluster_size, # The size constraint used
    first_pole_method = first_pole_method,  # Method used for first pole
    subsequent_pole_quantile = subsequent_pole_quantile, # Quantile used
    assignment_method = assignment_method,  # Assignment method used
    outlier_filter = outlier_filter,   # Whether outlier filtering was used
    outlier_threshold = if(outlier_filter) outlier_threshold else NULL,
    outliers_detected = outlier_units, # Units identified as outliers
    n_outliers = length(outlier_units), # Number of outliers
    call = match.call(),               # Function call for reference
    k = k,                             # Number of clusters
    n = n                              # Number of units
  )
  
  class(result) <- "polar_cluster"
  return(result)
}

# ========== Print method ==========
print.polar_cluster <- function(x, ...) {
  cat("Polar Ordination Clustering\n")
  cat("===========================\n\n")
  cat("Number of clusters (k):", x$k, "\n")
  cat("Number of units:", x$n, "\n")
  cat("Assignment method:", x$assignment_method, "\n")
  
  if (x$assignment_method == "nearest") {
    cat("Max cluster size:", 
        if (is.infinite(x$max_cluster_size)) "Inf (no limit)" else x$max_cluster_size, 
        "\n")
  }
  
  cat("First pole method:", x$first_pole_method, "\n")
  cat("Subsequent pole quantile:", x$subsequent_pole_quantile, "\n")
  cat("Outlier filter:", x$outlier_filter, "\n")
  
  if (x$outlier_filter) {
    cat("Outliers detected:", x$n_outliers, "units\n")
    if (x$n_outliers > 0 && x$n_outliers <= 5) {
      cat("  (", paste(x$outliers_detected, collapse = ", "), ")\n", sep = "")
    }
  }
  
  cat("Average silhouette width:", round(x$avg_silhouette, 3), "\n")
  
  if (x$reassignments > 0) {
    cat("Reassignments:", x$reassignments, "\n")
  }
  cat("\n")
  
  cat("Poles:\n")
  for (i in 1:x$k) {
    cat(sprintf("  Cluster %d: %s (n=%d, avg.dist=%.3f)\n", 
                i, x$poles[i], x$cluster_sizes[i], x$within_cluster_dist[i]))
  }
  
  cat("\nClustering vector:\n")
  print(head(x$clustering, 10))
  if (x$n > 10) cat("  ... (", x$n - 10, "more units)\n")
  
  invisible(x)
}

# ========== Plot method ==========
plot.polar_cluster <- function(x, ...) {
  par(mfrow = c(1, 2))
  
  # Plot 1: Cluster sizes
  bp <- barplot(x$cluster_sizes, 
                names.arg = paste0("C", 1:x$k),
                main = paste("Cluster Sizes (", x$assignment_method, ")", sep = ""),
                xlab = "Cluster",
                ylab = "Number of Units",
                col = rainbow(x$k),
                ylim = c(0, max(x$cluster_sizes) * 1.1))
  
  # Add max size reference line if using nearest method and limit is finite
  if (x$assignment_method == "nearest" && is.finite(x$max_cluster_size)) {
    abline(h = x$max_cluster_size, col = "red", lty = 2, lwd = 2)
    text(bp[x$k], x$max_cluster_size, 
         labels = paste("Max =", x$max_cluster_size), 
         pos = 3, col = "red", cex = 0.8)
  }
  
  # Add mean line if using balanced method
  if (x$assignment_method == "balanced") {
    abline(h = mean(x$cluster_sizes), col = "blue", lty = 2, lwd = 2)
    text(bp[x$k], mean(x$cluster_sizes), 
         labels = paste("Mean =", round(mean(x$cluster_sizes), 1)), 
         pos = 3, col = "blue", cex = 0.8)
  }
  
  # Plot 2: Silhouette plot
  colors <- rainbow(x$k)[x$clustering]
  barplot(sort(x$silhouette, decreasing = TRUE),
          main = paste0("Silhouette Plot (avg = ", 
                        round(x$avg_silhouette, 3), ")"),
          xlab = "Units",
          ylab = "Silhouette Width",
          col = colors[order(x$silhouette, decreasing = TRUE)],
          border = NA)
  abline(h = 0, lty = 2)
  abline(h = x$avg_silhouette, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

# ========== Summary method ==========
summary.polar_cluster <- function(object, ...) {
  cat("Polar Ordination Clustering Summary\n")
  cat("====================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Configuration:\n")
  cat("  Number of clusters (k):", object$k, "\n")
  cat("  Assignment method:", object$assignment_method, "\n")
  
  if (object$assignment_method == "nearest") {
    cat("  Max cluster size:", 
        if (is.infinite(object$max_cluster_size)) "Inf" else object$max_cluster_size, 
        "\n")
  }
  
  cat("  First pole method:", object$first_pole_method, "\n")
  cat("  Subsequent pole quantile:", object$subsequent_pole_quantile, "\n")
  cat("  Outlier filter:", object$outlier_filter, "\n")
  
  if (object$outlier_filter) {
    cat("  Outlier threshold (IQR):", object$outlier_threshold, "\n")
    cat("  Outliers detected:", object$n_outliers, "\n")
  }
  
  cat("  Reassignments:", object$reassignments, "\n\n")
  
  cat("Cluster Information:\n")
  cluster_info <- data.frame(
    Pole = object$poles,
    Size = object$cluster_sizes,
    Avg_Dist = round(object$within_cluster_dist, 3)
  )
  rownames(cluster_info) <- paste0("Cluster_", 1:object$k)
  print(cluster_info)
  cat("\n")
  
  cat("Silhouette Statistics:\n")
  cat("  Mean:  ", round(mean(object$silhouette), 3), "\n")
  cat("  Median:", round(median(object$silhouette), 3), "\n")
  cat("  Min:   ", round(min(object$silhouette), 3), "\n")
  cat("  Max:   ", round(max(object$silhouette), 3), "\n")
  
  if (object$n_outliers > 0) {
    cat("\nOutliers detected (", object$n_outliers, "):\n", sep = "")
    cat(" ", paste(object$outliers_detected, collapse = ", "), "\n")
  }
  
  invisible(object)
}


# Comprehensive k selection
select_optimal_k <- function(dist_matrix, k_range = 2:15,
                             first_pole_method = "median",
                             subsequent_pole_quantile = 0.75,
                             assignment_method = "balanced",
                             outlier_filter = TRUE,
                             min_acceptable_size = 5,
                             plot = TRUE) {
  #' Comprehensive selection of optimal k using multiple criteria
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @param min_acceptable_size Minimum acceptable cluster size
  #' @param plot If TRUE, generate diagnostic plots
  #' @return List with recommendations and diagnostic information
  
  cat("Testing k values from", min(k_range), "to", max(k_range), "...\n")
  
  # 1. Silhouette analysis
  cat("1. Running silhouette analysis...\n")
  sil_result <- find_optimal_k(dist_matrix, k_range, 
                               first_pole_method, subsequent_pole_quantile,
                               assignment_method, outlier_filter)
  
  # 2. Elbow method
  cat("2. Calculating within-cluster distances...\n")
  elbow_result <- calculate_wcsd(dist_matrix, k_range,
                                 first_pole_method, subsequent_pole_quantile,
                                 assignment_method, outlier_filter)
  
  # 3. Cluster size diagnostics
  cat("3. Analyzing cluster size distributions...\n")
  diagnostics <- lapply(sil_result$all_results, diagnose_k_quality, min_acceptable_size)
  
  # 4. Find k values with no small clusters
  valid_k <- k_range[sapply(diagnostics, function(d) d$n_small_clusters == 0)]
  
  # 5. Among valid k, find best silhouette
  if (length(valid_k) > 0) {
    valid_indices <- which(k_range %in% valid_k)
    best_valid_idx <- valid_indices[which.max(sil_result$silhouette_scores[valid_indices])]
    recommended_k <- k_range[best_valid_idx]
  } else {
    # If no k has zero small clusters, choose based on quality score
    quality_scores <- sapply(diagnostics, function(d) d$quality_score)
    recommended_k <- k_range[which.max(quality_scores)]
  }
  
  # Generate plots
  if (plot) {
    par(mfrow = c(2, 2))
    
    # Plot 1: Silhouette scores
    plot(k_range, sil_result$silhouette_scores, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Average Silhouette Width",
         main = "Silhouette Analysis", col = "blue", lwd = 2)
    abline(v = sil_result$optimal_k, col = "red", lty = 2)
    abline(h = 0.5, col = "gray", lty = 3)
    text(sil_result$optimal_k, max(sil_result$silhouette_scores), 
         labels = paste("k =", sil_result$optimal_k), pos = 4, col = "red")
    
    # Plot 2: Elbow plot
    plot(k_range, elbow_result$wcsd, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Distances",
         main = "Elbow Method", col = "darkgreen", lwd = 2)
    
    # Plot 3: Number of small clusters
    n_small_vec <- sapply(diagnostics, function(d) d$n_small_clusters)
    plot(k_range, n_small_vec, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = paste("# Clusters < ", min_acceptable_size),
         main = "Small Cluster Count", col = "orange", lwd = 2)
    abline(h = 0, col = "red", lty = 2)
    
    # Plot 4: Quality score
    quality_scores <- sapply(diagnostics, function(d) d$quality_score)
    plot(k_range, quality_scores, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Composite Quality Score",
         main = "Overall Quality", col = "purple", lwd = 2)
    abline(v = recommended_k, col = "red", lty = 2)
    text(recommended_k, max(quality_scores), 
         labels = paste("Recommended k =", recommended_k), pos = 4, col = "red")
    
    par(mfrow = c(1, 1))
  }
  
  # Print summary
  cat("\n=== K SELECTION SUMMARY ===\n")
  cat("Silhouette-optimal k:", sil_result$optimal_k, 
      "(silhouette =", round(sil_result$optimal_silhouette, 3), ")\n")
  cat("Valid k values (no small clusters):", 
      if(length(valid_k) > 0) paste(valid_k, collapse = ", ") else "None", "\n")
  cat("RECOMMENDED k:", recommended_k, "\n")
  
  cat("\nDetailed diagnostics by k:\n")
  for (i in seq_along(k_range)) {
    d <- diagnostics[[i]]
    cat(sprintf("k=%2d: Sil=%.3f, Small=%d/%d, CV=%.2f, Status=%s\n",
                d$k, d$avg_silhouette, d$n_small_clusters, d$k, 
                d$size_cv, d$recommendation))
  }
  
  return(list(
    recommended_k = recommended_k,
    silhouette_optimal_k = sil_result$optimal_k,
    valid_k_values = valid_k,
    silhouette_analysis = sil_result,
    elbow_analysis = elbow_result,
    diagnostics = diagnostics,
    all_results = sil_result$all_results
  ))
}

# Function to find optimal k using silhouette analysis
find_optimal_k <- function(dist_matrix, k_range = 2:10, 
                           first_pole_method = "median",
                           subsequent_pole_quantile = 0.75,
                           assignment_method = "balanced",
                           outlier_filter = TRUE) {
  #' Find optimal number of clusters using silhouette analysis
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test (default: 2:10)
  #' @param ... Other parameters passed to polar_cluster
  #' @return List with silhouette scores and recommended k
  
  results <- list()
  silhouette_scores <- numeric(length(k_range))
  cluster_size_stats <- list()
  
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    
    # Run clustering
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    
    results[[i]] <- result
    silhouette_scores[i] <- result$avg_silhouette
    
    # Track cluster size statistics
    cluster_size_stats[[i]] <- list(
      k = k_val,
      sizes = result$cluster_sizes,
      min_size = min(result$cluster_sizes),
      max_size = max(result$cluster_sizes),
      sd_size = sd(result$cluster_sizes),
      n_small = sum(result$cluster_sizes < 5)  # Number of small clusters
    )
  }
  
  # Find optimal k (maximum silhouette)
  optimal_idx <- which.max(silhouette_scores)
  optimal_k <- k_range[optimal_idx]
  
  return(list(
    k_range = k_range,
    silhouette_scores = silhouette_scores,
    optimal_k = optimal_k,
    optimal_silhouette = silhouette_scores[optimal_idx],
    all_results = results,
    cluster_size_stats = cluster_size_stats
  ))
}

# Function to calculate within-cluster sum of distances
calculate_wcsd <- function(dist_matrix, k_range = 2:10,
                           first_pole_method = "median",
                           subsequent_pole_quantile = 0.75,
                           assignment_method = "balanced",
                           outlier_filter = TRUE) {
  #' Calculate Within-Cluster Sum of Distances for different k
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @return List with WCSD values
  
  wcsd_values <- numeric(length(k_range))
  
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    
    # Calculate total within-cluster sum of distances
    wcsd_values[i] <- sum(result$dist_to_pole)
  }
  
  return(list(
    k_range = k_range,
    wcsd = wcsd_values
  ))
}

# Calculate gap statistic
calculate_gap_statistic <- function(dist_matrix, k_range = 2:10,
                                    B = 50,  # Number of bootstrap samples
                                    first_pole_method = "median",
                                    subsequent_pole_quantile = 0.75,
                                    assignment_method = "balanced",
                                    outlier_filter = TRUE) {
  #' Calculate Gap Statistic for optimal k selection
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @param B Number of bootstrap reference datasets
  #' @return List with gap statistics
  
  n <- nrow(as.matrix(dist_matrix))
  dist_mat <- as.matrix(dist_matrix)
  
  observed_wcsd <- numeric(length(k_range))
  expected_wcsd <- matrix(0, nrow = length(k_range), ncol = B)
  
  # Calculate observed WCSD
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    observed_wcsd[i] <- sum(result$dist_to_pole)
  }
  
  # Generate B reference datasets and calculate expected WCSD
  for (b in 1:B) {
    # Generate random distance matrix with similar properties
    # Using uniform distribution in the range of original distances
    rand_dist <- matrix(runif(n * n, min(dist_mat), max(dist_mat)), n, n)
    rand_dist <- (rand_dist + t(rand_dist)) / 2  # Make symmetric
    diag(rand_dist) <- 0
    rownames(rand_dist) <- colnames(rand_dist) <- rownames(dist_mat)
    
    for (i in seq_along(k_range)) {
      k_val <- k_range[i]
      result <- polar_cluster(
        rand_dist, 
        k = k_val,
        first_pole_method = first_pole_method,
        subsequent_pole_quantile = subsequent_pole_quantile,
        assignment_method = assignment_method,
        outlier_filter = FALSE  # Don't filter outliers in random data
      )
      expected_wcsd[i, b] <- sum(result$dist_to_pole)
    }
  }
  
  # Calculate gap statistic
  gap <- log(rowMeans(expected_wcsd)) - log(observed_wcsd)
  
  # Calculate standard deviation
  sdk <- apply(log(expected_wcsd), 1, sd)
  sk <- sdk * sqrt(1 + 1/B)
  
  # Find optimal k using 1-SE rule
  # Choose smallest k such that Gap(k) >= Gap(k+1) - s_{k+1}
  optimal_k <- k_range[1]
  for (i in 1:(length(k_range) - 1)) {
    if (gap[i] >= gap[i + 1] - sk[i + 1]) {
      optimal_k <- k_range[i]
      break
    }
  }
  
  return(list(
    k_range = k_range,
    gap = gap,
    gap_se = sk,
    optimal_k = optimal_k
  ))
}

# Diagnose cluster quality
diagnose_k_quality <- function(result, min_acceptable_size = 5) {
  #' Diagnose the quality of a clustering result
  #' 
  #' @param result A polar_cluster result object
  #' @param min_acceptable_size Minimum acceptable cluster size
  #' @return List with diagnostic information
  
  sizes <- result$cluster_sizes
  n_small <- sum(sizes < min_acceptable_size)
  prop_small <- n_small / result$k
  
  # Calculate size balance (coefficient of variation)
  cv_size <- sd(sizes) / mean(sizes)
  
  # Identify problematic clusters
  problematic <- which(sizes < min_acceptable_size)
  
  quality_score <- result$avg_silhouette * (1 - prop_small) * (1 / (1 + cv_size))
  
  diagnosis <- list(
    k = result$k,
    n_units = result$n,
    avg_silhouette = result$avg_silhouette,
    n_small_clusters = n_small,
    prop_small_clusters = prop_small,
    size_cv = cv_size,
    problematic_clusters = problematic,
    quality_score = quality_score,
    recommendation = ifelse(n_small > result$k * 0.2, 
                            "TOO MANY SMALL CLUSTERS - Consider reducing k",
                            ifelse(result$avg_silhouette < 0.25,
                                   "WEAK CLUSTERING - Try different k or parameters",
                                   "ACCEPTABLE"))
  )
  
  return(diagnosis)
}