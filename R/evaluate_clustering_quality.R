#' Evaluate Clustering Quality
#'
#' This function evaluates the quality of a clustering result based on a dissimilarity matrix
#' using several internal validation metrics: silhouette score, Davies-Bouldin Index,
#' Calinski-Harabasz Index, and Dunn Index. It uses k-medoids clustering (`pam`) to perform clustering.
#'
#' @param D A dissimilarity matrix or a `dist` object.
#' @param group A vector of true class labels for samples (used for label alignment).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{avg_silhouette}{Average silhouette width. Higher is better (ideal: > 0.7).}
#'   \item{dbi}{Davies-Bouldin Index. Lower is better (ideal: < 0.5).}
#'   \item{chi}{Calinski-Harabasz Index. Higher is better (ideal: > 500).}
#'   \item{dunn}{Dunn Index. Higher is better (ideal: > 1.5).}
#' }
#' @import cluster fpc clusterCrit clue
#' @export
evaluate_clustering_quality <- function(D, group) {
  if (!inherits(D, "dist")) {
    D <- as.dist(D)
  }
  
  n_clusters <- length(unique(group))
  pam_res <- cluster::pam(D, k = n_clusters, diss = TRUE)
  cluster_labels <- pam_res$clustering
  
  group_numeric <- as.numeric(as.factor(group))
  cluster_labels_matched <- .match_cluster_labels(cluster_labels, group_numeric)
  
  sil <- cluster::silhouette(cluster_labels, dist = D)
  avg_silhouette <- mean(sil[, 3])
  
  coords <- stats::cmdscale(D, k = 2)
  dbi <- .compute_dbi(coords, cluster_labels)
  
  ch_index <- clusterCrit::intCriteria(as.matrix(coords), as.integer(cluster_labels), "Calinski_Harabasz")
  chi <- ch_index$calinski_harabasz
  
  dunn_result <- clusterCrit::intCriteria(as.matrix(coords), as.integer(cluster_labels), "Dunn")
  
  return(list(
    avg_silhouette = avg_silhouette,
    dbi = dbi,
    chi = chi,
    dunn = dunn_result$dunn
  ))
}

#' @keywords internal
.match_cluster_labels <- function(pred_labels, true_labels) {
  contingency <- table(pred_labels, true_labels)
  cost_mat <- max(contingency) - contingency
  assignment <- clue::solve_LSAP(cost_mat)
  matched_labels <- assignment[pred_labels]
  return(matched_labels)
}

#' @keywords internal
.compute_dbi <- function(data, cluster_labels) {
  data <- as.matrix(data)
  crit <- clusterCrit::intCriteria(traj = data, part = as.integer(cluster_labels), crit = "Davies_Bouldin")
  return(crit$davies_bouldin)
}
