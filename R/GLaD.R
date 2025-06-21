#' Graph Laplacian Distance (GLaD) for Microbiome Samples
#'
#' Computes the GLaD dissimilarity matrix for samples in a [`phyloseq`] object,
#' incorporating phylogenetic relationships through a graph Laplacian. GLaD is
#' a flexible kernel-based approach that propagates tip-level abundances through
#' a phylogenetic tree, yielding distances that reflect both abundance and
#' topology.
#'
#' @param physeq A [`phyloseq`] object with a rooted phylogenetic tree.
#' @param rho A scalar shrinkage parameter in [0, 1); defaults to 0.5.
#' @param weighted Logical; if TRUE, use relative abundance, otherwise
#'   presence/absence.
#'
#' @return A symmetric matrix of class [`dist`] containing the pairwise GLaD
#'   dissimilarities between samples.
#'
#' @details The GLaD distance between samples $i$ and $j$ is computed as
#'   \deqn{ \sqrt{(p_i - p_j)^\top L^{-1} (p_i - p_j)} }
#'   where $p_i$ and $p_j$ are the full (tip + internal node) normalized
#'   abundance vectors and $L$ is the regularized graph Laplacian derived from
#'   the adjacency matrix of the phylogeny. Internal node abundances are
#'   calculated by summing descendants. The linear system is solved efficiently
#'   using sparse matrix techniques.
#'
#' @importFrom phyloseq phy_tree taxa_are_rows otu_table sample_names nsamples
#' @importFrom phangorn Descendants
#' @importFrom ape is.rooted Ntip Nnode
#' @importFrom Matrix Matrix solve diag
#' @importFrom methods as
#' @export

GLaD <- function(physeq, rho = 0.5, weighted = TRUE) {
  tree <- phy_tree(physeq)
  tip_labels <- tree$tip.label
  
  # Fix duplicated tip labels if needed
  if (any(duplicated(tip_labels))) {
    tip_labels_fixed <- make.unique(tip_labels)
    tree$tip.label <- tip_labels_fixed
    phy_tree(physeq) <- tree
    if (taxa_are_rows(physeq)) {
      rownames(otu_table(physeq)) <- tip_labels_fixed
    } else {
      colnames(otu_table(physeq)) <- tip_labels_fixed
    }
  }
  
  if (!is.rooted(tree)) stop("The phylogenetic tree must be rooted.")
  
  # Build adjacency and Laplacian
  A <- build_adjacency_matrix(tree)
  D <- diag(rowSums(A))
  L <- D - rho * A
  L_sparse <- Matrix(L + diag(1e-8, nrow(L)), sparse = TRUE)
  
  # Node labels and IDs
  ntips <- Ntip(tree)
  nnodes <- Nnode(tree)
  tip_names <- tree$tip.label
  internal_ids <- (ntips + 1):(ntips + nnodes)
  internal_names <- paste0("Node_", internal_ids)
  all_node_names <- c(tip_names, internal_names)
  
  # OTU table (samples x tips)
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  full_mat <- matrix(0, nrow = nsamples(physeq), ncol = length(all_node_names))
  rownames(full_mat) <- sample_names(physeq)
  colnames(full_mat) <- all_node_names
  full_mat[, tip_names] <- otu_mat[, tip_names, drop = FALSE]
  
  # Fill internal nodes
  for (i in seq_along(internal_ids)) {
    node_id <- internal_ids[i]
    node_name <- internal_names[i]
    desc <- Descendants(tree, node_id, type = "tips")[[1]]
    desc_names <- tip_names[desc]
    full_mat[, node_name] <- rowSums(otu_mat[, desc_names, drop = FALSE])
  }
  
  rel_abund <- if (weighted) {
    sweep(full_mat, 1, rowSums(full_mat), FUN = "/")
  } else {
    (full_mat > 0) * 1
  }
  
  nsamp <- nrow(rel_abund)
  dist_mat <- matrix(0, nsamp, nsamp)
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(rel_abund)
  
  for (i in 1:(nsamp - 1)) {
    for (j in (i + 1):nsamp) {
      diff <- rel_abund[i, ] - rel_abund[j, ]
      x <- solve(L_sparse, diff)
      dist <- sqrt(sum(diff * x))
      dist_mat[i, j] <- dist
      dist_mat[j, i] <- dist
    }
  }
  
  as.dist(dist_mat)
}

# Internal utility: build adjacency matrix from phylo tree
build_adjacency_matrix <- function(tree) {
  N <- max(tree$edge)
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]
    A[parent, child] <- 1
    A[child, parent] <- 1
  }
  A
}
