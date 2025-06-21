#' GLaD with Eigen-Decomposition (GLaD_eigen)
#'
#' Computes a GLaD-style phylogenetic distance matrix using eigen-decomposition
#' of the graph Laplacian. Compared to `GLaD()`, this variant uses spectral
#' embedding and avoids solving linear systems directly.
#'
#' @param physeq A [`phyloseq`] object containing a rooted phylogenetic tree.
#' @param rho Shrinkage parameter for Laplacian matrix; default is 1.
#' @param weighted Logical; if TRUE, use relative abundances, otherwise binary presence/absence.
#'
#' @return A symmetric [`dist`] object with pairwise distances between samples.
#'
#' @details The method forms the Laplacian matrix `L = D - \rho A`, then performs
#' an eigen-decomposition to obtain its eigenvectors and eigenvalues. Each sample's
#' abundance vector is projected into the eigenbasis and scaled by the inverse
#' square root of the eigenvalues, and Euclidean distance is computed in this space.
#'
#'
#' @importFrom phyloseq phy_tree taxa_are_rows otu_table sample_names nsamples
#' @importFrom phangorn Descendants
#' @importFrom ape is.rooted Ntip Nnode
#' @importFrom RSpectra eigs_sym
#' @importFrom Matrix Diagonal as
#' @export

GLaD_eigen <- function(physeq, rho = 1, weighted = TRUE) {
  tree <- phy_tree(physeq)
  tip_labels <- tree$tip.label
  
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
  
  A <- build_adjacency_matrix(tree)
  D <- Matrix::Diagonal(x = rowSums(A))
  L <- D - rho * A
  L <- L + Matrix::Diagonal(n = nrow(L), x = 1e-8)
  L_sparse <- as(L, "dgCMatrix")
  
  k <- min(100, nrow(L) - 1)
  eigs_res <- RSpectra::eigs_sym(L_sparse, k = k, which = "SM")
  keep <- which(abs(eigs_res$values) > 1e-8)
  lambda <- eigs_res$values[keep]
  U <- eigs_res$vectors[, keep]
  
  ntips <- Ntip(tree)
  nnodes <- Nnode(tree)
  tip_names <- tree$tip.label
  internal_ids <- (ntips + 1):(ntips + nnodes)
  internal_names <- paste0("Node_", internal_ids)
  all_node_names <- c(tip_names, internal_names)
  
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  full_mat <- matrix(0, nrow = nsamples(physeq), ncol = length(all_node_names))
  rownames(full_mat) <- sample_names(physeq)
  colnames(full_mat) <- all_node_names
  full_mat[, tip_names] <- otu_mat[, tip_names, drop = FALSE]
  
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
  
  proj <- rel_abund %*% U
  scaled_proj <- sweep(proj, 2, sqrt(lambda), FUN = "/")
  
  dist(scaled_proj, method = "euclidean")
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
