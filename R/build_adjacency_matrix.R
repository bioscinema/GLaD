#' Build an adjacency matrix for a phylogenetic tree
#'
#' Generates a symmetric \(N \times N\) adjacency matrix that encodes the
#' parent–child relationships of a phylogenetic tree stored as a
#' \code{\link[ape]{phylo}} object.  Tip and internal node labels (when
#' available) are used as row and column names.
#'
#' @param tree A rooted or unrooted tree of class \code{"phylo"} (from the
#'   \strong{ape} package).  The component \code{$edge} is assumed to index tips
#'   in \eqn{1,\dots,N_\text{tip}} and internal nodes in
#'   \eqn{N_\text{tip}+1,\dots,N_\text{tip}+N_\text{node}}.
#'
#' @return A square numeric matrix whose \eqn{(i,j)} entry equals 1 when nodes
#'   \(i\) and \(j\) are directly connected by an edge and 0 otherwise.
#'
#' @details
#' The routine loops over each edge of the tree and populates both \eqn{A[i,j]}
#' and \eqn{A[j,i]} so that the result is symmetric.  Internal nodes receive
#' automatic labels when \code{tree$node.label} is \code{NULL}.
#'
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   set.seed(1)
#'   toy_tree <- ape::rtree(5)          # random 5‑tip tree
#'   A <- build_adjacency_matrix(toy_tree)
#'   A
#' }
#'
#' @importFrom ape rtree
#' @export
build_adjacency_matrix <- function(tree) {
  
  Ntip  <- length(tree$tip.label)
  Nnode <- tree$Nnode
  N     <- Ntip + Nnode
  
  # Construct row/column names: tips first, then internal nodes
  node_names <- c(
    tree$tip.label,
    if (!is.null(tree$node.label))
      tree$node.label
    else
      as.character((Ntip + 1):N)
  )
  
  # Empty adjacency matrix
  A <- matrix(0, nrow = N, ncol = N,
              dimnames = list(node_names, node_names))
  
  # Populate symmetric adjacency matrix
  for (k in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[k, 1]
    child  <- tree$edge[k, 2]
    A[parent, child] <- 1
    A[child,  parent] <- 1
  }
  
  A
}
