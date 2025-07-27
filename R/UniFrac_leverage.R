#' UniFrac OTU‑level leverage
#'
#' Computes a branch‑length “leverage” table for every operational taxonomic unit
#' (OTU) in a \code{\link[phyloseq]{phyloseq}} object and reports whether the
#' phylogenetic tree is strictly binary.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object that **must** contain
#'   a rooted phylogenetic tree with branch lengths.  The function queries that
#'   tree via \code{\link[phyloseq]{phy_tree}()}.
#'
#' @return A list with two elements
#' \describe{
#'   \item{\code{leverage}}{A data frame with one row per OTU and four columns:
#'         \code{OTU}, internal \code{parent} node, terminal \code{child} node,
#'         and branch \code{length} (the edge length subtending the OTU).}
#'   \item{\code{binary.tree}}{Logical; \code{TRUE} if the tree is strictly
#'         bifurcating, \code{FALSE} otherwise (polytomies present).}
#' }
#'
#' @details
#' The function filters the tree’s \code{[ape]} matrix to rows whose
#' \code{child} index corresponds to a tip.  These rows are then joined with the
#' vector of branch lengths to form the leverage table.  No UniFrac distance is
#' actually computed here; the routine merely extracts the information needed to
#' gauge each OTU’s contribution (“leverage”) to any weighted or unweighted
#' UniFrac calculation performed downstream.
#'
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data(GlobalPatterns, package = "phyloseq")
#'   UniFrac_leverage(GlobalPatterns)
#' }
#'
#' @import phyloseq
#' @import dplyr
#' @export
UniFrac_leverage <- function(physeq) {
  ### phylogenetic leverage analysis of UniFrac
  tree <- phyloseq::phy_tree(physeq)
  if (is.null(tree$edge.length)) {
    stop("Tree does not contain branch lengths.")
  }
  
  otu_labels <- tree$tip.label
  otu_indices <- seq_along(otu_labels)
  otu_edges <- as.data.frame(tree$edge)
  colnames(otu_edges) <- c("parent", "child")
  otu_edges$length <- tree$edge.length
  
  otu_branch_lengths <- otu_edges %>%
    dplyr::filter(child %in% otu_indices) %>%
    dplyr::mutate(OTU = otu_labels[child]) %>%
    dplyr::select(OTU, parent, child, length)
  
  binary.tree <- ape::is.binary(tree)  # TRUE if no polytomies
  
  list(leverage = otu_branch_lengths, binary.tree = binary.tree)
}
