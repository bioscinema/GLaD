% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_structure_preservation.R
\name{evaluate_structure_preservation}
\alias{evaluate_structure_preservation}
\title{Evaluate Structure Preservation in Low-Dimensional Embeddings}
\usage{
evaluate_structure_preservation(D, ord_dim = 2)
}
\arguments{
\item{D}{A numeric symmetric dissimilarity matrix or object of class [`dist`].}

\item{ord_dim}{Target number of dimensions for the embedding (default: 2).}
}
\value{
A named list with components:
  * `FNI` – fraction of negative eigenvalues from PCoA; closer to 0 is better.
  * `trustworthiness` – local false neighbor rate from low- to high-dim; in [0, 1], higher is better.
  * `continuity` – local recall rate from high- to low-dim; in [0, 1], higher is better.
  * `mantel_r` – correlation between high- and low-dimensional distances.
  * `mantel_p` – significance of the Mantel test.
  * `stress` – Kruskal stress-1 from raw distances; lower is better.
}
\description{
Quantifies how well a distance or dissimilarity matrix `D` is preserved after
projection into a lower-dimensional Euclidean space using classical
multidimensional scaling (MDS). Includes metrics such as the Fraction of
Negative Eigenvalues (FNI), Trustworthiness, Continuity, Mantel correlation,
and Kruskal stress.
}
\details{
Classical MDS is applied via [`cmdscale`]. Trustworthiness and
continuity are computed by comparing nearest neighbors across the original
and embedded spaces. The number of neighbors `k` adapts to sample size.
}
\seealso{
[`cmdscale`], [`vegan::mantel`]
}
