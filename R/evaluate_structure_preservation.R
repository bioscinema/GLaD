#' Evaluate Structure Preservation in Low-Dimensional Embeddings
#'
#' Quantifies how well a distance or dissimilarity matrix `D` is preserved after
#' projection into a lower-dimensional Euclidean space using classical
#' multidimensional scaling (MDS). Includes metrics such as the Fraction of
#' Negative Eigenvalues (FNI), Trustworthiness, Continuity, Mantel correlation,
#' and Kruskal stress.
#'
#' @param D A numeric symmetric dissimilarity matrix or object of class [`dist`].
#' @param ord_dim Target number of dimensions for the embedding (default: 2).
#'
#' @return A named list with components:
#'   * `FNI` – fraction of negative eigenvalues from PCoA; closer to 0 is better.
#'   * `trustworthiness` – local false neighbor rate from low- to high-dim; in [0, 1], higher is better.
#'   * `continuity` – local recall rate from high- to low-dim; in [0, 1], higher is better.
#'   * `mantel_r` – correlation between high- and low-dimensional distances.
#'   * `mantel_p` – significance of the Mantel test.
#'   * `stress` – Kruskal stress-1 from raw distances; lower is better.
#'
#' @details Classical MDS is applied via [`cmdscale`]. Trustworthiness and
#' continuity are computed by comparing nearest neighbors across the original
#' and embedded spaces. The number of neighbors `k` adapts to sample size.
#'
#' @seealso [`cmdscale`], [`vegan::mantel`]
#'
#' @importFrom vegan mantel
#' @importFrom stats dist
#' @importFrom Matrix diag
#' @export

evaluate_structure_preservation <- function(D, ord_dim = 2) {
  D <- as.matrix(D)
  n <- nrow(D)
  
  # 1. FNI (Fraction of Negative Eigenvalues from Gram matrix)
  J <- diag(n) - matrix(1, n, n) / n
  Gram <- -0.5 * J %*% (D^2) %*% J
  eig <- eigen(Gram)
  lambda <- eig$values
  lambda[abs(lambda) < 1e-12] <- 0
  FNI <- sum(abs(lambda[lambda < 0])) / sum(abs(lambda))
  
  # 2. Trustworthiness and Continuity ------------------------------------------
  coords <- cmdscale(D, k = ord_dim)
  D_low  <- as.matrix(dist(coords))
  low_rank  <- apply(D_low, 1, order)[-1, ]   # Exclude self
  high_rank <- apply(D, 1, order)[-1, ]
  
  if (n <= 50) {
    k <- 5
  } else if (n <= 500) {
    k <- ceiling(n / 10)
  } else {
    k <- 50
  }
  
  trust <- 0
  cont  <- 0
  for (i in 1:n) {
    u <- low_rank[1:k, i]
    v <- high_rank[1:k, i]
    ranks_high <- match(u, high_rank[, i])
    ranks_low  <- match(v, low_rank[, i])
    trust <- trust + sum((ranks_high[ranks_high > k] - k))
    cont  <- cont  + sum((ranks_low[ranks_low > k] - k))
  }
  
  norm <- n * k * (2 * n - 3 * k - 1) / 2
  trustworthiness <- 1 - (2 / norm) * trust
  continuity      <- 1 - (2 / norm) * cont
  
  # 3. Mantel test -------------------------------------------------------------
  mantel_out <- vegan::mantel(as.dist(D), dist(coords), permutations = 999)
  mantel_r   <- mantel_out$statistic
  mantel_p   <- mantel_out$signif
  
  # 4. Stress ------------------------------------------------------------------
  D_high_vec <- as.vector(D)
  D_low_vec  <- as.vector(D_low)
  stress     <- sqrt(sum((D_high_vec - D_low_vec)^2) / sum(D_high_vec^2))
  
  list(
    FNI             = FNI,
    trustworthiness = trustworthiness,
    continuity      = continuity,
    mantel_r        = mantel_r,
    mantel_p        = mantel_p,
    stress          = stress
  )
}
