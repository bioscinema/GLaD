# GLaD

**Graph Laplacian-Based Distances and Diagnostics for Microbiome Data**

The `GLaD` R package provides a set of tools to compute novel graph Laplacian-based dissimilarities for microbiome compositional data and evaluate the structural fidelity of distance matrices and clustering results. It is designed for integration with phylogenetic information and compositional abundance profiles, especially in the context of microbiome beta-diversity.

---

## Key Functions

### `GLaD()`

Computes pairwise sample dissimilarities using a regularized graph Laplacian constructed from a rooted phylogenetic tree.

- **Input**: A `phyloseq` object with a rooted tree.
- **Distance Formula**:  
  \[
  d(i, j) = \sqrt{(p_i - p_j)^\top L^{-1} (p_i - p_j)}
  \]
- **Options**:
  - `rho`: Shrinkage parameter (default: 0.5)
  - `weighted`: Use relative abundances or presence/absence

---

### `GLaD_eigen()`

A spectral variant of GLaD that projects sample vectors onto the eigenbasis of the Laplacian matrix and computes Euclidean distances in the transformed space.

- **Faster** and more scalable for large trees
- Uses top `k` eigenvectors from the Laplacian spectrum

---

### `evaluate_clustering_quality()`

Quantitatively evaluates clustering results using:

- **Silhouette Width** (`avg_silhouette`)
- **Davies-Bouldin Index** (`dbi`)
- **Calinski-Harabasz Index** (`chi`)
- **Dunn Index** (`dunn`)

All metrics are based on dissimilarity input (e.g., from `GLaD()`).

---

### `evaluate_structure_preservation()`

Assesses how well a distance matrix is preserved under dimensionality reduction (e.g., MDS):

- **FNI**: Fraction of negative eigenvalues
- **Trustworthiness & Continuity**: Local neighborhood preservation
- **Mantel r & p-value**: Correlation with original distances
- **Kruskal Stress**: Distance embedding distortion

---

### `UniFrac_leverage()`

Extracts OTU‑level branch‑length “leverage” values from a
`phyloseq` object and flags whether the underlying tree is strictly
binary.

* **Returns**:  
  * `$leverage` – data frame of *OTU, parent, child, length*  
  * `$binary.tree` – logical (`TRUE` if bifurcating)

---

### `GLaD_leverage()`

Computes Laplacian‑based leverage scores for each OTU under the GLaD
framework.

* **Input**: `phyloseq` object, `rho` (0 ≤ ρ < 1)  
* **Returns**: list (`$leverage`, `$binary.tree`) as above

---

### `build_adjacency_matrix()`

Low‑level helper that turns any `phylo` tree into a symmetric
adjacency matrix. Most users will not call this directly, but it is
exported for advanced work with graph algorithms.

---

## Example Usage

```r
library(GLaD)
library(phyloseq)

# Compute GLaD distance
D_glad <- GLaD(physeq_object, rho = 0.5)

# Compute spectral variant
D_eig <- GLaD_eigen(physeq_object)

# Evaluate clustering result
clust_eval <- evaluate_clustering_quality(D_glad, group = sample_labels)

# Evaluate embedding structure
embed_eval <- evaluate_structure_preservation(D_glad)

# OTU leverage
lev_uni  <- UniFrac_leverage(physeq_object)
lev_glad <- GLaD_leverage(physeq_object, rho = 0.5)
```

## Installation

```r
# install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/GLaD")
library(GLaD)
```



