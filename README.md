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


### `benchmark_beta_methods()`

Runs a full benchmarking workflow over multiple beta-diversity distances, evaluates each using your
existing evaluators, computes per-metric Borda scores (ranked **within metric & dataset**), and
optionally makes a tournament plot.

- **Input**
  - `phyloseq` object
  - `methods`: any mix of  
    `Weighted Unifrac`, `Unweighted Unifrac`, `Generalized Unifrac`,  
    `Weighted Glad1`, `Unweighted Glad1`, `Weighted Glad`, `Unweighted Glad`,  
    `manhattan`, `euclidean`, `bray`, `jaccard`, `chisq`, `chord`, `hellinger`, `aitchison`
  - `group`: sample grouping (column name in `sample_data(ps)` or a vector)
  - `rho`: GLAD shrinkage (0 ≤ ρ < 1), `pseudocount` for Aitchison, `figure = TRUE/FALSE`,
    `dataset_name` label
- **Assumes these functions already exist in your environment**:
  - `evaluate_clustering_quality(D, group)`
  - `evaluate_structure_preservation(D, ord_dim)`
- **Returns**
  - `$table`: one row per method with metrics
  - `$plot`: ggplot object (or `NULL` if `figure = FALSE`)

**Requires**: `phyloseq`, `GUniFrac`, `vegan`, `dplyr`, `tidyr`, `ggplot2`, `ggbump`.


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

# Benchmark a panel of distances and draw the tournament plot
res <- benchmark_beta_methods(
  physeq_object,
  group = sample_data(physeq_object)$Group,   # or a vector the same length as samples
  methods = c("Weighted Unifrac","Unweighted Unifrac","Generalized Unifrac",
              "Weighted Glad1","Unweighted Glad1","Weighted Glad","Unweighted Glad",
              "manhattan","euclidean","bray","jaccard","chisq","chord","hellinger","aitchison"),
  rho = 0.5,
  pseudocount = 0.5,
  figure = TRUE,
  dataset_name = "Gut"
)
res$table
print(res$plot)
```

## Installation

```r
# install.packages("devtools")
library(devtools)
devtools::install_github("bioscinema/GLaD")
library(GLaD)
```



