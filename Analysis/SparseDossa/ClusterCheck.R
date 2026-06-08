suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(SparseDOSSA2)
  library(MiSPU)   # throat.otu.tab / throat.tree
})

## -------------------- user settings --------------------
FIT_RDS <- Sys.getenv("FIT_RDS", "SparseDOSSA2_fit_URT_lambda0.1.rds")
OUTDIR  <- Sys.getenv("OUTDIR", "pam_cluster_summary_output")

K_GRID <- c(2L, 3L, 4L, 8L)

if (!file.exists(FIT_RDS)) stop("FIT_RDS not found: ", FIT_RDS)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## -------------------- helpers --------------------
get_est_mean_abs <- function(fit_obj) {
  f <- fit_obj$EM_fit$fit
  if (is.null(f)) f <- fit_obj$fit
  if (is.null(f)) stop("Cannot find EM_fit$fit or fit in FIT_RDS object.")
  
  pi0 <- f$pi0
  mu  <- f$mu
  
  if (!is.null(f$sigma)) {
    sigma <- f$sigma
  } else if (!is.null(f$sigma2)) {
    sigma <- sqrt(f$sigma2)
  } else {
    stop("Cannot find sigma or sigma2 in fit object.")
  }
  
  if (is.null(names(mu)) || is.null(names(pi0)) || is.null(names(sigma))) {
    stop("mu/pi0/sigma must be named by OTU IDs.")
  }
  
  common_ids <- Reduce(intersect, list(names(mu), names(pi0), names(sigma)))
  if (!length(common_ids)) stop("No overlapping OTU IDs across mu/pi0/sigma.")
  
  mu    <- mu[common_ids]
  pi0   <- pi0[common_ids]
  sigma <- sigma[common_ids]
  
  (1 - pi0) * exp(mu + 0.5 * sigma^2)
}

upper_tri_values <- function(mat) {
  if (nrow(mat) <= 1L) return(numeric(0))
  mat[upper.tri(mat, diag = FALSE)]
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

summarize_within_cop <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(data.frame(
      within_cop_n_pair = 0L,
      within_cop_mean = NA_real_,
      within_cop_sd = NA_real_,
      within_cop_min = NA_real_,
      within_cop_q1 = NA_real_,
      within_cop_median = NA_real_,
      within_cop_q3 = NA_real_,
      within_cop_max = NA_real_,
      within_cop_iqr = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    within_cop_n_pair = length(x),
    within_cop_mean = mean(x),
    within_cop_sd = if (length(x) > 1L) sd(x) else NA_real_,
    within_cop_min = min(x),
    within_cop_q1 = safe_quantile(x, 0.25),
    within_cop_median = safe_quantile(x, 0.50),
    within_cop_q3 = safe_quantile(x, 0.75),
    within_cop_max = max(x),
    within_cop_iqr = IQR(x),
    stringsAsFactors = FALSE
  )
}

## -------------------- load template data --------------------
data("throat.otu.tab")
data("throat.tree")

counts <- as.matrix(throat.otu.tab)

common <- intersect(colnames(counts), throat.tree$tip.label)
counts <- counts[, common, drop = FALSE]

tree0 <- ape::drop.tip(
  throat.tree,
  setdiff(throat.tree$tip.label, colnames(counts))
)

## same OTU filtering as your pasted script
keep <- colSums(counts > 0) > 1
counts_f <- counts[, keep, drop = FALSE]

tree0 <- ape::drop.tip(
  tree0,
  setdiff(tree0$tip.label, colnames(counts_f))
)

## make count table and tree exactly aligned
counts_f <- counts_f[, tree0$tip.label, drop = FALSE]

## -------------------- get fitted mean absolute abundance --------------------
fit <- readRDS(FIT_RDS)

est_abs_mean <- get_est_mean_abs(fit)

missing_otus <- setdiff(colnames(counts_f), names(est_abs_mean))
if (length(missing_otus) > 0) {
  stop(
    "These OTUs are in counts_f but missing from est_abs_mean: ",
    paste(head(missing_otus, 20), collapse = ", ")
  )
}

est_abs_mean <- est_abs_mean[colnames(counts_f)]

## -------------------- COP distance matrix --------------------
cop <- ape::cophenetic.phylo(tree0)
cop <- cop[colnames(counts_f), colnames(counts_f)]

stopifnot(
  all(rownames(cop) == colnames(counts_f)),
  all(colnames(cop) == colnames(counts_f)),
  all(names(est_abs_mean) == colnames(counts_f))
)

## -------------------- PAM cluster summary over K --------------------
cluster_summary_by_k <- vector("list", length(K_GRID))
names(cluster_summary_by_k) <- paste0("K", K_GRID)

pam_by_k <- vector("list", length(K_GRID))
names(pam_by_k) <- paste0("K", K_GRID)

for (K_COMM in K_GRID) {
  
  cat(sprintf("[PAM] Running K = %d\n", K_COMM))
  
  pam_fit <- cluster::pam(
    x = as.dist(cop),
    k = K_COMM,
    diss = TRUE
  )
  
  pam_k <- pam_fit$clustering
  
  if (is.null(names(pam_k))) {
    names(pam_k) <- rownames(cop)
  }
  
  pam_k <- pam_k[rownames(cop)]
  
  stopifnot(
    length(pam_k) == nrow(cop),
    all(names(pam_k) == rownames(cop)),
    length(unique(pam_k)) == K_COMM,
    !anyNA(pam_k),
    all(table(pam_k) > 0)
  )
  
  pam_by_k[[paste0("K", K_COMM)]] <- pam_k
  
  ## cluster-level overall abundance and number of OTUs
  abund_sum <- tapply(est_abs_mean[names(pam_k)], pam_k, sum, na.rm = TRUE)
  n_otus    <- tapply(pam_k, pam_k, length)
  
  ## rank clusters by overall abundance
  ord_labels <- names(sort(abund_sum, decreasing = TRUE))
  rank_map   <- setNames(seq_along(ord_labels), ord_labels)
  
  cluster_summary_k <- vector("list", K_COMM)
  
  for (g in seq_len(K_COMM)) {
    
    otus_g <- names(pam_k)[pam_k == g]
    
    cop_g <- cop[otus_g, otus_g, drop = FALSE]
    within_cop_vals <- upper_tri_values(cop_g)
    within_cop_sum <- summarize_within_cop(within_cop_vals)
    
    cluster_summary_k[[g]] <- cbind(
      data.frame(
        k = K_COMM,
        cluster = g,
        cluster_ordered_by_abundance = as.integer(rank_map[as.character(g)]),
        n_otus = as.integer(n_otus[as.character(g)]),
        overall_abundance = as.numeric(abund_sum[as.character(g)]),
        overall_abundance_fraction = as.numeric(abund_sum[as.character(g)]) / sum(est_abs_mean),
        medoid_otu = as.character(pam_fit$medoids[g]),
        stringsAsFactors = FALSE
      ),
      within_cop_sum
    )
  }
  
  cluster_summary_k <- do.call(rbind, cluster_summary_k)
  rownames(cluster_summary_k) <- NULL
  
  cluster_summary_by_k[[paste0("K", K_COMM)]] <- cluster_summary_k
  
  out_k <- file.path(OUTDIR, sprintf("pam_cluster_summary_K%d.csv", K_COMM))
  write.csv(cluster_summary_k, out_k, row.names = FALSE)
  
  cat(sprintf("[PAM] Wrote K = %d summary to: %s\n", K_COMM, out_k))
}

cluster_summary_all <- do.call(rbind, cluster_summary_by_k)
rownames(cluster_summary_all) <- NULL

outfile_csv <- file.path(OUTDIR, "pam_cluster_summary_K2_K3_K4_K8.csv")
outfile_rds <- file.path(OUTDIR, "pam_cluster_summary_K2_K3_K4_K8.rds")

write.csv(cluster_summary_all, outfile_csv, row.names = FALSE)

saveRDS(
  list(
    cluster_summary_all = cluster_summary_all,
    cluster_summary_by_k = cluster_summary_by_k,
    pam_by_k = pam_by_k,
    K_GRID = K_GRID,
    cop = cop,
    est_abs_mean = est_abs_mean
  ),
  outfile_rds
)

cat("\n[DONE] PAM cluster summaries written to:\n")
cat("  ", outfile_csv, "\n")
cat("  ", outfile_rds, "\n\n")

print(cluster_summary_all)