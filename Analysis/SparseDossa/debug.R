args <- 1
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1L) stop("rep_id must be a positive integer.")

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(SparseDOSSA2)
  library(phyloseq)
  library(GUniFrac)
  library(vegan)
  library(MiSPU)   # for throat.otu.tab / throat.tree data in your original workflow
})

## -------------------- user / HPC settings --------------------
FIT_RDS      <- Sys.getenv("FIT_RDS", "SparseDOSSA2_fit_URT_lambda0.1.rds")
OUTDIR       <- Sys.getenv("OUTDIR", "beta_div_sparse12_out")
BASE_SEED    <- as.integer(Sys.getenv("BASE_SEED", "12345"))
N_SAMPLES    <- as.integer(Sys.getenv("N_SAMPLES", "100"))
ADONIS_PERM  <- as.integer(Sys.getenv("ADONIS_PERM", "999"))
K_COMM_GRID  <- c(2L, 4L, 6L, 8L, 10L, 12L)

if (!file.exists(FIT_RDS)) stop("FIT_RDS not found: ", FIT_RDS)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "raw"), recursive = TRUE, showWarnings = FALSE)

## beta settings: beta = 0 is null / type I setting
beta_grid <- c(0, 0.6, 0.8, 1.2, 1.6, 2.0)
signal_mode_grid <- c("abundance", "prevalence", "both")

## -------------------- helpers --------------------
row_tss <- function(x) {
  x <- as.matrix(x)
  rs <- rowSums(x)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  keep <- is.finite(rs) & rs > 0
  if (any(keep)) {
    out[keep, ] <- x[keep, , drop = FALSE] / rs[keep]
  }
  out
}

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

safe_distance <- function(physeq, method) {
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  out <- tryCatch({
    if (method == "GLaD_0.5_weighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD(physeq)
    } else if (method == "GLaD_0.5_unweighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD(physeq, weighted = FALSE)
    } else if (method == "GLaD_1_weighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD_eigen(physeq)
    } else if (method == "GLaD_1_unweighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD_eigen(physeq, weighted = FALSE)
    } else if (method == "GUniFrac") {
      tree <- phy_tree(physeq)
      unifrac_arr <- GUniFrac::GUniFrac(otu_mat, tree, alpha = c(0, 0.5, 1))$unifracs
      as.dist(unifrac_arr[, , "d_0.5"])
    } else if (method == "aitchison") {
      vegan::vegdist(otu_mat, method = "aitchison", pseudocount = 0.5)
    } else if (method == "unifrac_unweighted") {
      phyloseq::distance(physeq, method = "unifrac", weighted = FALSE)
    } else if (method == "unifrac_weighted_phyloseq") {
      phyloseq::distance(physeq, method = "unifrac", weighted = TRUE)
    } else {
      vegan::vegdist(otu_mat, method = method)
    }
  }, error = function(e) {
    NULL
  })
  
  if (is.matrix(out)) out <- as.dist(out)
  out
}

safe_adonis_p <- function(D, x, permutations = 999) {
  if (is.null(D)) return(NA_real_)
  if (any(!is.finite(as.numeric(D)))) return(NA_real_)
  
  out <- tryCatch({
    tab <- vegan::adonis2(D ~ x, permutations = permutations)
    as.numeric(tab$`Pr(>F)`[1])
  }, error = function(e) {
    NA_real_
  })
  out
}

run_permanova <- function(physeq, distance_methods, permutations = 999, seed = 1L) {
  x <- as.numeric(sample_data(physeq)$x)
  pvals <- setNames(rep(NA_real_, length(distance_methods)), distance_methods)
  
  for (j in seq_along(distance_methods)) {
    method <- distance_methods[j]
    set.seed(seed + j)
    D <- safe_distance(physeq, method)
    pvals[method] <- safe_adonis_p(D, x, permutations = permutations)
  }
  
  pvals
}

make_phyloseq_full <- function(counts_mat, tree, y_vec) {
  counts_mat <- as.matrix(counts_mat)
  stopifnot(nrow(counts_mat) == length(y_vec))
  
  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("sample_", seq_len(nrow(counts_mat)))
  }
  
  phyloseq(
    otu_table(counts_mat, taxa_are_rows = FALSE),
    sample_data(data.frame(
      x = as.numeric(y_vec),
      row.names = rownames(counts_mat)
    )),
    phy_tree(tree)
  )
}

make_metadata_matrix <- function(y_vec) {
  y_vec <- as.integer(y_vec)
  out <- matrix(
    y_vec,
    ncol = 1,
    dimnames = list(paste0("sample_", seq_along(y_vec)), "Y")
  )
  out
}

make_spike_config <- function(feature_ids, signal_mode, effect_size) {
  props <- switch(
    signal_mode,
    abundance = "abundance",
    prevalence = "prevalence",
    both = c("abundance", "prevalence"),
    stop("Unknown signal_mode: ", signal_mode)
  )
  
  do.call(rbind, lapply(props, function(prop) {
    data.frame(
      metadata_datum = 1,
      feature_spiked = feature_ids,
      associated_property = prop,
      effect_size = effect_size,
      stringsAsFactors = FALSE
    )
  }))
}

## -------------------- fixed manuscript-style setup --------------------
data("throat.otu.tab")
data("throat.tree")

counts <- as.matrix(throat.otu.tab)
common <- intersect(colnames(counts), throat.tree$tip.label)
counts <- counts[, common, drop = FALSE]
tree0  <- ape::drop.tip(throat.tree, setdiff(throat.tree$tip.label, colnames(counts)))

keep <- colSums(counts > 0) > 1
counts_f <- counts[, keep, drop = FALSE]
tree0    <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(counts_f)))

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

## -------------------- community distance matrix reused across K --------------------
cop <- ape::cophenetic.phylo(tree0)
cop <- cop[colnames(counts_f), colnames(counts_f)]

distance_methods <- c(
  "GLaD_0.5_weighted",
  "GLaD_0.5_unweighted",
  "GLaD_1_weighted",
  "GLaD_1_unweighted",
  "euclidean",
  "bray",
  "jaccard",
  "aitchison",
  "unifrac_unweighted",
  "unifrac_weighted_phyloseq",
  "GUniFrac"
)

total_settings <- sum(K_COMM_GRID) * length(beta_grid) * length(signal_mode_grid)

pam_list       <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
rank_map_list  <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
abund_sum_list <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
n_otus_list    <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
settings_list  <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))

cat(sprintf(
  "[START] rep_id=%d | total_settings=%d | K grid=%s | N=%d | adonis_perm=%d\n",
  rep_id, total_settings, paste(K_COMM_GRID, collapse = ","), N_SAMPLES, ADONIS_PERM
))

rows_out <- vector("list", total_settings)
row_idx <- 0L

for (K_COMM in K_COMM_GRID) {
  
  ## -------------------- community definition: same style as first script --------------------
#  pam_k <- cluster::pam(as.dist(cop), k = K_COMM, diss = TRUE)$clustering
  
  pc_k <- polar_cluster(dist_matrix = as.matrix(cop), 
                        k = 3, # K_COMM
                        max_cluster_size = Inf,
                        first_pole_method = "median", 
                        subsequent_pole_quantile = 1,
                        assignment_method = "nearest", # try nearest and balanced
                        outlier_filter = TRUE,
                        outlier_threshold = 1.5)$clustering
  
  abund_sum <- tapply(est_abs_mean[names(pam_k)], pam_k, sum, na.rm = TRUE)
  n_otus    <- tapply(pam_k, pam_k, length)
  
  ord_labels <- names(sort(abund_sum, decreasing = TRUE))
  rank_map   <- setNames(seq_along(ord_labels), ord_labels)
  
  ## settings: one signal community at a time, binary Y, no confounder
  settings <- expand.grid(
    cluster_idx = seq_len(K_COMM),
    beta = beta_grid,
    signal_mode = signal_mode_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  pam_list[[as.character(K_COMM)]]       <- pam_k
  rank_map_list[[as.character(K_COMM)]]  <- rank_map
  abund_sum_list[[as.character(K_COMM)]] <- abund_sum
  n_otus_list[[as.character(K_COMM)]]    <- n_otus
  settings_list[[as.character(K_COMM)]]  <- settings
  
  cat(sprintf("[K=%d] settings=%d\n", K_COMM, nrow(settings)))
  
  for (s in seq_len(nrow(settings))) {
    row_idx <- row_idx + 1L
    job <- settings[s, , drop = FALSE]
    
    cluster_idx_chr <- as.character(job$cluster_idx)
    C_otus <- names(pam_k)[pam_k == job$cluster_idx]
    
    if (!length(C_otus)) {
      rows_out[[row_idx]] <- data.frame(
        rep_id = rep_id,
        setting_id = s,
        k = K_COMM,
        signal_cluster = job$cluster_idx,
        cluster_ordered = NA_integer_,
        cluster_abund_sum_template = NA_real_,
        cluster_n_otus_template = NA_integer_,
        signal_mode = as.character(job$signal_mode),
        beta = job$beta,
        n_samples = NA_integer_,
        y_mean = NA_real_,
        seed_sim = NA_integer_,
        seed_test = NA_integer_,
        runtime_sec = NA_real_,
        matrix(
          NA_real_,
          nrow = 1,
          ncol = length(distance_methods),
          dimnames = list(NULL, paste0("p_", distance_methods))
        ),
        stringsAsFactors = FALSE
      )
      next
    }
    
    cluster_ordered <- as.integer(rank_map[cluster_idx_chr])
    cluster_abund_sum_template <- as.numeric(abund_sum[cluster_idx_chr])
    cluster_n_otus_template <- as.integer(n_otus[cluster_idx_chr])
    
    seed_sim  <- BASE_SEED + (rep_id - 1L) * 10000L + row_idx
    seed_test <- BASE_SEED + 2000000L + (rep_id - 1L) * 10000L + row_idx
    
    cat(sprintf(
      "[rep %03d | %3d/%3d] K=%02d cluster=%02d signal=%s beta=%.2f\n",
      rep_id, row_idx, total_settings, K_COMM, job$cluster_idx,
      as.character(job$signal_mode), job$beta
    ))
    
    ## ---------- create balanced-binomial Y first, then spike signal during simulation ----------
    set.seed(seed_sim + 500000L)
    Y <- rbinom(N_SAMPLES, 1, 0.5)
    metadata_matrix <- make_metadata_matrix(Y)
    spike_config <- make_spike_config(
      feature_ids = C_otus,
      signal_mode = as.character(job$signal_mode),
      effect_size = job$beta
    )
    
    set.seed(seed_sim)
    sim <- SparseDOSSA2::SparseDOSSA2(
      template = fit,
      n_sample = N_SAMPLES,
      new_features = FALSE,
      spike_metadata = spike_config,
      metadata_matrix = metadata_matrix,
      verbose = FALSE
    )
    
    A_abs      <- t(sim$simulated_matrices$a_spiked)
    counts_sim <- t(sim$simulated_data)
    
    tree <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(A_abs)))
    A_abs      <- A_abs[, tree$tip.label, drop = FALSE]
    counts_sim <- counts_sim[, tree$tip.label, drop = FALSE]
    
    libsize <- rowSums(counts_sim)
    keep_samp <- is.finite(libsize) & libsize > 0
    if (!all(keep_samp)) {
      A_abs      <- A_abs[keep_samp, , drop = FALSE]
      counts_sim <- counts_sim[keep_samp, , drop = FALSE]
      Y          <- Y[keep_samp]
    }
    
    counts_sim_tss_full <- row_tss(counts_sim)
    
    C_curr <- intersect(C_otus, colnames(A_abs))
    if (!length(C_curr)) {
      rows_out[[row_idx]] <- data.frame(
        rep_id = rep_id,
        setting_id = s,
        k = K_COMM,
        signal_cluster = job$cluster_idx,
        cluster_ordered = cluster_ordered,
        cluster_abund_sum_template = cluster_abund_sum_template,
        cluster_n_otus_template = cluster_n_otus_template,
        signal_mode = as.character(job$signal_mode),
        beta = job$beta,
        n_samples = 0L,
        y_mean = NA_real_,
        seed_sim = seed_sim,
        seed_test = seed_test,
        runtime_sec = NA_real_,
        matrix(
          NA_real_,
          nrow = 1,
          ncol = length(distance_methods),
          dimnames = list(NULL, paste0("p_", distance_methods))
        ),
        stringsAsFactors = FALSE
      )
      next
    }
    
    ## ---------- whole OTU table beta-diversity tests ----------
    physeq_full <- make_phyloseq_full(counts_mat = counts_sim, tree = tree, y_vec = Y)
    
    S <- rowSums(A_abs[, C_curr, drop = FALSE])
    
    t0 <- proc.time()[["elapsed"]]
    pvals <- run_permanova(
      physeq = physeq_full,
      distance_methods = distance_methods,
      permutations = ADONIS_PERM,
      seed = seed_test
    )
    t1 <- proc.time()[["elapsed"]]
    
    rows_out[[row_idx]] <- data.frame(
      rep_id = rep_id,
      setting_id = s,
      k = K_COMM,
      signal_cluster = job$cluster_idx,
      cluster_ordered = cluster_ordered,
      cluster_abund_sum_template = cluster_abund_sum_template,
      cluster_n_otus_template = cluster_n_otus_template,
      signal_mode = as.character(job$signal_mode),
      beta = job$beta,
      is_null = as.integer(job$beta == 0),
      n_samples = nrow(counts_sim),
      y_mean = mean(Y),
      signal_mean_abs = mean(S, na.rm = TRUE),
      signal_sd_abs = sd(S, na.rm = TRUE),
      seed_sim = seed_sim,
      seed_test = seed_test,
      runtime_sec = (t1 - t0),
      as.data.frame(as.list(setNames(as.numeric(pvals), paste0("p_", names(pvals))))),
      stringsAsFactors = FALSE
    )
  }
}

res_df <- do.call(rbind, rows_out[seq_len(row_idx)])

outfile_csv <- file.path(OUTDIR, "raw", sprintf("beta_div_sparse12_rep%03d.csv", rep_id))
outfile_rds <- file.path(OUTDIR, "raw", sprintf("beta_div_sparse12_rep%03d.rds", rep_id))

write.csv(res_df, outfile_csv, row.names = FALSE)
saveRDS(
  list(
    results = res_df,
    settings = settings_list[["12"]],
    settings_by_k = settings_list,
    pam12 = pam_list[["12"]],
    pam_by_k = pam_list,
    rank_map = rank_map_list[["12"]],
    rank_map_by_k = rank_map_list,
    abund_sum = abund_sum_list[["12"]],
    abund_sum_by_k = abund_sum_list,
    n_otus = n_otus_list[["12"]],
    n_otus_by_k = n_otus_list,
    K_COMM_grid = K_COMM_GRID,
    distance_methods = distance_methods,
    signal_mode_grid = signal_mode_grid,
    rep_id = rep_id,
    n_samples = N_SAMPLES,
    adonis_perm = ADONIS_PERM
  ),
  outfile_rds
)

cat(sprintf("[DONE] rep_id=%d | wrote:\n  %s\n  %s\n", rep_id, outfile_csv, outfile_rds))
