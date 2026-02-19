# ================== SimuDA_FixedTree: DA-only simulator using a given tree ==================
# - Uses a user-provided phylogenetic tree (phylo). No random/NJ tree construction.
# - Taxa are taken from tree$tip.label (or a user-chosen subset via taxa_use).
# - Supports x ∈ {binary, continuous, ordinal} to drive differential abundance.
# - Signal selection: sparsity (da_p) + tree concentration (da_conc).
# - Output: phyloseq with otu_table, sample_data (includes x), tax_table(DA_sign, DA_coef),
#           and the pruned tree.
# ===========================================================================================

suppressPackageStartupMessages({
  library(phyloseq)
})

# ---------- Safe Dirichlet sampler ----------
.rdirichlet_safe <- function(n, alpha){
  k <- length(alpha)
  out <- matrix(NA_real_, nrow = n, ncol = k)
  for(i in seq_len(n)){
    g <- rgamma(k, shape = alpha, rate = 1)
    s <- sum(g)
    out[i,] <- if (s == 0 || any(!is.finite(g))) rep(1/k, k) else g/s
  }
  out
}

# ---------- Robust local sampler (guarantees m positives) ----------
.local_select_indices <- function(K, m, dist_to_anchor, conc){
  if (m <= 0) return(integer(0))
  if (is.null(dist_to_anchor) || all(!is.finite(dist_to_anchor))) {
    return(sample.int(K, size = m, replace = FALSE))
  }
  d <- as.numeric(dist_to_anchor)
  tau <- max(quantile(d[d > 0], probs = 0.2, na.rm = TRUE), .Machine$double.eps)
  bw  <- max(tau * (1 - conc + 1e-3), .Machine$double.eps)
  w_local <- exp(-d / bw)
  w <- (1 - conc) * 1 + conc * w_local
  
  pos <- which(w > 0 & is.finite(w))
  if (length(pos) < m) {
    ord <- order(d, decreasing = FALSE)
    return(ord[seq_len(m)])
  }
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (s <= 0) {
    return(sample.int(K, size = m, replace = FALSE))
  } else {
    w <- w / s
    return(sample.int(K, size = m, replace = FALSE, prob = w))
  }
}

# ---------- Strict subtree selection (conc≈1); augments if too small ----------
.strict_subtree_pick <- function(tree, m){
  if (!inherits(tree, "phylo")) return(NULL)
  if (!requireNamespace("phangorn", quietly = TRUE)) return(NULL)
  K <- ape::Ntip(tree); if (m <= 0 || m > K) return(NULL)
  desc <- phangorn::Descendants(tree, node = (K+1):(K+tree$Nnode), type = "tips")
  sizes <- vapply(desc, length, 1L)
  j <- which.min(abs(sizes - m))
  if (length(j)) {
    tips <- unique(desc[[j]])
    tips <- tips[tips >= 1 & tips <= K]
    if (length(tips)) return(tips)
  }
  NULL
}

# ---------- Two-pole DA signal selection (positive & negative blocks) ----------
.pick_signal_sets_2pole <- function(tree, m_up, m_down, conc){
  if (!inherits(tree, "phylo")) stop("tree must be a phylo object.")
  if (!requireNamespace("ape", quietly = TRUE)) stop("Please install 'ape' for cophenetic distances.")
  
  K <- ape::Ntip(tree)
  tips <- seq_len(K)
  
  m_up   <- max(0L, min(as.integer(m_up),   K))
  m_down <- max(0L, min(as.integer(m_down), K - m_up))
  if (m_up + m_down == 0L) {
    return(list(up = integer(0), down = integer(0), anchor_up = NA_integer_, anchor_down = NA_integer_))
  }
  
  # cophenetic distance matrix among tips
  D <- ape::cophenetic.phylo(tree)
  
  # ---- Positive anchor & concentrated block ----
  anchor_up <- sample.int(K, 1)
  # order by distance to anchor_up (closest first), include anchor_up itself
  ord_up <- order(D[, anchor_up], decreasing = FALSE)
  
  n_up_conc <- min(as.integer(round(conc * m_up)), m_up)
  up_conc <- head(ord_up, n_up_conc)
  
  remain <- setdiff(tips, up_conc)
  n_up_rand <- m_up - length(up_conc)
  up_rand <- if (n_up_rand > 0L) sample(remain, n_up_rand, replace = FALSE) else integer(0)
  up_set <- c(up_conc, up_rand)
  
  # ---- Negative anchor (farthest from positive anchor) & concentrated block ----
  anchor_down <- NA_integer_
  down_set <- integer(0)
  if (m_down > 0L) {
    remain <- setdiff(tips, up_set)
    if (length(remain) == 0L) {
      # no space left; should not happen due to cap above
      return(list(up = up_set, down = integer(0), anchor_up = anchor_up, anchor_down = NA_integer_))
    }
    # farthest tip from anchor_up among remaining
    anchor_down <- remain[ which.max(D[remain, anchor_up]) ]
    
    # order remaining tips by distance to anchor_down (closest first)
    ord_down <- remain[ order(D[remain, anchor_down], decreasing = FALSE) ]
    
    n_down_conc <- min(as.integer(round(conc * m_down)), m_down)
    down_conc <- head(ord_down, n_down_conc)
    
    remain2 <- setdiff(remain, down_conc)
    n_down_rand <- m_down - length(down_conc)
    down_rand <- if (n_down_rand > 0L) sample(remain2, n_down_rand, replace = FALSE) else integer(0)
    
    down_set <- c(down_conc, down_rand)
  }
  
  list(up = up_set, down = down_set, anchor_up = anchor_up, anchor_down = anchor_down)
}


# ========================== Main: DA-only using input tree ==========================
SimuDA_FixedTree <- function(
    # required tree & optional taxa subset
  tree,                       # REQUIRED: phylo with tip.label
  taxa_use = NULL,            # optional character vector subset of tree$tip.label
  
  # cohort size / counts
  n,                          # total samples
  N_pts = 10000,              # per-sample depth (scalar or length n)
  
  # cohort-level background composition
  base_alpha = 1.0,           # Dirichlet center for cohort composition q (scalar or length K)
  sample_disp = 200,          # sample-level Dirichlet concentration for p_i around q
  dm_theta = Inf,             # Dirichlet–Multinomial concentration (Inf => Multinomial)
  
  # predictor x (binary / continuous / ordinal)
  x = NULL,                   # vector length n; if NULL, generated per x_type
  x_type = c("binary","continuous","ordinal"),
  p1 = 0.5,                   # if x is NULL & binary: P(x=1)
  cont_mean = 0, cont_sd = 1, # if x is NULL & continuous: Normal(mean, sd)
  ord_levels = 3, ord_probs = NULL, # if x is NULL & ordinal: levels 0..L-1
  
  # DA controls
  da_p = 0.10,                # sparsity: proportion in (0,1) or integer >=1
  da_conc = 0.0,              # tree concentration in [0,1]
  da_logFC = 1.0,             # effect size (see below)
  da_up_frac = 0.5,           # fraction of signals up (rest down)
  cont_standardize = TRUE,    # standardize continuous x before applying slope
  
  # optional filtering
  min_prev = 0,               # drop low-prevalence taxa (0=off; <1 proportion)
  min_var  = 0,               # drop near-constant taxa (0=off)
  
  # reproducibility
  seed = 1
){
  set.seed(seed)
  
  # ---- 0) Validate and set taxa from the provided tree ----
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  if (is.null(tree$tip.label) || anyNA(tree$tip.label))
    stop("The input tree must have non-NA tip labels.")
  
  all_tips <- as.character(tree$tip.label)
  if (!is.null(taxa_use)) {
    taxa_use <- as.character(taxa_use)
    bad <- setdiff(taxa_use, all_tips)
    if (length(bad)) stop("taxa_use contains names not in tree tips: ", paste(bad, collapse = ", "))
    # prune tree to taxa_use
    if (!requireNamespace("ape", quietly = TRUE)) stop("Please install 'ape' for tree pruning.")
    tree <- ape::drop.tip(tree, setdiff(all_tips, taxa_use))
    tax_names <- taxa_use
  } else {
    tax_names <- all_tips
  }
  K <- length(tax_names)
  if (K < 2) stop("Tree (after pruning) must have at least 2 tips/taxa.")
  # (re-)order tip labels to match tax_names
  if (!identical(tree$tip.label, tax_names)) {
    # reorder tree tip order to tax_names if needed
    perm <- match(tax_names, tree$tip.label)
    if (anyNA(perm)) stop("Internal mismatch: failed to match tax_names to tree$tip.label.")
    # ape doesn't provide a trivial reorder by tip; for simulation, matching names is enough
    tree$tip.label <- tax_names
  }
  
  # ---- 1) Resolve/generate x and n ----
  if (missing(n) || is.null(n) || n < 1) stop("'n' must be a positive integer.")
  if (!is.null(x)) {
    if (length(x) != n) stop("Length of 'x' must be equal to 'n'.")
    # infer x_type
    if (is.factor(x)) {
      x_type_res <- if (is.ordered(x)) "ordinal" else if (nlevels(x)==2) "binary" else "ordinal"
    } else if (is.numeric(x)) {
      ux <- sort(unique(x))
      x_type_res <- if (all(ux %in% c(0,1)) && length(ux)<=2) "binary" else "continuous"
    } else stop("Unsupported 'x' type; use numeric/factor/ordered, or leave NULL.")
  } else {
    x_type_res <- match.arg(x_type)
    if (x_type_res == "binary") {
      x <- rbinom(n, 1, prob = p1)
    } else if (x_type_res == "continuous") {
      x <- rnorm(n, mean = cont_mean, sd = cont_sd)
    } else {
      L <- ord_levels
      if (is.null(L) || L < 2) stop("'ord_levels' must be >= 2 for ordinal x.")
      if (is.null(ord_probs)) ord_probs <- rep(1/L, L)
      if (length(ord_probs) != L) stop("'ord_probs' length must equal 'ord_levels'.")
      ord_probs <- ord_probs / sum(ord_probs)
      x <- sample(0:(L-1), size = n, replace = TRUE, prob = ord_probs)
      x <- ordered(x, levels = 0:(L-1))
    }
  }
  
  # Normalize x for internal math
  x_bin <- NULL; x_cont <- NULL; x_ord <- NULL
  if (x_type_res == "binary") {
    if (is.factor(x)) x <- as.integer(x) - min(as.integer(x))
    x <- as.numeric(x)
    if (!all(x %in% c(0,1))) stop("Binary x must be in {0,1}.")
    x_bin <- x
  } else if (x_type_res == "continuous") {
    x_cont <- as.numeric(x)
    if (cont_standardize) {
      sdx <- sd(x_cont)
      if (is.finite(sdx) && sdx > 0) x_cont <- (x_cont - mean(x_cont)) / sdx
    }
  } else {
    if (is.factor(x)) {
      if (!is.ordered(x)) x <- ordered(x)
      x_ord <- as.integer(x) - 1L
    } else {
      x_ord <- as.integer(x)
      if (min(x_ord) < 0) x_ord <- x_ord - min(x_ord)
    }
  }
  
  # ---- 2) Depth vector ----
  if (length(N_pts) == 1) N_pts <- rep(N_pts, n)
  if (length(N_pts) != n) stop("Length of N_pts must be n or scalar.")
  
  # ---- 3) Cohort center composition q ----
  alpha0 <- if (length(base_alpha)==1) rep(base_alpha, K) else as.numeric(base_alpha)
  if (length(alpha0) != K || any(alpha0 <= 0)) stop("'base_alpha' must be positive (length K) or scalar.")
  q <- as.numeric(.rdirichlet_safe(1, alpha0))
  
  # ---- 4) Signal selection on the provided tree (two-pole + partial random) ----
  m_sig <- if (da_p < 1) round(da_p * K) else min(K, as.integer(da_p))
  m_sig <- max(0L, min(as.integer(m_sig), K))
  
  m_up   <- max(0L, min(as.integer(round(m_sig * da_up_frac)), m_sig))
  m_down <- m_sig - m_up
  
  sgn <- rep(0L, K)
  if (m_sig > 0L) {
    ss <- .pick_signal_sets_2pole(tree = tree, m_up = m_up, m_down = m_down, conc = da_conc)
    if (length(ss$up))   sgn[ss$up]   <- +1L
    if (length(ss$down)) sgn[ss$down] <- -1L
  }
  beta <- sgn * da_logFC  # per-taxon coefficient
  
  # ---- 5) Generate per-sample compositions and apply DA per x ----
  P <- matrix(NA_real_, nrow = n, ncol = K)
  for (i in seq_len(n)){
    p0 <- as.numeric(.rdirichlet_safe(1, q * sample_disp))
    if (any(sgn != 0L)) {
      mul <- rep(1.0, K)
      if (x_type_res == "binary") {
        if (x_bin[i] == 1) {
          idx <- (sgn != 0L)
          mul[idx] <- exp(da_logFC * sgn[idx])
        }
      } else if (x_type_res == "continuous") {
        idx <- (sgn != 0L)
        mul[idx] <- exp( (da_logFC * sgn[idx]) * x_cont[i] )
      } else { # ordinal
        idx <- (sgn != 0L)
        mul[idx] <- exp( (da_logFC * sgn[idx]) * x_ord[i] )
      }
      p <- p0 * mul
      s <- sum(p); if (!is.finite(s) || s <= 0) s <- 1
      P[i, ] <- p / s
    } else {
      P[i, ] <- p0
    }
  }
  
  # ---- 6) Draw counts ----
  OTU <- matrix(0L, nrow = n, ncol = K)
  if (is.finite(dm_theta)) {
    for (i in seq_len(n)){
      alpha <- P[i, ] * dm_theta
      Pstar <- as.numeric(.rdirichlet_safe(1, alpha))
      OTU[i, ] <- as.vector(rmultinom(1, size = N_pts[i], prob = Pstar))
    }
  } else {
    for (i in seq_len(n)){
      OTU[i, ] <- as.vector(rmultinom(1, size = N_pts[i], prob = P[i, ]))
    }
  }
  colnames(OTU) <- tax_names
  rownames(OTU) <- paste0("S", seq_len(n))
  
  # ---- 7) Optional filtering (tree & annotations pruned consistently) ----
  keep <- rep(TRUE, K)
  if (!is.null(min_prev) && min_prev > 0){
    thr <- if (min_prev < 1) ceiling(min_prev * n) else as.integer(min_prev)
    thr <- max(thr, 1L)
    keep <- keep & (colSums(OTU > 0) >= thr)
  }
  if (!is.null(min_var) && min_var > 0){
    keep <- keep & apply(OTU, 2, function(z) var(as.numeric(z)) > min_var)
  }
  if (!all(keep)) {
    kept_names <- colnames(OTU)[keep]
    OTU <- OTU[, keep, drop = FALSE]
    P   <- P[, keep, drop = FALSE]
    sgn <- sgn[keep]
    beta <- beta[keep]
    if (!requireNamespace("ape", quietly = TRUE)) stop("Please install 'ape' for tree pruning.")
    tree <- ape::drop.tip(tree, setdiff(tax_names, kept_names))
    tax_names <- kept_names
  }
  
  # ---- 8) Pack into phyloseq ----
  ps_otu <- otu_table(OTU, taxa_are_rows = FALSE)
  # keep x class (ordered for ordinal)
  x_for_sd <- x; if (is.ordered(x_for_sd)) x_for_sd <- x_for_sd  # noop; ensures ordered stays ordered
  ps_sd <- sample_data(data.frame(
    SampleID = rownames(OTU),
    x = x_for_sd,
    N = N_pts,
    x_type = if (is.factor(x)) if (is.ordered(x)) "ordinal" else if (nlevels(x)==2) "binary" else "ordinal"
    else if (is.numeric(x)) if (all(unique(x) %in% c(0,1))) "binary" else "continuous" else NA_character_,
    row.names = rownames(OTU),
    check.names = FALSE
  ))
  da_label <- ifelse(sgn > 0, "up", ifelse(sgn < 0, "down", "none"))
  ps_tax <- tax_table(cbind(
    DA_sign = as.character(da_label),
    DA_coef = as.character(beta)
  ))
  rownames(ps_tax) <- colnames(OTU)
  
  phyloseq(ps_otu, ps_sd, ps_tax, phy_tree(tree))
}
# ======================== End: SimuDA_FixedTree ==================================
