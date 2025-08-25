#' Benchmark a panel of beta-diversity distances from a phyloseq object
#'
#' Computes multiple beta-diversity distances from a \code{phyloseq} object,
#' evaluates each distance using user-supplied evaluators, aggregates metrics,
#' computes per-metric Borda scores (within-metric, within-dataset), and
#' optionally returns a tournament plot.
#'
#' @param ps A \code{phyloseq} object.
#' @param methods Character vector of method names (see \code{\link{compute_distance_by_name}}).
#'   Defaults to a mix of UniFrac/GLAD and \pkg{vegan} distances.
#' @param rho Numeric in \eqn{[0,1)} passed to GLAD distances.
#' @param pseudocount Numeric pseudocount used for the Aitchison (CLR+Euclidean) distance.
#' @param figure Logical; if \code{TRUE}, also return a ggplot object (see Details).
#' @param group Sample grouping: either a column name in \code{sample_data(ps)}
#'   or a vector of length \code{nsamples(ps)}.
#' @param dataset_name Character tag stored in the output table and used in the plot title.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{table}}{Data frame with one row per method and evaluation metrics
#'   (FNI, stress, trustworthiness, continuity, Mantel r/p, avg silhouette, DBI, CH, Dunn)
#'   plus a Borda score column.}
#'   \item{\code{plot}}{A \pkg{ggplot2} object from \code{\link{make_tournament_plot}},
#'   or \code{NULL} when \code{figure = FALSE} or the plotting function is unavailable.}
#' }
#'
#' @details
#' This function assumes the following evaluators are available in the calling
#' environment:
#' \itemize{
#'   \item \code{evaluate_clustering_quality(D, group)}
#'   \item \code{evaluate_structure_preservation(D, ord_dim)}
#' }
#' Distances are computed via \code{\link{compute_distance_by_name}}:
#' UniFrac/GUniFrac are computed with \pkg{phyloseq}/\pkg{GUniFrac}; GLAD variants
#' are delegated to user-provided \code{GLaD()} / \code{GLaD_eigen()} functions;
#' all other named methods are implemented using \pkg{vegan} \code{vegdist()}
#' (with the appropriate preprocessing for chord/hellinger/chisq and CLR for Aitchison).
#'
#' Borda scores are summed over metrics after ranking \emph{within each metric and dataset}
#' (lower-is-better for \code{FNI} and \code{stress}; higher-is-better otherwise).
#' If \code{figure = TRUE} and \code{\link{make_tournament_plot}} exists, a tournament
#' plot is produced for \code{dataset_name}.
#'
#' @examples
#' \dontrun{
#' res <- benchmark_beta_methods(ps, group = "GroupVar", figure = TRUE)
#' res$table
#' print(res$plot)
#' }
#'
#' @seealso \code{\link{compute_distance_by_name}}, \code{\link{make_tournament_plot}}
#'
#' @importFrom phyloseq sample_names sample_data nsamples taxa_are_rows otu_table UniFrac phy_tree
#' @importFrom GUniFrac GUniFrac
#' @importFrom vegan vegdist decostand
#' @importFrom dplyr select mutate group_by summarise ungroup left_join n_distinct bind_rows filter arrange pull
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point scale_y_reverse dup_axis theme_minimal labs theme element_text
#' @importFrom stats as.dist
#' @export
benchmark_beta_methods <- function(
    ps,
    methods = c("Weighted Unifrac", "Unweighted Unifrac", "Generalized Unifrac",
                "Weighted Glad1", "Unweighted Glad1", "Weighted Glad", "Unweighted Glad",
                "manhattan", "euclidean", "bray", "jaccard", "chisq", "chord", "hellinger", "aitchison"),
    rho = 0.5,
    pseudocount = 0.5,
    figure = TRUE,
    group,                         # sample_data column name OR vector of length nsamples(ps)
    dataset_name = "Dataset"
) {
  stopifnot(inherits(ps, "phyloseq"))
  
  # resolve grouping
  samps <- sample_names(ps)
  if (length(group) == 1 && is.character(group) && group %in% colnames(sample_data(ps))) {
    grp <- as.vector(sample_data(ps)[[group]])
  } else if (length(group) == nsamples(ps)) {
    grp <- as.vector(group)
  } else stop("`group` must be a sample_data() column name or a vector of length nsamples(ps).")
  names(grp) <- samps
  
  rows <- list()
  for (m in methods) {
    D <- tryCatch(
      compute_distance_by_name(ps, m, rho = rho, pseudocount = pseudocount),
      error = function(e) { warning(sprintf("Method '%s' failed: %s", m, e$message)); NULL }
    )
    if (is.null(D)) next
    
    # ensure distance order matches sample order
    if (!all(labels(D) == samps)) {
      M <- as.matrix(D); M <- M[samps, samps]; D <- as.dist(M)
    }
    
    cq <- evaluate_clustering_quality(D, grp)        # assumed to exist
    sp <- evaluate_structure_preservation(D, 2)      # assumed to exist
    
    rows[[length(rows) + 1]] <- data.frame(
      data_name       = dataset_name,
      method          = m,
      FNI             = sp$FNI,
      trustworthiness = sp$trustworthiness,
      continuity      = sp$continuity,
      mantel_r        = sp$mantel_r,
      mantel_p        = sp$mantel_p,
      stress          = sp$stress,
      avg_silhouette  = cq$avg_silhouette,
      dbi             = cq$dbi,
      chi             = cq$chi,
      dunn            = cq$dunn,
      stringsAsFactors = FALSE
    )
  }
  
  if (!length(rows)) stop("No distances were successfully computed.")
  df <- bind_rows(rows)
  
  # Borda score (higher is better). Metrics: lower-is-better = FNI, stress; higher-is-better = rest (except p-values).
  metrics_for_borda <- c("FNI","stress","trustworthiness","continuity","mantel_r","avg_silhouette")
  rank_df <- df %>%
    select(data_name, method, all_of(metrics_for_borda)) %>%
    pivot_longer(cols = all_of(metrics_for_borda), names_to = "metric", values_to = "value") %>%
    group_by(data_name, metric) %>%
    mutate(
      rank_metric = ifelse(metric %in% c("FNI","stress"),
                           rank(value, ties.method = "min"),
                           rank(-value, ties.method = "min")),
      N_metric = dplyr::n_distinct(method)
    ) %>%
    ungroup()
  
  borda <- rank_df %>%
    group_by(data_name, method) %>%
    summarise(borda = sum(N_metric - rank_metric + 1), .groups = "drop")
  
  df <- df %>% left_join(borda, by = c("data_name","method"))
  
  # Optional figure: uses your make_tournament_plot(df, ...) if it exists
  p <- NULL
  if (isTRUE(figure)) {
    if (exists("make_tournament_plot")) {
      p <- make_tournament_plot(df, selected_dataset = dataset_name, dataset_name = dataset_name,
                                metrics_to_remove = c("chi","dunn","mantel_p","dbi"))
    } else {
      message("Figure requested but make_tournament_plot() not found — returning table only.")
    }
  }
  
  list(table = df, plot = p)
}


# ---- Helpers -----------------------------------------------------------------
get_sample_matrix <- function(ps) {
  X <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) X <- t(X)
  X[sample_names(ps), , drop = FALSE]
}

# For methods other than UniFrac / GLAD, compute distances using vegdist only
compute_distance_with_vegdist <- function(ps, method, pseudocount = 1e-6) {
  X <- get_sample_matrix(ps)
  
  m <- tolower(method)
  if (m %in% c("manhattan","euclidean","bray", "jaccard", "chord", "hellinger", "chisq")) {
    return(vegdist(X, method = m))
  }
  
  
  if (m == "aitchison") {
    return(vegdist(X, method = "aitchison", pseudocount = pseudocount))
  }
  
  stop(sprintf("Unsupported vegdist method: %s", method))
}

# UniFrac / GUniFrac
dist_weighted_unifrac   <- function(ps) phyloseq::UniFrac(ps, weighted = TRUE,  normalized = TRUE)
dist_unweighted_unifrac <- function(ps) phyloseq::UniFrac(ps, weighted = FALSE, normalized = TRUE)
dist_gunifrac           <- function(ps, alpha = 0.5) {
  X  <- get_sample_matrix(ps)
  tr <- phyloseq::phy_tree(ps)
  out <- GUniFrac::GUniFrac(X, tr, alpha = c(alpha))$unifracs
  as.dist(out[, , paste0("d_", alpha)])
}

compute_distance_by_name <- function(ps, method, rho = 0.5, pseudocount = 1e-6) {
  m <- tolower(method)
  switch(m,
         "weighted unifrac"    = dist_weighted_unifrac(ps),
         "unweighted unifrac"  = dist_unweighted_unifrac(ps),
         "generalized unifrac" = dist_gunifrac(ps, alpha = 0.5),
         
         "weighted glad1"      = GLaD_eigen(ps, weighted = T),
         "unweighted glad1"    = GLaD_eigen(ps, weighted = F),
         "weighted glad"       = GLaD(ps,  rho = rho, weighted = T),
         "unweighted glad"     = GLaD(ps, rho = rho, weighted = F),
         
         # everything else strictly via vegdist (possibly after vegan::decostand)
         compute_distance_with_vegdist(ps, method, pseudocount = pseudocount)
  )
}

# ------------------------- MAIN ENTRY -----------------------------------------
benchmark_beta_methods <- function(
    ps,
    methods = c("Weighted Unifrac", "Unweighted Unifrac", "Generalized Unifrac",
                "Weighted Glad1", "Unweighted Glad1", "Weighted Glad", "Unweighted Glad",
                "manhattan", "euclidean", "bray", "jaccard", "chisq", "chord", "hellinger", "aitchison"),
    rho = 0.5,
    pseudocount = 0.5,
    figure = TRUE,
    group,                         # sample_data column name OR vector of length nsamples(ps)
    dataset_name = "Dataset"
) {
  stopifnot(inherits(ps, "phyloseq"))
  
  # resolve grouping
  samps <- sample_names(ps)
  if (length(group) == 1 && is.character(group) && group %in% colnames(sample_data(ps))) {
    grp <- as.vector(sample_data(ps)[[group]])
  } else if (length(group) == nsamples(ps)) {
    grp <- as.vector(group)
  } else stop("`group` must be a sample_data() column name or a vector of length nsamples(ps).")
  names(grp) <- samps
  
  rows <- list()
  for (m in methods) {
    D <- tryCatch(
      compute_distance_by_name(ps, m, rho = rho, pseudocount = pseudocount),
      error = function(e) { warning(sprintf("Method '%s' failed: %s", m, e$message)); NULL }
    )
    if (is.null(D)) next
    
    # ensure distance order matches sample order
    if (!all(labels(D) == samps)) {
      M <- as.matrix(D); M <- M[samps, samps]; D <- as.dist(M)
    }
    
    cq <- evaluate_clustering_quality(D, grp)        # assumed to exist
    sp <- evaluate_structure_preservation(D, 2)      # assumed to exist
    
    rows[[length(rows) + 1]] <- data.frame(
      data_name       = dataset_name,
      method          = m,
      FNI             = sp$FNI,
      trustworthiness = sp$trustworthiness,
      continuity      = sp$continuity,
      mantel_r        = sp$mantel_r,
      mantel_p        = sp$mantel_p,
      stress          = sp$stress,
      avg_silhouette  = cq$avg_silhouette,
      dbi             = cq$dbi,
      chi             = cq$chi,
      dunn            = cq$dunn,
      stringsAsFactors = FALSE
    )
  }
  
  if (!length(rows)) stop("No distances were successfully computed.")
  df <- bind_rows(rows)
  
  # Borda score (higher is better). Metrics: lower-is-better = FNI, stress; higher-is-better = rest (except p-values).
  metrics_for_borda <- c("FNI","stress","trustworthiness","continuity","mantel_r","avg_silhouette")
  rank_df <- df %>%
    select(data_name, method, all_of(metrics_for_borda)) %>%
    pivot_longer(cols = all_of(metrics_for_borda), names_to = "metric", values_to = "value") %>%
    group_by(data_name, metric) %>%
    mutate(
      rank_metric = ifelse(metric %in% c("FNI","stress"),
                           rank(value, ties.method = "min"),
                           rank(-value, ties.method = "min")),
      N_metric = dplyr::n_distinct(method)
    ) %>%
    ungroup()
  
  borda <- rank_df %>%
    group_by(data_name, method) %>%
    summarise(borda = sum(N_metric - rank_metric + 1), .groups = "drop")
  
  df <- df %>% left_join(borda, by = c("data_name","method"))
  
  # Optional figure: uses your make_tournament_plot(df, ...) if it exists
  p <- NULL
  if (isTRUE(figure)) {
    if (exists("make_tournament_plot")) {
      p <- make_tournament_plot(df, selected_dataset = dataset_name, dataset_name = dataset_name,
                                metrics_to_remove = c("chi","dunn","mantel_p","dbi"))
    } else {
      message("Figure requested but make_tournament_plot() not found — returning table only.")
    }
  }
  
  list(table = df, plot = p)
}


make_tournament_plot <- function(df,
                                 selected_dataset = NULL,
                                 dataset_name = NULL,
                                 metrics_to_remove = c("chi","dunn","mantel_p","dbi")) {
  
  # Choose dataset
  datasets <- unique(df$data_name)
  if (is.null(selected_dataset) || !(selected_dataset %in% datasets)) {
    selected_dataset <- datasets[1]
  }
  if (is.null(dataset_name)) dataset_name <- selected_dataset
  
  # Keep only columns we might plot; ignore any input 'borda' and recompute it here
  candidates <- c("FNI","trustworthiness","continuity","mantel_r","stress","avg_silhouette",
                  "mantel_p","dbi","chi","dunn")
  have <- intersect(candidates, names(df))
  metric_cols <- setdiff(have, metrics_to_remove)
  if (length(metric_cols) == 0L) stop("No metrics left to plot after filtering.")
  
  df_ds <- df %>% filter(data_name == selected_dataset)
  
  # Long form + ranking per metric
  df_long <- df_ds %>%
    select(data_name, method, all_of(metric_cols)) %>%
    pivot_longer(-c(data_name, method), names_to = "metric", values_to = "value")
  
  df_ranked <- df_long %>%
    group_by(data_name, metric) %>%
    mutate(
      rank_metric = if_else(metric %in% c("FNI","stress"),
                            rank(value, ties.method = "min"),
                            rank(-value, ties.method = "min"))
    ) %>%
    ungroup()
  
  # Borda based on the metrics actually plotted
  N <- n_distinct(df_ranked$method)
  borda_score <- df_ranked %>%
    mutate(score = N - rank_metric + 1) %>%
    group_by(data_name, method) %>%
    summarise(borda = sum(score), .groups = "drop")
  
  # Tie-broken final ranks within each metric
  df_ranked2 <- df_ranked %>%
    left_join(borda_score, by = c("data_name","method")) %>%
    group_by(data_name, metric) %>%
    arrange(rank_metric, desc(borda), method, .by_group = TRUE) %>%
    mutate(final_rank = row_number()) %>%
    ungroup()
  
  # Add Borda as a metric
  df_borda <- borda_score %>%
    mutate(metric = "BordaCount",
           value  = borda,
           rank_metric = rank(-borda, ties.method = "min")) %>%
    group_by(data_name, metric) %>%
    arrange(rank_metric, desc(borda), method, .by_group = TRUE) %>%
    mutate(final_rank = row_number()) %>%
    ungroup()
  
  df_plot <- bind_rows(
    df_ranked2 %>% select(data_name, method, metric, value, rank_metric, borda, final_rank),
    df_borda
  ) %>% filter(data_name == selected_dataset)
  
  # Metric order & labels
  preferred <- c("FNI","trustworthiness","continuity","mantel_r","stress","avg_silhouette","BordaCount")
  metric_order <- intersect(preferred, unique(df_plot$metric))
  metric_labels <- c(
    FNI = "Neg. Eigenvalue Frac.",
    trustworthiness = "Trustworthiness",
    continuity = "Continuity",
    mantel_r = "Mantel r",
    stress = "Stress-1",
    avg_silhouette = "Avg. Silhouette",
    BordaCount = "Borda Score"
  )
  df_plot$metric_label <- metric_labels[df_plot$metric]
  
  # Left (first metric) and right (Borda) y-axis label orders
  first_metric <- metric_order[1]
  last_metric  <- metric_order[length(metric_order)]
  
  left_order <- df_plot %>%
    filter(metric == first_metric) %>%
    arrange(final_rank) %>%
    pull(method)
  
  right_order <- df_plot %>%
    filter(metric == last_metric) %>%
    arrange(final_rank) %>%
    pull(method)
  
  # Ensure consistent breaks
  n_methods <- n_distinct(df_plot$method)
  y_breaks  <- seq_len(n_methods)
  
  ggplot(df_plot, aes(
    x = factor(metric_label, levels = metric_labels[metric_order]),
    y = final_rank, group = method, color = method
  )) +
    ggbump::geom_bump(size = 0.9, alpha = 0.85) +
    geom_point(size = 2, alpha = 0.95) +
    scale_y_reverse(
      breaks = y_breaks,
      labels = left_order,
      sec.axis = dup_axis(breaks = y_breaks, labels = right_order, name = NULL)
    ) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Tournament Plot -", dataset_name), x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
      axis.text.y.left  = element_text(size = 10, face = "bold"),
      axis.text.y.right = element_text(size = 10, face = "bold"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

