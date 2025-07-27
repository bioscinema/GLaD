# ── libraries ────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)      # for label_parse()

# ── read data ────────────────────────────────────────────────────────────────
result <- read.csv("CaseStudyResult/filtered_CaseStudyResult.csv")

# ── helper vectors ───────────────────────────────────────────────────────────

# 1. exact method order required on the y‑axis
ordered_methods <- c(
  "GLaD_0.5_unweighted", "GLaD_0.5_weighted",
  "GLaD_1_unweighted",   "GLaD_1_weighted",
  "unifrac_unweighted",  "unifrac_weighted_phyloseq",
  "GUniFrac", "aitchison", "bray", "euclidean", "jaccard"
)

# LaTeX‑like labels (no bold so `label_parse()` works out‑of‑box)
method_latex_labels <- c(
  "GLaD_0.5_unweighted"        = "Unweighted~GLaD~(rho==0.5)",
  "GLaD_0.5_weighted"          = "Weighted~GLaD~(rho==0.5)",
  "GLaD_1_unweighted"          = "Unweighted~GLaD~(rho==1)",
  "GLaD_1_weighted"            = "Weighted~GLaD~(rho==1)",
  "unifrac_unweighted"         = "Unweighted~UniFrac",
  "unifrac_weighted_phyloseq"  = "Weighted~UniFrac",
  "GUniFrac"                   = "GUniFrac",
  "aitchison"                  = "Aitchison",
  "bray"                       = "Bray-Curtis",
  "euclidean"                  = "Euclidean",
  "jaccard"                    = "Jaccard"
)

# dataset labels for the x‑axis
dataset_labels <- c(
  "Concordance_oral_16S" = "ORIGINS Oral Microbiome",
  "covid"                = "COVID-Gut Microbiome",
  "HIV_gut"              = "HIV-Gut Microbiome",
  "IBD_16s"              = "IBD Dataset (16S)",
  "IBD_wgs"              = "IBD Dataset (WGS)",
  "Pig_gut"              = "Pig Gut Microbiome"
)

# ── wrangle & rank ───────────────────────────────────────────────────────────
plot_df <- result %>%
  rename(
    Dataset     = source,
    Method      = Dissimilarity,
    PERMANOVA_F = pseudo_F,
    MiRKAT_R2   = MiRKAT_Rsquared
  ) %>%
  filter(Method %in% ordered_methods) %>%
  mutate(Method = factor(Method, levels = ordered_methods)) %>%     # strict order
  group_by(Dataset) %>%
  mutate(
    F_rank  = rank(-PERMANOVA_F),   # larger F  ⇒ better (rank 1)
    R2_rank = rank(-MiRKAT_R2)      # larger R² ⇒ better
  ) %>%
  pivot_longer(c(F_rank, R2_rank),
               names_to  = "Metric",
               values_to = "Rank") %>%
  ungroup() %>%
  mutate(
    Metric        = recode(Metric,
                           F_rank  = "PERMANOVA~F",
                           R2_rank = "MiRKAT~R^2"),
    Dataset_label = recode(Dataset, !!!dataset_labels),
    Method_latex  = recode(as.character(Method), !!!method_latex_labels),
    Method_latex  = factor(Method_latex,                        # keep y‑axis order
                           levels = method_latex_labels[ordered_methods]),
    Rank_f        = factor(Rank, levels = sort(unique(Rank)))   # ordered factor
  )


method_order <- plot_df %>%
  group_by(Method) %>%
  summarise(mean_rank = mean(Rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>%
  pull(Method)


plot_df <- result %>%
  rename(
    Dataset     = source,
    Method      = Dissimilarity,
    PERMANOVA_F = pseudo_F,
    MiRKAT_R2   = MiRKAT_Rsquared
  ) %>%
  filter(Method %in% method_order) %>%
  mutate(Method = factor(Method, levels = method_order)) %>%     # strict order
  group_by(Dataset) %>%
  mutate(
    F_rank  = rank(-PERMANOVA_F),   # larger F  ⇒ better (rank 1)
    R2_rank = rank(-MiRKAT_R2)      # larger R² ⇒ better
  ) %>%
  pivot_longer(c(F_rank, R2_rank),
               names_to  = "Metric",
               values_to = "Rank") %>%
  ungroup() %>%
  mutate(
    Metric        = recode(Metric,
                           F_rank  = "PERMANOVA~F",
                           R2_rank = "MiRKAT~R^2"),
    Dataset_label = recode(Dataset, !!!dataset_labels),
    Method_latex  = recode(as.character(Method), !!!method_latex_labels),
    Method_latex  = factor(Method_latex,                        # keep y‑axis order
                           levels = method_latex_labels[method_order]),
    Rank_f        = factor(Rank, levels = sort(unique(Rank)))   # ordered factor
  )

# ── palette & size: same length as #ranks ────────────────────────────────────
n_ranks   <- nlevels(plot_df$Rank_f)                        # e.g. 11
rank_pal  <- colorRampPalette(c("#d7312d", "#f2724d","#fee395", "#fef9b7","#acd2e5", "#6090c1"))(n_ranks)
size_vals <- rescale(rev(seq_len(n_ranks)), to = c(2, 8))   # rank 1 = largest

# ── plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(plot_df,
            aes(x    = Dataset_label,
                y    = Method_latex,
                fill = Rank_f,
                size = Rank_f)) +
  geom_point(shape = 21, colour = "black") +
  facet_wrap(~ Metric, scales = "free_x", labeller = label_parsed) +
  scale_fill_manual(values = rank_pal,
                    name   = "Rank\n(1 = best)") +
  scale_size_manual(values = size_vals,
                    name   = "Rank\n(1 = best)") +
  guides(
    fill = guide_legend(
      override.aes = list(size   = size_vals,
                          shape  = 21,
                          colour = "black")),
    size = "none"
  ) +
  scale_y_discrete(
    limits = rev(levels(plot_df$Method_latex)),   # ← reverse the order
    labels = label_parse()
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1.2, "lines"),
    strip.text    = element_text(size = 14)
  ) +
  labs(x = "", y = "")

print(p)

ggsave("CaseStudyResult/bubble_plot.png", plot = p, width = 10, height = 6, dpi = 1000)


########################################## Radar plot #########################################

# ── packages ───────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)

# ── ordering, labels, colours ──────────────────────────────────────────────────
ordered_methods <- c(
  "jaccard",
  "GLaD_0.5_unweighted", "GLaD_0.5_weighted",
  "GLaD_1_unweighted",   "GLaD_1_weighted",
  "unifrac_unweighted",  "unifrac_weighted_phyloseq",
  "GUniFrac", "aitchison", "bray", "euclidean"
)

method_latex_labels <- c(
  "GLaD_0.5_unweighted"      = "bold(Unweighted~GLaD~(rho==0.5))",
  "GLaD_0.5_weighted"        = "bold(Weighted~GLaD~(rho==0.5))",
  "GLaD_1_unweighted"        = "bold(Unweighted~GLaD~(rho==1))",
  "GLaD_1_weighted"          = "bold(Weighted~GLaD~(rho==1))",
  "unifrac_unweighted"       = "bold(Unweighted~UniFrac)",
  "unifrac_weighted_phyloseq"= "bold(Weighted~UniFrac)",
  "GUniFrac"                 = "bold(GUniFrac)",
  "aitchison"                = "bold(Aitchison)",
  "bray"                     = "bold(Bray-Curtis)",
  "euclidean"                = "bold(Euclidean)",
  "jaccard"                  = "bold(Jaccard)"
)

method_colors <- c(
  "GLaD_0.5_unweighted"      = "mediumorchid4",
  "GLaD_0.5_weighted"        = "mediumorchid1",
  "GLaD_1_unweighted"        = "mediumpurple1",
  "GLaD_1_weighted"          = "mediumpurple3",
  "unifrac_unweighted"       = "royalblue3",
  "unifrac_weighted_phyloseq"= "steelblue1",
  "GUniFrac"                 = "steelblue3",
  "aitchison"                = "gray48",
  "bray"                     = "coral3",
  "euclidean"                = "mediumseagreen",
  "jaccard"                  = "goldenrod2"
)

# ── build summary table ────────────────────────────────────────────────────────
fni_summary <- result %>%
  filter(Dissimilarity %in% ordered_methods) %>%
  mutate(IsEuclideanC = if_else(IsEuclidean == "YES", 1, 0)) %>%
  group_by(Dissimilarity) %>%
  summarise(FNI_zero_count = sum(IsEuclideanC), .groups = "drop") %>%
  mutate(method = factor(Dissimilarity, levels = ordered_methods))

# helper: parse labels on demand
label_fun <- function(x) parse(text = method_latex_labels[x])

# ── plot ───────────────────────────────────────────────────────────────────────
r <- ggplot(fni_summary,
            aes(x = method, y = FNI_zero_count, fill = method)) +
  geom_bar(stat = "identity", width = 1, colour = "black", show.legend = FALSE) +
  coord_polar(start = 0, clip = "off") +                    # don’t clip labels
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(breaks = pretty(fni_summary$FNI_zero_count)) +
  scale_x_discrete(labels = label_fun) +                    # parse LaTeX‑style
  theme_minimal(base_size = 13) +
  theme(
    axis.title      = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1, colour = "black"),
    panel.grid      = element_line(colour = "grey80", linetype = "dotted"),
    panel.grid.major.y = element_line(colour = "grey60"),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.margin     = margin(t = 10, r = 40, b = 10, l = 40) # extra room
  ) 

print(r)
ggsave("CaseStudyResult/radar_plot.png", plot = r, width = 9, height = 6, dpi = 1000)
