library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(scales)
result <- read.csv("CaseStudyResult/filtered_CaseStudyResult.csv")
# Define custom order of distance methods
ordered_methods <- c(
  "GLaD_0.5_unweighted",
  "GLaD_0.5_weighted",
  "GLaD_1_unweighted",
  "GLaD_1_weighted",
  "unifrac_unweighted",
  "unifrac_weighted_phyloseq",
  "GUniFrac",
  "aitchison",
  "bray",
  "euclidean",
  "jaccard"
)

# Define LaTeX-friendly method labels
method_latex_labels <- c(
  "GLaD_0.5_unweighted" = "Unweighted~GLaD~(rho==0.5)",
  "GLaD_0.5_weighted" = "Weighted~GLaD~(rho==0.5)",
  "GLaD_1_unweighted" = "Unweighted~GLaD~(rho==1)",
  "GLaD_1_weighted" = "Weighted~GLaD~(rho==1)",
  "unifrac_unweighted" = "Unweighted~UniFrac",
  "unifrac_weighted_phyloseq" = "Weighted~UniFrac",
  "GUniFrac" = "GUniFrac",
  "aitchison" = "Aitchison",
  "bray" = "Bray-Curtis",
  "euclidean" = "Euclidean",
  "jaccard" = "Jaccard"
)

# Define dataset labels
dataset_labels <- c(
  "Concordance_oral_16S" = "ORIGINS Oral Microbiome",
  "covid" = "COVID-Gut Microbiome",
  "HIV_gut" = "HIV-Gut Microbiome",
  "IBD_16s" = "IBD Dataset (16S)",
  "IBD_wgs" = "IBD Dataset (WGS)",
  "Pig_gut" = "Pig Gut Microbiome"
)

# Process and reshape data
plot_df <- result %>%
  rename(
    Dataset = source,
    Method = Dissimilarity,
    PERMANOVA_F = pseudo_F,
    MiRKAT_R2 = MiRKAT_Rsquared
  ) %>%
  filter(Method %in% ordered_methods) %>%
  mutate(Method = factor(Method, levels = ordered_methods)) %>%
  group_by(Dataset) %>%
  mutate(
    F_rank = rank(-PERMANOVA_F),
    R2_rank = rank(-MiRKAT_R2)
  ) %>%
  pivot_longer(cols = c(F_rank, R2_rank), names_to = "Metric", values_to = "Rank") %>%
  ungroup() %>%
  mutate(
    Metric = recode(Metric,
                    F_rank = "PERMANOVA~F",
                    R2_rank = "MiRKAT~R^2"),
    Dataset_label = recode(Dataset, !!!dataset_labels),
    Method_latex = recode(as.character(Method), !!!method_latex_labels)
  )

# Plot
p <- ggplot(plot_df, aes(x = Dataset_label, y = Method_latex, size = Rank)) +
  geom_point(aes(fill = Rank), shape = 21, color = "black") +
  facet_wrap(~ Metric, scales = "free_x", labeller = label_parsed) +
  scale_size_continuous(range = c(2, 8), trans = "reverse") +
  scale_fill_gradientn(
    colours = c("#d73027", "#fee08b", "#1a9850"),
    values = scales::rescale(c(1, 5.5, 11)),  # 中间值可调
    name = "Rank\n(1 = best)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1.2, "lines"),
    strip.text = element_text(size = 14)
  ) +
  labs(x = "Dataset", y = "Method", size = "Rank", fill = "Rank") +
  scale_y_discrete(limits = rev(method_latex_labels[ordered_methods]), labels = scales::label_parse())

p
ggsave("CaseStudyResult/bubble_plot.png", plot = p, width = 10, height = 6, dpi = 1000)

