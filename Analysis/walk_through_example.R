source("Eva_Diss_Structure.R")
source("Eva_Cluster_Structure.R")
source("GLaD.R") # GLaD use internal nodes

################################

### step 1: filter the phyloseq and get a group variable
load("IBD_16s_data_V4.RData")
sample_data <- sample_data(phy1)
physeq_sub <- subset_samples(phy1, diagnosis %in% c("cd", "uc"))
otu_mat <- otu_table(physeq_sub)
if (!taxa_are_rows(physeq_sub)) {
  otu_mat <- t(otu_mat)
}
prevalence <- apply(otu_mat, 1, function(x) sum(x > 0))
sample_size <- nsamples(physeq_sub)
threshold <- 1 / sample_size
keep_otus <- prevalence >= threshold * sample_size  
physeq_filtered <- prune_taxa(keep_otus, physeq_sub)

### step 2: calculate dissimilarities (GLaD, weighted unifrac, unweighted unifrac, generalized unifrac, Euclidean, Aitchison, Bray-Curtis, Jaccard)
wunifrac_dist <- UniFrac(physeq_filtered, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE)
unifrac_dist <- UniFrac(physeq_filtered, weighted = FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE)
wglad_dist <- GLaD(physeq_filtered, weighted = T)
uglad_dist <- GLaD(physeq_filtered, weighted = F)
wglad1_dist <- GLaD_eigen(physeq_filtered, weighted = T)
uglad1_dist <- GLaD_eigen(physeq_filtered, weighted = F)
euc_dist <- distance(physeq_filtered,method = "euclidean")

### step 3: evaluation based on dissimilarity structure
wunifrac_eva1 <- evaluate_structure_preservation(wunifrac_dist)
unifrac_eva1 <- evaluate_structure_preservation(unifrac_dist)
wglad_eva1 <- evaluate_structure_preservation(wglad_dist)
uglad_eva1 <- evaluate_structure_preservation(uglad_dist)
wglad1_eva1 <- evaluate_structure_preservation(wglad1_dist)
uglad1_eva1 <- evaluate_structure_preservation(uglad1_dist)
euc_eva1 <- evaluate_structure_preservation(euc_dist)

### step 4: evaluation based on cluster structure
group_var = as.factor(physeq_filtered@sam_data$diagnosis)
wunifrac_eva2 <- evaluate_clustering_quality(wunifrac_dist, group_var)
unifrac_eva2 <- evaluate_clustering_quality(unifrac_dist, group_var)
wglad_eva2 <- evaluate_clustering_quality(wglad_dist, group_var)
uglad_eva2 <- evaluate_clustering_quality(uglad_dist, group_var)
wglad1_eva2 <- evaluate_clustering_quality(wglad1_dist, group_var)
uglad1_eva2 <- evaluate_clustering_quality(uglad1_dist, group_var)
euc_eva2 <- evaluate_clustering_quality(euc_dist, group_var)


