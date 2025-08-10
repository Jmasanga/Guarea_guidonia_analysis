setwd("~/Fungal_composition")

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(ANCOMBC)
library(DT)
library(tidyr)

Otu_table <- read.csv("otu_table.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)
Taxa <- read.csv("Otu_taxa.csv", header=T, dec=".",sep=",", skipNul = T)
Metadata <- read.csv("Guarea_metadata.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)

OTU <- data.matrix(Otu_table)
TAX <- as.matrix(Taxa)
Metadata <- Guarea_Metadata
all(rownames(OTU) == rownames(TAX))
ps <- phyloseq(otu_table(OTU, taxa_are_rows=TRUE), tax_table(TAX), sample_data(Metadata))
ps

###Prune
#Subset taxa
#remove samples with NAs
nrow(ps_otu1@tax_table)
ps.prune = subset_taxa(ps_otu1, Phylum=="Phylum")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps_otu1, Phylum !="NA")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, Class != "NA")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, Family != "NA")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, Genus != "NA")
nrow(ps.prune@tax_table)

# prune OTUs that are not present in at least one sample
nrow(ps.prune@otu_table)
ps.prune <- prune_taxa(taxa_sums(ps.prune) > 0, ps.prune)
nrow(ps.prune@otu_table)

# prune samples with no OTUS
ncol(ps.prune@otu_table)
ps.prune <- prune_samples(sample_sums(ps.prune) > 0, ps.prune)
ncol(ps.prune@otu_table)

#### Rarefy ####
min(colSums(ps.prune@otu_table))
max(colSums(ps.prune@otu_table))

# rarecurve(t(otu_table(ps.prune)), sample = ,step=50, cex=0.5, label = F) 
ps.rare <- rarefy_even_depth(ps.prune, sample.size = 1000, rngseed = 336, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare <- prune_samples(sample_sums(ps.rare) > 0, ps.rare)

min(colSums(ps.rare@otu_table))
max(colSums(ps.rare@otu_table))

ps
ps.prune
ps.rare

sample_data(ps.even) <- sample_data(ps.prune)[sample_names(ps.even), ]

Phylum = get_taxa_unique(ps.rare, "Phylum")
classes = get_taxa_unique(ps.rare, "Class")
orders = get_taxa_unique(ps.rare, "Order")
family = get_taxa_unique(ps.rare, "Family")
genus = get_taxa_unique(ps.rare, "Genus")
species = get_taxa_unique(ps.rare, "Species")

# Relative Abundance Transformation
ps.even = transform_sample_counts(ps_otu, function(otu) otu/sum(otu))

ps.even
ps.rare

length(get_taxa_unique(ps.rare, "Phylum")) 
length(get_taxa_unique(ps.rare, "Class"))
length(get_taxa_unique(ps.rare, "Order"))
length(get_taxa_unique(ps.rare, "Family")) 
length(get_taxa_unique(ps.rare, "Genus"))
length(get_taxa_unique(ps.rare, "Species"))


##relative abundance of the phylum
physeq_phylum <- tax_glom(ps.rel, taxrank = "Phylum")
phylum_abundance <- psmelt(physeq_phylum) %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
write.csv(phylum_abundance, "results/Fungal_Phylum_relative_abundance.csv")
p_phylum <- ggplot(phylum_abundance, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative Abundance", title = "Fungal Phylum Composition")
ggsave("results/Fungal_Phylum_Composition.png", p_phylum, width = 10, height = 6)

## Alpha diversity
alpha_div <- estimate_richness(ps.rare, measures = c("Shannon", "Chao1"))
alpha_div <- cbind(alpha_div, sample_data(ps.rare))
write.csv(alpha_div, "results/Fungal_alpha_diversity.csv")

## Beta diversity
ps.pcoa <- ordinate(ps.rare, method = "PCoA", distance = "bray")
p_pcoa <- plot_ordination(ps.rare, ps.pcoa, type = "samples",
                          color = "Microbiome_treatment",
                          shape = "Soil_inoculum") +
  geom_point(size = 3) +
  theme_minimal()
ggsave("results/Fungal_PCoA.png", p_pcoa, width = 8, height = 6)

## PERMANOVA
dist_mat <- phyloseq::distance(ps.rare, method = "bray")
metadata_df <- data.frame(sample_data(ps.rare))
permanova <- adonis2(dist_mat ~ Microbiome_treatment + Soil_inoculum +
                       Density + Moisture,
                     data = metadata_df, permutations = 999)
write.csv(as.data.frame(permanova), "results/Fungal_PERMANOVA_results.csv")


## Differential Abundance using ANCOM-BC2 
library(ANCOMBC)
if (requireNamespace("microbiome", quietly = TRUE)) {
  data(atlas1006, package = "microbiome")
  # subset to baseline
  pseq = phyloseq::subset_samples(atlas1006, time == 0)
  
  # run ancombc function
  set.seed(123)
  out = ancombc(data = pseq, tax_level = "Family",
                formula = "age + nationality + bmi_group",
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
                tol = 1e-5, max_iter = 100, conserve = TRUE,
                alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
} else {
  message("The 'microbiome' package is not installed. Please install it to use this example.")
}


output <- ancombc2(
  data = ps.rare, tax_level = "Family",
  fix_formula = "Microbiome_treatment + Soil_inoculum + Density + Moisture",
  p_adj_method = "holm", pseudo_sens = TRUE,
  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
  group = "Microbiome_treatment", struc_zero = TRUE, neg_lb = TRUE,
  alpha = 0.05, n_cl = 2, verbose = TRUE,
  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE)

res_prim <- output$res
write.csv(res_prim, "results/Fungal_ANCOMBC2_results.csv")

## Structural zeros table
tab_zero <- output$zero_ind
datatable(tab_zero, caption = "Structural Zeros")
sig_taxa <- res_prim %>% filter(q_val < 0.05) %>% pull(taxon)
ps_sig <- prune_taxa(sig_taxa, ps.rel)
mat_sig <- otu_table(ps_sig)
heatmap(as.matrix(mat_sig), scale = "row", Colv = NA, Rowv = NA,
        main = "Significant Taxa Heatmap")
