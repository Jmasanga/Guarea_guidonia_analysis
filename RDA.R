
setwd("~/Metabolomics_analyses")

library(vegan)       # RDA, PCA, ecological stats
library(ggplot2)     # Plotting
library(ggrepel)     # Labeling
library(tibble)      # Data wrangling
library(dplyr)
library(tidyr)
library(patchwork)   # Combine plots

## Load data
otu_raw <- read.csv("Rarefied_otu.csv", row.names = 1, check.names = FALSE)
metabolite_df <- read.csv("metabolite_df.csv", row.names = 1, check.names = FALSE)
otu_log <- decostand(otu_raw, method = "log")
pca_result <- prcomp(metabolite_df_t, scale. = TRUE)
pc_data <- as.data.frame(pca_result$x[, 1:10])
rda_model <- rda(otu_log_t ~ ., data = pc_data)

anova(rda_model, permutations = 999)    
anova(rda_model, by = "axis", permutations = 999)
anova(rda_model, by = "term", permutations = 999)
RsquareAdj(rda_model)                             

scores_sites <- scores(rda_model, display = "sites", scaling = 2)
scores_env <- scores(rda_model, display = "bp", scaling = 2)

p_rda <- ggplot() +
  geom_point(aes(x = scores_sites[,1], y = scores_sites[,2]),
             color = "blue", size = 3) +
  geom_segment(aes(x = 0, y = 0, xend = scores_env[,1], yend = scores_env[,2]),
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  labs(x = "RDA1", y = "RDA2",
       title = "RDA: OTUs vs Metabolite PCs") +
  theme_minimal(base_size = 14)
ggsave("results/RDA_OTUs_vs_Metabolite_PCs.png", p_rda, width = 8, height = 6)

##Metabolite loadings 
loadings_df <- as.data.frame(pca_result$rotation[, 1:10]) %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "PC", values_to = "loading") %>%
  mutate(abs_loading = abs(loading),
         Direction = ifelse(loading > 0, "Up", "Down")) %>%
  group_by(PC) %>%
  arrange(desc(abs_loading), .by_group = TRUE) %>%
  mutate(Top10 = row_number() <= 10) %>%
  ungroup()

write.csv(loadings_df,
          "results/Supplementary_Table_PC_Loadings.csv",
          row.names = FALSE)

##Top 10 contributors per PC 
plot_pc_loadings <- function(pc_name) {
  pc_data <- loadings_df %>%
    filter(PC == pc_name, Top10) %>%
    mutate(Metabolite = factor(Metabolite, levels = Metabolite[order(loading)]))
  
  ggplot(pc_data, aes(x = loading, y = Metabolite, fill = Direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(title = pc_name, x = "Loading value", y = NULL) +
    theme_classic(base_size = 14) +
    theme(panel.border = element_rect(color = "black", fill = NA))
}

pc_plots <- lapply(paste0("PC", 1:10), plot_pc_loadings)
combined_pcs <- wrap_plots(pc_plots, ncol = 2)

ggsave("results/PC_Loading_Top10_All.png", combined_pcs,
       width = 12, height = 12)