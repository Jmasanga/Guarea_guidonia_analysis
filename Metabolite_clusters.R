setwd("~/Metabolomics_analyses")

library(cluster)
library(tidyverse)
library(Hmisc)
library(robustbase)

###Helper functions
test_k <- function(dist_matrix, k_init = 3, k_final = 10) {
  for (k in k_init:k_final) {
    clustering <- pam(dist_matrix, k = k)
    avs <- mean(silhouette(clustering)[, 3])
    message("k = ", k, "; AVS = ", round(avs, 3))
  }
}

cluster_rep_features <- function(dist_matrix, k) {
  final_clusters <- pam(dist_matrix, k = k)
  cluster_df <- tibble(
    featureid = as.character(names(final_clusters$clustering)),
    cluster = factor(final_clusters$clustering))
  sil_df <- tibble(
    featureid = as.character(rownames(silhouette(final_clusters))),
    cluster = factor(silhouette(final_clusters)[, 1]),
    sil_width = silhouette(final_clusters)[, 3])
  reps <- sil_df %>% group_by(cluster) %>% filter(sil_width > mean(sil_width)) %>% ungroup()
  list(all_features = cluster_df, representatives = reps)
}

cluster_medians <- function(df, cluster, env_df) {
  sample_ids <- intersect(env_df$SampleID, colnames(df))
  df_tem <- df %>% select(-cluster) %>% column_to_rownames("featureid") %>% select(all_of(sample_ids))
  medians <- robustbase::colMedians(as.matrix(df_tem), na.rm = TRUE)
  tibble(cluster = paste0("Cluster", cluster), !!!as.list(medians))
}

env_cor <- function(env_var, cluster_df, env_df) {
  cluster_df <- cluster_df %>% column_to_rownames("cluster") %>% select(all_of(env_df$SampleID))
  corr <- rcorr(t(cluster_df), env_df[[env_var]], type = "spearman")
  tibble(
    cluster = rownames(cluster_df),
    env_variable = env_var,
    rho = corr$r[1:nrow(cluster_df), "y"],
    pvalue = corr$P[1:nrow(cluster_df), "y"],
    pval_adj = p.adjust(corr$P[1:nrow(cluster_df), "y"]))
}

matrix_cor <- function(cluster_df, other_df, type) {
  cluster_df <- cluster_df %>% column_to_rownames("cluster") %>% select(all_of(colnames(other_df)))
  if (type == "Traits") {
    sel_rows <- rowSums(!is.na(other_df)) > ceiling(ncol(other_df) / 2)
  } else if (type == "abundance") {
    sel_rows <- rowSums(other_df > 0) > ceiling(ncol(other_df) / 2)
  }
  other_df_filt <- other_df[sel_rows, ]
  corr <- rcorr(t(cluster_df), t(other_df_filt), type = "spearman")
  corr_r <- as.data.frame(corr$r) %>% rownames_to_column("Feature") %>% pivot_longer(-Feature, names_to = "cluster", values_to = "rho") %>% drop_na()
  corr_p <- as.data.frame(corr$P) %>% rownames_to_column("Feature") %>% pivot_longer(-Feature, names_to = "cluster", values_to = "pvalue") %>% drop_na()
  corr_p_adj <- p.adjust(corr$P, method = "fdr") %>% matrix(nrow = nrow(corr$P), dimnames = dimnames(corr$P)) %>% as.data.frame() %>% rownames_to_column("Feature") %>% pivot_longer(-Feature, names_to = "cluster", values_to = "p_adj") %>% drop_na()
  corr_r %>% inner_join(corr_p, by = c("Feature", "cluster")) %>% inner_join(corr_p_adj, by = c("Feature", "cluster"))
}

###Data filtering
filter_half_presence <- function(df) {
  df[rowSums(!is.na(df)) > 0.5 * ncol(df), ]
}

Soil_fungicide_filt <- filter_half_presence(Soil_fungicide_df1)
Heat_treatment_filt <- filter_half_presence(Heat_treatment_df1)
Control_filt <- filter_half_presence(Control_df)
Foliar_fungicide_filt <- filter_half_presence(Foliar_fungicide_df)

##Distance matrices
Soil_fungicide_dist <- daisy(Soil_fungicide_filt, metric = "manhattan")
Heat_treatment_dist <- daisy(Heat_treatment_filt, metric = "manhattan")
Control_dist <- daisy(Control_filt, metric = "manhattan")
Foliar_fungicide_dist <- daisy(Foliar_fungicide_filt, metric = "manhattan")

##Optimal k
test_k(Soil_fungicide_dist)
test_k(Control_dist)
test_k(Heat_treatment_dist)
test_k(Foliar_fungicide_dist)

##Clustering & representatives
k_opt <- 6
Soil_fungicide_clusters <- cluster_rep_features(Soil_fungicide_dist, k_opt)
rep_features <- Soil_fungicide_clusters$representatives %>% left_join(Soil_fungicide_meta, by = c("featureid" = "Sample"))
write_csv(rep_features, "Soil_fungicide_cluster_repFeatures.csv")

rep_features_Soil_fungicide <- Soil_fungicide_filt %>% rownames_to_column("featureid") %>% right_join(rep_features, by = "featureid")
write_csv(rep_features_Soil_fungicide, "Cluster_analysis_rep_features_Soil_fungicide.csv")

##Consensus cluster values
rep_values <- rep_features_Soil_fungicide %>% split(.$cluster)
consensus_list <- imap(rep_values, ~ cluster_medians(.x, .y, Traits))
consensus_df <- bind_rows(consensus_list)

##Correlations with traits
env_vars <- colnames(Traits)[2:8]
corr_env_df <- map_dfr(env_vars, ~ env_cor(.x, consensus_df, Traits))
write_csv(corr_env_df, "MetaboliteClustersTrait_correlations.csv")

##Optimized consensus
consensus_clusters <- Cluster_analysis_rep_features_Traits %>%
  pivot_longer(-c(featureid, cluster), names_to = "SampleID", values_to = "Abundance") %>%
  group_by(cluster, SampleID) %>%
  summarise(ConsensusValue = median(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = ConsensusValue)
