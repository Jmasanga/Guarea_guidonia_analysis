Legacy of warming and microbial treatments shape root exudates and rhizosphere fungal communities of a tropical tree

This repository contains 4 R scripts used for analysis of data in the article "Legacy of warming and microbial treatments shape root exudates and rhizosphere fungal communities of a tropical tree".
Fungal_diversity.R contains scripts used to analyse Fungal diversity, community composition and differential abundance

Metabolite_diversity.R contains R code used to analyze diversity and composition of root exudate metabolite features
The analyses integrate:
**Chemodiversity analysis** using the `chemodiv` package
**Linear Discriminant Analysis (LDA)** for treatment classification
**LIMMA** for metabolite enrichment

Metabolite_clusters.R contains R code used to categorize metabolite features into clusters.  
**Cluster analysis** of metabolite features uses the 'cluster' package adopted from Lin & Peddada 2020 (https://www.nature.com/articles/s41467-020-17041-7)

Fungal_diversity.R contains R code for analysis of fungal diversity, community composition and differential abundance
Analysis adopts `phyloseq` and `ANCOMBC2` packages

RDA.R contains R code used for linking fungal OTUs with metabolite principal components

