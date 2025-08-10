# Chemodiversity, LDA, and Metabolite Enrichment Analysis – Guarea guidonia

**Author:** Joel Masanga  
**Description:** This repository contains R scripts for chemodiversity analysis, 
Linear Discriminant Analysis (LDA), and LIMMA-based metabolite enrichment for G. guidonia samples, 
along with metabolite clustering, fungal community analysis, and metabolite–OTU correlation.

---

###Overview

This analysis integrates:
- **Chemodiversity analysis** using the `chemodiv` package
- **Linear Discriminant Analysis (LDA)** for treatment classification
- **LIMMA** for metabolite enrichment
- **Cluster analysis** of metabolite features
- **Fungal community and enrichment analyses** via `phyloseq` and `ANCOMBC2`
- **Redundancy Analysis (RDA)** linking fungal OTUs with metabolite principal components

Results are written to the `/results` directory.

---

##Repository Structure

```
.
├── chemodiversity_LDA_limma_analysis.R  # Main analysis script
├── data/                                # Input CSV data files
│   ├── AnnotationMetab.csv
│   ├── GuaguiCompData.csv
│   ├── metadata.csv
│   ├── Guagui_Normalized_Metabolites.csv
│   ├── otu_table_reheader.csv
│   ├── taxa.csv
│   ├── merge.csv
│   ├── Rarefied_otu.csv
│   └── metabolite_transposed.csv
├── results/                             # Output CSV tables, plots, and figures
└── README.md                            # This file
```

---

## Requirements

**R version:** >= 4.0.0  

**Required packages:**
```r
install.packages(c(
  "chemodiv", "vegan", "dplyr", "MASS", "limma", "cluster",
  "tidyverse", "Hmisc", "robustbase", "phyloseq", "ggplot2",
  "ANCOMBC", "DT", "tidyr", "ggrepel", "tibble", "patchwork"
))
```

---

##How to Run

1. Clone the repository:
```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
```

2. Ensure the `data/` directory contains all required input files.

3. Open R or RStudio and run:
```r
source("chemodiversity_LDA_limma_analysis.R")
```

4. All outputs (tables, plots) will be saved in the `results/` directory.

---

##Output Files

**Example outputs include:**
- `limma_results_Soil_fungicide_vs_Control.csv`  
- `Fungal_alpha_diversity.csv`  
- `Fungal_Phylum_Composition.png`  
- `MetaboliteClustersTrait_correlations.csv`  
- `RDA_OTUs_vs_Metabolite_PCs.png`  
- `PC_Loading_Top10_All.png`

