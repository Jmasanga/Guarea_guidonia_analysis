setwd("~/Metabolomics_analyses")

library(chemodiv)
library(vegan)
library(dplyr)
library(MASS)
library(limma)

AnnotationMetab <- read.csv("data/AnnotationMetab.csv")
GuaguiCompData <- read.csv("data/GuaguiCompData.csv", row.names = 1)
metadata <- read.csv("data/Guarea_metadata.csv")


chemoDivCheck(compoundData = AnnotationMetab, sampleData = GuaguiCompData)

GuaguiNPC <- NPCTable(compoundData = AnnotationMetab)
print(GuaguiNPC[1, ])

GuaguiCompData[GuaguiCompData == 0] <- 1e-6

GuaguiCompDis <- compDis(compoundData = AnnotationMetab, type = "PubChemFingerprint")
compDisMat <- GuaguiCompDis$fingerDisMat[[1]]

GuaguiSampDis <- sampDis(sampleData = GuaguiCompData,
                         compDisMat = compDisMat,
                         type = "GenUniFrac")
adonis2(GuaguiSampDis$GenUniFrac ~ Density, data = metadata,
        method = "euclidean", permutations = 999)

GuaguiDivrichness <- calcDiv(sampleData = GuaguiCompData,
                             compDisMat = compDisMat,
                             type = "FuncHillDiv", q = 0)
GuaguiDiv <- calcDiv(sampleData = GuaguiCompData,
                     compDisMat = compDisMat,
                     type = "FuncHillDiv", q = 1)

GuaguiDivProf <- calcDivProf(sampleData = GuaguiCompData,
                             compDisMat = compDisMat,
                             type = "FuncHillDiv")
print(head(GuaguiDivProf$divProf)[, 1:5])

GuaguiBetaDiv <- calcBetaDiv(sampleData = GuaguiCompData,
                             compDisMat = compDisMat,
                             type = "FuncHillDiv")

# LDA
Guagui_df <- read.csv("data/Guagui_Normalized_Metabolites.csv")
Guarea_Metadata <- metadata

Guagui_df$combined_class <- interaction(Guagui_df$Soil_inoculum, 
                                        Guagui_df$Microbiome_treatment, 
                                        Guagui_df$Density, 
                                        Guagui_df$Moisture)

lda_result <- lda(combined_class ~ ., data = Guagui_df)
lda_pred <- predict(lda_result)
lda_scores <- as.data.frame(lda_pred$x)

lda_scores$Microbiome_treatment <- Guarea_Metadata$Microbiome_treatment
lda_scores$Soil_inoculum <- Guarea_Metadata$Soil_inoculum
lda_scores$Density <- Guarea_Metadata$Density
lda_scores$Moisture <- Guarea_Metadata$Moisture

lda_scores$Shape <- factor(
  ifelse(lda_scores$Soil_inoculum == "Ambient" & lda_scores$Moisture == "WS", "Ambient_WS",
         ifelse(lda_scores$Soil_inoculum == "Ambient" & lda_scores$Moisture == "WW", "Ambient_WW",
                ifelse(lda_scores$Soil_inoculum == "Warmed" & lda_scores$Moisture == "WS", "Warmed_WS",
                       ifelse(lda_scores$Soil_inoculum == "Warmed" & lda_scores$Moisture == "WW", "Warmed_WW", NA)))))

# Metabolite Enrichment Analysis
metab_matrix <- read.csv("data/Guagui_Normalized_Metabolites.csv", row.names = 1)
metadata$Microbiome_treatment <- factor(metadata$Microbiome_treatment,
                                        levels = c("Control", "Soil_fungicide", "Foliar_fungicide", "Heat_treatment"))

design <- model.matrix(~ Microbiome_treatment, data = metadata)
colnames(design) <- gsub("Microbiome_treatment", "", colnames(design))

fit <- lmFit(t(metab_matrix), design)
fit <- eBayes(fit)

results_soil <- topTable(fit, coef = "Soil_fungicide", number = Inf, adjust = "fdr")
results_foliar <- topTable(fit, coef = "Foliar_fungicide", number = Inf, adjust = "fdr")
results_heat <- topTable(fit, coef = "Heat_treatment", number = Inf, adjust = "fdr")

dir.create("results", showWarnings = FALSE)
write.csv(results_soil, "results/limma_results_Soil_fungicide_vs_Control.csv")
write.csv(results_foliar, "results/limma_results_Foliar_fungicide_vs_Control.csv")
write.csv(results_heat, "results/limma_results_Heat_treatment_vs_Control.csv")

volcanoplot(fit, coef = "Soil_fungicide", highlight = 10,
            main = "Soil fungicide vs Control", names = rownames(metab_matrix))
