# 03_QC_EDA.R
# QC = Quality Control
# EDA = Exploratory Data Analysis
# Goal: Examine data distributions, identify outliers, and perform PCA.

library(tidyverse)
library(limma)  # for normalization and visualization

# Load processed expression data and metadata
expr_frontal <- readRDS("data_processed/expr_frontal.rds")
pheno_frontal <- readRDS("data_processed/pheno_frontal.rds")

# -------------------------------
# 1. Check basic dimensions
# -------------------------------
dim(expr_frontal)
table(pheno_frontal$diagnosis)

# -------------------------------
# 2. Boxplot of raw expression values
# -------------------------------
pdf("results/figures/boxplot_raw_expression.pdf", width = 10, height = 6)
boxplot(expr_frontal,
        main = "Expression Distributions (Raw)",
        las = 2, outline = FALSE,
        col = "lightblue")
dev.off()

# -------------------------------
# 3. Normalize expression (important!)
# -------------------------------
expr_norm <- normalizeBetweenArrays(expr_frontal, method = "quantile")

# Save normalized expression
saveRDS(expr_norm, "data_processed/expr_frontal_normalized.rds")

# Boxplot after normalization
pdf("results/figures/boxplot_normalized_expression.pdf", width = 10, height = 6)
boxplot(expr_norm,
        main = "Expression Distributions (Normalized)",
        las = 2, outline = FALSE,
        col = "lightgreen")
dev.off()

# -------------------------------
# 4. PCA (Principal Component Analysis)
# -------------------------------
expr_t <- t(expr_norm)  # PCA needs samples as rows

pca <- prcomp(expr_t, scale. = TRUE)

pca_df <- data.frame(
  Sample = pheno_frontal$sample_id,
  Diagnosis = pheno_frontal$diagnosis,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# PCA Plot
pdf("results/figures/pca_frontal_ASD_vs_Control.pdf", width = 6, height = 5)
ggplot(pca_df, aes(PC1, PC2, color = Diagnosis)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Frontal Cortex Samples") +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "% variance)"))
dev.off()

