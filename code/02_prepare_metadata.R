# 02_prepare_metadata.R
# Goal:
# - Load the raw ExpressionSet for GSE28521
# - Separate expression values and sample metadata
# - Create clean variables: diagnosis, brain region
# - Keep only frontal cortex samples
# - Save processed objects for later steps

# --- Packages ---
library(GEOquery)   # understands GEO data structures
library(Biobase)    # ExpressionSet class: exprs(), pData()
library(tidyverse)  # data wrangling: mutate, filter, etc.

# --- 1. Load the raw ExpressionSet we saved earlier ---
eset <- readRDS("data_raw/GSE28521_eset.rds")

# Expression matrix: genes x samples
expr <- exprs(eset)

# Sample metadata (phenotype data): one row per sample
pheno <- pData(eset)

# Quick checks
dim(expr)   # how many genes, how many samples?
dim(pheno)  # should have the same number of columns/rows as expr samples

# Look at what metadata columns we have
colnames(pheno)
head(pheno$title)
head(pheno$source_name_ch1)
head(pheno$characteristics_ch1)


# --- 2. Build clean diagnosis and region variables ---
# We will scan text columns for keywords like "autism", "control", "frontal", etc.

# (a) First, combine all character columns into one long text field per sample
pheno_df <- as.data.frame(pheno)

char_cols <- pheno_df %>%
  select(where(is.character))

text_all <- apply(char_cols, 1, paste, collapse = " ")

pheno_clean <- pheno_df %>%
  mutate(
    sample_id = rownames(pheno_df),
    text_all = tolower(text_all),
    
    # Diagnosis: ASD vs Control
    diagnosis = case_when(
      str_detect(text_all, "autism") ~ "ASD",
      str_detect(text_all, "control") ~ "Control",
      TRUE ~ NA_character_
    ),
    
    # Brain region: frontal / temporal / cerebellum
    region = case_when(
      str_detect(text_all, "frontal") ~ "Frontal",
      str_detect(text_all, "temporal") ~ "Temporal",
      str_detect(text_all, "cerebell") ~ "Cerebellum",
      TRUE ~ NA_character_
    )
  )

# See how many samples we have by diagnosis and region
table(pheno_clean$diagnosis, pheno_clean$region, useNA = "ifany")


# --- 3. Keep only frontal cortex samples with clear diagnosis ---
pheno_frontal <- pheno_clean %>%
  filter(region == "Frontal", !is.na(diagnosis))

# Subset the expression matrix to these frontal samples
expr_frontal <- expr[, pheno_frontal$sample_id]

# Double-check dimensions
dim(expr_frontal)
nrow(pheno_frontal)


# --- 4. Save processed data objects for later steps ---
if (!dir.exists("data_processed")) dir.create("data_processed")

saveRDS(pheno_frontal, "data_processed/pheno_frontal.rds")
saveRDS(expr_frontal,  "data_processed/expr_frontal.rds")