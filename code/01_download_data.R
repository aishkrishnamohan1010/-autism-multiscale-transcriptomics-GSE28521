# 01_download_data.R
# Download GSE28521 dataset from GEO and save as an R object

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages (run once)
BiocManager::install(c("GEOquery", "Biobase"), ask = FALSE)

# Load libraries
library(GEOquery)
library(Biobase)

# Create data_raw folder if it doesn't exist
if (!dir.exists("data_raw")) dir.create("data_raw")

# Download GSE28521 (this needs internet access)
gse <- getGEO("GSE28521", GSEMatrix = TRUE)

# Usually it returns a list; we take the first ExpressionSet
length(gse)
eset <- gse[[1]]

# Save ExpressionSet as an .rds file for later use
saveRDS(eset, file = "data_raw/GSE28521_eset.rds")
