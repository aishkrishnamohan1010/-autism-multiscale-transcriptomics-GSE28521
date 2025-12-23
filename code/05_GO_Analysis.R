#STEP 5.1 — Prepare gene list for enrichment
#extra significant genes
deg <- read.csv(
  "results/tables/DEG_ASD_vs_Control_Frontal.csv",
  row.names = 1
)

# Significant genes
deg_sig <- deg %>%
  filter(adj.P.Val < 0.05)

nrow(deg_sig)

#Separate up- and down-regulated genes (important!)
genes_up <- rownames(deg_sig)[deg_sig$logFC > 0]
genes_down <- rownames(deg_sig)[deg_sig$logFC < 0]

length(genes_up)
length(genes_down)

#STEP 5.2 — Install & load enrichment packages
#We’ll use clusterProfiler, a gold-standard package.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Step 5.3 We need the RIGHT annotation database
head(rownames(deg))

BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)

# Step 5.4 Map PROBE IDs → ENTREZ IDs (correct way)
library(AnnotationDbi)

probe_ids <- rownames(deg)

probe2entrez <- mapIds(
  illuminaHumanv4.db,
  keys = probe_ids,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

head(probe2entrez)


# Step 5.5 Add Entrez IDs to DEG table

deg$ENTREZID <- probe2entrez

#Remove genes without Entrez IDs:
  
deg_clean <- deg %>% filter(!is.na(ENTREZID))


# Step 5.6 Recreate up/down gene lists (now correct)
deg_sig <- deg_clean %>% filter(adj.P.Val < 0.05)

up_entrez <- deg_sig$ENTREZID[deg_sig$logFC > 0]
down_entrez <- deg_sig$ENTREZID[deg_sig$logFC < 0]

length(up_entrez)
length(down_entrez)


# Step 5.7 NOW run GO enrichment

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)

library(clusterProfiler)
ego_up <- enrichGO(
  gene          = up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_down <- enrichGO(
  gene          = down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Step 5.8 Plot

library(enrichplot)
library(ggplot2)

dotplot(ego_up, showCategory = 10) + ggtitle("GO BP: Upregulated in ASD")
ggsave("results/figures/GO_BP_ASD_Up.pdf", width = 7, height = 5)

dotplot(ego_down, showCategory = 10) + ggtitle("GO BP: Downregulated in ASD")
ggsave("results/figures/GO_BP_ASD_Down.pdf", width = 7, height = 5)


# Step 5.9 How to SEE the actual genes (very important)
ego_down@result[, c("Description", "Count", "geneID")]
ego_up@result[, c("Description", "Count", "geneID")]

# Step 5.10 To View the Full GO result table

go_table <- ego_down@result
head(go_table)


# Step 5.11 View ONLY the useful columns (recommended)

go_table_simple <- go_table[, c(
  "Description",
  "Count",
  "GeneRatio",
  "p.adjust",
  "geneID"
)]

go_table_simple


# Step 5.12 Same thing but even clearer

library(tidyr)

go_long <- go_table_simple %>%
  separate_rows(geneID, sep = "/")

go_long