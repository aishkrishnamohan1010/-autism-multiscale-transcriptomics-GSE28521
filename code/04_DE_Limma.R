# 04_DE_limma.R
# Differential expression: ASD vs Control (Frontal cortex)
# STEP 4.1 — Load data and packages
library(limma)
library(tidyverse)

# Load processed data
expr_frontal  <- readRDS("data_processed/expr_frontal.rds")
pheno_frontal <- readRDS("data_processed/pheno_frontal.rds")


#STEP 4.2 — Make sure diagnosis is a factor (VERY important)

# Ensure correct group ordering
pheno_frontal$diagnosis <- factor(
  pheno_frontal$diagnosis,
  levels = c("Control", "ASD")
)

table(pheno_frontal$diagnosis)

# STEP 4.3 — Build the design matrix
design <- model.matrix(~ diagnosis, data = pheno_frontal)
design

#STEP 4.4 — Normalize (between arrays)

expr_norm <- normalizeBetweenArrays(expr_frontal)

#STEP 4.5 — Fit the linear model

fit <- lmFit(expr_norm, design)
fit <- eBayes(fit)



#STEP 4.6 — Extract results (this is the key output)

deg <- topTable(
  fit,
  coef = "diagnosisASD",
  number = Inf,
  adjust.method = "BH"
)

head(deg)

# STEP 4.7 — Save DEG table (portfolio-critical)

if (!dir.exists("results/tables")) {
  dir.create("results/tables", recursive = TRUE)
}

write.csv(
  deg,
  "results/tables/DEG_ASD_vs_Control_Frontal.csv",
  row.names = TRUE
)

# STEP 4.8 — How many DE genes do we have?

deg_sig <- deg %>%
  filter(adj.P.Val < 0.05)

nrow(deg_sig)


# Also check direction:
#TRUE → higher in ASD

#FALSE → lower in ASD

table(deg_sig$logFC > 0)


# STEP 4.9 — Volcano plot (visual summary)

deg <- deg %>%
   mutate(significant = ifelse(adj.P.Val < 0.05, "Significant", "Not significant"))


ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: ASD vs Control (Frontal Cortex)",
    x = "log2 Fold Change (ASD / Control)",
    y = "-log10(FDR)",
    color = ""
  )

#Save it:
ggsave(
  "results/figures/Volcano_ASD_vs_Control_Frontal.pdf",
  width = 6,
  height = 5
)

