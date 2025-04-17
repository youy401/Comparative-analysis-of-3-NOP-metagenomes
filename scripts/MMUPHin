# 📘 MMUPHin + MaAsLin2: Meta-analysis and differential abundance modeling

# Description:
# This script performs multi-study normalization and meta-analysis using MMUPHin,
# and compares results with a pooled differential abundance analysis using MaAsLin2.

# Input:
# - your_metadata.txt: Sample metadata (tab-delimited, includes studyID, dose, period, etc.)
# - your_abundance_data.txt: Normalized feature abundance data (e.g., CPM)
#   Includes features such as hydrogenase, reductase, or other target gene counts

# Analysis:
# - Batch effect correction across studies using `adjust_batch()` from MMUPHin
# - ADONIS (PERMANOVA) to assess batch structure pre/post adjustment
# - Meta-analysis using `lm_meta()` (CPLM method)
# - Forest plots and export of meta-analysis fits
# - Pooled MaAsLin2 analysis using adjusted abundance data

# Notes:
# - `dose` is treated as the exposure variable
# - Batch is set to `studyID`
# - No normalization or transformation is applied (data already normalized to CPM)

# Output:
# - "meta_fits_output.xlsx": Meta-analysis results
# - "forest_plot.pdf": Forest plot of meta-analysis effect sizes
# - "maaslin2_output/": Pooled MaAsLin2 output directory

# ============================================
# 📦 Load libraries
# ============================================
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(Maaslin2)
library(openxlsx)
library(vegan)

# ============================================
# 📁 Set working directory and load data
# ============================================
# Set your working directory as needed
setwd("your/path/to/project")

# Input files
metadata <- read.table("your_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
abundance <- read.table("your_abundance_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
rownames(abundance) <- abundance[,1]
abundance <- abundance[,-1]

metadata$studyID <- as.factor(metadata$studyID)

# ============================================
# ⚙️ Batch correction using MMUPHin
# ============================================
fit_adjust_batch <- adjust_batch(
  feature_abd = abundance,
  batch = "studyID",
  covariates = c("study_condition"),
  data = metadata,
  control = list(verbose = FALSE)
)

abundance_adj <- fit_adjust_batch$feature_abd_adj

# ============================================
# 📊 ADONIS (PERMANOVA) before and after batch correction
# ============================================
D_before <- vegdist(t(abundance))
D_after <- vegdist(t(abundance_adj))

set.seed(1)
fit_adonis_before <- adonis(D_before ~ studyID, data = metadata, permutations = 999)
fit_adonis_after <- adonis(D_after ~ studyID, data = metadata, permutations = 999)

# ============================================
# 📈 Meta-analysis across studies (study-wise modeling)
# ============================================
fit_lm_meta <- lm_meta(
  feature_abd = abundance_adj,
  batch = "studyID",
  exposure = c("dose"),
  covariates = c("breed", "period", "forage", "concentrate"),
  data = metadata,
  control = list(
    normalization = "NONE",
    transform = "NONE",
    analysis_method = "CPLM",
    output = "meta_analysis_output",
    forest_plot = "forest_plot.pdf",
    verbose = TRUE
  )
)

write.xlsx(fit_lm_meta, "meta_fits_output.xlsx", rowNames = FALSE)

# Plot significant features (qval.fdr < 0.05)
fit_lm_meta$meta_fits %>% 
  filter(qval.fdr < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()

# ============================================
# 📉 Maaslin2 analysis (global pooled model)
# ============================================
fit_maaslin2 <- Maaslin2(
  input_data = abundance_adj,
  input_metadata = metadata,
  output = "maaslin2_output",
  fixed_effects = c("dose"),
  random_effects = c("breed", "period"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "CPLM",
  min_prevalence = 0.1,
  min_abundance = 0.001
)
