# 📘 ANCOM-BC2: Differential abundance testing at the Genus level

# Description:
# This script performs differential abundance analysis using ANCOM-BC2 at the Genus level.

# Input:
# - Taxonomic abundance (absolute count) data
# - Sample metadata
# - Taxonomic annotations
# - All data are formatted into a phyloseq object

# Study details:
# - The 'Dose' variable (levels: 0, 53, 161, 345) indicates 3-NOP supplementation levels.
# - This dataset corresponds to the Beef1 study:
#   Romero-Perez et al. (2014), *The potential of 3-nitrooxypropanol to lower enteric methane emissions from beef cattle*
# - 'Dose' is treated as a categorical variable with 0 (control) as the reference level.

# Analysis:
# - ANCOM-BC2 is applied using the formula: Dose + Period + time
# - P-values are adjusted using the Benjamini-Hochberg (BH) method
# - Results are exported to an Excel file


# Required libraries
library(ANCOMBC)
library(phyloseq)
library(readxl)
library(writexl)

# Set working directory
setwd("your/working/directory")

# Load metadata and convert dose to factor with 0 as reference level
metadata <- read_excel("your_input_file.xlsx", sheet = "metadata")
metadata$Dose <- factor(metadata$Dose, levels = c(0, 53, 161, 345))

# Load taxonomy and count data (absolute abundance)
taxonomy <- read_excel("your_input_file.xlsx", sheet = "Taxonomy")
counts <- read_excel("your_input_file.xlsx", sheet = "Sample")

# Convert count data to matrix format
otu_mat <- as.matrix(counts[,-1])
rownames(otu_mat) <- counts[[1]]

# Convert taxonomy to matrix format
tax_mat <- as.matrix(taxonomy[,-1])
rownames(tax_mat) <- taxonomy[[1]]

# Format metadata
sample_df <- data.frame(metadata)
rownames(sample_df) <- sample_df$ID

# Construct phyloseq object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(sample_df)
ps <- phyloseq(OTU, TAX, samples)

# Run ANCOM-BC2
ancombc_result <- ancombc2(
  data = ps,
  fix_formula = "Dose + Period + time",  # adjust covariates as needed
  p_adj_method = "BH",
  struc_zero = TRUE,
  tax_level = "Genus",
  group = "Dose"
)

# Save result
res <- ancombc_result$res
write_xlsx(res, "Dose_BH_ancombc2_results.xlsx")
