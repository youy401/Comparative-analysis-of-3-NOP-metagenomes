# 📘 MaAsLin2 analysis for archaeal taxa

# Description:
# This script performs differential abundance analysis using MaAsLin2

# Input:
# - feature_abundance: Archaeal taxonomic abundance data
# - metadata: Sample metadata with variables such as Dose, Period, time
# - Input files should be Excel format with sheets named appropriately (e.g., "metadata", "abundance")

# Model specification:
# - Fixed effects: Dose, Period, time
# - Normalization: NONE
# - Transformation: NONE
# - Analysis method: CPLM

# Optional diagnostic:
# - Zero inflation: Proportion of zeros across taxa
# - Overdispersion: Variance-to-mean ratio per taxon

# Output:
# - MaAsLin2 results saved to "your_output_folder/"
# - Histograms saved as PNG:
#     - "abundance_histogram.png"
#     - "feature_histograms.png"

# === Load packages ===
library(dplyr)
library(ggplot2)
library(Maaslin2)
library(openxlsx)
library(tidyverse)
library(readxl)

# === Set working directory ===
setwd("your_directory_path")

# === Load data ===
file_path <- "your_input_file.xlsx"

metadata <- read_excel(file_path, sheet = "metadata")
abundance <- read_excel(file_path, sheet = "abundance")

metadata_df <- as.data.frame(metadata)
abundance_df <- as.data.frame(abundance)

rownames(abundance_df) <- abundance_df$ID
abundance_df$ID <- NULL

rownames(metadata_df) <- metadata_df$ID
metadata_df$ID <- NULL

# === Optional: Check zero inflation and overdispersion ===
zero_proportions <- colMeans(abundance_df == 0)
mean(zero_proportions)  # Mean proportion of zeros

variance_mean_ratios <- apply(abundance_df, 2, function(x) var(x)/mean(x))
mean(variance_mean_ratios)  # Mean variance-to-mean ratio

# === Visualize abundance distribution ===
png("abundance_histogram.png", width=800, height=600)
hist(as.matrix(abundance_df), breaks=50, main="Distribution of Abundance Values")
dev.off()

png("feature_histograms.png", width=1200, height=900, res=150)
par(mfrow=c(2,2))
for(i in 1:ncol(abundance_df)) {
  hist(abundance_df[,i],
       main=colnames(abundance_df)[i],
       xlab="Abundance",
       breaks=20)
}
dev.off()

# === Run MaAsLin2 analysis ===
fit_data <- Maaslin2(
  input_data = abundance_df,
  input_metadata = metadata_df,
  output = "your_output_folder",
  fixed_effects = c("Dose", "Period", "time"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "CPLM"
)
