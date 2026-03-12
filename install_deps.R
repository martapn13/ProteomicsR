# Run this script once to install all dependencies before using ProteomicsR

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c(
  "limma",
  "edgeR",
  "NOISeq",
  "RUVSeq",
  "SummarizedExperiment",
  "Biobase",
  "BiocSingular",
  "ComplexHeatmap",
  "BiocParallel",
  "DelayedArray"
), ask = FALSE)

# CRAN packages
install.packages(c(
  "missForest",
  "caret",
  "doParallel",
  "ggplot2",
  "ggfortify",
  "pheatmap",
  "Rdimtools",
  "ruv",
  "remotes"
))

# GitHub-only packages
remotes::install_github("AbhishekSinha28/tcgaCleaneR")

message("All dependencies installed successfully!")
