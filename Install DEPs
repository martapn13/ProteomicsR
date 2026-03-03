# install_deps.R
# Run this script once to install all dependencies before using myProteomicsR

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c(
  "limma",
  "edgeR",
  "NOISeq",
  "RUVSeq",
  "SummarizedExperiment",
  "Biobase"
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
  "ruv"
))

# GitHub-only packages
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("martapn13/tcgaCleaneR")

message("All dependencies installed successfully!")
