# ProteomicsR
Multi-organ proteomics analysis pipeline
# ProteomicsR

An R pipeline for multi-organ proteomics analysis, including imputation, batch correction (RUViii-PRPS), differential expression (limma/voom), and QC visualizations.

## Installation

```r
# Install dependencies first
source("install_deps.R")

# Then install this package
remotes::install_github("martapn13/ProteomicsR")
```

## Usage

```r
library(myProteomicsR)

run_tissue_analysis(
  pheno_file  = "phenoAdrenals2.txt",
  exprs_file  = "exprsAdrenals2.txt",
  tissue_name = "Adrenals",
  use_halfmin = FALSE,
  best_ncomp  = 5
)
```

## Pipeline Steps

1. **Data loading** — phenotype and expression matrices
2. **NZV filtering** — removes near-zero variance features
3. **Imputation** — missForest (default) or half-minimum
4. **SPECU ranking** — identifies negative control features
5. **RUViii-PRPS** — batch correction using replicate samples
6. **ARSyNseq** — additional noise removal
7. **limma/voom** — differential expression (KO vs WT, DAPA vs WT, KO vs DAPA)
8. **Outputs** — volcano plots, heatmaps, and CSV result tables

## Dependencies

| Package | Source |
|---|---|
| limma, edgeR, NOISeq, RUVSeq, SummarizedExperiment | Bioconductor |
| missForest, caret, doParallel, ggplot2, ggfortify, pheatmap, Rdimtools, ruv | CRAN |
| tcgaCleaneR | GitHub |

## License
MIT
