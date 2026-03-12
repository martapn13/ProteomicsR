# ProteomicsR

An R pipeline for multi-organ proteomics analysis, including imputation, batch correction (RUViii-PRPS), differential expression (limma/voom), and QC visualizations.

Original pipeline by Leong Ng. Modified and packaged by Marta Nascimento.

## Installation

### 1. Install Bioconductor packages
```r
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "NOISeq", "RUVSeq", 
                       "SummarizedExperiment", "Biobase",
                       "BiocSingular", "ComplexHeatmap", 
                       "BiocParallel", "DelayedArray"))
```

### 2. Install CRAN packages
```r
install.packages(c("missForest", "caret", "doParallel", "ggplot2", 
                   "ggfortify", "pheatmap", "Rdimtools", "ruv"))
```

### 3. Install tcgaCleaneR from GitHub
```r
install.packages("remotes")
remotes::install_github("AbhishekSinha28/tcgaCleaneR")
```

### 4. Install ProteomicsR
```r
remotes::install_github("martapn13/ProteomicsR")
```

## Usage

### Tissue analysis
```r
library(ProteomicsR)

setwd("path/to/your/data")

run_tissue_analysis(
  pheno_file  = "phenoAdrenals.txt",
  exprs_file  = "exprsAdrenals.txt",
  tissue_name = "Adrenals"
)
```

### Plasma analysis
```r
run_plasma_analysis(
  pheno_file = "phenoPlasma.txt",
  exprs_file = "exprsPlasma.txt",
  raw_file   = "Plasma_Raw.txt"
)
```

## Pipeline Steps

1. **Data loading** — phenotype and expression matrices
2. **NZV filtering** — removes near-zero variance features
3. **Imputation** — half-minimum of zero values
4. **SPECU ranking** — identifies negative control features
5. **RUViii-PRPS** — batch correction using replicate samples (tissue only)
6. **ARSyNseq** — additional noise removal
7. **limma/voom** — differential expression (KO vs WT, DAPA vs WT, KO vs DAPA)
8. **Outputs** — volcano plots, heatmaps, and CSV result tables

## Output files

- `ExprsprpsNoiseq_<tissue>.csv` — corrected expression matrix
- `Res_noiseq_<tissue>.csv` — full limma results
- `<tissue>_sig_KO_vs_WT.csv`, `<tissue>_sig_DAPA_vs_WT.csv`, `<tissue>_sig_KO_vs_DAPA.csv` — significant DE proteins
- `volcano_<tissue>_<comparison>.png` — volcano plots

Plasma additionally saves:
- `sig_*_Plasma_withID.csv` — significant DE proteins with protein IDs mapped

## Dependencies

| Package | Source |
|---|---|
| limma, edgeR, NOISeq, RUVSeq, SummarizedExperiment, Biobase, BiocSingular, ComplexHeatmap | Bioconductor |
| missForest, caret, doParallel, ggplot2, ggfortify, pheatmap, Rdimtools, ruv | CRAN |
| tcgaCleaneR | GitHub (AbhishekSinha28/tcgaCleaneR) |

## License

MIT
## Dependencies

| Package | Source |
|---|---|
| limma, edgeR, NOISeq, RUVSeq, SummarizedExperiment | Bioconductor |
| missForest, caret, doParallel, ggplot2, ggfortify, pheatmap, Rdimtools, ruv | CRAN |
| tcgaCleaneR | GitHub |

## License
MIT
