# ProteomicsR

An R pipeline for multi-organ proteomics analysis, including half-minimum imputation, batch correction (RUViii-PRPS), differential expression (limma/voom), UniProt-to-gene symbol conversion, and QC visualizations.

Original pipeline by Leong Ng. Modified and packaged by Marta Nascimento.

---

## Installation

### 1. Install Bioconductor packages

```r
install.packages("BiocManager")
BiocManager::install(c(
  "limma", "edgeR", "NOISeq", "RUVSeq",
  "SummarizedExperiment", "Biobase",
  "AnnotationDbi", "org.Mm.eg.db", "UniProt.ws"
))
```

### 2. Install CRAN packages

```r
install.packages(c(
  "caret", "doParallel", "ggplot2",
  "ggfortify", "pheatmap", "Rdimtools", "ruv"
))
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

---

## Usage

### Tissue analysis (Liver, Adrenals, BAT, WAT, Kidney, Muscle, Spleen)

```r
library(ProteomicsR)

run_tissue_analysis(
  pheno_file  = "phenoAdrenals.txt",
  exprs_file  = "exprsAdrenals.txt",
  tissue_name = "Adrenals",
  best_ncomp  = 5
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

---

## Pipeline Steps

### Tissue pipeline (`run_tissue_analysis`)
1. **Data loading** — phenotype and expression matrices
2. **NZV filtering** — removes near-zero variance features
3. **Half-minimum imputation** — replaces zero values with half the minimum detected value
4. **SPECU ranking** — identifies negative control features for RUViii
5. **RUViii-PRPS** — batch correction using replicate samples
6. **ARSyNseq/TMM normalisation** — additional noise removal
7. **limma/voom** — differential expression (KO vs WT, DAPA vs WT, KO vs DAPA)
8. **Gene symbol conversion** — UniProt accessions mapped to gene symbols via `org.Mm.eg.db` (with `UniProt.ws` fallback)
9. **Outputs** — volcano plots, heatmaps with gene labels, and CSV result tables

### Plasma pipeline (`run_plasma_analysis`)
Steps 1–3 and 6–9 as above. RUViii-PRPS and SPECU are not applied to plasma.
Protein IDs are already MGI symbols in the raw file — no additional conversion required.

---

## Output files

### Tissue
- `ExprsprpsNoiseq_<tissue>.csv` — corrected expression matrix
- `Res_noiseq_<tissue>.csv` — full limma results
- `<tissue>_sig_KO_vs_WT.csv`, `<tissue>_sig_DAPA_vs_WT.csv`, `<tissue>_sig_KO_vs_DAPA.csv` — significant DEPs (UniProt row names)
- `<tissue>_sig_KO_vs_WT_MGI.csv`, `<tissue>_sig_DAPA_vs_WT_MGI.csv`, `<tissue>_sig_KO_vs_DAPA_MGI.csv` — significant DEPs with `GeneName` and `UniProt` columns
- `volcano_<tissue>_<comparison>.png` — volcano plots

### Plasma
- `Res2_noiseq2_plasma.csv` — full limma results
- `sig_*_Plasma.csv` — significant DEPs
- `sig_*_Plasma_withID.csv` — significant DEPs with MGI gene symbols mapped
- `volcano_*_full.png` — volcano plots

---

## Dependencies

| Package | Source |
|---|---|
| limma, edgeR, NOISeq, RUVSeq, SummarizedExperiment, Biobase, AnnotationDbi, org.Mm.eg.db, UniProt.ws | Bioconductor |
| caret, doParallel, ggplot2, ggfortify, pheatmap, Rdimtools, ruv | CRAN |
| tcgaCleaneR | GitHub (AbhishekSinha28/tcgaCleaneR) |

---

## License

MIT
