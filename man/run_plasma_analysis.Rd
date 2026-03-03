\name{run_plasma_analysis}
\alias{run_plasma_analysis}
\title{Plasma Proteomics Analysis Pipeline}
\usage{
run_plasma_analysis(
  pheno_file,
  exprs_file,
  raw_file
)
}
\arguments{
\item{pheno_file}{Path to the phenotype .txt file (tab-separated). Must contain
one row per sample with the following columns:
\itemize{
  \item \code{phenostate} - experimental group (e.g. "WT", "KO", "DAPA", "QC")
  \item \code{sampleName} - unique sample identifier (used as row names)
  \item \code{sampleQC} - sample type, either "Sample" or "QC"
  \item \code{Replicates} - replicate label; matches sampleName for biological
    samples and "QC" for quality control samples
}
QC samples (where \code{Replicates == "QC"}) are automatically removed before analysis.
}
\item{exprs_file}{Path to the expression .txt file (tab-separated). Must be
formatted as:
\itemize{
  \item First column = protein identifiers as row names (e.g. UniProt IDs such as "Q9JHZ2" or protein groups separated by semicolons such as "O88746;O88746-2")
  \item Remaining columns = samples (column names must match \code{sampleName} in the pheno file)
  \item Values = raw protein intensity (non-log scale, non-negative)
  \item Missing values should be coded as 0
  \item Includes both biological samples and quality control samples (e.g. "Quality_control_1", "Quality_control_2")
}
Example column structure: Quality_control_1, Quality_control_2, KO_single_sample_4, WT_single_sample_1, DAPA_single_sample_7, etc.
}
\item{raw_file}{Path to the raw plasma .txt file for ID mapping.}
}
\value{A named list with df_norm, efit, and deg_list.}
\description{Runs a full proteomics analysis pipeline for plasma samples,
using half-minimum imputation and ARSyNseq normalisation (no RUViii-PRPS).}
\examples{
\dontrun{
run_plasma_analysis(
  pheno_file = "phenoPlasma.txt",
  exprs_file = "exprsPlasma.txt",
  raw_file   = "Plasma_Raw.txt"
)
}
}
