\name{run_tissue_analysis}
\alias{run_tissue_analysis}
\title{Multi-Organ Proteomics Analysis Pipeline}
\usage{
run_tissue_analysis(
  pheno_file,
  exprs_file,
  tissue_name,
  use_halfmin = FALSE,
  best_ncomp = 5
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
\item{tissue_name}{Short label used in all output file names (e.g. "Adrenals").}
\item{use_halfmin}{Logical. If TRUE, uses half-minimum imputation. Default FALSE.}
\item{best_ncomp}{Integer. Number of components for RUViii. Default 5.}
}
\value{A named list with df122C, efit, and deg_list.}
\description{Runs a full proteomics analysis pipeline for a single tissue.}
\examples{
\dontrun{
run_tissue_analysis(
  pheno_file  = "phenoAdrenals.txt",
  exprs_file  = "exprsAdrenals.txt",
  tissue_name = "Adrenals"
)
}
}
