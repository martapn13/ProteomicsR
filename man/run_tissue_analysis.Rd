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
\item{pheno_file}{Path to the phenotype .txt file (tab-separated).}
\item{exprs_file}{Path to the expression .txt file (tab-separated).}
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
