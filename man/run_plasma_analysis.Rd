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
\item{pheno_file}{Path to the phenotype .txt file (tab-separated).}
\item{exprs_file}{Path to the expression .txt file (tab-separated).}
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
