# ============================================================
#  Run ProteomicsR pipeline for all tissues
#  Edit file paths below to match your data location
#  Set use_halfmin = TRUE for tissues where missForest fails
# ============================================================

run_tissue_analysis(
  pheno_file  = "phenoAdrenals.txt",
  exprs_file  = "exprsAdrenals.txt",
  tissue_name = "Adrenals",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoKidney.txt",
  exprs_file  = "exprsKidney.txt",
  tissue_name = "Kidney",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoLiver.txt",
  exprs_file  = "exprsLiver.txt",
  tissue_name = "Liver",
  use_halfmin = TRUE,   # half-min imputation for Liver
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoBAT.txt",
  exprs_file  = "exprsBAT.txt",
  tissue_name = "BAT",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoWAT.txt",
  exprs_file  = "exprsWAT.txt",
  tissue_name = "WAT",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoSpleen.txt",
  exprs_file  = "exprsSpleen.txt",
  tissue_name = "Spleen",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoMuscle.txt",
  exprs_file  = "exprsMuscle.txt",
  tissue_name = "Muscle",
  use_halfmin = FALSE,
  best_ncomp  = 5
)

run_tissue_analysis(
  pheno_file  = "phenoHeart.txt",
  exprs_file  = "exprsHeart.txt",
  tissue_name = "Heart",
  use_halfmin = FALSE,
  best_ncomp  = 5
)
