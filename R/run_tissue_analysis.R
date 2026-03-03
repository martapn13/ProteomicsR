#' @title Multi-Organ Proteomics Analysis Pipeline
#'
#' @description
#' Runs a full proteomics analysis pipeline for a single tissue, including
#' filtering, imputation, batch correction (RUViii-PRPS), differential
#' expression (limma/voom), and QC visualizations (RLE, PCA, volcano plots,
#' heatmaps).
#'
#' @param pheno_file Path to the phenotype \code{.txt} file (tab-separated). Must contain
#'   one row per sample with the following columns:
#'   \itemize{
#'     \item \code{phenostate} - experimental group (e.g. "WT", "KO", "DAPA", "QC")
#'     \item \code{sampleName} - unique sample identifier (used as row names)
#'     \item \code{sampleQC} - sample type, either "Sample" or "QC"
#'     \item \code{Replicates} - replicate label; matches sampleName for biological
#'       samples and "QC" for quality control samples
#'   }
#'   QC samples (where \code{Replicates == "QC"}) are automatically removed before analysis.
#' @param exprs_file Path to the expression \code{.txt} file (tab-separated). Must be
#'   formatted as:
#'   \itemize{
#'     \item First column = protein identifiers as row names (e.g. UniProt IDs such as "Q9JHZ2"
#'       or protein groups separated by semicolons such as "O88746;O88746-2")
#'     \item Remaining columns = samples (column names must match \code{sampleName} in the pheno file)
#'     \item Values = raw protein intensity (non-log scale, non-negative)
#'     \item Missing values should be coded as 0
#'     \item Includes both biological samples and quality control samples
#'       (e.g. "Quality_control_1", "Quality_control_2")
#'   }
#'   Example column structure: Quality_control_1, Quality_control_2, KO_single_sample_4,
#'   WT_single_sample_1, DAPA_single_sample_7, etc.
#' @param tissue_name Short label used in all output file names (e.g. \code{"Adrenals"}, \code{"Liver"}).
#' @param use_halfmin Logical. If \code{TRUE}, skips missForest and uses only half-minimum
#'   imputation. Recommended for tissues where \code{missForest} fails (e.g. Liver).
#'   Default is \code{FALSE}.
#' @param best_ncomp Integer. Number of surrogate variables (components) for RUViii.
#'   Default is \code{5}.
#'
#' @return Invisibly returns a named list with:
#' \describe{
#'   \item{df122C}{Corrected expression matrix (proteins x samples).}
#'   \item{efit}{The \code{eBayes} fitted model object from limma.}
#'   \item{deg_list}{Named list of DE results: \code{KO_vs_WT}, \code{DAPA_vs_WT}, \code{KO_vs_DAPA}.}
#' }
#'
#' @author Leong Ng \email{lln1@@leicester.ac.uk} (original pipeline)
#' @author Marta Nascimento \email{marta.parente.nascimento@@outlook.com} (modifications)
#'
#' @examples
#' \dontrun{
#' run_tissue_analysis(
#'   pheno_file  = "phenoAdrenals.txt",
#'   exprs_file  = "exprsAdrenals.txt",
#'   tissue_name = "Adrenals",
#'   use_halfmin = FALSE,
#'   best_ncomp  = 5
#' )
#' }
#'
#' @export
run_tissue_analysis <- function(pheno_file,
                                exprs_file,
                                tissue_name,
                                use_halfmin = FALSE,
                                best_ncomp  = 5) {

  requireNamespace("ggfortify", quietly = TRUE)

  message("\n========== Starting analysis: ", tissue_name, " ==========\n")

  ColorBatch <- c("darkgreen", "blue", "red4", "gold1", "yellow4", "brown", "lightblue",
                  "coral4", "maroon4", "chartreuse4",
                  "#660099", "#CC0066", "#FF9999", "#FF9900")

  # ── 1. Read data ─────────────────────────────────────────────────────────────
  pheno <- read.table(pheno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(pheno) <- pheno$sampleName
  pheno1 <- pheno[!(pheno$Replicates == "QC"), ]

  exprs1 <- read.delim(exprs_file, row.names = 1, check.names = FALSE)
  exprs1 <- as.matrix(exprs1)

  # ── 2. Remove low-variance features (caret NZV) ──────────────────────────────
  myfilt1 <- t(exprs1)
  nzv     <- caret::preProcess(myfilt1, method = "nzv", uniqueCut = 10)
  myfilt2 <- t(predict(nzv, myfilt1))
  message("Dimensions after NZV filter: ", nrow(myfilt2), " x ", ncol(myfilt2))
  message("Min value after NZV filter: ", min(myfilt2))

  # ── 3. Imputation ─────────────────────────────────────────────────────────────
  # Half-minimum always runs first
  df12  <- as.matrix(myfilt2)
  df12T <- as.data.frame(t(df12))
  cols  <- 1:ncol(df12T)
  df12T[cols] <- lapply(df12T[cols], function(x)
    replace(x, x == 0, min(x[x > 0], na.rm = TRUE) / 2))
  df122 <- t(df12T)
  message("Min value after half-min: ", min(df122))
  hist(log2(df122), main = paste("log2 intensity -", tissue_name))

  if (!use_halfmin) {
    message("Using missForest imputation for: ", tissue_name)
    df122 <- t(missForest::missForest(t(df122))$ximp)
    message("Min value after missForest: ", min(df122))
  } else {
    message("Using half-minimum imputation only for: ", tissue_name)
  }

  # ── 4. SPECU feature ranking ──────────────────────────────────────────────────
  X    <- t(df122)
  out1 <- Rdimtools::do.specu(X, ndim = ncol(X) * 1, sigma = 5,
                              preprocess = "cscale", ranking = "method1")
  yyyy    <- out1[["featidx"]]
  ssss    <- colnames(X)
  llll    <- ssss[c(yyyy)]
  qqq     <- seq(max(length(llll)))
  Ordered <- data.frame(llll[qqq])
  colnames(Ordered) <- "Peakorder"
  rownames(Ordered) <- Ordered$Peakorder
  Ordered$ordrow    <- 1:nrow(Ordered)

  NumAll  <- nrow(Ordered)
  NumPeak <- as.integer(nrow(Ordered) * 0.75)
  NegC    <- Ordered[c(NumPeak:NumAll), ]

  # ── 5. RUViii-PRPS batch correction ──────────────────────────────────────────
  df120 <- log2(df122[, rownames(pheno)])
  M100  <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(raw = df120),
    colData = pheno
  )

  raw.data1  <- as.matrix(SummarizedExperiment::assay(M100, "raw"))
  rep_groups <- unique(pheno$Replicates)

  ReplicateMatrix1 <- matrix(0, nrow = nrow(pheno), ncol = length(rep_groups))
  colnames(ReplicateMatrix1) <- paste0("mouse", rep_groups)
  rownames(ReplicateMatrix1) <- pheno$sampleName
  for (i in seq_along(rep_groups)) {
    ReplicateMatrix1[, i] <- as.numeric(pheno$Replicates == rep_groups[i])
  }

  ncg.set1 <- rownames(raw.data1) %in% NegC$Peakorder
  message("Negative control features: ", sum(ncg.set1))

  df10 <- tcgaCleaneR::runRUV_III_PRPS(
    ruv.data    = t(raw.data1),
    ruv.rep     = ReplicateMatrix1,
    ncg.set     = ncg.set1,
    k           = best_ncomp,
    return.info = TRUE
  )
  ruv.iii.10 <- t(df10$new.ruv.data)

  # ── 6. QC plots: RLE and PCA ──────────────────────────────────────────────────
  print(
    ruv::ruv_rle(Y = t(raw.data1), ylim = c(-3, 3), rowinfo = pheno) +
      ggplot2::geom_point(ggplot2::aes(x = rle.x.factor, y = middle, colour = state)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(colour = "state") +
      ggplot2::scale_color_manual(values = ColorBatch) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "cyan")
  )
  print(
    ruv::ruv_rle(Y = t(ruv.iii.10), ylim = c(-3, 3), rowinfo = pheno) +
      ggplot2::geom_point(ggplot2::aes(x = rle.x.factor, y = middle, colour = state)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(colour = "state") +
      ggplot2::scale_color_manual(values = ColorBatch) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "cyan")
  )

  pca_res <- prcomp(t(as.matrix(raw.data1)), scale = TRUE)
  print(autoplot(pca_res, data = pheno, colour = "state") +
          ggplot2::scale_colour_manual(values = rainbow(4)))

  pca_res1 <- prcomp(t(as.matrix(ruv.iii.10)), scale = TRUE)
  print(autoplot(pca_res1, data = pheno, colour = "state") +
          ggplot2::scale_colour_manual(values = rainbow(4)))

  # ── 7. NOISeq / ARSyNseq ─────────────────────────────────────────────────────
  metadataB  <- data.frame(labelDescription = colnames(pheno))
  phenoDataB <- new("AnnotatedDataFrame", data = pheno, varMetadata = metadataB)
  mydata2    <- Biobase::ExpressionSet(assayData = ruv.iii.10, phenoData = phenoDataB)

  mydatacorr32 <- NOISeq::ARSyNseq(mydata2, factor = "state", batch = FALSE, norm = "tmm",
                                   Variability = 0.75, beta = 2, logtransf = TRUE)
  df122C <- 2^mydatacorr32@assayData[["exprs"]]

  write.csv(df122C, file = paste0("ExprsprpsNoiseq_", tissue_name, ".csv"))

  pca_res2 <- prcomp(t(log2(df122C)), scale = TRUE)
  print(autoplot(pca_res2, data = pheno, colour = "state") +
          ggplot2::scale_colour_manual(values = rainbow(4)))

  print(
    ruv::ruv_rle(Y = t(log2(df122C)), ylim = c(-3, 3), rowinfo = pheno) +
      ggplot2::geom_point(ggplot2::aes(x = rle.x.factor, y = middle, colour = state)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(colour = "state") +
      ggplot2::scale_color_manual(values = ColorBatch) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "cyan")
  )

  # ── 8. limma differential expression ─────────────────────────────────────────
  state  <- factor(pheno1$state, levels = c("WT", "KO", "DAPA"))
  mouse  <- factor(pheno1$Replicates)

  unpaired.design <- model.matrix(~ 0 + state)
  colnames(unpaired.design) <- c("WT", "KO", "DAPA")
  rownames(unpaired.design) <- rownames(pheno1)

  d  <- edgeR::DGEList(counts = df122C[, rownames(unpaired.design)])
  d0 <- edgeR::calcNormFactors(d, method = "TMM")
  v  <- limma::voom(counts = d0, design = unpaired.design,
                    lib.size = NULL, normalize.method = "none", plot = TRUE)

  contrast.matrix <- limma::makeContrasts(KO - WT, DAPA - WT, KO - DAPA, levels = unpaired.design)
  corfit <- limma::duplicateCorrelation(v, design = unpaired.design, ndups = 1, block = pheno1$mouse)
  vfit   <- limma::lmFit(v, design = unpaired.design, block = mouse, cor = corfit$consensus)
  vfit   <- limma::contrasts.fit(vfit, contrasts = contrast.matrix)
  efit   <- limma::eBayes(vfit, robust = TRUE)

  limma::plotSA(vfit, main = " model: Mean-variance trend")
  limma::plotSA(efit, main = "Final model: Mean-variance trend")

  dt5 <- limma::decideTests(efit, p.value = 0.05, lfc = log2(1.5))
  summary(dt5)
  limma::write.fit(efit, results = dt5,
                   file     = paste0("Res_noiseq_", tissue_name, ".csv"),
                   digits   = NULL, adjust = "BH", method = "separate",
                   F.adjust = "BH", quote = TRUE, sep = ",", row.names = TRUE)

  deg_KO_vs_WT   <- limma::topTable(efit, coef = "KO - WT",   adjust.method = "fdr", number = Inf)
  deg_DAPA_vs_WT <- limma::topTable(efit, coef = "DAPA - WT", adjust.method = "fdr", number = Inf)
  deg_KO_vs_DAPA <- limma::topTable(efit, coef = "KO - DAPA", adjust.method = "fdr", number = Inf)

  logFC_threshold <- log2(1.5)
  FDR_threshold   <- 0.05

  sig_KO_vs_WT   <- deg_KO_vs_WT[abs(deg_KO_vs_WT$logFC) > logFC_threshold &
                                    deg_KO_vs_WT$adj.P.Val < FDR_threshold, ]
  sig_DAPA_vs_WT <- deg_DAPA_vs_WT[abs(deg_DAPA_vs_WT$logFC) > logFC_threshold &
                                      deg_DAPA_vs_WT$adj.P.Val < FDR_threshold, ]
  sig_KO_vs_DAPA <- deg_KO_vs_DAPA[abs(deg_KO_vs_DAPA$logFC) > logFC_threshold &
                                      deg_KO_vs_DAPA$adj.P.Val < FDR_threshold, ]

  write.csv(sig_KO_vs_WT,   paste0(tissue_name, "_sig_KO_vs_WT.csv"),   row.names = TRUE)
  write.csv(sig_DAPA_vs_WT, paste0(tissue_name, "_sig_DAPA_vs_WT.csv"), row.names = TRUE)
  write.csv(sig_KO_vs_DAPA, paste0(tissue_name, "_sig_KO_vs_DAPA.csv"), row.names = TRUE)

  # ── 9. Volcano plots ──────────────────────────────────────────────────────────
  deg_list <- list(KO_vs_WT   = deg_KO_vs_WT,
                   DAPA_vs_WT = deg_DAPA_vs_WT,
                   KO_vs_DAPA = deg_KO_vs_DAPA)

  for (name in names(deg_list)) {
    df <- deg_list[[name]]
    df$category <- "Not significant"
    df$category[df$adj.P.Val < FDR_threshold & df$logFC >  logFC_threshold] <- "Upregulated"
    df$category[df$adj.P.Val < FDR_threshold & df$logFC < -logFC_threshold] <- "Downregulated"
    df$category <- factor(df$category, levels = c("Upregulated", "Downregulated", "Not significant"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = category)) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::scale_color_manual(values = c("Upregulated"     = "blue",
                                             "Downregulated"   = "red",
                                             "Not significant" = "grey")) +
      ggplot2::geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed") +
      ggplot2::labs(title = paste("Volcano Plot:", tissue_name, "-", name),
                    x     = expression(Log[2] ~ Fold ~ Change),
                    y     = expression(-Log[10] ~ adj.P ~ value)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     legend.title = ggplot2::element_blank())

    ggplot2::ggsave(paste0("volcano_", tissue_name, "_", name, ".png"),
                    plot = p, width = 10, height = 8, dpi = 300)
    print(p)

    message(tissue_name, " | ", name, ": ",
            sum(df$adj.P.Val < FDR_threshold & df$logFC >  logFC_threshold), " up, ",
            sum(df$adj.P.Val < FDR_threshold & df$logFC < -logFC_threshold), " down")
  }

  # ── 10. Heatmaps ──────────────────────────────────────────────────────────────
  expr_mat <- df122C
  qc_cols  <- grep("Quality_control", colnames(expr_mat), value = TRUE)
  if (length(qc_cols) > 0) expr_mat <- expr_mat[, !colnames(expr_mat) %in% qc_cols, drop = FALSE]

  pheno1_h       <- pheno1[colnames(expr_mat), , drop = FALSE]
  annotation_col <- data.frame(state = pheno1_h$state)
  rownames(annotation_col) <- rownames(pheno1_h)

  for (comp in names(deg_list)) {
    df   <- deg_list[[comp]]
    hits <- subset(df, adj.P.Val < FDR_threshold & abs(logFC) > logFC_threshold)

    if (nrow(hits) == 0) {
      message(tissue_name, " | ", comp, ": No significant DE proteins.")
      next
    }

    topP      <- head(rownames(hits[order(hits$adj.P.Val), ]), 10)
    mat       <- expr_mat[topP, , drop = FALSE]
    col_order <- rownames(pheno1_h)[order(pheno1_h$state)]

    pheatmap::pheatmap(
      t(scale(t(mat[, col_order, drop = FALSE]))),
      annotation_col = annotation_col[col_order, , drop = FALSE],
      show_rownames  = TRUE,
      show_colnames  = TRUE,
      cluster_cols   = FALSE,
      fontsize_row   = 8,
      main           = paste("Top 10 DE proteins:", tissue_name, "-", comp)
    )
  }

  message("\n========== Finished: ", tissue_name, " ==========\n")
  invisible(list(df122C = df122C, efit = efit, deg_list = deg_list))
}

  message("\n========== Finished: ", tissue_name, " ==========\n")
  invisible(list(df122C = df122C, efit = efit, deg_list = deg_list))
}
