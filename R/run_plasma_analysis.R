#' @title Plasma Proteomics Analysis Pipeline
#'
#' @description
#' Runs a full proteomics analysis pipeline for plasma samples, including
#' filtering, half-minimum imputation, ARSyNseq normalisation, differential
#' expression (limma/voom), ID mapping, volcano plots, and heatmaps.
#' Note: plasma uses half-minimum imputation only (no missForest, no RUViii-PRPS).
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
#' @param raw_file Path to the raw plasma \code{.txt} file (tab-separated) used for ID
#'   mapping. The first column must contain protein identifiers and will be renamed to \code{ID}.
#'   Row order must match the expression matrix.
#'
#' @return Invisibly returns a named list with:
#' \describe{
#'   \item{df_norm}{Corrected expression matrix (proteins x samples).}
#'   \item{efit}{The \code{eBayes} fitted model object from limma.}
#'   \item{deg_list}{Named list of DE results: \code{KO_vs_WT}, \code{DAPA_vs_WT}, \code{KO_vs_DAPA}.}
#' }
#'
#' @author Leong Ng \email{lln1@@leicester.ac.uk} (original pipeline)
#' @author Marta Nascimento \email{marta.parente.nascimento@@outlook.com} (modifications)
#'
#' @examples
#' \dontrun{
#' run_plasma_analysis(
#'   pheno_file = "phenoPlasma.txt",
#'   exprs_file = "exprsPlasma.txt",
#'   raw_file   = "Plasma_Raw.txt"
#' )
#' }
#'
#' @export
run_plasma_analysis <- function(pheno_file,
                                exprs_file,
                                raw_file) {

  requireNamespace("ggfortify", quietly = TRUE)

  message("\n========== Starting plasma analysis ==========\n")

  ColorBatch <- c("darkgreen", "blue", "red4", "gold1", "yellow4", "brown", "lightblue",
                  "coral4", "maroon4", "chartreuse4", "#660099", "#CC0066",
                  "#FF9999", "#FF9900")

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

  # ── 3. Half-minimum imputation ────────────────────────────────────────────────
  df12  <- as.matrix(myfilt2)
  df12T <- as.data.frame(t(df12))
  cols  <- 1:ncol(df12T)
  df12T[cols] <- lapply(df12T[cols], function(x)
    replace(x, x == 0, min(x[x > 0], na.rm = TRUE) / 2))
  df122 <- t(df12T)
  message("Min value after half-min: ", min(df122))
  hist(log2(df122), main = "log2 intensity - Plasma")

  # ── 4. Log2 transform (non-QC samples only) ──────────────────────────────────
  df_log <- log2(df122[, rownames(pheno1)])

  # ── 5. ARSyNseq normalisation ─────────────────────────────────────────────────
  metadataB  <- data.frame(labelDescription = colnames(pheno1))
  phenoDataB <- new("AnnotatedDataFrame", data = pheno1, varMetadata = metadataB)
  mydata2    <- Biobase::ExpressionSet(assayData = df_log, phenoData = phenoDataB)

  mydatacorr <- NOISeq::ARSyNseq(mydata2, factor = "state", batch = FALSE, norm = "tmm",
                                 Variability = 0.75, beta = 2, logtransf = TRUE)
  df_norm <- 2^mydatacorr@assayData[["exprs"]]

  # ── 6. QC plots: RLE and PCA ──────────────────────────────────────────────────
  print(
    ruv::ruv_rle(Y = t(df_log), ylim = c(-3, 3), rowinfo = pheno1) +
      ggplot2::geom_point(ggplot2::aes(x = rle.x.factor, y = middle, colour = state)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(colour = "state") +
      ggplot2::scale_color_manual(values = ColorBatch) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "cyan")
  )

  pca_res <- prcomp(t(as.matrix(df_log)), scale = TRUE)
  print(autoplot(pca_res, data = pheno1, colour = "state") +
          ggplot2::scale_colour_manual(values = rainbow(4)))

  pca_res2 <- prcomp(t(log2(df_norm)), scale = TRUE)
  print(autoplot(pca_res2, data = pheno1, colour = "state") +
          ggplot2::scale_colour_manual(values = rainbow(4)))

  print(
    ruv::ruv_rle(Y = t(log2(df_norm)), ylim = c(-3, 3), rowinfo = pheno1) +
      ggplot2::geom_point(ggplot2::aes(x = rle.x.factor, y = middle, colour = state)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(colour = "state") +
      ggplot2::scale_color_manual(values = ColorBatch) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "cyan")
  )

  # ── 7. limma differential expression ─────────────────────────────────────────
  state  <- factor(pheno1$state, levels = c("WT", "KO", "DAPA"))
  mouse  <- factor(pheno1$Replicates)

  unpaired.design <- model.matrix(~ 0 + state)
  colnames(unpaired.design) <- c("WT", "KO", "DAPA")
  rownames(unpaired.design) <- rownames(pheno1)

  d  <- edgeR::DGEList(counts = df_norm)
  d0 <- edgeR::calcNormFactors(d, method = "TMM")
  v  <- limma::voom(counts = d0, design = unpaired.design,
                    lib.size = NULL, normalize.method = "none", plot = TRUE)

  contrast.matrix <- limma::makeContrasts(KO - WT, DAPA - WT, KO - DAPA, levels = unpaired.design)
  corfit <- limma::duplicateCorrelation(v, design = unpaired.design, ndups = 1, block = pheno1$Replicates)
  vfit   <- limma::lmFit(v, design = unpaired.design, block = mouse, cor = corfit$consensus)
  vfit   <- limma::contrasts.fit(vfit, contrasts = contrast.matrix)
  efit   <- limma::eBayes(vfit, robust = TRUE)

  limma::plotSA(vfit, main = "Model: Mean-variance trend")
  limma::plotSA(efit, main = "Final model: Mean-variance trend")

  dt5 <- limma::decideTests(efit, p.value = 0.05, lfc = log2(1.5))
  summary(dt5)
  limma::write.fit(efit, results = dt5,
                   file     = "Res2_noiseq2_plasma.csv",
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

  write.csv(sig_KO_vs_WT,   "sig_KO_vs_WT_Plasma.csv",   row.names = TRUE)
  write.csv(sig_DAPA_vs_WT, "sig_DAPA_vs_WT_Plasma.csv", row.names = TRUE)
  write.csv(sig_KO_vs_DAPA, "sig_KO_vs_DAPA_Plasma.csv", row.names = TRUE)

  # ── 8. ID mapping ─────────────────────────────────────────────────────────────
  Plasma_raw <- read.table(raw_file, header = TRUE, sep = "\t",
                           check.names = FALSE, stringsAsFactors = FALSE)
  colnames(Plasma_raw)[1] <- "ID"

  add_IDs_exact <- function(sig_df, raw_df, out_prefix = "sig") {
    rn   <- rownames(sig_df)
    nums <- suppressWarnings(as.numeric(gsub("^row", "", rn, ignore.case = TRUE)))
    if (any(is.na(nums))) {
      warning("Some rownames did not match pattern '^row' and will produce NA indices.")
    }
    n_raw  <- nrow(raw_df)
    ID_vec <- rep(NA_character_, length(nums))
    valid  <- !is.na(nums) & nums >= 1 & nums <= n_raw
    if (any(valid)) ID_vec[valid] <- raw_df$ID[nums[valid]]

    sig_out <- sig_df
    sig_out$OrigRowName  <- rn
    sig_out$RowNumber    <- nums
    sig_out$RowIndexUsed <- nums
    sig_out$ID           <- ID_vec
    cols    <- colnames(sig_out)
    cols    <- c("ID", setdiff(cols, "ID"))
    sig_out <- sig_out[, cols, drop = FALSE]

    matched <- sum(!is.na(ID_vec))
    total   <- length(ID_vec)
    message(sprintf("%s: matched %d / %d (%.1f%%)", out_prefix, matched, total, 100 * matched / total))

    if (matched < total) {
      unmatched_df <- data.frame(OrigRowName = rn[is.na(ID_vec)],
                                 RowNumber   = nums[is.na(ID_vec)],
                                 stringsAsFactors = FALSE)
      write.csv(unmatched_df, paste0("unmatched_rows_", out_prefix, ".csv"), row.names = FALSE)
      message(sprintf("Wrote unmatched rows to: unmatched_rows_%s.csv", out_prefix))
    }
    return(sig_out)
  }

  sig_KO_vs_WT_ID   <- add_IDs_exact(sig_KO_vs_WT,   Plasma_raw, out_prefix = "KO_vs_WT")
  sig_DAPA_vs_WT_ID <- add_IDs_exact(sig_DAPA_vs_WT, Plasma_raw, out_prefix = "DAPA_vs_WT")
  sig_KO_vs_DAPA_ID <- add_IDs_exact(sig_KO_vs_DAPA, Plasma_raw, out_prefix = "KO_vs_DAPA")

  write.csv(sig_KO_vs_WT_ID,   "sig_KO_vs_WT_Plasma_withID.csv",   row.names = FALSE)
  write.csv(sig_DAPA_vs_WT_ID, "sig_DAPA_vs_WT_Plasma_withID.csv", row.names = FALSE)
  write.csv(sig_KO_vs_DAPA_ID, "sig_KO_vs_DAPA_Plasma_withID.csv", row.names = FALSE)

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
      ggplot2::labs(title = paste("Volcano Plot: Plasma -", name),
                    x     = expression(Log[2] ~ Fold ~ Change),
                    y     = expression(-Log[10] ~ adj.P ~ value)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     legend.title = ggplot2::element_blank())

    ggplot2::ggsave(paste0("volcano_", name, "_full.png"),
                    plot = p, width = 10, height = 8, dpi = 300)
    print(p)

    message("Plasma | ", name, ": ",
            sum(df$adj.P.Val < FDR_threshold & df$logFC >  logFC_threshold), " up, ",
            sum(df$adj.P.Val < FDR_threshold & df$logFC < -logFC_threshold), " down")
  }

  # ── 10. Heatmaps ──────────────────────────────────────────────────────────────
  expr_mat <- df_norm
  qc_cols  <- grep("Quality_control", colnames(expr_mat), value = TRUE)
  if (length(qc_cols) > 0) {
    message("Removing QC columns: ", paste(qc_cols, collapse = ", "))
    expr_mat <- expr_mat[, !colnames(expr_mat) %in% qc_cols, drop = FALSE]
  }

  expr_samples   <- colnames(expr_mat)
  pheno1         <- pheno1[expr_samples, , drop = FALSE]
  annotation_col <- data.frame(state = pheno1$state)
  rownames(annotation_col) <- rownames(pheno1)

  for (comp in names(deg_list)) {
    df   <- deg_list[[comp]]
    hits <- subset(df, adj.P.Val < FDR_threshold & abs(logFC) > logFC_threshold)

    if (nrow(hits) == 0) {
      message("Plasma | ", comp, ": No significant DE proteins.")
      next
    }

    topP      <- head(rownames(hits[order(hits$adj.P.Val), ]), 10)
    mat       <- expr_mat[topP, , drop = FALSE]
    col_order <- rownames(pheno1)[order(pheno1$state)]

    pheatmap::pheatmap(
      t(scale(t(mat[, col_order, drop = FALSE]))),
      annotation_col = annotation_col[col_order, , drop = FALSE],
      show_rownames  = TRUE,
      show_colnames  = TRUE,
      cluster_cols   = FALSE,
      fontsize_row   = 8,
      main           = paste("Top 10 DE proteins: Plasma -", comp)
    )
  }

  message("\n========== Finished: Plasma ==========\n")
  invisible(list(df_norm = df_norm, efit = efit, deg_list = deg_list))
}
