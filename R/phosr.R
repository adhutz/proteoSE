# PhosR signalome and PTM normalisation
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' Prepares Phosphorylation Data for Analysis
#'
#' This function takes a SummarizedExperiment object and prepares the data for
#' phosphorylation analysis using the PhosR package. It also generates a color palette
#' for plotting.
#'
#' @param se A SummarizedExperiment object.
#' @param species A character string specifying the species (default is "human").
#'        Must be one of "human", "mouse", or "rat".
#' @param numMotif An integer specifying the number of motifs (default is 5).
#' @param numSub An integer specifying the number of substrates (default is 1).
#' @param top An integer specifying the top number of substrates to select (default is 30).
#' @return A list containing the following elements:
#'         - kinaseSubstrateScore: A list of kinase substrate scores.
#'         - kinaseSubstratePred: A matrix of kinase substrate predictions.
#'         - mat: A matrix of collapsed phosphorylation data.
#'         - kinase_all_color: A color palette for plotting.
#'         - kinase_signalome_color" = A color palette for plotting.
#'         - seq: A vector of sequences.
#' @importFrom PhosR phosCollapse kinaseSubstrateScore kinaseSubstratePred
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
prep_phosR <- function(se, species = "human", numMotif = 5, numSub = 1, top = 30){
  
  # Load PhosphoSitePlus into the function frame, not globalenv (R CMD check NOTEs on the latter).
  # KinaseFamily is unused here and kinaseSubstrateScore() loads KinaseMotifs
  # itself, so PhosphoSitePlus is the only set we need.
  utils::data("PhosphoSitePlus", package = "PhosR", envir = environment())
  
  mat <- assay(se)
  rownames(mat) <- paste0(gsub("_", ";", rownames(mat)), ";",gsub("(.*?);.*", "\\1", rowData(se)$sequence_window))
  
  mat <- phosCollapse(assay(se), id=gsub("(.*;.*;)[1-9]*;(.*)", "\\1\\2", rownames(mat), perl = TRUE), 
                      stat=apply(abs(assay(se)), 1, max), by = "max")
  
  seq <- gsub(".*;.*;(.*)", "\\1", rownames(mat))
  
  rownames(mat) <- gsub("(.*;.*;).*","\\1",rownames(mat))
  colnames(mat) <- se$ID
    
  if(species == "human"){
    psite = PhosphoSite.human
  } else if(species == "mouse"){
    psite = PhosphoSite.mouse
  } else if(species == "rat"){
    psite = PhosphoSite.mouse
  } else{
    return("Species must be one of human, mouse or rat")
  }
  
  kss <- kinaseSubstrateScore(psite, mat, seq,numMotif = numMotif, numSub = numSub, species = "human")
  set.seed(1)
  ksp <- kinaseSubstratePred(kss, top=top)
  
  
  # Color palette
  my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  kinase_all_color <- my_color_palette(ncol(kss$combinedScoreMatrix))
  names(kinase_all_color) <- colnames(kss$combinedScoreMatrix)
  kinase_signalome_color <- kinase_all_color[colnames(ksp)]
  
  
  return(list("kinaseSubstrateScore" = kss, "kinaseSubstratePred" = ksp, "mat" = mat, "kinase_all_color" = kinase_all_color, "kinase_signalome_color" = kinase_signalome_color, "seq" = seq))
}




#' Plot Signalome Map
#'
#' This function creates a Signalome Map plot using the PhosR::plotSignalomeMap function
#' and customizes the appearance with additional ggplot2 layers and a custom theme.
#'
#' @param signalome_res A data frame containing the Signalome results.
#' @param kinase_signalome_color A vector of colors for the kinases in the Signalome Map.
#'
#' @return A ggplot object representing the Signalome Map plot.
#' @import ggplot2
#' @importFrom PhosR plotSignalomeMap
#' 
#' @export
# Signalome Map
plot_signalome_map <- function(signalome_res, kinase_signalome_color){
  
  p <- PhosR::plotSignalomeMap(signalome_res, kinase_signalome_color) +
    scale_size_continuous(range = c(1,5)) + 
    scale_y_continuous(labels = paste0("Cluster ", levels(as.factor(signalome_res$proteinModules))), breaks = 1:length(unique(signalome_res$proteinModules))) +
    theme_proteoSE()
  
  return(p)
}



#' Create a bar plot with error bars for the given data
#'
#' This function takes a matrix and a signalome result, processes the data, and creates a bar plot with error bars.
#' It returns a list of ggplot objects, one for each module in the signalome result.
#'
#' @param mat A matrix containing the data to be plotted
#' @param signalome_res A list containing the results of the signalome analysis
#' @return A list of ggplot objects, one for each module in the signalome result
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @export

module_barplot <- function(mat, signalome_res){
  
  
  # Get data
  df <- mat %>% as.data.frame() %>% mutate(gene_names = gsub("(.*);(.*);", "\\1", rownames(mat))) %>% tidyr::pivot_longer(-gene_names, names_to = "ID", values_to = "scaled_intensity") %>%
    mutate(condition = gsub("(.*)_.*", "\\1", ID))
  
  mod_list <- split(names(signalome_res$proteinModules), signalome_res$proteinModules)
  names(mod_list) <- paste0("Cluster_", names(mod_list))
  
  mod_list <- lapply(mod_list, function(x) filter(df, gene_names %in% x))
  
  l <- do.call(cbind, lapply(mod_list, function(x) group_by(x, condition) %>% summarize(m = mean(scaled_intensity))))
  
  # Function to calculate mean and standard deviation
  mean_sd <- function(x) {
    y = mean(x, na.rm = TRUE)
    return(c(y = y,
             ymin = y - sd(x, na.rm = TRUE),
             ymax = y + sd(x, na.rm = TRUE)))
  }
  
  # Create the bar plot with error bars
  res <- lapply(seq_along(mod_list), function(z) ggplot(mod_list[[z]], aes(x = factor(condition), y = scaled_intensity, fill = factor(condition))) +
                  geom_boxplot()+
                  xlab("Condition") +
                  ylab("Scaled intensity") +
                  labs(title = names(mod_list)[[z]]) +
                  theme_proteoSE())
  
  return(res)
}



#' Prepares phosphoproteomics data for PhosR analysis
#'
#' This function processes phosphoproteomics data, calculates kinase-substrate scores, and generates predictions for the given data.
#' It also generates color palettes for visualization.
#'
#' @param ppe Phosphoproteomics data object
#' @param species Species of the data, one of "human", "mouse", or "rat" (default: "human")
#' @param numMotif Number of motifs to use (default: 5)
#' @param numSub Number of substrates to use (default: 1)
#' @param top Number of top kinases to consider (default: 30)
#' @param assay Type of assay data to use (default: "z_scored")
#'
#' @return A list containing kinaseSubstrateScore, kinaseSubstratePred, mat, kinase_all_color, kinase_signalome_color, and seq
#' @importFrom PhosR PhosphoExperiment kinaseSubstrateScore kinaseSubstratePred
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
prep_phosR_from_ppe <- function(ppe, species = "human", numMotif = 5, numSub = 1, top = 30, assay = "z_scored"){
  
  # Load PhosphoSitePlus into the function frame, not globalenv (R CMD check NOTEs on the latter).
  # KinaseFamily is unused here and kinaseSubstrateScore() loads KinaseMotifs
  # itself, so PhosphoSitePlus is the only set we need.
  utils::data("PhosphoSitePlus", package = "PhosR", envir = environment())
  
  # Get data from provided object
  mat <- ppe@assays@data[[assay]]
  
  seq <- gsub(".*;.*;(.*)", "\\1", rownames(mat))
  
  rownames(mat) <- gsub(".*;(.*;.*;).*","\\1",rownames(mat))
  colnames(mat) <- ppe$ID
  
  if(species == "human"){
    psite = PhosphoSite.human
  } else if(species == "mouse"){
    psite = PhosphoSite.mouse
  } else if(species == "rat"){
    psite = PhosphoSite.rat
  } else{
    return("Species must be one of human, mouse or rat")
  }
  
  kss <- kinaseSubstrateScore(psite, mat, seq, numMotif = numMotif, numSub = numSub, species = "human")
  set.seed(1)
  ksp <- kinaseSubstratePred(kss, top=top)
  
  
  # Color palette
  my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  kinase_all_color <- my_color_palette(ncol(kss$combinedScoreMatrix))
  names(kinase_all_color) <- colnames(kss$combinedScoreMatrix)
  kinase_signalome_color <- kinase_all_color[colnames(ksp)]
  
  
  return(list("kinaseSubstrateScore" = kss, "kinaseSubstratePred" = ksp, "mat" = mat, "kinase_all_color" = kinase_all_color, "kinase_signalome_color" = kinase_signalome_color, "seq" = seq))
}


#' AOV Scale PhosR Function
#'
#' This function takes a ppe object and applies a series of transformations,
#' including mean abundance calculation, ANOVA filtering, subsetting, standardization,
#' and z-scoring. The transformed data is then added as a new assay to the ppe object.
#'
#' @param ppe A ppe object.
#' @param p_cut A numeric threshold for selecting rows based on ANOVA p-value (default: 0.05).
#' @param fc_cut A numeric threshold for selecting rows based on fold change (default: 0.5).
#' @param assay A character specifying the assay to use from ppe (default: "normalised").
#' @return A ppe object with the transformed data added as a new assay named "z_scored".
#' @importFrom PhosR PhosphoExperiment meanAbundance matANOVA
#' @export
aov_scale_phosR <- function(ppe, p_cut = 0.05, fc_cut = 0.5, assay = "normalised"){
  mat <- ppe@assays@data[[assay]]
  mat.mean <- meanAbundance(mat, grps = colData(ppe)$condition)
  
  aov <- matANOVA(mat=mat, grps = colData(ppe)$condition)
  idx <- (aov < p_cut) & (rowSums(mat.mean > fc_cut) > 0)
  
  mat.reg <- mat[idx, ,drop = FALSE]
  
  mat.std <- PhosR::standardise(mat.reg)
  colnames(mat.std) <- ppe$ID
  
  ppe <- ppe[idx,]
  ppe <- add_assay(ppe, mat.std, "z_scored")
  
  return(ppe)
}


#' Test Differential Phosphorylation
#'
#' This function tests differential phosphorylation using a linear model fit and an empirical Bayes approach.
#'
#' @param ppe A PhosR object containing the phosphorylation data.
#' @param contrast A character vector specifying the contrasts to be tested. Default is "treatment_vs_control".
#' @param test_all A logical value indicating whether to test all possible contrasts or just the specified ones. Default is FALSE.
#' @param assay A character string specifying the assay to use for the analysis. Default is "scaled".
#'
#' @return A modified PhosR object with the differential phosphorylation test results added to the rowData.
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom S4Vectors expand.grid
#' @importFrom BiocGenerics cbind
#' @import dplyr
#' @export
test_diff_phosR <- function(ppe, contrast = c("treatment_vs_control"), test_all = FALSE, assay = "scaled"){
  
  # Create design matrix
  design <- model.matrix(~ ppe$condition - 1)
  colnames(design) <- gsub("ppe\\$condition", "", colnames(design))
  
  fit <- lmFit(ppe@assays@data[[assay]], design)
  
  # Create contrasts
  if(test_all == FALSE){
    contrast_ <- gsub("_vs_", "-", contrast)
  } else{
    contrast_ <- expand.grid(colnames(design), colnames(design)) %>% 
      filter(!Var1 == Var2) %>% mutate(contrast_ = paste0(Var1, "-", Var2)) %>%
      select(contrast_) %>%
      unlist()
  }
  
  # Create contrast matrix
  contrast.matrix <- makeContrasts(contrasts = contrast_, levels=design)
  # topTable(coef = ) below indexes by the contrast string, but limma >= 3.68
  # names these columns after names(contrast_) -- which unlist() above sets to
  # "contrast_1", "contrast_2", ... See test_diff_limma() for the same fix.
  colnames(contrast.matrix) <- contrast_
  
  # fit model
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  fit2 <- eBayes(fit2)
  
  # Create a named list of contrasts
  names(contrast_) <- contrast_
  
  # Get a list of all tested contrasts and rename the columns accordingly
  l <- lapply(contrast_, function(x) topTable(fit2, coef = x, number = Inf, sort.by = "none", confint = TRUE))
  l <- lapply(seq_along(l), function(i) dplyr::rename_with(l[[i]], ~ paste0(names(l)[i], "_", .x)))
  
  # Merge the results and add the F statistic 
  l_com <- cbind(do.call(cbind, l), topTable(fit2, number = Inf, sort.by = "none", confint = TRUE))
  colnames(l_com) <- gsub("-", "_vs_", 
                          gsub("logFC", "diff", 
                               gsub("adj.P.Val", "p.adj", 
                                    gsub("P.Value", "p.val", colnames(l_com)))))
  
  
  rowData(ppe) <- rowData(ppe)[, !colnames(rowData(ppe)) %in% c(colnames(l_com), paste0(contrast, "_significant"), "significant")]
  rowData(ppe) <- cbind(rowData(ppe), l_com)
  
  return(ppe)
}




###########
#
############
#' Normalize PTM Intensities Against Background Protein Abundance
#'
#' @description
#' Normalizes post-translational modification (PTM) intensities against their
#' corresponding background protein abundances using a ratio-median approach.
#' For each PTM feature \eqn{i} and sample \eqn{j}, the normalized value is:
#'
#' \deqn{
#'   \hat{x}_{ij} = \frac{PTM_{ij} / BKG_{ij}}{median(PTM_i) / median(BKG_i)} \times median(PTM_i)
#' }
#'
#' which simplifies to \eqn{PTM_{ij} / BKG_{ij} \times median(BKG_i)},
#' correcting each PTM measurement for the abundance of its parent protein
#' while rescaling to a biologically interpretable range.
#'
#' Column alignment between PTM and background assays is performed on the
#' alphabetically sorted intersection of shared sample names, then restored
#' to the original PTM sample order before output.
#'
#' @param bkg_se A \code{SummarizedExperiment} containing the background
#'   (protein-level) abundance data. Rows are proteins; columns are samples.
#' @param ptm_se A \code{SummarizedExperiment} containing PTM-level abundance
#'   data. Rows are PTM features; columns are samples.
#' @param id_bkg Character scalar. Name of the \code{rowData} column in
#'   \code{bkg_se} that holds the gene/protein identifiers used for matching
#'   (e.g., \code{"gene_name"}).
#' @param id_ptm Character scalar. Name of the \code{rowData} column in
#'   \code{ptm_se} that holds the gene/protein identifiers used to look up
#'   each PTM's parent protein in \code{bkg_se} (e.g., \code{"gene_name"}).
#' @param assay_bkg Character scalar. Name of the assay in \code{bkg_se} to
#'   use as background (e.g., \code{"spectronaut_assay_prot"}).
#' @param assay_ptm Character scalar. Name of the assay in \code{ptm_se} to
#'   normalize (e.g., \code{"spectronaut_assay_ptm"}).
#' @param bkg_na Character scalar. Strategy for PTM rows whose parent protein
#'   has no match in \code{bkg_se}. One of:
#'   \describe{
#'     \item{\code{"set_as_na"}}{(default) Leave unmatched rows as \code{NA}
#'       in the normalized assay. Safe and explicit -- downstream analyses will
#'       see these rows as missing.}
#'     \item{\code{"use_ptm"}}{Copy the raw (uncorrected) PTM values into the
#'       normalized assay for unmatched rows. Use with caution: these values
#'       are not background-corrected and will be on a different scale.}
#'   }
#' @param return_unmatched Logical. If \code{TRUE} (default), the function
#'   returns a named list with components \code{norm_se} (the normalized
#'   \code{SummarizedExperiment}) and \code{unmatched_ptm} (a
#'   \code{SummarizedExperiment} of PTM rows for which no background protein
#'   was found). If \code{FALSE}, returns only \code{norm_se} directly.
#'
#' @return
#' If \code{return_unmatched = TRUE}: a named \code{list} with:
#' \describe{
#'   \item{\code{norm_se}}{A \code{SummarizedExperiment} with two assays:
#'     \code{norm_spectronaut_assay_ptm} (normalized values) and
#'     \code{spectronaut_assay_ptm} (original PTM values), plus a
#'     \code{bkg_match_found} logical column appended to \code{rowData}.}
#'   \item{\code{unmatched_ptm}}{A \code{SummarizedExperiment} containing the
#'     subset of \code{ptm_se} rows with no background match.}
#' }
#' If \code{return_unmatched = FALSE}: the \code{norm_se}
#' \code{SummarizedExperiment} directly.
#'
#' @details
#' **Matching logic.** Each PTM row is matched to its parent protein by
#' comparing \code{rowData(ptm_se)[[id_ptm]]} against
#' \code{rowData(bkg_se)[[id_bkg]]} using \code{match()}, which returns the
#' first hit. If a gene appears multiple times in \code{bkg_se}, only the first
#' row is used as the background reference.
#'
#' **Column alignment.** Only columns (samples) present in both \code{bkg_se}
#' and \code{ptm_se} are used. Samples exclusive to either input are dropped
#' with a warning. Alignment is performed on alphabetically sorted column
#' names; the original PTM column order is restored in the output.
#'
#' **Background medians of zero.** If the median background abundance for any
#' row is zero, division will produce \code{Inf} values. A warning is issued,
#' but computation proceeds. Consider filtering such rows upstream.
#'
#' **\code{bkg_match_found} flag.** A logical column is appended to
#' \code{rowData} of the output \code{SummarizedExperiment} indicating whether
#' a background protein was found for each PTM row (\code{TRUE}) or not
#' (\code{FALSE}).
#'
#' @section Normalization formula:
#' For matched PTM row \eqn{i} and sample column \eqn{j}:
#' \deqn{
#'   \hat{x}_{ij} = \frac{PTM_{ij}}{BKG_{ij}} \div
#'   \frac{median_i(PTM)}{median_i(BKG)} \times median_i(PTM)
#' }
#' Medians are computed across samples ignoring \code{NA}s.
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}},
#' \code{\link[SummarizedExperiment]{assay}},
#' \code{\link[SummarizedExperiment]{rowData}}
#'
#' @references
#' Tuechler N, Burtscher ML, et al. (2025). Dynamic multi-omics and mechanistic
#' modeling approach uncovers novel mechanisms of kidney fibrosis progression.
#' \emph{Molecular Systems Biology}, 21(8):1030-1065.
#' \doi{10.1038/s44320-025-00116-2}
#'
#' GitHub: \url{https://github.com/saezlab/kidneyfibrosis_multiomicsmodel_paper}
#'
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' result <- norm_ptm(
#'   bkg_se    = my_protein_se,
#'   ptm_se    = my_ptm_se,
#'   id_bkg    = "gene_name",
#'   id_ptm    = "gene_name",
#'   assay_bkg = "spectronaut_assay_prot",
#'   assay_ptm = "spectronaut_assay_ptm",
#'   bkg_na    = "set_as_na",
#'   return_unmatched = TRUE
#' )
#'
#' norm_se       <- result$norm_se
#' unmatched_se  <- result$unmatched_ptm
#' }
#'
#' @export
norm_ptm <- function(
    bkg_se, ptm_se,
    id_bkg, id_ptm,
    assay_bkg, assay_ptm,
    bkg_na = c("set_as_na", "use_ptm"),
    return_unmatched = TRUE
) {
  # -- Argument validation ----------------------------------------------------
  
  # Normalize bkg_na to a single valid choice
  bkg_na <- match.arg(bkg_na)
  
  # Guard: rownames of ptm_se must be unique for safe alignment in Step 11
  if (anyDuplicated(rownames(ptm_se))) {
    stop(
      "Duplicate rownames detected in ptm_se. ",
      "Cannot safely align normalized matrix to PTM SE. ",
      "Please ensure rownames are unique before calling norm_ptm()."
    )
  }
  
  # Guard: requested assay names must exist in their respective SEs
  if (!assay_bkg %in% SummarizedExperiment::assayNames(bkg_se)) {
    stop(
      "Assay '", assay_bkg, "' not found in bkg_se. ",
      "Available assays: ", paste(SummarizedExperiment::assayNames(bkg_se), collapse = ", ")
    )
  }
  if (!assay_ptm %in% SummarizedExperiment::assayNames(ptm_se)) {
    stop(
      "Assay '", assay_ptm, "' not found in ptm_se. ",
      "Available assays: ", paste(SummarizedExperiment::assayNames(ptm_se), collapse = ", ")
    )
  }
  
  # Guard: id columns must exist in rowData
  if (!id_bkg %in% colnames(SummarizedExperiment::rowData(bkg_se))) {
    stop("Column '", id_bkg, "' not found in rowData(bkg_se).")
  }
  if (!id_ptm %in% colnames(SummarizedExperiment::rowData(ptm_se))) {
    stop("Column '", id_ptm, "' not found in rowData(ptm_se).")
  }
  
  # -- Step 1: Extract assay matrices ----------------------------------------
  mat_bkg <- SummarizedExperiment::assay(bkg_se, assay_bkg)
  mat_ptm <- SummarizedExperiment::assay(ptm_se, assay_ptm)
  
  # -- Step 2: Extract gene IDs for row alignment ----------------------------
  gene_bkg <- SummarizedExperiment::rowData(bkg_se)[[id_bkg]]
  gene_ptm <- SummarizedExperiment::rowData(ptm_se)[[id_ptm]]
  
  # Warn if background gene IDs are not unique -- only the first match is used
  if (anyDuplicated(gene_bkg)) {
    warning(
      sum(duplicated(gene_bkg)), " duplicate gene ID(s) found in rowData(bkg_se)[[id_bkg]]. ",
      "Only the first occurrence of each ID will be used as background reference."
    )
  }
  
  # -- Step 3: Set background rownames for fast lookup ------------------------
  rownames(mat_bkg) <- gene_bkg
  
  # -- Step 4: Match PTM rows to background protein rows --------------------
  # match() returns NA for unmatched rows -- recorded for downstream use
  idx <- match(gene_ptm, gene_bkg)
  bkg_match_found <- !is.na(idx)
  
  if (!any(bkg_match_found)) {
    stop(
      "No PTM rows could be matched to any background protein. ",
      "Check that id_ptm and id_bkg refer to the same identifier type ",
      "(e.g., both gene symbols, both UniProt IDs)."
    )
  }
  
  # Prepare background matrix aligned row-wise with PTM (NA rows for no match)
  mat_bkg2 <- mat_bkg[idx, , drop = FALSE]
  rownames(mat_bkg2) <- rownames(mat_ptm)
  
  # -- Step 5: Column alignment ----------------------------------------------
  # Store original PTM column order for restoration after normalization
  orig_order <- colnames(mat_ptm)
  
  # Warn about columns present in PTM but absent from background (or vice versa)
  dropped_ptm <- setdiff(colnames(mat_ptm), colnames(mat_bkg2))
  dropped_bkg <- setdiff(colnames(mat_bkg2), colnames(mat_ptm))
  if (length(dropped_ptm) > 0) {
    warning(
      length(dropped_ptm), " PTM column(s) have no matching background column and will be dropped: ",
      paste(dropped_ptm, collapse = ", ")
    )
  }
  if (length(dropped_bkg) > 0) {
    warning(
      length(dropped_bkg), " background column(s) have no matching PTM column and will be ignored: ",
      paste(dropped_bkg, collapse = ", ")
    )
  }
  
  # Identify shared columns
  common_cols <- intersect(colnames(mat_ptm), colnames(mat_bkg2))
  if (length(common_cols) == 0) {
    stop("No matching columns between PTM and background assays.")
  }
  
  # Sort alphabetically for normalization (temporary alignment)
  sorted_cols <- sort(common_cols)
  mat_ptm_sorted <- mat_ptm[, sorted_cols, drop = FALSE]
  mat_bkg_sorted <- mat_bkg2[, sorted_cols, drop = FALSE]
  message("Columns matched alphabetically for normalization: ", paste(sorted_cols, collapse = ", "))
  
  # -- Step 6: Compute per-row medians ---------------------------------------
  med_ptm <- apply(mat_ptm_sorted, 1, median, na.rm = TRUE)
  med_bkg <- apply(mat_bkg_sorted, 1, median, na.rm = TRUE)
  
  # -- Step 7: Compute ratio_median normalization factor ---------------------
  ratio_median <- med_ptm / med_bkg
  
  # Warn if any background medians are zero -- will produce Inf in the normalized output
  zero_bkg_med <- bkg_match_found & !is.na(med_bkg) & (med_bkg == 0)
  if (any(zero_bkg_med, na.rm = TRUE)) {
    warning(
      sum(zero_bkg_med, na.rm = TRUE), " background row median(s) are exactly zero. ",
      "Division by zero will produce Inf values in the normalized assay. ",
      "Consider filtering these features before normalization."
    )
  }
  
  # -- Step 8: Initialize normalized output matrix ---------------------------
  norm_mat_sorted <- matrix(NA_real_, nrow = nrow(mat_ptm_sorted), ncol = ncol(mat_ptm_sorted))
  colnames(norm_mat_sorted) <- sorted_cols
  rownames(norm_mat_sorted) <- rownames(mat_ptm_sorted)
  
  # -- Step 9: Apply normalization formula (core -- do not modify) ------------
  # Formula: ((PTM_ij / BKG_ij) / (median_PTM_i / median_BKG_i)) * median_PTM_i
  # Note: for rows with no background match, med_bkg[i] is NA, so ratio_median[i]
  # is also NA, and the result for those rows remains NA regardless -- the
  # bkg_na handling below then decides what to put in those rows.
  for (i in seq_along(sorted_cols)) {
    norm_mat_sorted[, i] <- ((mat_ptm_sorted[, i] / mat_bkg_sorted[, i]) / ratio_median) * med_ptm
  }
  
  # -- Step 9b: Handle PTMs without background -------------------------------
  if (any(!bkg_match_found)) {
    n_missing <- sum(!bkg_match_found)
    message(n_missing, " PTM row(s) have no background match.")
    if (bkg_na == "use_ptm") {
      warning(
        "Using uncorrected PTM values for ", n_missing, " row(s) without background ",
        "(bkg_na = 'use_ptm'). These values are NOT background-corrected and will ",
        "be on a different scale than the normalized rows."
      )
      norm_mat_sorted[!bkg_match_found, ] <- mat_ptm_sorted[!bkg_match_found, ]
    } else {
      warning(
        "Leaving ", n_missing, " row(s) without background as NA ",
        "(bkg_na = 'set_as_na'). Check the 'bkg_match_found' column in ",
        "rowData() of the output to identify these features."
      )
    }
  }
  
  # -- Step 10: Restore original sample order --------------------------------
  available_cols <- intersect(orig_order, colnames(norm_mat_sorted))
  norm_mat <- norm_mat_sorted[, available_cols, drop = FALSE]
  message("Restored original sample order: ", paste(available_cols, collapse = ", "))
  
  # -- Step 11: Align rows to PTM SE -----------------------------------------
  if (!identical(rownames(norm_mat), rownames(ptm_se))) {
    warning("Row names differ -- aligning to PTM SE shared rows only.")
    shared_rows <- intersect(rownames(norm_mat), rownames(ptm_se))
    norm_mat <- norm_mat[shared_rows, , drop = FALSE]
    row_data_out <- SummarizedExperiment::rowData(ptm_se)[shared_rows, , drop = FALSE]
    bkg_match_found <- bkg_match_found[match(shared_rows, rownames(ptm_se))]
  } else {
    row_data_out <- SummarizedExperiment::rowData(ptm_se)
  }
  
  # Attach the match flag to rowData so downstream users can filter easily
  row_data_out$bkg_match_found <- bkg_match_found
  
  # -- Step 12: Construct the primary SummarizedExperiment output -------------
  # Include normalized assay (first) and original PTM assay (second)
  orig_assay <- SummarizedExperiment::assay(ptm_se, assay_ptm)
  orig_assay <- orig_assay[rownames(norm_mat), colnames(norm_mat), drop = FALSE]
  
  norm_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      norm_spectronaut_assay_ptm = norm_mat,
      spectronaut_assay_ptm      = orig_assay
    ),
    rowData = row_data_out,
    colData = SummarizedExperiment::colData(ptm_se),
    metadata = list(
      matched_columns     = common_cols,
      sorted_for_matching = sorted_cols,
      restored_order      = available_cols,
      bkg_na_handling     = bkg_na,
      n_ptm_matched       = sum(bkg_match_found),
      n_ptm_unmatched     = sum(!bkg_match_found),
      note = paste0(
        "Normalization performed on alphabetically sorted columns, then restored ",
        "to original order. Formula: ((PTM_ij / BKG_ij) / ratio_median_i) * med_ptm_i."
      )
    )
  )
  
  # -- Step 13: Optionally return unmatched PTM SE ----------------------------
  if (return_unmatched) {
    no_match_idx <- which(!bkg_match_found)
    unmatched_ptm_se <- if (length(no_match_idx)) ptm_se[no_match_idx, ] else ptm_se[FALSE, ]
    return(list(norm_se = norm_se, unmatched_ptm = unmatched_ptm_se))
  } else {
    return(norm_se)
  }
}

