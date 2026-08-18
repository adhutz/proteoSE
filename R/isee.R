# iSEE app launchers
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.

#' Create iSEE app with all custom panels
#'
#' Includes GeneAnnoTable, GOTable and RowDataTable. 
#'
#' @param se summarized experiment object
#'
#' @return an iSEE app
#' @importFrom iSEE iSEE RowDataTable
#' @importFrom iSEEu VolcanoPlot FeatureSetTable
#' @export
isee_mini <- function(se) {
  .assert_se(se)
  gat <- GeneAnnoTable(PanelWidth = 8L)
  got <- GOTable(PanelWidth = 8L)
  rst <- RowDataTable(PanelWidth = 12L)
  iSEE::iSEE(se, initial = list(got, rst, VolcanoPlot(PanelWidth = 6L), gat, FeatureSetTable(PanelWidth = 6L)))
}



#' Minimize results
#' 
#' Reduces the columns of rowData to the minimum. Names and IDs of genes and
#' proteins are retained together with t-test results including fold change,
#' p-value, adjusted p-value and if the result is significant (adj. p-value <0.05).
#'
#' @param se summarized Experiment
#'
#' @return se summarized Experiment with reduced rowData
#' @importFrom dplyr ends_with
#' @export
make_mini <- function(se) {
  rowData(se) <- rowData(se) %>%
    as.data.frame() %>%
    select(gene_names,
           protein_ids,
           majority_protein_ids,
           protein_names,
           dplyr::ends_with("_diff|_p.val|p.adj|significant"))

  return(se)
}



#' Make se compatible to iSEE 
#' 
#' Transform summarized experiment. The result is a summarized experiment
#' of the type "SingleCellExperiment" to allow dimensional reduction plots.
#' GO-enrichment is registered from the metadata as a feature set collection and
#' p-values and differences are registered to be plotted via volcano plots.
#' #'
#' @param se summarized experiment
#' @param PValuePatterns Character vector of column-name patterns identifying
#'   p-value columns to register for volcano plots (default \code{"p.val"}).
#' @param LogFCPatterns Character vector of column-name patterns identifying
#'   log fold-change columns to register for volcano plots (default \code{"_diff"}).
#'
#' @return se summarized experiment with registered features
#' @export
#' @importFrom iSEEu registerFeatureSetCollections registerPValuePatterns registerLogFCPatterns
#' @importFrom SingleCellExperiment reducedDims reducedDims<-
#' @importFrom stats prcomp
se_to_isee <- function(se, PValuePatterns = "p.val", LogFCPatterns = "_diff"){
  .assert_se(se)
  se <- as(se, "SingleCellExperiment")

  if("GO_enrichment" %in% names(metadata(se))){
  # Add GO-term enrichment as FeatureSet to allow plotting
  se <- registerFeatureSetCollections(se, metadata(se)$GO_enrichment)
  }
  
  if(!any(is.na(assay(se)))){
  # Calculate pca data
  pca_data <- prcomp(t(assay(se)), rank=50)

  # Add PCA data to the experiment object
  reducedDims(se)<-list(PCA=pca_data$x)
  }
  
  #Register columns that contain p-values and log differences (for volcano mainly)
  se<-registerPValuePatterns(se, PValuePatterns)
  se<-registerLogFCPatterns(se, LogFCPatterns)

}

