#' sev: The SummarizedExperiment Viewer for Proteomics Analysis
#'
#' A utility toolkit for (phospho)proteomics data analysis organised around the
#' Bioconductor [SummarizedExperiment::SummarizedExperiment] class. It provides
#' convenience functions for data import, missing-value imputation,
#' differential-abundance testing, GO-term and gene-set enrichment, PhosR-based
#' signalome analysis, STRING network retrieval and publication-ready plotting,
#' built on top of DEP2, iSEE/iSEEu, PhosR, limma and clusterProfiler.
#'
#' The package source is grouped by theme (see the file names in `R/`). For an
#' overview of the refactoring history and outstanding cleanup tasks, see
#' `AUDIT.md` and `CHANGES_phase0-1.md` in the package root.
#'
#' @keywords internal
#' @importFrom grDevices dev.off pdf
#' @importFrom stats aggregate formula lm median na.omit reorder rnorm runif sd
#'   setNames terms.formula
#' @importFrom utils data install.packages installed.packages read.delim
#'   type.convert write.table
"_PACKAGE"
