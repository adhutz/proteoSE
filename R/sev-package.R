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
#' @section Namespace conflict policy:
#' **dplyr wins every name collision.** Where another imported package exports a
#' name dplyr also exports, dplyr's version is the one bound in the sev2
#' namespace, and the other package's version must be called fully qualified
#' (e.g. `AnnotationDbi::select()`, `MSnbase::combine()`). Two consequences for
#' anyone editing `R/`:
#'
#' - Never `@importFrom` a name that dplyr also exports. Qualify the call site
#'   instead. Relying on NAMESPACE ordering to break the tie is fragile — the
#'   file is alphabetically generated, so the winner can flip when an import is
#'   added.
#' - The same applies inside strings that are later `parse()`d and `eval()`ed
#'   (see the iSEE panel command blocks in `R/isee-panels.R`): they resolve
#'   through this namespace, so a bare `select()` there means `dplyr::select()`.
#'
#' @keywords internal
#' @importFrom grDevices dev.off pdf
#' @importFrom stats aggregate formula lm median na.omit reorder rnorm runif sd
#'   setNames terms.formula
#' @importFrom utils data install.packages installed.packages read.delim
#'   type.convert write.table
"_PACKAGE"
