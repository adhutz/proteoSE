#' proteoSE: Proteomics and Phosphoproteomics Analysis on SummarizedExperiment
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
#' @section Caching:
#' The four functions that query an external API -- [get_network()] (STRING),
#' [genes_from_kegg()] (KEGG), [find_kws()] and [fetch_kw_accessions()]
#' (UniProt) -- are memoised on disk, so a repeated call with identical
#' arguments is served from the cache instead of the network. Entries live in
#' `tools::R_user_dir("proteoSE", "cache")`, expire after 30 days, and are
#' pruned once the cache exceeds 512 MB. Use [proteoSE_cache_info()] to see the
#' location and current size, and [proteoSE_cache_clear()] to empty it.
#'
#' Caching requires the suggested packages \pkg{memoise} and \pkg{cachem}; if
#' either is missing the functions simply fetch every time. Set
#' `options(proteoSE.cache = FALSE)` before loading the package to opt out.
#'
#' @section Progress messages:
#' Import, enrichment and network functions report what they are doing while
#' they run. The messages are on at an interactive prompt and off everywhere
#' else -- scripts, `knitr`, `R CMD check` -- so nothing prints where it would
#' only be noise. Override with [proteoSE_verbose()], or set
#' `options(proteoSE.verbose = TRUE/FALSE)` directly.
#'
#' Because the four cached functions above are memoised, their bodies do not
#' run on a cache hit; the absence of a "Fetching ..." message is how you tell
#' a cached answer from a live one.
#'
#' @section Namespace conflict policy:
#' **dplyr wins every name collision.** Where another imported package exports a
#' name dplyr also exports, dplyr's version is the one bound in the proteoSE
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
