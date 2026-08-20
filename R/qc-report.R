# One-call QC report
#
# Nothing here plots anything new. The report is a parameterised R Markdown
# template that calls the DEP2 panels the worked example already walks through;
# this file only locates it, feeds it an SE and copies the result out.

#' Put the requested assay first
#'
#' Every DEP2 plot function reads `assay(se)`, i.e. the *first* assay, so
#' "report on the imputed values" means "move that assay to the front".
#'
#' @param se A `SummarizedExperiment`.
#' @param assay Assay name or index.
#' @return `se` with `assay` first.
#' @noRd
.qc_assay_first <- function(se, assay) {
  a <- SummarizedExperiment::assays(se, withDimnames = FALSE)
  i <- if (is.character(assay)) match(assay, names(a)) else suppressWarnings(as.integer(assay))
  if (length(i) != 1L || is.na(i) || i < 1L || i > length(a)) {
    stop("qc_report() has no assay ", dQuote(assay), ".\n",
         "  Available assays: ", .or_none(names(a)), ".", call. = FALSE)
  }
  SummarizedExperiment::assays(se, withDimnames = FALSE) <-
    a[c(i, setdiff(seq_along(a), i))]
  se
}

#' Render a quality-control report to a self-contained HTML file
#'
#' The standard opening sequence of a proteomics analysis -- identification
#' counts, missingness pattern, intensity distributions, replicate variability
#' and sample clustering -- rendered in one call so it can be emailed to
#' collaborators who do not run R.
#'
#' Every section degrades rather than fails: a panel whose prerequisite is
#' missing (no missing values to plot, no batch column, too few features for a
#' PCA) prints a short note and the report carries on.
#'
#' @param se A `SummarizedExperiment` with a `condition` column in `colData`.
#' @param group `colData` column used to colour and annotate samples in the PCA
#'   and correlation panels.
#' @param batch Optional second `colData` column, annotated alongside `group` --
#'   the quickest way to see whether the samples cluster by batch instead of by
#'   condition.
#' @param assay Assay to report on, by name or index. Defaults to the first,
#'   which for a proteoSE object is the raw intensity matrix; pass e.g.
#'   `"imputed_perseus"` to QC the imputed values instead.
#' @param output Path of the HTML file to write.
#' @param open Open the finished report in a browser.
#'
#' @return The path of the written report, invisibly.
#'
#' @examples
#' \dontrun{
#' se <- optimized_spectronaut_to_se(
#'   system.file("extdata", "example_project_report.tsv", package = "proteoSE")
#' )
#' qc_report(se, output = "qc.html")
#' }
#'
#' @export
qc_report <- function(se,
                      group  = "condition",
                      batch  = NULL,
                      assay  = 1,
                      output = "qc_report.html",
                      open   = interactive()) {

  .require_pkg("rmarkdown", "render the QC report")
  .assert_se(se, require_coldata = c(group, batch))
  se <- .qc_assay_first(se, assay)

  template <- system.file("rmarkdown", "qc_report.Rmd", package = "proteoSE")
  if (!nzchar(template)) {
    stop("qc_report() cannot find its template; reinstall proteoSE.", call. = FALSE)
  }

  # Render in a temp copy so knitr never writes into the installed package.
  wd <- tempfile("qc_report")
  dir.create(wd)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  file.copy(template, file.path(wd, "qc_report.Rmd"))

  out <- normalizePath(output, winslash = "/", mustWork = FALSE)
  rmarkdown::render(
    file.path(wd, "qc_report.Rmd"),
    output_file = out,
    params      = list(se = se, group = group, batch = batch),
    envir       = new.env(parent = globalenv()),
    quiet       = TRUE
  )

  if (isTRUE(open)) utils::browseURL(out)
  invisible(out)
}
