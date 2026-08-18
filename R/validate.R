# Input validation
#
# One checker with two faces: validate_se() for users, .assert_se() for the
# package's own entry points. The value is entirely in the error message --
# it names the function, what it wanted, and what the object actually has.

#' Contrast prefixes present in a SummarizedExperiment
#'
#' Differential-testing results follow the `<contrast>_<suffix>` convention
#' (`_diff`, `_p.val`, `_p.adj`, `_significant`). The `_diff` columns are the
#' authoritative list, since every workflow that writes results writes one.
#'
#' @param se A `SummarizedExperiment`.
#' @return Character vector of contrast prefixes, possibly empty.
#' @noRd
.se_contrasts <- function(se) {
  cols <- colnames(SummarizedExperiment::rowData(se))
  sub("_diff$", "", grep("_diff$", cols, value = TRUE))
}

# "a, b and c", or "none" for an empty vector.
.or_none <- function(x) {
  if (length(x) == 0) "none" else paste(x, collapse = ", ")
}

#' Check that a SummarizedExperiment carries what a function needs
#'
#' Fails early with a message that names the missing column, assay or contrast
#' *and* lists what the object does have — instead of the
#' `undefined columns selected` error three frames deep that you get otherwise.
#'
#' Requirements are additive and all optional; an argument left at its default
#' checks nothing. Only ask for what the calling function genuinely needs: an
#' object that is missing `replicate` still plots fine, so do not require it.
#'
#' @param se Object to check.
#' @param require_rowdata Character vector of column names that must be in
#'   `rowData(se)`.
#' @param require_coldata Character vector of column names that must be in
#'   `colData(se)`.
#' @param require_assays Character vector of names that must be in
#'   `assayNames(se)`.
#' @param require_contrast Character vector of contrast prefixes (e.g.
#'   `"treat_vs_ctrl"`) that must have differential-testing results in
#'   `rowData(se)`.
#' @param class Class `se` must inherit from. Use `"PhosphoExperiment"` for the
#'   PhosR-based functions; it inherits from `SummarizedExperiment`, so the same
#'   checks apply.
#' @param .caller Name of the function to blame in the error message. Internal;
#'   the package's own `.assert_se()` fills it in.
#'
#' @return `se`, invisibly, so the call composes:
#'   `se <- validate_se(se, require_coldata = "condition")`.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(intensity = matrix(1:6, nrow = 3)),
#'   rowData = S4Vectors::DataFrame(gene_names = c("A", "B", "C")),
#'   colData = S4Vectors::DataFrame(condition = c("ctrl", "treat"))
#' )
#' validate_se(se, require_rowdata = "gene_names", require_coldata = "condition")
#'
#' # Tells you what is there, not just what is missing:
#' try(validate_se(se, require_assays = "imputed"))
#'
#' @export
validate_se <- function(se,
                        require_rowdata  = character(),
                        require_coldata  = character(),
                        require_assays   = character(),
                        require_contrast = NULL,
                        class            = "SummarizedExperiment",
                        .caller          = NULL) {

  who <- if (is.null(.caller)) "This object " else paste0(.caller, "() ")
  bail <- function(...) stop(who, ..., call. = FALSE)

  if (!methods::is(se, class)) {
    bail("needs a ", class, " object; got ",
         paste(class(se), collapse = "/"), ".")
  }

  missing_row <- setdiff(require_rowdata, colnames(SummarizedExperiment::rowData(se)))
  if (length(missing_row)) {
    bail("needs rowData column(s) ", .or_none(missing_row), ".\n",
         "  rowData has: ", .or_none(colnames(SummarizedExperiment::rowData(se))), ".")
  }

  missing_col <- setdiff(require_coldata, colnames(SummarizedExperiment::colData(se)))
  if (length(missing_col)) {
    bail("needs colData column(s) ", .or_none(missing_col), ".\n",
         "  colData has: ", .or_none(colnames(SummarizedExperiment::colData(se))), ".")
  }

  missing_assay <- setdiff(require_assays, SummarizedExperiment::assayNames(se))
  if (length(missing_assay)) {
    bail("needs assay(s) ", .or_none(missing_assay), ".\n",
         "  Available assays: ", .or_none(SummarizedExperiment::assayNames(se)), ".")
  }

  if (!is.null(require_contrast)) {
    have <- .se_contrasts(se)
    missing_contrast <- setdiff(require_contrast, have)
    if (length(missing_contrast)) {
      bail("needs a contrast ", .or_none(missing_contrast), ".\n",
           "  rowData has _diff columns for: ", .or_none(have), ".")
    }
  }

  invisible(se)
}

#' Internal front-end to validate_se() that blames the calling function
#'
#' @inheritParams validate_se
#' @return `se`, invisibly.
#' @noRd
.assert_se <- function(se, ...) {
  caller <- tryCatch(deparse(sys.call(-1)[[1]]), error = function(e) NULL)
  validate_se(se, ..., .caller = caller)
}
