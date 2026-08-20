# Progress and status messages.
#
# One toggle, no logging framework: everything chatty goes through .msg() or
# .done(), both of which are no-ops unless .verbose() is TRUE. The default is
# interactive(), so examples, tests and R CMD check stay silent without anyone
# having to remember to switch them off.
#
# Both helpers hand the caller's frame to cli, so message strings interpolate
# variables from the calling function: .msg("read {n} rows") works, and
# paste()/sprintf() at the call site is unnecessary.
#
# There is deliberately no .warn(): cli::cli_warn() already emits a real
# condition that tryCatch() and R CMD check can see, which a styled-alert
# wrapper would not.

.verbose <- function() {
  isTRUE(getOption("proteoSE.verbose", interactive()))
}

.msg <- function(...) {
  if (.verbose()) cli::cli_alert_info(..., .envir = parent.frame())
  invisible(NULL)
}

.done <- function(...) {
  if (.verbose()) cli::cli_alert_success(..., .envir = parent.frame())
  invisible(NULL)
}

#' Turn proteoSE progress messages on or off
#'
#' Import, enrichment and network functions report what they are doing while
#' they run. The messages are on at an interactive prompt and off everywhere
#' else (scripts, `R CMD check`, `knitr`); this function overrides that.
#'
#' It is a thin wrapper over `options(proteoSE.verbose = )`, which you can set
#' directly -- in `.Rprofile`, for example.
#'
#' @param on `TRUE` to print progress messages, `FALSE` to silence them, or
#'   `NULL` to restore the default (chatty when interactive, quiet otherwise).
#'
#' @return The previous value of the option, invisibly.
#' @export
#'
#' @examples
#' old <- proteoSE_verbose(FALSE)   # silence
#' proteoSE_verbose(old$proteoSE.verbose)
proteoSE_verbose <- function(on = TRUE) {
  stopifnot(is.null(on) || (is.logical(on) && length(on) == 1L))
  invisible(options(proteoSE.verbose = on))
}
