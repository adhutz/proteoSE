# Shiny application launchers

#' Launch the DEP2 LFQ Shiny application
#'
#' Starts the interactive Shiny app for label-free quantification (LFQ) analysis
#' with DEP2. The application lives in the installed package under
#' `inst/shiny/dep_lfq/` and is launched with [shiny::runApp()].
#'
#' @param ... Additional arguments passed on to [shiny::runApp()] (for example
#'   `launch.browser` or `port`).
#'
#' @return Invisibly, the value returned by [shiny::runApp()]. Called for its
#'   side effect of launching the app.
#'
#' @examplesIf interactive()
#' run_dep_lfq_app()
#'
#' @export
run_dep_lfq_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app. ",
         "Install it with install.packages('shiny').", call. = FALSE)
  }
  app_dir <- system.file("shiny", "dep_lfq", package = "sev")
  if (!nzchar(app_dir)) {
    stop("Could not locate the 'dep_lfq' Shiny app. Try re-installing 'sev'.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
