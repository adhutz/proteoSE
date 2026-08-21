# Reproducible pipeline object
#
# A pipeline is a recipe, not a scheduler: a list of
# `list(fn = "filter_perseus", args = list(perc_na = 0.33))` plus the input
# object. run_pipeline() folds do.call() over it. That is enough for the three
# things the object is for -- provenance, replay, and an auto-generated methods
# paragraph -- and it serialises with saveRDS() for free.
#
# ponytail: no `targets` backend, no dependency graph, no step store. This is a
# linear chain of 4-6 steps. Revisit only if intermediate caching measurably
# matters, i.e. if someone complains about re-running an import.
#
# ponytail: one generic add_step() rather than a step_filter()/step_impute()/
# step_test()/step_enrich_go() family. The wrappers were ~40 lines of
# boilerplate that also needed extending for every new pipeline-able function;
# write_methods() gets its per-function prose from .methods_templates instead,
# and falls back to a generic sentence for anything not listed there.

#' Build a reproducible analysis pipeline
#'
#' Records the functions and parameters of an analysis instead of running them,
#' so that the recipe can be inspected, saved, replayed, and turned into a
#' methods paragraph. A `SummarizedExperiment` already carries the *data* a
#' pipeline produced; what it does not carry is which functions with which
#' parameters produced it.
#'
#' Steps are added with [add_step()] and nothing executes until
#' [run_pipeline()] is called. The object is a plain list with an S3 class, so
#' `saveRDS()` / `readRDS()` round-trips it.
#'
#' @param input The object the pipeline starts from -- typically the
#'   `SummarizedExperiment` returned by [se_read_in()]. Not validated: whatever
#'   the first step accepts is fair game.
#'
#' @return An object of class `proteoSE_pipeline`.
#' @seealso [add_step()], [run_pipeline()], [write_methods()]
#' @export
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(intensity = matrix(1:6, nrow = 3)),
#'   rowData = S4Vectors::DataFrame(gene_names = c("A", "B", "C")),
#'   colData = S4Vectors::DataFrame(condition = c("ctrl", "treat"))
#' )
#'
#' pipe <- proteoSE_pipeline(se) |>
#'   add_step("filter_perseus", perc_na = 0.33) |>
#'   add_step("test_diff_limma")
#' pipe
#'
#' write_methods(pipe)
proteoSE_pipeline <- function(input) {
  structure(
    list(input = input, steps = list(),
         result = NULL, run_at = NULL, session = NULL),
    class = "proteoSE_pipeline"
  )
}

#' Add a step to a pipeline
#'
#' Appends a function call to the recipe and returns the pipeline, so steps
#' chain with `|>`. Nothing is executed here -- see [run_pipeline()].
#'
#' Arguments are evaluated and stored by value, which is what makes the
#' pipeline a provenance record: `add_step("filter_perseus", perc_na = cutoff)`
#' stores the number `cutoff` held when the step was added, not the name
#' `cutoff`.
#'
#' The function is looked up immediately, so a typo fails while you are
#' building the pipeline rather than an hour into a run.
#'
#' @param pipe A `proteoSE_pipeline` from [proteoSE_pipeline()].
#' @param fn Name of the function to call, as a string, e.g.
#'   `"filter_perseus"`. It receives the current object as its first argument.
#' @param ... Further arguments to `fn`. Named, and evaluated now.
#'
#' @return The pipeline, with one more step.
#' @export
#'
#' @examples
#' proteoSE_pipeline(1:10) |>
#'   add_step("head", n = 3) |>
#'   add_step("rev")
add_step <- function(pipe, fn, ...) {
  .assert_pipeline(pipe, "add_step")
  if (!is.character(fn) || length(fn) != 1L) {
    stop("add_step() needs the function *name* as a string, e.g. ",
         'add_step(pipe, "filter_perseus").', call. = FALSE)
  }
  match.fun(fn)   # fail now, not an hour into run_pipeline()
  pipe$steps <- c(pipe$steps, list(list(fn = fn, args = list(...))))
  pipe
}

#' Run a pipeline
#'
#' Folds the recorded steps over the input object, passing the result of each
#' step as the first argument of the next. Only the final object is kept;
#' intermediates are discarded.
#'
#' The result, the time of the run and `utils::sessionInfo()` are stored on the
#' returned pipeline, so a saved pipeline is a complete record of what was run
#' and where.
#'
#' @param pipe A `proteoSE_pipeline` with at least one step.
#'
#' @return The pipeline, with `result`, `run_at` and `session` filled in. Use
#'   [pipeline_result()] to get the analysed object itself.
#' @export
#'
#' @examples
#' pipe <- proteoSE_pipeline(1:10) |>
#'   add_step("head", n = 3) |>
#'   add_step("rev")
#' pipeline_result(run_pipeline(pipe))
run_pipeline <- function(pipe) {
  .assert_pipeline(pipe, "run_pipeline")
  if (length(pipe$steps) == 0) {
    stop("run_pipeline() needs at least one step; add one with add_step().",
         call. = FALSE)
  }

  obj <- pipe$input
  for (i in seq_along(pipe$steps)) {
    s <- pipe$steps[[i]]
    .msg("Step {i}/{length(pipe$steps)}: {s$fn}()")
    f <- match.fun(s$fn)
    obj <- tryCatch(
      do.call(f, c(list(obj), s$args)),
      error = function(e) {
        stop("pipeline step ", i, " (", s$fn, "()) failed: ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }

  pipe$result  <- obj
  pipe$run_at  <- Sys.time()
  pipe$session <- utils::sessionInfo()
  .done("Pipeline finished ({length(pipe$steps)} steps).")
  pipe
}

#' The object a pipeline produced
#'
#' @param pipe A `proteoSE_pipeline` that has been run.
#' @return The result of the final step.
#' @export
#' @examples
#' pipe <- run_pipeline(add_step(proteoSE_pipeline(1:10), "rev"))
#' pipeline_result(pipe)
pipeline_result <- function(pipe) {
  .assert_pipeline(pipe, "pipeline_result")
  if (is.null(pipe$run_at)) {
    stop("pipeline_result(): this pipeline has not been run yet. ",
         "Call run_pipeline() first.", call. = FALSE)
  }
  pipe$result
}

#' @param x A `proteoSE_pipeline`.
#' @param ... Ignored.
#' @rdname proteoSE_pipeline
#' @export
print.proteoSE_pipeline <- function(x, ...) {
  cat("<proteoSE pipeline>\n")
  cat("  input:  ", .describe(x$input), "\n", sep = "")
  if (length(x$steps) == 0) {
    cat("  steps:  none yet -- add one with add_step()\n")
  } else {
    cat("  steps:\n")
    for (i in seq_along(x$steps)) {
      s <- x$steps[[i]]
      cat("   ", i, ". ", s$fn, "(", .fmt_args(s$args), ")\n", sep = "")
    }
  }
  if (is.null(x$run_at)) {
    cat("  status: not run\n")
  } else {
    cat("  status: run ", format(x$run_at, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("  result: ", .describe(x$result), "\n", sep = "")
  }
  invisible(x)
}

#' Draft the methods paragraph for a pipeline
#'
#' Turns the recipe into prose: one sentence per step with the parameters
#' filled in, followed by the software versions the steps actually came from.
#' It is a template fill, not a substitute for reading it -- check the wording
#' and the numbers before they go into a manuscript.
#'
#' Parameters left at their defaults are reported too, because a methods
#' section has to state the threshold that was used, not the one that was
#' typed. The pipeline does not need to have been run.
#'
#' @param pipe A `proteoSE_pipeline`.
#' @param file Optional path to write the text to. When `NULL` (default) the
#'   text is printed.
#'
#' @return The text, as a character vector of paragraphs, invisibly.
#' @export
#'
#' @examples
#' proteoSE_pipeline(1:10) |>
#'   add_step("head", n = 3) |>
#'   write_methods()
write_methods <- function(pipe, file = NULL) {
  .assert_pipeline(pipe, "write_methods")

  sentences <- vapply(pipe$steps, function(s) {
    args <- .with_defaults(s$fn, s$args)
    tmpl <- .methods_templates[[s$fn]]
    if (is.null(tmpl)) {
      paste0(s$fn, "() was applied",
             if (length(args)) paste0(" with ", .fmt_args(args)) else "", ".")
    } else {
      tmpl(args)
    }
  }, character(1))

  txt <- c(
    if (length(sentences)) paste(sentences, collapse = " ") else
      "This pipeline has no steps.",
    .versions_sentence(pipe)
  )

  if (is.null(file)) {
    cat(txt, sep = "\n\n")
    cat("\n")
  } else {
    writeLines(txt, file)
    .done("Methods text written to {file}.")
  }
  invisible(txt)
}


# -- internals ---------------------------------------------------------------

.assert_pipeline <- function(pipe, caller) {
  if (!inherits(pipe, "proteoSE_pipeline")) {
    stop(caller, "() needs a proteoSE_pipeline object; got ",
         paste(class(pipe), collapse = "/"),
         ". Start one with proteoSE_pipeline().", call. = FALSE)
  }
  invisible(pipe)
}

# "SummarizedExperiment [1808 x 12]" / "integer [10]"
.describe <- function(x) {
  dims <- if (!is.null(dim(x))) paste(dim(x), collapse = " x ") else length(x)
  paste0(class(x)[1], " [", dims, "]")
}

# "perc_na = 0.33, filter_mode = \"each_group\""
.fmt_args <- function(args) {
  if (length(args) == 0) return("")
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  vals <- vapply(args, function(v) paste(deparse(v, width.cutoff = 200L),
                                         collapse = " "), character(1))
  paste0(ifelse(nzchar(nms), paste0(nms, " = "), ""), vals, collapse = ", ")
}

# Supplied arguments on top of the function's own defaults, minus the first
# argument (the object being piped). Defaults that are missing, or that are a
# match.arg() vector, are reduced to the value the function would actually use.
.with_defaults <- function(fn, args) {
  f    <- match.fun(fn)
  fmls <- formals(f)
  fmls <- fmls[setdiff(names(fmls)[-1], "...")]   # drop the piped object and ...
  env  <- environment(f)
  if (is.null(env)) env <- baseenv()
  defaults <- lapply(names(fmls), function(nm) {
    d <- fmls[[nm]]
    if (identical(d, quote(expr = ))) return(NULL)          # no default
    v <- tryCatch(eval(d, env), error = function(e) d)
    if (length(v) > 1 && is.character(v)) v[[1]] else v     # match.arg() vector
  })
  names(defaults) <- names(fmls)
  defaults <- defaults[!vapply(defaults, is.null, logical(1))]
  utils::modifyList(defaults, args)
}

# One sentence per pipeline-able subsystem. Anything not listed here gets the
# generic "f() was applied with a = 1" fallback, which is why there is no
# pressure to keep this exhaustive.
.methods_templates <- list(
  filter_perseus = function(a) paste0(
    "Features were filtered to those quantified in at least ",
    round(100 * (1 - a$perc_na)), "% of samples (filter mode: ",
    a$filter_mode, ")."),
  impute_perseus = function(a) paste0(
    "Missing values were imputed from a downshifted normal distribution ",
    "(width ", a$width, " SD, downshift ", a$downshift, " SD), as ",
    "implemented in ",
    "Perseus."),
  impute_DEP2 = function(a) paste0(
    "Missing values were imputed using the ", a$fun, " method."),
  test_diff_limma = function(a) paste0(
    "Differential abundance was tested with limma using the design ",
    paste(deparse(a$design_formula), collapse = " "), "."),
  enrich_go_se = function(a) paste0(
    "Gene Ontology enrichment (", a$ont, ") was performed against ", a$OrgDb,
    " at an adjusted p-value cutoff of ", a$pvalueCutoff, ".")
)

# "Analysis used R 4.5.3, proteoSE 0.3.0, DEP2 1.2.0."
.versions_sentence <- function(pipe) {
  pkgs <- vapply(pipe$steps, function(s) {
    environmentName(environment(match.fun(s$fn)))
  }, character(1))
  pkgs <- setdiff(unique(pkgs), c("R_GlobalEnv", "base", "", "proteoSE"))
  vers <- vapply(pkgs, function(p) {
    tryCatch(paste(p, utils::packageVersion(p)), error = function(e) p)
  }, character(1))
  # tryCatch: during devtools::load_all() the package is not installed under
  # this name and packageVersion() fails.
  self <- tryCatch(paste("proteoSE", utils::packageVersion("proteoSE")),
                   error = function(e) "proteoSE")
  paste0("Analysis was performed with ", R.version.string, ", ",
         paste(c(self, vers), collapse = ", "), ".")
}
