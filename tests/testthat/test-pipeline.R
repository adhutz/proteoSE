# The pipeline's mechanics are proteomics-free -- it folds do.call() over a
# list -- so these use base functions as steps. Whether filter_perseus() works
# is filter_perseus()'s tests' problem.

test_that("building a pipeline executes nothing", {
  fired <- FALSE
  side_effect <- function(x) {
    fired <<- TRUE
    x
  }

  pipe <- proteoSE_pipeline(1:10) |>
    add_step("side_effect") |>
    add_step("head", n = 3)

  expect_false(fired)
  expect_length(pipe$steps, 2)
  expect_null(pipe$result)
})

test_that("a step defined in the caller's script is found at run time", {
  double_it <- function(x) x * 2

  pipe <- proteoSE_pipeline(1:3) |> add_step("double_it")
  expect_equal(pipeline_result(run_pipeline(pipe)), c(2, 4, 6))
})

test_that("add_step() rejects a typo immediately", {
  expect_error(add_step(proteoSE_pipeline(1:10), "no_such_functon"),
               "'no_such_functon' of mode 'function' was not found")
  expect_error(add_step(proteoSE_pipeline(1:10), head),
               "function \\*name\\* as a string")
  expect_error(add_step(1:10, "head"), "needs a proteoSE_pipeline")
})

test_that("run_pipeline() gives the same answer as calling the steps by hand", {
  pipe <- proteoSE_pipeline(1:10) |>
    add_step("head", n = 3) |>
    add_step("rev") |>
    add_step("sum")

  expect_equal(pipeline_result(run_pipeline(pipe)), sum(rev(head(1:10, 3))))
})

test_that("run_pipeline() records when it ran and in what session", {
  pipe <- run_pipeline(add_step(proteoSE_pipeline(1:10), "rev"))
  expect_s3_class(pipe$run_at, "POSIXct")
  expect_s3_class(pipe$session, "sessionInfo")
})

test_that("a pipeline round-trips through saveRDS()", {
  pipe <- proteoSE_pipeline(1:10) |>
    add_step("head", n = 3) |>
    add_step("rev")

  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(pipe, f)
  restored <- readRDS(f)

  expect_equal(restored, pipe)
  expect_equal(pipeline_result(run_pipeline(restored)),
               pipeline_result(run_pipeline(pipe)))
})

test_that("an unrun pipeline has no result to hand out", {
  expect_error(pipeline_result(proteoSE_pipeline(1:10)), "has not been run")
  expect_error(run_pipeline(proteoSE_pipeline(1:10)), "at least one step")
})

test_that("a failing step is named in the error", {
  pipe <- proteoSE_pipeline(1:10) |>
    add_step("rev") |>
    add_step("log", base = "not a number")

  expect_error(run_pipeline(pipe), "step 2 \\(log\\(\\)\\)")
})

test_that("print() works on a pipeline in either state", {
  empty <- proteoSE_pipeline(1:10)
  expect_output(print(empty), "none yet")
  expect_output(print(empty), "not run")

  pipe <- add_step(empty, "head", n = 3)
  expect_output(print(pipe), "head\\(n = 3\\)")
  expect_output(print(run_pipeline(pipe)), "status: run")
})

test_that("write_methods() reports every step and its parameters", {
  # filter_perseus()'s filter_mode is left at its default, and must still be
  # reported: a methods section states the value used, not the value typed.
  txt <- paste(
    write_methods(
      proteoSE_pipeline(1:10) |>
        add_step("filter_perseus", perc_na = 0.25) |>
        add_step("head", n = 3)
    ),
    collapse = " "
  )

  expect_match(txt, "75% of samples")
  expect_match(txt, "each_group", fixed = TRUE)
  expect_match(txt, "head\\(\\) was applied with n = 3")
  expect_match(txt, R.version.string, fixed = TRUE)
})

# A template that names an argument the function does not have silently
# paste0()s a NULL, i.e. nothing -- "downshift  SD". Catch the gap, not the
# prose.
test_that("every methods template fills in all its parameters", {
  for (fn in names(.methods_templates)) {
    sentence <- .methods_templates[[fn]](.with_defaults(fn, list()))
    expect_no_match(sentence, "  |NULL|NA\b", label = fn)
  }
})

test_that("write_methods() writes to a file when asked", {
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  txt <- write_methods(add_step(proteoSE_pipeline(1:10), "rev"), file = f)
  expect_equal(readLines(f), unlist(strsplit(txt, "\n")))
})
