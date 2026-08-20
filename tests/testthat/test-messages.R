test_that(".msg() and .done() respect the verbose option", {
  withr_opts <- options(proteoSE.verbose = TRUE)
  on.exit(options(withr_opts), add = TRUE)

  expect_message(.msg("hello"), "hello")
  expect_message(.done("done"), "done")

  options(proteoSE.verbose = FALSE)
  expect_silent(.msg("hello"))
  expect_silent(.done("done"))
})

test_that(".msg() interpolates from the calling function's frame", {
  old <- options(proteoSE.verbose = TRUE)
  on.exit(options(old), add = TRUE)

  f <- function() {
    n <- 7
    .msg("{n} rows")
  }
  expect_message(f(), "7 rows")
})

test_that("proteoSE_verbose() sets and restores the option", {
  old <- proteoSE_verbose(FALSE)
  on.exit(options(old), add = TRUE)
  expect_false(.verbose())

  proteoSE_verbose(TRUE)
  expect_true(.verbose())

  proteoSE_verbose(NULL)          # back to the interactive() default
  expect_identical(.verbose(), interactive())

  expect_error(proteoSE_verbose("yes"))
})

test_that("instrumented functions are silent when verbose is off", {
  # The property that keeps R CMD check and knitr quiet. genes_from_kegg()
  # is instrumented and hits the network, so poke the guard directly instead.
  old <- options(proteoSE.verbose = FALSE)
  on.exit(options(old), add = TRUE)
  expect_silent(.msg("Fetching KEGG pathway {.val hsa00190}"))
})

test_that("find_kws() says so when nothing matched", {
  # Real UniProt call: useful locally, not worth a flaky CI failure.
  skip_on_ci()
  skip_if_offline()
  expect_error(
    find_kws("zzzzznotakeyword"),
    "No UniProt keywords matched"
  )
})
