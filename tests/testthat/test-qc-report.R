test_that(".qc_assay_first moves the requested assay to the front", {
  se <- make_test_se()
  SummarizedExperiment::assays(se)$imputed <- SummarizedExperiment::assay(se)

  expect_equal(SummarizedExperiment::assayNames(.qc_assay_first(se, "imputed")),
               c("imputed", "intensity"))
  expect_equal(SummarizedExperiment::assayNames(.qc_assay_first(se, 2)),
               c("imputed", "intensity"))
  expect_equal(SummarizedExperiment::assayNames(.qc_assay_first(se, 1)),
               c("intensity", "imputed"))
})

test_that(".qc_assay_first names the available assays when asked for a missing one", {
  se <- make_test_se()
  expect_error(.qc_assay_first(se, "imputed"), "Available assays: intensity")
  expect_error(.qc_assay_first(se, 3), "Available assays: intensity")
})

test_that("qc_report() rejects a missing grouping column before rendering", {
  expect_error(qc_report(make_test_se(), group = "batch_id"),
               "colData column\\(s\\) batch_id")
})

test_that("qc_report() renders a self-contained HTML file", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_on_cran()

  out <- tempfile(fileext = ".html")
  on.exit(unlink(out), add = TRUE)

  path <- qc_report(make_test_se(), output = out, open = FALSE)
  expect_true(file.exists(path))
  expect_gt(file.size(path), 10000)

  # The fixture has no missing values, so the missingness panels must degrade
  # rather than abort the render.
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "Not available")
  expect_match(html, "Sample relationships")
})
