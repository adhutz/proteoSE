# plot_antigen() had never been called with its own defaults, and all three of
# its defect were in that path: `min.segment.length` was used in the body but
# was not a formal (only a globalVariables() entry), the `targets` default was
# an unevaluable data.frame(), and `if (additional_sets == "all")` errored for
# any vector of length != 1 -- the documented multi-value usage.

# plot_antigen() matches samples with "^<condition>_[0-9]+$".
antigen_se <- function() {
  se <- make_test_se(20)
  colnames(se) <- paste(SummarizedExperiment::colData(se)$condition,
                        SummarizedExperiment::colData(se)$replicate, sep = "_")
  SummarizedExperiment::colData(se)$sample <- colnames(se)
  SummarizedExperiment::colData(se)$ID     <- colnames(se)
  rd <- SummarizedExperiment::rowData(se)
  rd$treat_vs_ctrl_diff  <- seq_len(nrow(se)) / 10
  rd$treat_vs_ctrl_p.adj <- 0.01
  SummarizedExperiment::rowData(se) <- rd
  se
}

test_that("plot_antigen() runs on its own defaults", {
  p <- suppressWarnings(suppressMessages(
    plot_antigen(antigen_se(), "treat_vs_ctrl")))
  expect_s3_class(p, "ggplot")
})

test_that("plot_antigen() accepts additional_sets as a vector, not just a keyword", {
  se <- antigen_se()
  expect_s3_class(suppressWarnings(suppressMessages(
    plot_antigen(se, "treat_vs_ctrl", additional_sets = "all"))), "ggplot")
  expect_s3_class(suppressWarnings(suppressMessages(
    plot_antigen(se, "treat_vs_ctrl", additional_sets = c("ctrl", "treat")))), "ggplot")
})

test_that("plot_antigen()'s targets default is evaluable", {
  expect_silent(eval(formals(plot_antigen)$targets))
})

test_that("plot_antigen() rejects an unknown contrast readably", {
  expect_error(plot_antigen(antigen_se(), "nonexistent"),
               "plot_antigen\\(\\) needs a contrast nonexistent")
})
