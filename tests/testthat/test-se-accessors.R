test_that("add_assay puts the new assay first and retains the old one", {
  se <- make_test_se(n_features = 5)
  original <- SummarizedExperiment::assay(se, "intensity")
  new_mat <- log2(abs(original) + 1)

  se2 <- add_assay(se, new_mat, "log2")

  # New assay is first and named correctly; old assay retained under its name
  expect_equal(SummarizedExperiment::assayNames(se2)[1], "log2")
  expect_true("intensity" %in% SummarizedExperiment::assayNames(se2))
  expect_equal(length(SummarizedExperiment::assays(se2)), 2)

  # Values are preserved
  expect_equal(SummarizedExperiment::assay(se2, "log2"), new_mat)
  expect_equal(SummarizedExperiment::assay(se2, "intensity"), original)
})

test_that("add_assay does not change the dimensions of the object", {
  se <- make_test_se(n_features = 8)
  new_mat <- SummarizedExperiment::assay(se) * 2
  se2 <- add_assay(se, new_mat, "scaled")
  expect_equal(dim(se2), dim(se))
})
