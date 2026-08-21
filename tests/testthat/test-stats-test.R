# test_diff_limma() had three defects that only a default call reached:
# fdr.type was never match.arg()ed, a missing contrast fell through to
# makeContrasts(NULL), and an unnamed contrast deleted the comparison column.

test_that("test_diff_limma() says what it needs when given no contrast", {
  expect_error(test_diff_limma(make_test_se()), "needs at least one contrast")
  expect_error(test_diff_limma(make_test_se()), "ctrl, treat")
})

test_that("test_diff_limma() names an unnamed contrast <a>_vs_<b>", {
  res <- test_diff_limma(make_test_se(n_features = 300),
                         advanced_contrast = "treat - ctrl",
                         fdr.type = "BH")
  cols <- colnames(SummarizedExperiment::rowData(res))

  expect_true(all(paste0("treat_vs_ctrl", c("_diff", "_p.val", "_p.adj")) %in%
                    cols))
  # ... and an explicit name still wins
  named <- test_diff_limma(make_test_se(n_features = 300),
                           advanced_contrast = c(mine = "treat - ctrl"),
                           fdr.type = "BH")
  expect_true("mine_diff" %in% colnames(SummarizedExperiment::rowData(named)))
})

test_that("test_diff_limma() runs on its default fdr.type", {
  res <- suppressWarnings(
    test_diff_limma(make_test_se(n_features = 300),
                    advanced_contrast = "treat - ctrl")
  )
  expect_true("treat_vs_ctrl_p.adj" %in%
                colnames(SummarizedExperiment::rowData(res)))
})
