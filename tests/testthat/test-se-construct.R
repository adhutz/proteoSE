# A minimal SE carrying a "<contrast>_diff" rowData column, as produced by the
# DEP2/limma test step that log2fc_to_fc() consumes.
make_diff_se <- function(diffs = c(1, -1, NA, 2)) {
  n <- length(diffs)
  mat <- matrix(stats::rnorm(n * 4), nrow = n)
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = mat),
    rowData = S4Vectors::DataFrame(
      name              = paste0("G", seq_len(n)),
      cond_vs_ctrl_diff = diffs
    )
  )
}

test_that("log2fc_to_fc adds 2^diff fold-change columns and preserves NA", {
  se  <- make_diff_se(c(1, -1, NA, 2))
  se2 <- suppressMessages(log2fc_to_fc(se))

  rd <- SummarizedExperiment::rowData(se2)
  expect_true("cond_vs_ctrl_foldchange" %in% colnames(rd))
  expect_equal(rd$cond_vs_ctrl_foldchange, c(2, 0.5, NA, 4))
  # original diff column is untouched
  expect_equal(rd$cond_vs_ctrl_diff, c(1, -1, NA, 2))
})

test_that("log2fc_to_fc validates its input", {
  expect_error(log2fc_to_fc(data.frame(x = 1)), "SummarizedExperiment")

  no_diff <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = matrix(1:4, nrow = 2)),
    rowData = S4Vectors::DataFrame(name = c("a", "b"))
  )
  expect_error(suppressMessages(log2fc_to_fc(no_diff)), "_diff")
})

# A minimal phospho-style SE carrying the per-contrast statistics columns that
# add_significance() pivots over: site_id_mult plus <contrast>_{diff,p.adj,
# p.val,CI.L,CI.R}. The pivot splits each name on its last underscore, so the
# contrast label is everything up to "_<stat>".
make_sig_se <- function() {
  n <- 3L
  mat <- matrix(stats::rnorm(n * 4), nrow = n)
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = mat),
    rowData = S4Vectors::DataFrame(
      site_id_mult       = c("s1", "s2", "s3"),
      cond_vs_ctrl_diff  = c(2, 0.5, -3),   # |diff| > 1 for s1 and s3 only
      cond_vs_ctrl_p.adj = c(0.01, 0.20, 0.001),
      cond_vs_ctrl_p.val = c(0.005, 0.10, 0.0005),
      cond_vs_ctrl_CI.L  = c(1, -1, -4),
      cond_vs_ctrl_CI.R  = c(3, 2, -2)
    )
  )
}

test_that("add_significance flags rows passing both p and diff thresholds", {
  se  <- make_sig_se()
  se2 <- suppressWarnings(add_significance(se, p_thr = 0.05, diff_thr = 1))
  rd  <- SummarizedExperiment::rowData(se2)

  # Per-contrast and overall significance columns are added
  expect_true("cond_vs_ctrl_significant" %in% colnames(rd))
  expect_true("significant" %in% colnames(rd))

  # s1: p.adj 0.01 < 0.05 & diff 2 > 1            -> TRUE
  # s2: p.adj 0.20 not < 0.05                      -> FALSE
  # s3: p.adj 0.001 < 0.05 & diff -3 < -1          -> TRUE
  expect_equal(rd$cond_vs_ctrl_significant, c(TRUE, FALSE, TRUE))
  expect_equal(rd$significant, c(TRUE, FALSE, TRUE))
})

test_that("add_significance respects a stricter diff threshold", {
  se  <- make_sig_se()
  # diff_thr = 2.5 now excludes s1 (|2| < 2.5); s3 (|-3|) still passes
  se2 <- suppressWarnings(add_significance(se, p_thr = 0.05, diff_thr = 2.5))
  rd  <- SummarizedExperiment::rowData(se2)
  expect_equal(rd$cond_vs_ctrl_significant, c(FALSE, FALSE, TRUE))
})

test_that("add_randna() works without fdrtool attached", {
  # HybridMTest reaches fdrtool::gcmlcm() through the search path only, so
  # add_randna() used to fail with "could not find function 'gcmlcm'" in any
  # session that had not run library(fdrtool).
  skip_if_not_installed("HybridMTest")
  skip_if("package:fdrtool" %in% search(), "fdrtool already attached")

  se <- make_test_se(n_features = 300, replicates = 4)
  a <- SummarizedExperiment::assay(se)
  set.seed(2)
  a[sample(length(a), 600)] <- NA
  SummarizedExperiment::assay(se) <- a

  # HybridMTest's own grenander step warns on small/degenerate p-value sets;
  # that is not what this test is about.
  out <- suppressWarnings(add_randna(se))
  expect_type(SummarizedExperiment::rowData(out)$randna, "logical")
  expect_length(SummarizedExperiment::rowData(out)$randna, nrow(se))
  # and it puts the search path back the way it found it
  expect_false("package:fdrtool" %in% search())
})
