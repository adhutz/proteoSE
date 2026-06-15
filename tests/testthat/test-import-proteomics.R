# Characterization tests for the protein-level importers.
#
# These pin the CURRENT output of each importer on the real example data so the
# Phase 5 unification (folding se_read_in / spectronaut_read_in / fragpipe_read_in
# into read_proteomics()) can be proven behaviour-preserving. They are golden
# tests: if the refactor changes a number here, that is a regression to explain,
# not a test to "fix".
#
# Only MaxQuant example data ships in inst/extdata, so only se_read_in is pinned
# here for now. Spectronaut/FragPipe fixtures + their characterization tests are
# added once sample files exist (then this file gains read_proteomics(source=)
# equivalents).

test_that("se_read_in (MaxQuant) reproduces the reference SummarizedExperiment", {
  skip_on_cran()
  f <- system.file("extdata", "example_1_proteinGroups.txt", package = "sev")
  skip_if(f == "", "example_1_proteinGroups.txt not installed")

  se <- suppressWarnings(suppressMessages(se_read_in(f)))

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(1808L, 30L))
  expect_equal(SummarizedExperiment::assayNames(se), "lfq_raw")

  cd <- SummarizedExperiment::colData(se)
  expect_true(all(c("label", "sample", "condition", "replicate") %in% colnames(cd)))
  expect_equal(nrow(cd), 30L)
  expect_setequal(
    unique(cd$condition),
    c("griesser_heat", "griesser_no_heat", "no_glom_heat",
      "scarred_glom_heat", "sds_heat", "sds_no_heat")
  )

  expect_equal(head(rownames(se), 5),
               c("IGLV3-9", "IGKV2D-28", "IGKV3D-11", "IGHV3-74", "P0DPI2"))

  m <- SummarizedExperiment::assay(se, "lfq_raw")
  expect_equal(dim(m), c(1808L, 30L))
  expect_equal(sum(is.na(m)), 33642L)
  expect_equal(sum(m[, 1], na.rm = TRUE), 18637.33, tolerance = 1e-2)
})
