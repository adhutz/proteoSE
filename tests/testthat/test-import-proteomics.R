# Characterization tests for the protein-level importers.
#
# These pin the CURRENT output of each importer on the bundled synthetic example
# data (data-raw/make_proteome_example.R) so the
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
  f <- system.file("extdata", "example_1_proteinGroups.txt", package = "proteoSE")
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
    c("cond_a", "cond_b", "cond_c", "cond_d", "cond_e", "cond_f")
  )

  expect_equal(head(rownames(se), 5),
               c("MIR1263", "SUN5", "MTND5P13", "DEPTOR", "MIR4731"))

  m <- SummarizedExperiment::assay(se, "lfq_raw")
  expect_equal(dim(m), c(1808L, 30L))
  expect_equal(sum(is.na(m)), 6603L)
  expect_equal(sum(m[, 1], na.rm = TRUE), 42583.81, tolerance = 1e-2)
})
