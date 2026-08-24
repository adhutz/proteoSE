# The point of validate_se() is the error text, so the tests assert on it.

# make_test_se() has no results columns; add two contrasts' worth.
with_results <- function(se, contrasts = c("treat_vs_ctrl", "drugA_vs_dmso")) {
  rd <- SummarizedExperiment::rowData(se)
  for (cn in contrasts) {
    rd[[paste0(cn, "_diff")]]  <- seq_len(nrow(se)) / 10
    rd[[paste0(cn, "_p.adj")]] <- 0.01
  }
  SummarizedExperiment::rowData(se) <- rd
  se
}

test_that("validate_se passes silently and returns the object invisibly", {
  se <- make_test_se()
  expect_silent(validate_se(se,
                            require_rowdata = c("gene_names", "protein_ids"),
                            require_coldata = c("condition", "replicate"),
                            require_assays  = "intensity"))
  expect_invisible(validate_se(se))
  expect_identical(validate_se(se), se)
})

test_that("a missing rowData column is named alongside the ones present", {
  se <- make_test_se()
  expect_error(validate_se(se, require_rowdata = "sequence_window"),
               "sequence_window")
  expect_error(validate_se(se, require_rowdata = "sequence_window"),
               "rowData has: gene_names, protein_ids")
})

test_that("a missing colData column is named alongside the ones present", {
  se <- make_test_se()
  expect_error(validate_se(se, require_coldata = "batch"),
               "colData has: sample, condition, replicate")
})

test_that("a missing assay lists assayNames()", {
  se <- make_test_se()
  expect_error(validate_se(se, require_assays = "imputed"),
               "Available assays: intensity")
})

test_that("a missing contrast lists the contrasts that do exist", {
  se <- with_results(make_test_se())
  expect_silent(validate_se(se, require_contrast = "treat_vs_ctrl"))
  expect_error(validate_se(se, require_contrast = "nonexistent"),
               "rowData has _diff columns for: treat_vs_ctrl, drugA_vs_dmso")
})

test_that("an SE with no results at all still reports readably", {
  expect_error(validate_se(make_test_se(), require_contrast = "any"),
               "_diff columns for: none")
})

test_that("non-SE input is rejected by its actual class", {
  expect_error(validate_se(data.frame(a = 1)),
               "needs a SummarizedExperiment object; got data.frame")
})

test_that(".assert_se blames the calling function", {
  caller <- function(se) .assert_se(se, require_assays = "imputed")
  expect_error(caller(make_test_se()), "^caller\\(\\) needs assay")
})

test_that("plot_volcano rejects an unknown contrast readably", {
  se <- with_results(make_test_se())
  expect_error(plot_volcano(se, contrast = "nonexistent"),
               "plot_volcano\\(\\) needs a contrast nonexistent")
  # ...instead of the undefined-columns error from three frames deep
  expect_error(plot_volcano(se, contrast = "nonexistent"),
               "treat_vs_ctrl")
})

test_that("the guarded entry points reject a non-SE by name", {
  df <- data.frame(a = 1)
  expect_error(filter_perseus(df), "filter_perseus\\(\\) needs a SummarizedExperiment")
  expect_error(impute_perseus(df), "impute_perseus\\(\\) needs a SummarizedExperiment")
  expect_error(write_prot(df),     "write_prot\\(\\) needs a SummarizedExperiment")
})

# ---------------------------------------------------------------------------
# Chunk I: the seven entry points that used to fail with a raw dplyr error
# twenty frames deep, plus the two the required-argument probe reached.

test_that("the Chunk I guards name the missing column and what is present", {
  se <- with_results(make_test_se())

  expect_error(add_significance(se),
               "add_significance\\(\\) needs rowData column\\(s\\) site_id_mult")
  expect_error(make_mini(se),
               "make_mini\\(\\) needs rowData column\\(s\\) majority_protein_ids")
  expect_error(write_prot(se),
               "write_prot\\(\\) needs rowData column\\(s\\) protein_descriptions")
  expect_error(write_phos(se),
               "write_phos\\(\\) needs rowData column\\(s\\) protein")
  expect_error(test_diff_long(se),
               "test_diff_long\\(\\) needs rowData column\\(s\\) significant")
  expect_error(phospho_ora(se),
               "phospho_ora\\(\\) needs rowData column\\(s\\) significant")
  expect_error(get_df_wide_append_assay(se),
               "get_df_wide_append_assay\\(\\) needs rowData column\\(s\\) description")
  expect_error(prep_ksea(se, "treat_vs_ctrl"),
               "prep_ksea\\(\\) needs rowData column\\(s\\) proteins")
  expect_error(plot_antigen_missing(se, "treat", "ctrl"),
               "plot_antigen_missing\\(\\) needs colData column\\(s\\) label")

  # the whole point: the message lists what the object does have
  expect_error(make_mini(se), "rowData has: gene_names, protein_ids")
})

test_that("prep_ksea reports an unknown contrast rather than an empty selection", {
  se <- with_results(make_test_se())
  SummarizedExperiment::rowData(se)$proteins   <- "P1"
  SummarizedExperiment::rowData(se)$amino_acid <- "S"
  SummarizedExperiment::rowData(se)$position   <- 1
  expect_error(prep_ksea(se, "nonexistent"),
               "prep_ksea\\(\\) needs a contrast nonexistent")
})
