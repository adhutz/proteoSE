test_that("shorten_terms preserves NA and vector length and type", {
  x <- c("antigen processing", "positive regulation of process", NA)
  out <- shorten_terms(x)
  expect_type(out, "character")
  expect_length(out, length(x))
  expect_true(is.na(out[3]))
})

test_that("shorten_terms validates its arguments", {
  expect_error(shorten_terms(123), "character")
  expect_error(
    shorten_terms("antigen processing", extra_words = c("antigen", "x")),
    "named"
  )
})

test_that("shorten_terms applies user replacements", {
  out <- shorten_terms(
    "antigen processing",
    extra_words = c("antigen" = "antig."),
    extra_first = TRUE
  )
  expect_false(grepl("antigen", out, fixed = TRUE))
})

test_that(".parse_ratio turns clusterProfiler k/n strings into numbers", {
  # Replaces DOSE::parse_ratio(), which current DOSE no longer exports.
  expect_equal(.parse_ratio("12/345"), 12 / 345)
  expect_equal(.parse_ratio(c("1/2", "3/4")), c(0.5, 0.75))
  expect_equal(.parse_ratio(character(0)), numeric(0))
  expect_true(is.na(.parse_ratio(NA_character_)))
})

test_that("enrich_go_se() says so when nothing has a gene symbol", {
  # gseGO() aborts with "NAs in names(stats) are not allowed" if any feature
  # lacks a symbol -- real MaxQuant output always has a few. enrich_go_se()
  # drops those; when that leaves nothing, it must say why.
  se <- make_test_se(n_features = 6)
  SummarizedExperiment::rowData(se)$gene_names <- NA_character_
  SummarizedExperiment::rowData(se)$treat_vs_ctrl_diff <- seq_len(nrow(se)) / 10

  expect_error(
    enrich_go_se(se, col_names = "treat_vs_ctrl_diff"),
    "No usable gene symbols"
  )
})

# Both defaults are bare package objects. They resolve only through a real
# @importFrom: listing them in globalVariables() silences the check NOTE
# without making the lookup work, which is how `ontology_db = GO.db` stayed
# broken until a vignette render hit it.
test_that("get_go_terms() default annotation databases resolve", {
  ns <- environment(get_go_terms)
  expect_s4_class(eval(formals(get_go_terms)$organism_db, ns), "OrgDb")
  expect_s4_class(eval(formals(get_go_terms)$ontology_db, ns), "GODb")
})
