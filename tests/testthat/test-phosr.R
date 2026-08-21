# prep_phosR()/prep_phosR_from_ppe() need real phospho data, so the species
# argument is all that is tested here -- it used to silently hand rat data to
# the mouse table and hard-code "human" downstream.
test_that("prep_phosR rejects an unknown species before touching the data", {
  expect_error(prep_phosR(NULL, species = "chicken"), "should be one of")
  expect_error(prep_phosR_from_ppe(NULL, species = "chicken"), "should be one of")
})

test_that("each species maps to its own PhosphoSitePlus table", {
  skip_if_not_installed("PhosR")
  e <- new.env()
  utils::data("PhosphoSitePlus", package = "PhosR", envir = e)
  for (sp in c("human", "mouse", "rat")) {
    expect_true(exists(paste0("PhosphoSite.", sp), envir = e))
  }
  expect_false(identical(e$PhosphoSite.rat, e$PhosphoSite.mouse))
})
