test_that("the hard-coded iSEE multi-select flag still matches iSEE", {
  # R/isee-panels.R inlines "INTERNAL_multi_select" rather than calling
  # iSEE:::.flagMultiSelect, which R CMD check NOTEs on. If iSEE ever renames
  # the constant, the panel silently stops seeing row selections -- so assert
  # the two agree here instead.
  skip_if_not_installed("iSEE")
  expect_identical(
    utils::getFromNamespace(".flagMultiSelect", "iSEE"),
    "INTERNAL_multi_select"
  )
})
