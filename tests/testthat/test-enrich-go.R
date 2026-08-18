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
