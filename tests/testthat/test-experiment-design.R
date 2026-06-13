test_that("block_randomize is reproducible for a fixed seed", {
  skip_if_not_installed("anticlust")

  df <- data.frame(
    batch         = 1,
    condition     = rep(c("a", "b"), 4),
    bio_replicate = 1:8,
    stringsAsFactors = FALSE
  )

  r1 <- block_randomize(df, seed = 1)
  r2 <- block_randomize(df, seed = 1)

  expect_false(is.null(r1))
  expect_equal(r1, r2)
  expect_equal(nrow(r1), nrow(df))
  expect_true(all(c("block", "new_sample_number", "condition") %in% colnames(r1)))
  expect_equal(r1$seed_used[1], "1")
})
