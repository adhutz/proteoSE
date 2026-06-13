test_that("impute_seeded is reproducible for identical input", {
  noisy <- function(x) x + stats::rnorm(length(x))

  r1 <- impute_seeded(1:10, noisy)
  r2 <- impute_seeded(1:10, noisy)
  expect_identical(r1, r2)
})

test_that("impute_seeded derives the seed from the input, not global RNG state", {
  noisy <- function(x) x + stats::rnorm(length(x))

  set.seed(1)
  a <- impute_seeded(1:10, noisy)
  set.seed(999)            # different global state ...
  b <- impute_seeded(1:10, noisy)
  expect_identical(a, b)   # ... still the same result for the same input

  # Different input generally yields a different draw
  c <- impute_seeded(11:20, noisy)
  expect_false(isTRUE(all.equal(unname(a), unname(c))))
})

test_that("impute_seeded passes extra arguments through to impute_fun", {
  add_const <- function(x, k = 0) x + k
  expect_equal(impute_seeded(rep(0, 3), add_const, k = 5), rep(5, 3))
})
