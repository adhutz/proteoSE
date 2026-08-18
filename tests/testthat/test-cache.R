test_that(".onLoad() memoises every function named in .cached_fns", {
  skip_if_not_installed("memoise")
  skip_if_not_installed("cachem")
  skip_if_not(isTRUE(getOption("proteoSE.cache", TRUE)))

  for (fn in .cached_fns) {
    expect_true(memoise::is.memoised(get(fn, envir = asNamespace("proteoSE"))),
                info = fn)
  }
})

test_that("cache_info() and cache_clear() report and empty the cache dir", {
  dir <- file.path(tempdir(), "proteoSE-cache-test")
  dir.create(dir, showWarnings = FALSE)
  old <- options(proteoSE.cache_dir = dir)
  on.exit({ options(old); unlink(dir, recursive = TRUE) }, add = TRUE)

  writeBin(raw(100), file.path(dir, "a.rds"))
  writeBin(raw(200), file.path(dir, "b.rds"))

  info <- suppressMessages(proteoSE_cache_info())
  expect_equal(info$dir, dir)
  expect_equal(info$n, 2)
  expect_equal(info$size, 300)

  freed <- suppressMessages(proteoSE_cache_clear())
  expect_equal(freed, 300)
  expect_equal(suppressMessages(proteoSE_cache_info())$n, 0)
})
