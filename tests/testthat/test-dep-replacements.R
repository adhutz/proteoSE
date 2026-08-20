# Guards for the three internal helpers that replaced dependencies
# (DescTools::Closest, gdata::cbindX, GOSemSim:::getOffsprings/getAncestors).
#
# These started out as equivalence tests against the originals. The expected
# values below are the ones those packages produced, hard-coded once the
# migration was proven, so the test suite no longer depends on packages the
# package itself dropped.

test_that(".closest returns the nearest value, smallest on a tie", {
  x <- c(1, 4, 4.5, 9)
  expect_equal(.closest(x, 4.2), 4)
  expect_equal(.closest(x, 0), 1)
  expect_equal(.closest(c(6, 4), 5), 4)  # tie: smallest of the two

  # min(DescTools::Closest(x, a)) for each of these
  expect_equal(.closest(x, -3), 1)
  expect_equal(.closest(x, 0.5), 1)
  expect_equal(.closest(x, 4.25), 4)     # tie between 4 and 4.5
  expect_equal(.closest(x, 100), 9)
})

test_that(".cbind_pad pads short data frames with NA", {
  dfs <- list(data.frame(a = 1:3), data.frame(b = 1:5), data.frame(c = 1))
  out <- .cbind_pad(dfs)
  expect_equal(dim(out), c(5, 3))
  expect_equal(out$a, c(1:3, NA, NA))
  expect_equal(out$b, 1:5)
  expect_equal(out$c, c(1, NA, NA, NA, NA))
})

test_that(".go_relatives returns the GO.db mappings GOSemSim used internally", {
  skip_if_not_installed("GO.db")
  off <- .go_relatives("BP", "OFFSPRING")
  anc <- .go_relatives("BP", "ANCESTOR")
  expect_type(off, "list")
  expect_true(length(off) > 1000)
  # the BP root is its own ancestor set boundary and has many offspring
  expect_true("GO:0008150" %in% names(off))
  expect_true("GO:0008150" %in% names(anc))
  expect_error(.go_relatives("XX", "OFFSPRING"), "Unknown ontology")
})
