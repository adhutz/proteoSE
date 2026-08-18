# Guards for the three internal helpers that replaced dependencies
# (DescTools::Closest, gdata::cbindX, GOSemSim:::getOffsprings/getAncestors).
# Where the original package is still installed, the check is equivalence.

test_that(".closest matches min(DescTools::Closest())", {
  x <- c(1, 4, 4.5, 9)
  expect_equal(.closest(x, 4.2), 4)
  expect_equal(.closest(x, 0), 1)
  expect_equal(.closest(c(6, 4), 5), 4)  # tie: smallest of the two, as min(Closest())

  skip_if_not_installed("DescTools")
  for (a in c(-3, 0.5, 4.25, 100)) {
    expect_equal(.closest(x, a), min(DescTools::Closest(x = x, a = a)))
  }
})

test_that(".cbind_pad matches gdata::cbindX", {
  dfs <- list(data.frame(a = 1:3), data.frame(b = 1:5), data.frame(c = 1))
  out <- .cbind_pad(dfs)
  expect_equal(dim(out), c(5, 3))
  expect_equal(out$a, c(1:3, NA, NA))
  expect_equal(out$c, c(1, NA, NA, NA, NA))

  skip_if_not_installed("gdata")
  ref <- do.call(gdata::cbindX, dfs)
  expect_equal(unname(as.list(out)), unname(as.list(ref)))
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
