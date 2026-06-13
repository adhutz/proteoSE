test_that("center_substring extracts the centered substring", {
  expect_equal(center_substring("hello", 3), "ell")
  expect_equal(center_substring("abcdef", 2), "cd")
})

test_that("center_substring returns the whole string when n exceeds its length", {
  expect_equal(center_substring("hi", 5), "hi")
  expect_equal(center_substring("abc", 3), "abc")
})

test_that("split_genes keeps only the first entry by default", {
  df <- data.frame(gene_names = c("A;B;C", "D", "E;F"),
                   stringsAsFactors = FALSE)
  out <- split_genes(df, colname = "gene_names", keep_all = FALSE)
  expect_equal(nrow(out), 3)
  expect_equal(out$gene_names, c("A", "D", "E"))
})

test_that("split_genes expands all entries to long format when keep_all = TRUE", {
  df <- data.frame(gene_names = c("A;B;C", "D"),
                   value = c(1, 2),
                   stringsAsFactors = FALSE)
  # cSplit() emits an upstream type.convert() warning that is not our concern
  out <- suppressWarnings(as.data.frame(split_genes(df, colname = "gene_names", keep_all = TRUE)))
  # A;B;C -> 3 rows, D -> 1 row
  expect_equal(nrow(out), 4)
  expect_setequal(unique(out$gene_names), c("A", "B", "C", "D"))
})

test_that("split_genes honours a custom separator", {
  df <- data.frame(protein_ids = c("P1|P2", "P3"), stringsAsFactors = FALSE)
  out <- split_genes(df, colname = "protein_ids", keep_all = FALSE, sep = "|")
  expect_equal(out$protein_ids, c("P1", "P3"))
})

test_that("compare_dataframes_full flags identical frames as fully matching", {
  df1 <- data.frame(id = 1:3, val = c("a", NA, "c"), stringsAsFactors = FALSE)
  df2 <- data.frame(id = 1:3, val = c("a", NA, "c"), stringsAsFactors = FALSE)
  invisible(capture.output(res <- compare_dataframes_full(df1, df2)))
  expect_true(all(res$all_match))
  expect_true(all(res$match_pct == 100))
  # NA-vs-NA is counted as a match
  expect_equal(res$na_both[res$column == "val"], 1)
})

test_that("compare_dataframes_full detects a differing column", {
  df1 <- data.frame(id = 1:3, val = c(1, 2, 3))
  df2 <- data.frame(id = 1:3, val = c(1, 2, 99))
  invisible(capture.output(res <- compare_dataframes_full(df1, df2)))
  expect_true(res$all_match[res$column == "id"])
  expect_false(res$all_match[res$column == "val"])
  expect_equal(res$matching_rows[res$column == "val"], 2)
})
