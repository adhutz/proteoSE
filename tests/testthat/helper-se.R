# Test helpers for sev.
#
# `make_test_se()` builds a tiny SummarizedExperiment with the column/assay
# conventions sev expects (gene_names, condition, replicate), for use as a
# fixture in unit tests. Real tests are added in Phase 4 of the refactor
# (see AUDIT.md, section 8).

make_test_se <- function(n_features = 10, conditions = c("ctrl", "treat"),
                         replicates = 3, seed = 1) {
  set.seed(seed)
  samples <- expand.grid(condition = conditions,
                         replicate = seq_len(replicates),
                         stringsAsFactors = FALSE)
  n_samples <- nrow(samples)
  mat <- matrix(stats::rnorm(n_features * n_samples, mean = 20),
                nrow = n_features,
                dimnames = list(paste0("P", seq_len(n_features)), NULL))
  col_data <- S4Vectors::DataFrame(
    sample    = paste(samples$condition, samples$replicate, sep = "_rep_"),
    condition = samples$condition,
    replicate = samples$replicate
  )
  row_data <- S4Vectors::DataFrame(
    gene_names  = paste0("G", seq_len(n_features)),
    protein_ids = rownames(mat)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = mat),
    colData = col_data,
    rowData = row_data
  )
}
