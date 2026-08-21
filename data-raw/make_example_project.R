# Generate a fictional Spectronaut example project for proteoSE.
#
# Human U2OS DIA-MS, 3 conditions x 5 replicates = 15 runs, ~4000 protein
# groups drawn from org.Hs.eg.db (real symbols + UniProt accessions +
# descriptions). Built-in structure: proteins regulated only in TreatA, only in
# TreatB, in both, and (majority) unchanged; missingness is a mix of MAR
# (random) and MNAR (low-abundance dropout).
#
# Writes two tab-separated files to inst/extdata/ that round-trip through
# optimized_spectronaut_to_se(candidates = NULL). Run from the package root:
#   Rscript -e ".libPaths('C:/R/R-4.5.3/library'); source('data-raw/make_example_project.R')"

suppressMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
})

set.seed(2026)

n_proteins  <- 4000
conditions  <- c("Ctrl", "TreatA", "TreatB")
n_rep       <- 5

## ---- 1. Real human protein annotation ------------------------------------
ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  keytype = "SYMBOL",
  columns = c("UNIPROT", "GENENAME")
) %>%
  filter(!is.na(UNIPROT), !is.na(GENENAME)) %>%
  distinct(SYMBOL, .keep_all = TRUE)          # one accession per gene symbol

stopifnot(nrow(ann) >= n_proteins)
ann <- ann[sample(nrow(ann), n_proteins), ]

## ---- 2. Run / condition grid ---------------------------------------------
grid   <- expand.grid(rep = seq_len(n_rep), condition = conditions,
                      stringsAsFactors = FALSE)
runs   <- paste0(grid$condition, "_", grid$rep)          # e.g. "Ctrl_1"
n_runs <- length(runs)                                   # 15

## ---- 3. Base abundances + per-run offsets --------------------------------
base    <- rnorm(n_proteins, mean = 16, sd = 2)          # log2 scale
run_off <- rnorm(n_runs, mean = 0, sd = 0.15)            # loading differences

## ---- 4. Regulation classes -----------------------------------------------
# ~85% unchanged; remainder split across only-A, only-B, and both.
cls <- sample(
  c("none", "A_only", "B_only", "both"),
  size = n_proteins, replace = TRUE,
  prob = c(0.85, 0.06, 0.06, 0.03)
)

# per-protein log2 fold change vs Ctrl, drawn to mix strong and subtle effects
fc <- function(k) sample(c(-1, 1), k, TRUE) * runif(k, 0.6, 3)   # |log2FC| 0.6-3
shift_A <- numeric(n_proteins)
shift_B <- numeric(n_proteins)
shift_A[cls %in% c("A_only", "both")] <- fc(sum(cls %in% c("A_only", "both")))
shift_B[cls %in% c("B_only", "both")] <- fc(sum(cls %in% c("B_only", "both")))

## ---- 5. Assemble the intensity matrix ------------------------------------
mat <- matrix(NA_real_, nrow = n_proteins, ncol = n_runs,
              dimnames = list(NULL, runs))
for (j in seq_len(n_runs)) {
  cond <- grid$condition[j]
  eff  <- switch(cond, Ctrl = 0, TreatA = shift_A, TreatB = shift_B)
  mat[, j] <- base + eff + run_off[j] + rnorm(n_proteins, sd = 0.3)
}

## ---- 6. Missingness: MAR + MNAR ------------------------------------------
# MAR: ~4% of values dropped at random.
mar <- matrix(runif(length(mat)) < 0.04, nrow = n_proteins)
# MNAR: low-abundance values more likely missing (logistic in intensity).
p_mnar <- plogis(-(mat - quantile(mat, 0.15, na.rm = TRUE)) * 1.1)
mnar   <- matrix(runif(length(mat)) < p_mnar, nrow = n_proteins)
mat[mar | mnar] <- NA
# Guarantee no fully-missing row (keep at least one Ctrl value per protein).
allna <- rowSums(!is.na(mat)) == 0
if (any(allna)) mat[allna, 1] <- base[allna]

## ---- 7. Build the wide report -------------------------------------------
qcols <- paste0("[", seq_len(n_runs), "] ", runs, ".raw.PG.Log2Quantity")
report <- data.frame(
  PG.ProteinGroups       = ann$UNIPROT,
  PG.Genes               = ann$SYMBOL,
  PG.ProteinDescriptions = ann$GENENAME,
  PG.Organisms           = "Homo sapiens",
  check.names = FALSE
)
report <- cbind(report, setNames(as.data.frame(round(mat, 3)), qcols))

## ---- 8. Build the condition setup ----------------------------------------
# Run.Label carries the ".raw" suffix and Label the bare run name, matching the
# layout optimized_spectronaut_to_se() expects (it drops Label, renames from
# Run.Label). See data-raw round-trip validation.
conditions_tbl <- data.frame(
  Run.Label = paste0(runs, ".raw"),
  Label     = runs,
  Condition = grid$condition,
  Replicate = grid$rep,
  File.Name = paste0(runs, ".raw"),
  check.names = FALSE
)

## ---- 9. Write ------------------------------------------------------------
out_dir <- "inst/extdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
report_path <- file.path(out_dir, "example_project_report.tsv")
cond_path   <- file.path(out_dir, "example_project_conditions.tsv")
write.table(report,         report_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(conditions_tbl, cond_path,   sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote ", report_path, " (", nrow(report), " x ", ncol(report), ")")
message("Wrote ", cond_path,   " (", nrow(conditions_tbl), " x ", ncol(conditions_tbl), ")")
message("Missing values in report: ",
        sum(is.na(mat)), " / ", length(mat),
        " (", round(100 * mean(is.na(mat)), 1), "%)")
