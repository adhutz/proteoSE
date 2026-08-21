# Generate inst/extdata/example_1_proteinGroups.txt and example_1_design.txt
#
# A SYNTHETIC MaxQuant "proteinGroups.txt" file for vignettes/proteoSE.Rmd, the
# examples and the import tests. Every value in it is simulated -- intensities,
# LFQ values, accessions, protein names, peptide counts and scores are all
# generated here. No measured data goes into the repository.
#
# Only the columns the proteoSE proteome pipeline actually reads are written, so
# the file stays a few hundred KB (a real MaxQuant proteinGroups table has ~300
# columns, most of them per-sample identification metadata nothing here touches).
#
# Gene symbols are sampled from org.Hs.eg.db so that GO enrichment in the
# vignette has something real to map against; that is public annotation, not
# experimental data.
#
# Run from the package root:
#   Rscript data-raw/make_proteome_example.R

.libPaths("C:/R/R-4.5.3/library")
suppressMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(20260821)

n_groups   <- 1957                                    # rows before filtering
n_flagged  <- 149                                     # -> 1808 usable rows
# Deliberately generic, and already in the snake_case that janitor produces on
# import, so the design file below can name the columns without guessing.
conditions <- c("cond_a", "cond_b", "cond_c", "cond_d", "cond_e", "cond_f")
reps       <- 1:5
samples    <- as.vector(outer(paste0("_rep_", reps), conditions,
                              function(r, c) paste0(c, r)))

# ---------------------------------------------------------------------------
# 1. Identities: real symbols, invented accessions
# ---------------------------------------------------------------------------
symbols <- sort(unique(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")))
symbols <- grep("^[A-Z][A-Z0-9]{1,7}$", symbols, value = TRUE)

# Four GO biological processes seed the gene pool. Their members become the
# responders in cond_b below, two terms up and two down, so the vignette's GSEA
# step recovers real terms instead of plotting an empty enrichment. Random
# filler makes up the rest of the pool.
go_up   <- c("GO:0006954", "GO:0006281")   # inflammatory response, DNA repair
go_down <- c("GO:0007059", "GO:0006914")   # chromosome segregation, autophagy

go_members <- function(ids) {
  tab <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = ids, keytype = "GOALL", columns = "SYMBOL"))
  intersect(unique(tab$SYMBOL[!is.na(tab$SYMBOL)]), symbols)
}
# The four terms together have more members than there are protein groups, so
# take a slice of each: enough for GSEA to see the term, not the whole ontology.
up_genes   <- sample(go_members(go_up), 350)
down_genes <- sample(setdiff(go_members(go_down), up_genes), 250)

seeded <- c(up_genes, down_genes)
filler <- sample(setdiff(symbols, seeded), n_groups - length(seeded))
gene   <- sample(c(seeded, filler))

# A few protein groups sharing a symbol, as isoform groups and groups split on
# shared peptides do -- gseGO() rejects duplicate names, so enrich_go_se() has
# to collapse them.
gene[sample(n_groups, 12)] <- sample(gene, 12)

# A handful of groups with no gene symbol, as MaxQuant produces for unannotated
# entries -- enrich_go_se() has to survive them (see R/enrich-go.R).
gene[sample(n_groups, 6)] <- ""

# Accession-shaped, but deliberately not real UniProt accessions. A protein
# group holds one or more of them, semicolon separated, the way MaxQuant writes
# it; split_genes() in se_read_in() takes the first.
acc <- sprintf("SYN%05d", seq_len(n_groups))
extra <- vapply(seq_len(n_groups), function(i) {
  if (runif(1) < 0.25) paste0(";SYN", sprintf("%05d", sample(n_groups, 1))) else ""
}, character(1))
protein_ids <- paste0(acc, extra)

# ---------------------------------------------------------------------------
# 2. Intensities
#
#    Log-normal baseline per protein, a per-sample loading offset, and a spiked
#    responder set per condition (cond_a is the reference) so the vignette's
#    differential testing and GO enrichment have something to find. Missing
#    values are more likely at low intensity, the way real data behaves.
# ---------------------------------------------------------------------------
base_log <- rnorm(n_groups, mean = 26, sd = 2.4)
loading  <- setNames(rnorm(length(samples), 0, 0.12), samples)

effects <- lapply(conditions[-1], function(cond) {
  e <- numeric(n_groups)
  hit <- sample(n_groups, 120)
  e[hit] <- sample(c(-1, 1), 120, TRUE) * runif(120, 1.2, 3.0)
  e
})
names(effects) <- conditions[-1]

# cond_b instead moves whole GO terms, coherently, so the enrichment has a
# signal to find.
e_b <- numeric(n_groups)
e_b[gene %in% up_genes]   <-  runif(sum(gene %in% up_genes), 0.8, 2.2)
e_b[gene %in% down_genes] <- -runif(sum(gene %in% down_genes), 0.8, 2.2)
effects[["cond_b"]] <- e_b

lfq <- matrix(NA_real_, n_groups, length(samples), dimnames = list(NULL, samples))
for (s in samples) {
  cond <- sub("_rep_[0-9]+$", "", s)
  eff  <- if (cond %in% names(effects)) effects[[cond]] else numeric(n_groups)
  lg   <- base_log + eff + loading[[s]] + rnorm(n_groups, 0, 0.35)
  # low-abundance groups drop out more often; a small flat rate on top of that
  # gives add_randna() both MAR and MNAR rows to classify
  p_missing <- stats::plogis(-(lg - 23.5) * 1.3) * 0.6 + 0.01
  lg[stats::runif(n_groups) < p_missing] <- NA
  lfq[, s] <- lg
}

# MaxQuant reports raw intensities alongside LFQ; they track each other loosely.
intensity <- lfq + matrix(rnorm(length(lfq), 0.4, 0.25), nrow(lfq))

as_int <- function(m) {
  out <- round(2^m)
  out[is.na(out)] <- 0          # MaxQuant writes 0, not NA, for a missing value
  out
}
lfq_int <- as_int(lfq)
int_int <- as_int(intensity)

# ---------------------------------------------------------------------------
# 3. Assemble, in MaxQuant's column naming
# ---------------------------------------------------------------------------
n_pep <- pmax(1, rpois(n_groups, 6))

out <- data.frame(
  `Protein IDs`          = protein_ids,
  `Majority protein IDs` = acc,
  `Protein names`        = paste0("Synthetic protein ", seq_len(n_groups)),
  `Gene names`           = gene,
  `Fasta headers`        = paste0(">sp|", acc, "|SYN", seq_len(n_groups),
                                  "_HUMAN Synthetic protein ", seq_len(n_groups),
                                  " OS=Homo sapiens OX=9606"),
  `Number of proteins`   = 1L + (extra != ""),
  Peptides               = n_pep,
  `Razor + unique peptides` = n_pep,
  `Unique peptides`      = pmax(1, n_pep - rbinom(n_groups, 2, 0.3)),
  `Sequence coverage [%]` = round(runif(n_groups, 2, 80), 1),
  `Mol. weight [kDa]`    = round(runif(n_groups, 8, 250), 3),
  `Sequence length`      = sample(70:2200, n_groups, replace = TRUE),
  `Q-value`              = 0,
  Score                  = round(runif(n_groups, 5, 320), 3),
  Intensity              = rowSums(int_int),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

for (s in samples) out[[paste0("Intensity ", s)]]     <- int_int[, s]
for (s in samples) out[[paste0("LFQ intensity ", s)]] <- lfq_int[, s]
for (s in samples) out[[paste0("MS/MS count ", s)]]   <- rpois(n_groups, 8)

# ---------------------------------------------------------------------------
# 4. The three flag columns se_read_in() filters on. Exactly n_flagged rows
#    carry a flag, so the importer returns 1808 protein groups.
# ---------------------------------------------------------------------------
out$`Only identified by site` <- ""
out$Reverse                   <- ""
out$`Potential contaminant`   <- ""

flagged <- sample(n_groups, n_flagged)
which_flag <- sample(1:3, n_flagged, replace = TRUE, prob = c(0.4, 0.35, 0.25))
out$`Only identified by site`[flagged[which_flag == 1]] <- "+"
out$Reverse[flagged[which_flag == 2]]                   <- "+"
out$`Potential contaminant`[flagged[which_flag == 3]]   <- "+"

out$id <- seq_len(n_groups) - 1L

path <- file.path("inst", "extdata", "example_1_proteinGroups.txt")
utils::write.table(out, path, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "")
cat("wrote", path, "-", nrow(out), "groups,", ncol(out), "columns,",
    round(file.size(path) / 1024), "KB\n")

# ---------------------------------------------------------------------------
# 5. The matching experimental design, for the vignette's "read design from
#    file" path. The sample names are already snake_case, so the labels match
#    what janitor produces from the column headers on import.
# ---------------------------------------------------------------------------
design <- data.frame(
  label     = paste0("lfq_intensity_", samples),
  sample    = samples,
  ID        = sub("_rep_", "_", samples),
  condition = sub("_rep_[0-9]+$", "", samples),
  replicate = as.integer(sub(".*_rep_", "", samples)),
  stringsAsFactors = FALSE
)

path_design <- file.path("inst", "extdata", "example_1_design.txt")
utils::write.csv(design, path_design, row.names = FALSE, quote = FALSE)
cat("wrote", path_design, "-", nrow(design), "samples in",
    length(unique(design$condition)), "conditions\n")
