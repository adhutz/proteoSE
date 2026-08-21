# Generate inst/extdata/example_phospho_sites.txt
#
# A SYNTHETIC MaxQuant "Phospho (STY)Sites.txt" file for
# vignettes/Phosphoproteome_analysis.Rmd. Every value in it is simulated --
# intensities, occupancies, localisation probabilities, sequence windows and
# accessions are all generated here. No measured data goes into the repository.
#
# Only the columns the proteoSE phospho pipeline actually reads are written, so
# the file stays small (a real MaxQuant site table has ~196 columns).
#
# Gene symbols are sampled from org.Hs.eg.db so that GO over-representation in
# the vignette has something real to map against; that is public annotation, not
# experimental data.
#
# Run from the package root:
#   Rscript data-raw/make_phospho_example.R

.libPaths("C:/R/R-4.5.3/library")
suppressMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(20260821)

n_sites     <- 1200
conditions  <- c("CondA", "CondB", "CondC")   # deliberately generic
reps        <- 1:4
multiplicity <- 1:3
samples     <- as.vector(outer(paste0("_r", reps), conditions,
                               function(r, c) paste0(c, r)))

# ---------------------------------------------------------------------------
# 1. Identities: real symbols, invented accessions and positions
# ---------------------------------------------------------------------------
symbols <- sort(unique(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")))
symbols <- grep("^[A-Z][A-Z0-9]{1,7}$", symbols, value = TRUE)
gene <- sample(symbols, n_sites, replace = TRUE)

# Accession-shaped, but deliberately not real UniProt accessions.
acc <- sprintf("SYN%05d", seq_len(n_sites))

aa       <- sample(c("S", "T", "Y"), n_sites, replace = TRUE, prob = c(.85, .12, .03))
position <- sample(20:900, n_sites, replace = TRUE)

# 31-mer window centred on the modified residue, as MaxQuant writes it.
aa_alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
seq_window <- vapply(seq_len(n_sites), function(i) {
  flank <- function(n) paste(sample(aa_alphabet, n, replace = TRUE), collapse = "")
  paste0(flank(15), aa[i], flank(15))
}, character(1))

# ---------------------------------------------------------------------------
# 2. Intensities
#
#    Log-normal baseline per site, a modest per-sample loading offset, and a
#    spiked-in set of sites that respond in CondB and/or CondC so the vignette's
#    differential testing has something to find. Missing values are more likely
#    at low intensity, the way real data behaves.
# ---------------------------------------------------------------------------
base_log <- rnorm(n_sites, mean = 24, sd = 2.2)
loading  <- setNames(rnorm(length(samples), 0, 0.15), samples)

# CondA is the reference; CondB and CondC each move a disjoint set of sites.
responder_b <- sample(n_sites, 90)
responder_c <- sample(setdiff(seq_len(n_sites), responder_b), 60)
effect_b <- numeric(n_sites)
effect_c <- numeric(n_sites)
effect_b[responder_b] <- sample(c(-1, 1), 90, TRUE) * runif(90, 1.5, 3.5)
effect_c[responder_c] <- sample(c(-1, 1), 60, TRUE) * runif(60, 1.5, 3.5)

intensity <- matrix(NA_real_, n_sites, length(samples), dimnames = list(NULL, samples))
for (s in samples) {
  cond <- sub("_r[0-9]+$", "", s)
  eff <- switch(cond, CondB = effect_b, CondC = effect_c, numeric(n_sites))
  lg <- base_log + eff + loading[[s]] + rnorm(n_sites, 0, 0.55)
  # low-abundance sites drop out more often
  p_missing <- stats::plogis(-(lg - 22) * 1.1) * 0.55
  lg[stats::runif(n_sites) < p_missing] <- NA
  intensity[, s] <- round(2^lg)
}

# Per-multiplicity split: multiplicity 1 carries most of the signal.
mult_share <- c(0.72, 0.21, 0.07)
mult_intensity <- lapply(multiplicity, function(m) {
  x <- round(intensity * mult_share[m] * matrix(runif(length(intensity), 0.8, 1.2),
                                                nrow(intensity)))
  # a site is not necessarily observed at every multiplicity
  x[matrix(runif(length(x)) < c(0.02, 0.35, 0.6)[m], nrow(x))] <- NA
  x
})

# ---------------------------------------------------------------------------
# 3. Occupancy: a bounded fraction, loosely tracking intensity
# ---------------------------------------------------------------------------
occupancy <- round(
  stats::plogis((log2(pmax(intensity, 1)) - 24) / 3 + rnorm(length(intensity), 0, 0.4)) * 100,
  2
)
occupancy[is.na(intensity)] <- NA

# ---------------------------------------------------------------------------
# 4. Assemble, in MaxQuant's column naming
# ---------------------------------------------------------------------------
blank <- function(x) ifelse(is.na(x), "", as.character(x))

out <- data.frame(
  Proteins                  = acc,
  `Positions within proteins` = position,
  `Leading proteins`        = acc,
  Protein                   = acc,
  `Protein names`           = paste0("Synthetic protein ", seq_len(n_sites)),
  `Gene names`              = gene,
  `Localization prob`       = round(runif(n_sites, 0.5, 1), 3),
  `Amino acid`              = aa,
  `Sequence window`         = seq_window,
  Position                  = position,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

for (s in samples) out[[paste0("Intensity ", s)]]  <- blank(intensity[, s])
for (m in multiplicity) {
  for (s in samples) {
    out[[paste0("Intensity ", s, "___", m)]] <- blank(mult_intensity[[m]][, s])
  }
}
for (s in samples) out[[paste0("Occupancy ", s)]] <- blank(occupancy[, s])

# Kept rows must be empty in both; a handful of decoys/contaminants exercise the
# filter the importer applies.
out$Reverse <- ""
out$`Potential contaminant` <- ""
out$Reverse[sample(n_sites, 15)] <- "+"
out$`Potential contaminant`[sample(n_sites, 12)] <- "+"

path <- file.path("inst", "extdata", "example_phospho_sites.txt")
utils::write.table(out, path, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "")
cat("wrote", path, "-", nrow(out), "sites,", ncol(out), "columns,",
    round(file.size(path) / 1024), "KB\n")
