# future.md — High-impact functionality to add to `sev`

Forward-looking proposals for the *next* phase of `sev`, after the cleanup in `AUDIT.md` lands. Each item states **the problem**, **the proposal**, **why it's high-impact**, and a **sketch of the API with an example**. Ordered roughly by impact-to-effort.

These are intentionally designed to fit the existing `SummarizedExperiment`-centric design, so they compose with the functions already in the package.

---

## 1. One-call QC report (`qc_report()`)

**Problem.** Every analysis starts with the same manual checks — missingness, per-sample CVs, intensity distributions, PCA, sample-correlation heatmap, ID counts. Right now users stitch these together by hand from scattered functions.

**Proposal.** A single function that takes a `SummarizedExperiment` and renders a self-contained HTML QC report (parameterised R Markdown shipped in `inst/rmarkdown/`).

**Why high-impact.** It's the first thing every user does, it catches bad samples before hours of downstream work, and it makes results shareable with non-R collaborators. Highest ROI of anything here.

```r
qc_report(se,
          group = "condition",
          output = "qc_report.html")

# Returns the path; opens in browser. Sections:
#  - ID counts per sample (barplot)
#  - missing-value heatmap + % missing per sample/condition
#  - intensity boxplots / density (pre/post normalisation)
#  - coefficient-of-variation per condition
#  - PCA coloured by condition + batch
#  - sample-sample correlation heatmap (flags swaps/outliers)
```

Implementation reuses existing pieces (`clustered_heatmap`, `my_theme`, `add_randna`) inside an `.Rmd` template.

---

## 2. Cache layer for all external API calls (`with_cache()`)

**Problem.** `get_network` (STRING), `genes_from_kegg` (KEGG), `find_kws`/`fetch_kw_accessions` (UniProt), `pubmed_query`, `scopus_query` all hit the network on every call. That makes them slow, flaky, rate-limited, and **untestable**.

**Proposal.** A small on-disk memoisation layer (built on `memoise` + a `tools::R_user_dir("sev","cache")` directory). Wrap each network function once.

**Why high-impact.** Speeds up real analyses (re-runs are instant), removes the #1 source of CI flakiness, and is a prerequisite for testing the API functions properly (§8 of the audit).

```r
# user-facing: identical call, now cached for 30 days
net <- get_network(c("TP53","MDM2","CDKN1A"))   # first call: hits STRING
net <- get_network(c("TP53","MDM2","CDKN1A"))   # second call: from cache

sev_cache_clear()                 # housekeeping
sev_cache_info()                  # size, entries, age
```

Internally:
```r
get_network <- with_cache(get_network_impl, ttl = 30 * 24 * 3600)
```

---

## 3. SE structural validation (`validate_se()` / `.assert_se()`)

**Problem.** Most functions assume specific `rowData` columns (`gene_names`, `_diff`, `p.val`, …) or assay names and fail deep inside with cryptic errors when those are missing.

**Proposal.** A public `validate_se(se, require = c(...))` plus an internal `.assert_se()` used as a guard at the top of functions. Clear, early, actionable error messages.

**Why high-impact.** Turns "Error in `[.data.frame`: undefined columns" into "`se` is missing rowData column `gene_names` (required by `plot_volcano`). Available: ...". Saves users (and you, in support) enormous time.

```r
validate_se(se, require_rowdata = c("gene_names"),
                require_assays  = c("centered"),
                require_contrast = "treat_vs_ctrl")
#> Error: plot_volcano() needs a contrast 'treat_vs_ctrl'.
#>   rowData has _diff columns for: ctrl_vs_dmso, drugA_vs_dmso.
```

---

## 4. Unified import dispatcher (`read_proteomics()`)

**Problem.** Seven+ near-identical importers (`se_read_in`, `spectronaut_read_in`, `fragpipe_read_in`, `spectronaut_to_se`, `optimized_spectronaut_to_se`, `phos_read_in_*`). Users must know which to call; the cleaning logic is copy-pasted and drifts.

**Proposal.** One entry point that detects or is told the source, delegating to a shared internal `.table_to_se()` core (this is audit §4c, surfaced here as the user-facing payoff).

**Why high-impact.** One thing to learn, one place to fix bugs, consistent output regardless of search engine.

```r
se <- read_proteomics("report.tsv",
                      source = "spectronaut",   # | "maxquant" | "fragpipe" | "diann"
                      type   = "protein",       # | "phospho"
                      design = "design.txt")

# auto-detect from file contents when source is omitted:
se <- read_proteomics("proteinGroups.txt")
```

Adding DIA-NN support (currently absent, increasingly common) is a natural extension of the same dispatcher.

---

## 5. Standards-compliant export: SDRF / mzTab (`export_sdrf()`)

**Problem.** Researchers increasingly must deposit data to PRIDE/ProteomeXchange, which expects SDRF sample metadata and mzTab results. `sev` already holds all of this in `colData`/`rowData`.

**Proposal.** `export_sdrf(se, file)` and `export_mztab(se, file)` that serialise the SE metadata into the community formats.

**Why high-impact.** Turns a painful, error-prone, manual deposition step into one line — a concrete reason for other labs to adopt the package.

```r
export_sdrf(se, "experiment.sdrf.tsv")     # sample/condition/replicate/factor columns
export_mztab(se, "results.mztab")          # quant + differential results
```

---

## 6. Reproducible pipeline object (`sev_pipeline()`)

**Problem.** Analyses are bespoke scripts; re-running with the same parameters and capturing provenance is manual. `save_session_report()` (in `mlasse`) is a start but isn't wired into a workflow.

**Proposal.** A lightweight, chainable pipeline that records each step and its parameters, runnable end-to-end and serialisable.

**Why high-impact.** Reproducibility is the single biggest pain point in proteomics labs; a recorded, replayable pipeline plus an auto-generated methods paragraph is a major differentiator.

```r
pipe <- sev_pipeline(read_proteomics("report.tsv", source = "spectronaut")) |>
  step_filter(perc_na = 0.33) |>
  step_impute(method = "MinProb") |>
  step_test(contrast = "treat_vs_ctrl") |>
  step_enrich_go()

run(pipe)                          # executes, caches intermediates
write_methods(pipe, "methods.md")  # auto-generated methods text + package versions
saveRDS(pipe, "analysis.rds")      # full provenance, re-runnable
```

Could be backed by `targets` for the heavy lifting.

---

## 7. Design-formula & contrast helpers (`make_design()` / `list_contrasts()`)

**Problem.** `advanced_test`/`limma` require users to hand-build `model.matrix` and `makeContrasts` strings — a frequent error source (typos in contrast names produce silent nonsense).

**Proposal.** Helpers that build the design from `colData` and enumerate/validate all pairwise (or specified) contrasts, with name checking.

**Why high-impact.** Removes the most common statistical-setup mistake and makes the powerful `advanced_test` accessible to non-statisticians.

```r
design   <- make_design(se, ~ 0 + condition + batch)
contrasts <- list_contrasts(se, "condition")          # all pairwise, validated
res      <- test_diff_limma(se, design, contrasts["treat_vs_ctrl"])
```

---

## 8. Reactome / MSigDB enrichment + cached gene-set backends

**Problem.** Enrichment is GO/KEGG-only (`se_GOE`, `genes_from_kegg`, `phospho_ora`). Reactome and MSigDB (hallmark, etc.) are heavily requested and not covered.

**Proposal.** `enrich_reactome(se)` and `enrich_msigdb(se, collection = "H")` mirroring the existing `se_GOE` interface (results into metadata, same downstream plotting), built on `ReactomePA`/`msigdbr` and sharing the §2 cache.

**Why high-impact.** Broadens biological interpretation with near-zero new user-facing concepts (same pattern as existing GO enrichment).

```r
se <- enrich_reactome(se, contrast = "all")
se <- enrich_msigdb(se, collection = "H")   # hallmark gene sets
# reuse existing plotting on the new metadata
```

---

## 9. Power / sample-size planning (`estimate_power()`)

**Problem.** `block_randomize` and `plan_experiment` help *schedule* an experiment but not *size* it. "How many replicates do I need?" is asked constantly.

**Proposal.** Given a pilot SE (or assumed CV and effect size), estimate detectable fold-change vs. replicate count, with a plot.

**Why high-impact.** Moves the package upstream into experimental design — reviewers and grant applications increasingly demand power justification.

```r
estimate_power(pilot_se,
               fold_change = 1.5,
               alpha = 0.05,
               n = 2:10)        # -> table + ggplot of power vs n
```

---

## 10. Quality-of-life: progress, logging, and friendlier errors

**Problem.** Long-running functions (imports, enrichment, network) are silent; failures are raw.

**Proposal.** Adopt `cli` for consistent progress bars, status messages, and structured errors across the package; a single `sev_verbose()` toggle.

**Why high-impact.** Cheap to add, dramatically improves perceived quality and debuggability — the difference between "feels like a tool" and "feels like a script dump."

```r
sev_options(verbose = TRUE)
se <- read_proteomics("big_report.tsv", source = "spectronaut")
#> ℹ Reading Spectronaut report (1.2M rows)...
#> ✔ Parsed 8,432 protein groups across 24 samples [3.1s]
#> ⚠ 312 rows dropped (contaminants/reverse)
```

---

### Suggested order

1, 2, and 3 first — they are foundational, low-risk, and make everything else (and the tests in `AUDIT.md`) easier. 4 and 8 build directly on existing code. 5, 6, 7, 9 are larger features to schedule once the core is clean and tested. 10 can be woven in continuously.
