# future.md — High-impact functionality to add to `proteoSE`

Forward-looking proposals for the *next* phase of `proteoSE`. Each item states **the
problem**, **the proposal**, **why it's high-impact**, and a **sketch of the API with an
example**. Ordered roughly by impact-to-effort.

These are intentionally designed to fit the existing `SummarizedExperiment`-centric
design, so they compose with the functions already in the package.

**Baseline as of 2026-08-18** (measured, not remembered). The refactor is done: the
package is `proteoSE` 0.1.0 (formerly `sev`, then `sev2`), 23 themed `R/*.R` files, ~90
exports, renamed functions carry `.Deprecated()` aliases in `R/deprecated.R`, and
`devtools::test()` is **60 passing / 0 failing**. `R CMD check` on the *tracked* tree with
the CI workflow's own arguments gives **0 ERRORS / 2 WARNINGS / 1 NOTE** — both warnings
are the un-built `inst/doc` vignettes; the note covers the `:::` calls into
`GOSemSim`/`iSEE` plus the two org annotation packages, which are declared in
`Imports` so they get installed but are referenced by name at runtime, never
imported from.
Three vignettes exist plus a worked example (`docs/worked_example.Rmd`) driven by bundled
Spectronaut example data (`inst/extdata/example_project_report.tsv` + `_conditions.tsv`).

**GitHub Actions had never once passed** — 13 runs, 13 failures, all on a stale
`importFrom(DOSE, parse_ratio)` that made the package uninstallable anywhere a current
DOSE is resolved. Fixed on the `v0.2` branch (`60c99fe`); see `plans/deps-and-install.md`
for the evidence trail and the technique for reading CI errors without `gh`. That failure
is the strongest possible argument for item 11 below.

None of items 1-10 are implemented yet.

Items marked **[plan]** have a dedicated implementation plan in `plans/` (local-only,
gitignored). Items without one are proposals at the paragraph level.

---

## 1. One-call QC report (`qc_report()`) **[plan: `plans/qc-report.md`]**

**Problem.** Every analysis starts with the same manual checks — missingness, per-sample
CVs, intensity distributions, PCA, sample-correlation heatmap, ID counts. Right now users
stitch these together by hand from scattered functions; `docs/worked_example.Rmd` walks
through exactly that sequence, one call at a time.

**Proposal.** A single function that takes a `SummarizedExperiment` and renders a
self-contained HTML QC report (parameterised R Markdown shipped in `inst/rmarkdown/`).

**Why high-impact.** It's the first thing every user does, it catches bad samples before
hours of downstream work, and it makes results shareable with non-R collaborators.
Highest ROI of anything here. It is also now cheap to build and demo: the bundled example
project is a complete 4000 × 15 Spectronaut dataset with realistic MNAR missingness.

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

Implementation is assembly, not new plotting: reuse `plot_heatmap_clustered()`,
`theme_proteoSE()` and `add_randna()`, plus the DEP2 QC plots the worked example already
exercises (`plot_numbers`, `plot_frequency`, `plot_missval`, `plot_detect`, `plot_pca`,
`plot_cor`) inside one `.Rmd` template.

---

## 2. Cache layer for all external API calls **[plan: `plans/api-cache.md`]**

**Problem.** Six functions hit the network on every call:

| Function | File | Service |
|---|---|---|
| `get_network` | `R/network.R` | STRING (via `rbioapi`) |
| `genes_from_kegg` | `R/enrich-go.R` | KEGG (`KEGGREST`) |
| `find_kws`, `fetch_kw_accessions` | `R/uniprot.R` | UniProt (`uniprotREST`/`httr2`) |
| `pubmed_query` | `R/literature.R` | NCBI E-utilities (`rentrez`) |
| `scopus_query` | `R/literature.R` | Scopus |

That makes them slow, flaky, rate-limited, and **untestable** — which is why none of them
have a test today.

**Proposal.** A small on-disk memoisation layer (`memoise` + `cachem::cache_disk()` in
`tools::R_user_dir("proteoSE", "cache")`). Wrap each network function once.

**Why high-impact.** Speeds up real analyses (re-runs are instant), removes the #1 source
of CI flakiness, and is the prerequisite for testing these six functions at all.

```r
# user-facing: identical call, now cached for 30 days
net <- get_network(c("TP53","MDM2","CDKN1A"))   # first call: hits STRING
net <- get_network(c("TP53","MDM2","CDKN1A"))   # second call: from cache

proteoSE_cache_clear()            # housekeeping
proteoSE_cache_info()             # size, entries, age
```

Internally:
```r
get_network <- with_cache(get_network_impl, ttl = 30 * 24 * 3600)
```

---

## 3. SE structural validation (`validate_se()` / `.assert_se()`) **[plan: `plans/validate-se.md`]**

**Problem.** Most functions assume specific `rowData` columns (`gene_names`, `_diff`,
`p.val`, …) or assay names and fail deep inside with cryptic errors when those are
missing. This was flagged during the audit and never done.

**Proposal.** A public `validate_se(se, require = c(...))` plus an internal `.assert_se()`
used as a guard at the top of functions. Clear, early, actionable error messages.

**Why high-impact.** Turns "Error in `[.data.frame`: undefined columns" into "`se` is
missing rowData column `gene_names` (required by `plot_volcano`). Available: ...". Saves
users (and you, in support) enormous time.

```r
validate_se(se, require_rowdata = c("gene_names"),
                require_assays  = c("centered"),
                require_contrast = "treat_vs_ctrl")
#> Error: plot_volcano() needs a contrast 'treat_vs_ctrl'.
#>   rowData has _diff columns for: ctrl_vs_dmso, drugA_vs_dmso.
```

---

## 4. Unified import dispatcher (`read_proteomics()`) — **ON HOLD**

**Status: deferred by the maintainer (2026-06-15), and still deferred.** The importers
(`se_read_in`, `spectronaut_read_in`, `fragpipe_read_in`, `spectronaut_to_se`,
`optimized_spectronaut_to_se`, `phos_read_in_int/occ`) are to be optimised by hand first;
a dispatcher on top of code that is about to change is wasted work. Left here as a record
of the idea, not as a queued task.

**The original problem still stands.** Seven+ near-identical importers, users must know
which to call, cleaning logic is copy-pasted and drifts. The eventual shape:

```r
se <- read_proteomics("report.tsv",
                      source = "spectronaut",   # | "maxquant" | "fragpipe" | "diann"
                      type   = "protein",       # | "phospho"
                      design = "design.txt")
```

**What changed in the meantime, and matters when this is picked up:**

- `tests/testthat/test-import-proteomics.R` pins `se_read_in()`'s output on the real
  MaxQuant example data — the golden characterisation test that makes the rewrite safe.
- A Spectronaut fixture now ships too (`inst/extdata/example_project_report.tsv` +
  `_conditions.tsv`, exercised end-to-end by `docs/worked_example.Rmd` via
  `optimized_spectronaut_to_se()`), so a second characterisation test is now writable
  without new sample files. FragPipe still has no fixture.
- A WIP peptide-level QFeatures importer sits in `data-raw/spectronaut_to_qfeatures.R`
  (local-only). It calls three helpers that were never written, so it errors if called.
  Whoever unifies the importers should decide whether protein+peptide `QFeatures` output
  belongs in the same dispatcher.
- DIA-NN support is still absent and still increasingly common.

---

## 5. Standards-compliant export: SDRF / mzTab (`export_sdrf()`)

**Problem.** Researchers increasingly must deposit data to PRIDE/ProteomeXchange, which
expects SDRF sample metadata and mzTab results. `proteoSE` already holds all of this in
`colData`/`rowData`.

**Proposal.** `export_sdrf(se, file)` and `export_mztab(se, file)` that serialise the SE
metadata into the community formats.

**Why high-impact.** Turns a painful, error-prone, manual deposition step into one line —
a concrete reason for other labs to adopt the package.

```r
export_sdrf(se, "experiment.sdrf.tsv")     # sample/condition/replicate/factor columns
export_mztab(se, "results.mztab")          # quant + differential results
```

---

## 6. Reproducible pipeline object (`proteoSE_pipeline()`) **[plan: `plans/pipeline-object.md`]**

**Problem.** Analyses are bespoke scripts; re-running with the same parameters and
capturing provenance is manual. `save_session_report()` (in `R/utils.R`) is a start but
isn't wired into a workflow.

**Proposal.** A lightweight, chainable pipeline that records each step and its parameters,
runnable end-to-end and serialisable.

**Why high-impact.** Reproducibility is the single biggest pain point in proteomics labs;
a recorded, replayable pipeline plus an auto-generated methods paragraph is a major
differentiator.

```r
pipe <- proteoSE_pipeline(optimized_spectronaut_to_se(report, conditions)) |>
  step_filter(perc_na = 0.33) |>
  step_impute(method = "perseus") |>
  step_test(contrast = "TreatA_vs_Ctrl") |>
  step_enrich_go()

run(pipe)                          # executes, caches intermediates
write_methods(pipe, "methods.md")  # auto-generated methods text + package versions
saveRDS(pipe, "analysis.rds")      # full provenance, re-runnable
```

**Sequencing.** This is the largest item and should come *after* 2 and 3 — its steps
should wrap functions that already validate their inputs and cache their network calls.
See the plan for the argument that the first version is a recorded call list, not
`targets`.

---

## 7. Design-formula & contrast helpers (`make_design()` / `list_contrasts()`)

**Problem.** `test_diff_limma()` (formerly `advanced_test`) and `limma` require users to
hand-build `model.matrix` and `makeContrasts` strings — a frequent error source (typos in
contrast names produce silent nonsense).

**Proposal.** Helpers that build the design from `colData` and enumerate/validate all
pairwise (or specified) contrasts, with name checking.

**Why high-impact.** Removes the most common statistical-setup mistake and makes the
powerful `test_diff_limma()` accessible to non-statisticians. Pairs naturally with item 3
(same "tell the user what's actually available" error style).

```r
design    <- make_design(se, ~ 0 + condition + batch)
contrasts <- list_contrasts(se, "condition")          # all pairwise, validated
res       <- test_diff_limma(se, design, contrasts["treat_vs_ctrl"])
```

---

## 8. Reactome / MSigDB enrichment + cached gene-set backends

**Problem.** Enrichment is GO/KEGG-only (`enrich_go_se`, `genes_from_kegg`,
`phospho_ora`). Reactome and MSigDB (hallmark, etc.) are heavily requested and not
covered.

**Proposal.** `enrich_reactome(se)` and `enrich_msigdb(se, collection = "H")` mirroring
the existing `enrich_go_se()` interface (results into `metadata()`, same downstream
plotting), built on `ReactomePA`/`msigdbr` and sharing the item-2 cache.

**Why high-impact.** Broadens biological interpretation with near-zero new user-facing
concepts (same pattern as existing GO enrichment).

```r
se <- enrich_reactome(se, contrast = "all")
se <- enrich_msigdb(se, collection = "H")   # hallmark gene sets
# reuse existing plotting on the new metadata
```

Note the dependency cost: `ReactomePA` and `msigdbr` are two more packages on an already
heavy tree — declare them in **Suggests** with `requireNamespace()` guards (see item 11).

---

## 9. Power / sample-size planning (`estimate_power()`)

**Problem.** `block_randomize()` and `plan_experiment()` help *schedule* an experiment but
not *size* it. "How many replicates do I need?" is asked constantly.

**Proposal.** Given a pilot SE (or assumed CV and effect size), estimate detectable
fold-change vs. replicate count, with a plot.

**Why high-impact.** Moves the package upstream into experimental design — reviewers and
grant applications increasingly demand power justification.

```r
estimate_power(pilot_se,
               fold_change = 1.5,
               alpha = 0.05,
               n = 2:10)        # -> table + ggplot of power vs n
```

---

## 10. Quality-of-life: progress, logging, and friendlier errors **[plan: `plans/cli-progress.md`]**

**Problem.** Long-running functions (imports, enrichment, network) are silent; failures
are raw.

**Proposal.** Adopt `cli` for consistent progress bars, status messages, and structured
errors across the package; a single verbosity toggle.

**Why high-impact.** Cheap to add, dramatically improves perceived quality and
debuggability — the difference between "feels like a tool" and "feels like a script dump."

```r
proteoSE_options(verbose = TRUE)
se <- optimized_spectronaut_to_se(report, conditions)
#> i Reading Spectronaut report (1.2M rows)...
#> v Parsed 4,000 protein groups across 15 samples [3.1s]
#> ! 312 rows dropped (contaminants/reverse)
```

---

## 11. Dependency slimming + install reliability **[plan: `plans/deps-and-install.md`]**

**Problem.** `DESCRIPTION` declares **55 hard `Imports`** (against 18 `Suggests`),
spanning CRAN, Bioconductor and two GitHub-only packages (`DEP2`, `uniprotREST`, declared
in `Remotes:`). A new user's install has 55 ways to fail, and none of them can be
reproduced on the maintainer's machine, where every package is already present.

The empirical pass on 2026-08-18 (evidence trail in `plans/deps-and-install.md`) settled
most of what this item was guessing at:

- **The real failure mode was version drift, not unreachable repos.** A stale
  `importFrom(DOSE, parse_ratio)` made the package uninstallable wherever a current DOSE
  is resolved -- 13 CI runs, 13 failures, invisible locally because this machine holds
  DOSE 4.2.0. Fixed in `60c99fe`, which removed `DOSE` from `Imports` entirely. **The
  audit still worth doing is the rest of the `importFrom()` targets**, for the same class
  of bug.
- **Dependency *resolution* is healthy.** pak resolves CRAN + Bioconductor + both GitHub
  remotes on a bare runner every time, and `HybridMTest` is not archived.
- **`org.Hs.eg.db` / `org.Mm.eg.db` are now hard `Imports`** (`ac92684`). Every upstream
  package that uses an org db declares it in *Suggests* only, so a fresh install reliably
  ended up without one -- while clusterProfiler-based enrichment, which is most of what
  users come here for, needs it. Cost: they are large, and `R CMD check` gains a cosmetic
  "Namespaces in Imports field not imported from" NOTE because the OrgDb is passed as a
  string. Accepted deliberately.
- **`DO.db` is required by nothing** -- not one package in a fresh ~350-package tree
  declares it (DOSE moved to `HDO.db`), so the README's manual install step is stale.
- **Still true:** several `Imports` are used exactly once (`DescTools`, `IRanges`,
  `KEGGREST`, `SingleCellExperiment`, `ggraph`, `tidyselect`) -- a whole package pulled in
  for one call. That is what step 3 of the plan is for.

**Proposal.** (a) Verify the install empirically — CI already does a from-scratch
dependency resolve, and a clean-library install locally reproduces a new user exactly;
(b) declare the missing annotation packages; (c) move the optional/single-use ones to
`Suggests` behind `requireNamespace()` guards, following the pattern already used in
`R/utils.R` and `R/experiment-design.R`.

**Why high-impact.** This is the difference between "I installed it" and "I gave up".
Every other item on this list is worthless to a user who can't install the package.

---

## 12. pkgdown site **[plan: `plans/pkgdown-site.md`]**

**Problem.** The documentation exists — three vignettes, a full worked example with 12
real figures, an intro deck — but the only way to read it is to clone the repo or install
the package. GitHub renders none of it usefully.

**Proposal.** `_pkgdown.yml` (reference index grouped along the themed `R/` files,
articles from the vignettes) plus a GitHub Pages workflow.

**Why high-impact.** Cheap, and it turns work already done into something linkable from a
paper, a lab wiki, or an email. Note the conflict to resolve first: `docs/` currently
holds hand-made assets and pkgdown wants to own that directory.

---

## 13. Docs hygiene (small, do-now)

Not "future" work so much as loose ends left by the rename and the repo cleanup:

- `README.md` links to `docs/worked_example.md`, which was deleted from git and is now
  gitignored — a broken link on the GitHub landing page.
- `vignettes/worked_example.Rmd` is an untracked, pre-rename copy of
  `docs/worked_example.Rmd`: it still says `library(sev2)` / `sev2::plot_volcano()`, calls
  `get_df_wide()` (DEP2's, unqualified) and `impute_DEP()` (the export is `impute_DEP2`),
  and has no `%\VignetteIndexEntry{}` header, so it is not buildable as a vignette.
  Either fix and promote it, or delete it in favour of the `docs/` copy.

See `plans/docs-hygiene.md`.

---

### Suggested order

1. **11 (deps/install)** and **13 (docs hygiene)** first — they are small, and they gate
   whether anyone else can use the package at all.
2. **1, 2, 3** next: the foundational trio. `qc_report()` is the biggest visible win, the
   cache unblocks testing six untested functions, and `validate_se()` improves every
   error message in the package.
3. **12 (pkgdown)** once the vignettes are settled — it consumes them.
4. **8** builds directly on `enrich_go_se()`; **7** on `test_diff_limma()`. Both are
   contained, single-file additions.
5. **5, 9** are standalone features, schedule on demand.
6. **6 (pipeline)** last of the planned items — it should sit on top of 2 and 3.
7. **10 (cli)** is woven in continuously rather than scheduled.
8. **4 (importers)** stays on hold until the hand-optimisation happens.
