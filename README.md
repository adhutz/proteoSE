
# **proteoSE - Proteomics and Phosphoproteomics Analysis on SummarizedExperiment**

## Overview
`proteoSE` is a utility toolkit for (phospho)proteomics data analysis, organised around the
Bioconductor `SummarizedExperiment` class. One object carries your data from raw search-engine
output all the way to results: import (Spectronaut, MaxQuant, FragPipe), experiment design,
filtering, missing-value imputation, differential-abundance testing, GO-term and gene-set
enrichment, kinase-substrate and signalome analysis, STRING networks, and publication-ready
figures. It builds on *DEP2*, *limma*, *PhosR*, *clusterProfiler* and *iSEE* / *iSEEu*.

> **Note:** `proteoSE` was previously released as `sev` and then `sev2`. It installs under a
> distinct package name, so older installs keep working side by side
> (`library(sev)` vs `library(proteoSE)`). Renamed exports keep deprecated aliases for one
> release cycle — e.g. `theme_sev()` still works and points you at `theme_proteoSE()`.

## Installation
Open up R and run the following code:

```{r eval = FALSE, echo=T}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
# All dependencies -- CRAN, Bioconductor, the org.Hs.eg.db / org.Mm.eg.db
# annotation packages, and the two GitHub-only packages DEP2 and uniprotREST --
# are resolved automatically from DESCRIPTION.
devtools::install_github("adhutz/proteoSE")
```

A few optional features sit behind `Suggests` and are only needed if you use them:
`iSEEu` (volcano and feature-set panels), `ggraph` (STRING network plots), `cmapR`
(GCT / ssGSEA export) and `waldo` (detailed data-frame diffs). The functions that
need them say so and name the install command. To pull everything in up front, add
`dependencies = TRUE` to the call above.

## Worked example
[`docs/worked_example.Rmd`](docs/worked_example.Rmd) walks one fictional project all the
way from a raw Spectronaut report to an interactive `iSEE` browser, with real figures
generated from the bundled example data. The companion slide deck
`docs/proteoSE_pipeline_intro.pptx` introduces the `SummarizedExperiment` structure that
every step reads from and writes to.

## View vignette
An examplary analysis is included as a vignette.
```{r eval = FALSE, echo=T}
vignette("proteoSE")
```

## Run example
The examplary analysis is also included in the installation folder as an Rmarkdown document together with the proteinGroups.txt file. To open the corresponding folder directly, run the following code.
```{r eval = FALSE, echo=T}
utils::browseURL(system.file("extdata", package = "proteoSE"))
```

