
# **sev2 - The SummarizedExperiment Viewer**

## Overview
The sev2 package is a utility package that relies completely on the bioconductor packages *DEP2* and *iSEE* / *iSEEu* for the analysis, adding convenience functions, imputation methods, GO-term enrichment via *clusterProfiler* and additional *iSEE* containers.

> **Note:** `sev2` is the refactored successor to the original `sev`. It installs under a
> distinct package name, so you can keep `sev` and `sev2` installed side by side
> (`library(sev)` vs `library(sev2)`).

## Installation
Open up R and run the following code:

```{r eval = FALSE, echo=T}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
# DO.db is a Bioconductor annotation package, not a CRAN one.
if (!require("DO.db", quietly = TRUE))
    BiocManager::install("DO.db")

# All other dependencies (CRAN, Bioconductor and the two GitHub-only packages
# DEP2 and uniprotREST) are resolved automatically from DESCRIPTION.
devtools::install_github("adhutz/sev2")
```

## Worked example
[`docs/worked_example.md`](docs/worked_example.md) walks one fictional project all the
way from a raw Spectronaut report to an interactive `iSEE` browser, with real figures
generated from the bundled example data. The companion slide deck
`docs/sev2_pipeline_intro.pptx` introduces the `SummarizedExperiment` structure that
every step reads from and writes to.

## View vignette
An examplary analysis is included as a vignette.
```{r eval = FALSE, echo=T}
vignette("sev")
```

## Run example
The examplary analysis is also included in the installation folder as an Rmarkdown document together with the proteinGroups.txt file. To open the corresponding folder directly, run the following code.
```{r eval = FALSE, echo=T}
utils::browseURL(system.file("extdata", package = "sev2"))
```

