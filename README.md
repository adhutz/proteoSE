
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
    
if (!require("DO.db", quietly = TRUE))
    install.packages("DO.db")

#All other dependencies are installed automatically with sev2. 
devtools::install_github("adhutz/sev2")
```

## View vignette
An examplary analysis is included as a vignette. 
```{r eval = FALSE, echo=T}
vignette("sev")
```

## Run example
The examplary analysis is also included in the installation folder as an Rmarkdown document together with the proteinGroups.txt file. To open the corresponding folder directly, run the following code.
```{r eval = FALSE, echo=T}
utils::browseURL(system.file("extdata", package = "sev"))
```

