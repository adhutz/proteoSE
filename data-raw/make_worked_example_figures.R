# Regenerate the figures embedded in vignettes/worked_example.md.
#
# Runs the full example pipeline and writes PNGs to vignettes/articles/figures/.
# Run from the package root:
#   Rscript -e ".libPaths('C:/R/R-4.5.3/library'); source('data-raw/make_worked_example_figures.R')"
#
# Note: DEP2 also exports plot_volcano(), which masks proteoSE's when DEP2 is
# attached — hence the explicit proteoSE:: below.

suppressMessages({
  library(proteoSE)
  library(SummarizedExperiment)
  library(DEP2)
  library(ggplot2)
})

figdir <- "vignettes/articles/figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(name, expr, w = 7, h = 5, dpi = 110) {
  png(file.path(figdir, paste0(name, ".png")), width = w, height = h,
      units = "in", res = dpi)
  on.exit(dev.off())
  obj <- force(expr)
  if (inherits(obj, c("ggplot", "Heatmap", "HeatmapList"))) print(obj)
}

## ---- pipeline ------------------------------------------------------------
se <- optimized_spectronaut_to_se(
  report         = "inst/extdata/example_project_report.tsv",
  conditionSetup = "inst/extdata/example_project_conditions.tsv",
  candidates     = NULL)
assayNames(se)[1] <- "lfq_raw"
se_filt <- filter_perseus(se, perc_na = 0.33, filter_mode = "each_group")
se_imp  <- impute_perseus(se_filt)
se_diff <- add_rejections(
  test_diff(se_imp, type = "control", control = "Ctrl"),
  alpha = 0.05, lfc = 1)

## ---- dataset-description QC (DEP2) ---------------------------------------
save_fig("01_numbers",    plot_numbers(se))
save_fig("02_frequency",  plot_frequency(se))
save_fig("03_missval",    plot_missval(se_filt), h = 6)
save_fig("04_detect",     plot_detect(se_filt))
save_fig("05_coverage",   plot_coverage(se))

## ---- imputation / normalization (DEP2) ----------------------------------
save_fig("06_imputation", plot_imputation(se_filt, se_imp))
save_fig("07_norm",       plot_normalization(se_filt, se_imp))

## ---- multivariate QC (DEP2) ---------------------------------------------
save_fig("08_pca", plot_pca(se_imp, indicate = "condition"))
save_fig("09_cor", plot_cor(se_diff, indicate = "condition"), w = 6)

## ---- results (proteoSE volcano; DEP2 heatmap) -------------------------------
save_fig("10_volcano_A",
         proteoSE::plot_volcano(se_diff, contrast = "TreatA_vs_Ctrl",
                            sig_from_column = TRUE, top_n_sign = 12)$plot_out,
         h = 6)
save_fig("11_volcano_B",
         proteoSE::plot_volcano(se_diff, contrast = "TreatB_vs_Ctrl",
                            sig_from_column = TRUE, top_n_sign = 12)$plot_out,
         h = 6)
sig_se <- se_diff[which(rowData(se_diff)$significant), ]
save_fig("12_heatmap",
         DEP2::plot_heatmap(sig_se, type = "centered", kmeans = TRUE, k = 4,
                            indicate = "condition"), h = 7)

message("Wrote 12 figures to ", figdir)
