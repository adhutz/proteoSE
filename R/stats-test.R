# Differential-abundance testing and batch correction
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.


#' Long Test
#'
#' This function returns test results for an se object in long format.
#'
#' @param se An se object.
#' @return A data frame.
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect ends_with
#' @import dplyr
#' @examples
#' # Example usage:
#' # result <- test_diff_long(se)
#' @export
test_diff_long <- function(se){
  
  res <- se %>% 
    get_rowdata() %>% 
    select(-significant) %>% 
    tidyr::pivot_longer(c(ends_with("_p.adj"), 
                          ends_with("_p.val"), 
                          ends_with("_diff"), 
                          ends_with("_significant"), 
                          ends_with("_CI.L"), 
                          ends_with("_CI.R")), 
                        values_to = "val", 
                        names_to = c("contrast", "pn"), 
                        names_pattern = "(.*)_(.*)$") %>%
    dplyr::select(-ends_with("p.val"), -ends_with("p.adj")) %>%
    tidyr::pivot_wider(names_from = "pn", 
                       values_from = "val") %>% 
    mutate(significant = ifelse(significant == 0, 
                                FALSE, 
                                TRUE))
  return(res)
}


#' test_diff_limma
#'
#' @param se SummarizedExperiment with groups as condition column of colData
#' @param design_formula formula object specifying the design of the linear model. Best left unchanged
#' @param advanced_contrast character vector specifying the contrasts to be tested. e.g. c("condA_vs_condB" = "condA - condB", "condB_vs_sctrl" = "condB - (condA + condC)/2")
#' @param fdr.type character string specifying the method to use for multiple testing correction. Default is "Strimmer's qvalue(t)".
#'
#' @return SummarizedExperiment with the results of the test added to the rowData
#' @export
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom tidyr gather unite pivot_wider
#' @importFrom tibble rownames_to_column
#' @importFrom qvalue qvalue
#' @importFrom purrr map_df
#' @importFrom stats p.adjust
#' @importFrom fdrtool fdrtool
#' @importFrom SummarizedExperiment rowData colData assay 
#' @importFrom stats model.matrix 
#' @import dplyr

test_diff_limma <- function (se, design_formula = formula(~0 + condition), advanced_contrast = NULL, 
                       fdr.type = c("Strimmer's qvalue(t)", "Strimmer's qvalue(p)", 
                                    "BH", "Storey's qvalue")) 
{
  
  col_data <- colData(se)
  raw <- assay(se)
  
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>% 
    as.character() %>% .[-1]
  
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  
  cntrst = advanced_contrast
  
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst), 
                                      collapse = ", "))
  
  fit <- lmFit(raw, design = design)
  
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if (any(is.na(raw))) {
    for (i in cntrst) {
      
      covariates <- strsplit(i, " - |\\(|\\)| \\+ ") %>% unlist
      covariates <- covariates[!covariates == "" & !grepl("/", covariates)]
      
      single_contrast <- makeContrasts(contrasts = i, levels = covariates)
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 
                                                                             1]
      
    }
  }
  eB_fit <- eBayes(contrast_fit, trend = FALSE)
  retrieve_fun <- function(comp, fit = eB_fit, fdr.type) {
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf, 
                    confint = TRUE)
    res <- res[!is.na(res$t), ]
    if (fdr.type == "Strimmer's qvalue(t)") {
      fdr_res <- fdrtool(res$t, plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if (fdr.type == "Strimmer's qvalue(p)") {
      fdr_res <- fdrtool(res$P.Value, statistic = "pvalue", 
                         plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
    }
    if (fdr.type == "BH") {
      padj <- p.adjust(res$P.Value, method = "BH")
      res$qval = padj
    }
    if (fdr.type == "Storey's qvalue") {
      qval_res = qvalue::qvalue(res$P.Value)
      res$qval = qval_res$qvalues
    }
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  message(fdr.type)
  limma_res <- purrr::map_df(cntrst, retrieve_fun, fdr.type = fdr.type)
  
  limma_res$comparison = names(cntrst)[match(limma_res$comparison, 
                                             cntrst)]
  
  table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, 
                                       CI.R, t, P.Value, qval, comparison) %>% tidyr::gather(variable, 
                                                                                             value, -c(rowname, comparison)) %>% mutate(variable = recode(variable, 
                                                                                                                                                          logFC = "diff", t = "t.stastic", P.Value = "p.val", qval = "p.adj")) %>% 
    tidyr::unite(temp, comparison, variable) %>% tidyr::pivot_wider(names_from = temp, values_from = value)
  
  rowData(se) <- left_join(as.data.frame(rowData(se)), table, 
                       by = c("name" = "rowname"))
  return(se)
}

########################################################################
########################################################################
########################################################################
#' Extract and Combine LIMMA Results from a List of SummarizedExperiments
#'
#' Iterates over a named list of \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' objects produced by \pkg{DEP2}, extracts the LIMMA differential expression
#' results from \code{rowData()}, and returns a single tidy long-format data
#' frame with one row per protein x comparison x dataset combination.
#'
#' @param x A \strong{named} list of
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} objects.  Names
#'   are used to populate the \code{dataset} column in the output and are
#'   therefore required.  A common naming convention is to prefix names with
#'   \code{"imp_"} or \code{"nonimp_"} to distinguish imputed from
#'   non-imputed results (see \code{imputation_prefix}).
#' @param id_cols Character vector.  Non-metric \code{rowData} columns to
#'   retain as identifiers in the output.  Default:
#'   \code{c("name", "ID", "gene_names", "protein_ids")}.
#' @param imputation_prefix Named character vector mapping dataset name
#'   prefixes to imputation status labels, or \code{NULL} to skip adding
#'   an \code{imputation_status} column.  Default:
#'   \code{c(imp_ = "imputed", nonimp_ = "non-imputed")}.
#'
#' @details
#' ## Expected \code{rowData} structure
#' The function selects columns matching the regex
#' \code{(_p\\.adj|_diff|_p\\.val)$}, which is the column naming convention
#' used by \pkg{DEP2} / \pkg{limma} for adjusted p-values, log2 fold-changes,
#' and raw p-values respectively.  All other \code{rowData} columns are
#' dropped except those listed in \code{id_cols}.
#'
#' ## Pivoting logic
#' Metric columns are first pivoted to long format (one row per metric), the
#' comparison prefix and metric type are parsed from the column name, and the
#' result is pivoted back to wide format so that each row represents one
#' protein x comparison combination with dedicated \code{log2fc},
#' \code{pval}, and \code{adj_pval} columns.
#'
#' ## Imputation status
#' If \code{imputation_prefix} is non-\code{NULL}, an \code{imputation_status}
#' column is derived by matching the start of each \code{dataset} name against
#' the names of \code{imputation_prefix}.  Dataset names that match no prefix
#' receive \code{NA}.  Set \code{imputation_prefix = NULL} to suppress this
#' column entirely.
#'
#' ## Downstream convenience columns
#' \code{abs_log2fc} (absolute log2 fold-change) is always added to the
#' output as it is almost universally needed for volcano plots and ranking.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per protein x
#'   comparison x dataset, containing:
#' \describe{
#'   \item{\code{dataset}}{Name of the source SE object from \code{names(x)}.}
#'   \item{\code{row_id}}{Original \code{rownames} of the SE, useful for
#'     re-joining to the source object.}
#'   \item{id columns}{All columns listed in \code{id_cols} (e.g.
#'     \code{name}, \code{ID}, \code{gene_names}, \code{protein_ids}).}
#'   \item{\code{comparison}}{The contrast name as defined in \pkg{DEP2}
#'     (e.g. \code{"Treatment_vs_Control"}).}
#'   \item{\code{log2fc}}{Log2 fold-change (\code{_diff} column).}
#'   \item{\code{pval}}{Raw p-value (\code{_p.val} column).}
#'   \item{\code{adj_pval}}{Adjusted p-value (\code{_p.adj} column).}
#'   \item{\code{abs_log2fc}}{Absolute log2 fold-change.}
#'   \item{\code{imputation_status}}{Imputation label derived from the
#'     dataset name prefix, or \code{NA} if no prefix matched.  Omitted
#'     when \code{imputation_prefix = NULL}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' limma_df <- combine_limma_outputs(diff_final)
#'
#' # Custom identifier columns
#' limma_df <- combine_limma_outputs(
#'   diff_final,
#'   id_cols = c("name", "ID", "gene_names", "protein_ids", "gene_primary")
#' )
#'
#' # Custom imputation prefixes
#' limma_df <- combine_limma_outputs(
#'   diff_final,
#'   imputation_prefix = c(imp_ = "imputed", raw_ = "raw")
#' )
#'
#' # Suppress imputation_status column
#' limma_df <- combine_limma_outputs(diff_final, imputation_prefix = NULL)
#' }
#'
#' @importFrom purrr imap_dfr
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate case_when
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_remove str_detect str_starts
#'
#' @export
combine_limma_outputs <- function(
    x,
    id_cols           = c("name", "ID", "gene_names", "protein_ids"),
    imputation_prefix = c(imp_ = "imputed", nonimp_ = "non-imputed")
) {
  
  # -------------------------------------------------------------------------
  # 1.  Validate input: names are required to populate the dataset column
  # -------------------------------------------------------------------------
  if (is.null(names(x))) {
    stop(
      "The list must be named. ",
      "Names are used to identify each dataset in the output.",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 2.  Validate that id_cols exist in every SE before iterating
  #     Catches typos early rather than producing a cryptic pivot error later
  # -------------------------------------------------------------------------
  purrr::iwalk(x, function(se, nm) {
    rd_cols     <- colnames(as.data.frame(SummarizedExperiment::rowData(se)))
    missing_ids <- setdiff(id_cols, rd_cols)
    if (length(missing_ids) > 0) {
      stop(
        "SE '", nm, "': the following id_cols were not found in rowData: ",
        paste(missing_ids, collapse = ", "), ".",
        call. = FALSE
      )
    }
  })
  
  # -------------------------------------------------------------------------
  # 3.  Extract, reshape, and combine results across all SE objects
  # -------------------------------------------------------------------------
  all_results <- purrr::imap_dfr(x, function(se, dataset_name) {
    
    res <- as.data.frame(SummarizedExperiment::rowData(se))
    
    # --- 3a. Select identifier columns and LIMMA metric columns ---
    res_clean <- res |>
      dplyr::select(
        dplyr::all_of(id_cols),
        dplyr::matches("(_p\\.adj|_diff|_p\\.val)$")
      ) |>
      dplyr::mutate(
        dataset = dataset_name,
        row_id  = rownames(res)
      )
    
    # --- 3b. Pivot metric columns to long format ---
    #         Each metric column name encodes two pieces of information:
    #         the comparison name (prefix) and the metric type (suffix).
    #         We split these here before pivoting back to wide format.
    res_long <- res_clean |>
      tidyr::pivot_longer(
        cols      = dplyr::matches("(_p\\.adj|_diff|_p\\.val)$"),
        names_to  = "metric_full",
        values_to = "value"
      ) |>
      dplyr::mutate(
        # Strip the metric suffix to recover the comparison name
        comparison  = stringr::str_remove(metric_full, "(_p\\.adj|_diff|_p\\.val)$"),
        metric_type = dplyr::case_when(
          stringr::str_detect(metric_full, "_p\\.adj$") ~ "adj_pval",
          stringr::str_detect(metric_full, "_p\\.val$") ~ "pval",
          stringr::str_detect(metric_full, "_diff$")    ~ "log2fc"
        )
      ) |>
      dplyr::select(-metric_full)
    
    # --- 3c. Pivot back to wide: one row per protein x comparison ---
    tidyr::pivot_wider(
      res_long,
      names_from  = metric_type,
      values_from = value
    )
  })
  
  # -------------------------------------------------------------------------
  # 4.  Add convenience columns used in almost every downstream step
  # -------------------------------------------------------------------------
  all_results <- dplyr::mutate(
    all_results,
    abs_log2fc = abs(log2fc)
  )
  
  # -------------------------------------------------------------------------
  # 5.  Optionally derive imputation status from dataset name prefix
  #     Matches the start of each dataset name against the names() of
  #     imputation_prefix; first match wins, unmatched -> NA.
  # -------------------------------------------------------------------------
  if (!is.null(imputation_prefix)) {
    all_results <- dplyr::mutate(
      all_results,
      imputation_status = {
        # For each dataset name, find the first prefix that matches
        status <- rep(NA_character_, nrow(all_results))
        for (prefix in names(imputation_prefix)) {
          matched <- stringr::str_starts(all_results$dataset, prefix) & is.na(status)
          status[matched] <- imputation_prefix[[prefix]]
        }
        status
      }
    )
  }
  
  all_results
}

########################################################################
########################################################################
########################################################################
#' Remove Batch Effects for Visualisation
#'
#' Applies \code{\link[limma]{removeBatchEffect}} to the primary assay of a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} and stores the
#' result as a new assay named \code{"batch_corrected"}.  The original assay,
#' \code{rowData()}, and \code{colData()} are left entirely unchanged.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'   Must have a \code{batch} and a \code{condition} column in
#'   \code{colData(se)}, as these are passed directly to
#'   \code{limma::removeBatchEffect(batch = ...)} and
#'   \code{limma::removeBatchEffect(group = ...)} respectively.
#' @param batch_col Character scalar.  Name of the \code{colData} column
#'   containing the batch factor.  Default: \code{"batch"}.
#' @param group_col Character scalar.  Name of the \code{colData} column
#'   containing the biological group / condition factor.  This is passed as
#'   the \code{group} argument to \code{limma::removeBatchEffect} to protect
#'   biological variation from being removed along with the batch effect.
#'   Default: \code{"condition"}.
#' @param assay_name Character scalar.  Name of the new assay that will be
#'   appended to \code{se}.  Default: \code{"batch_corrected"}.
#'
#' @details
#' ## What this function does and does not do
#' \itemize{
#'   \item \strong{Does}: appends a new assay (\code{assay_name}) containing
#'     the batch-corrected intensity matrix.
#'   \item \strong{Does not}: modify the primary assay, \code{rowData()}, or
#'     \code{colData()}.  All original data remain accessible in the returned
#'     object.
#' }
#'
#' ## Visualisation only
#' Batch correction via \code{limma::removeBatchEffect} is intended for
#' \emph{visualisation} (PCA, heatmaps, volcano plots) and \strong{not} for
#' differential expression analysis.  For statistical testing, include
#' \code{batch} as a covariate in the design matrix passed to \pkg{limma}
#' or \pkg{DEP2} directly.
#'
#' ## Applying to a list of SE objects
#' Use \code{\link{lapply}} to apply this function to each element of a list:
#' \preformatted{
#' list_batch_corrected <- lapply(list_se, rm_batch)
#' }
#'
#' @return The input \code{SummarizedExperiment} with an additional assay
#'   named \code{assay_name} containing the batch-corrected matrix.  All
#'   other slots are identical to the input.
#'
#' @examples
#' \dontrun{
#' # Single SE object
#' se_corrected <- rm_batch(se)
#'
#' # Custom column names
#' se_corrected <- rm_batch(
#'   se,
#'   batch_col  = "plate",
#'   group_col  = "treatment",
#'   assay_name = "plate_corrected"
#' )
#'
#' # Apply to a list of SE objects
#' list_corrected <- lapply(list_se, rm_batch)
#'
#' # Inspect available assays after correction
#' SummarizedExperiment::assayNames(se_corrected)
#' }
#'
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importFrom limma removeBatchEffect
#'
#' @export
rm_batch <- function(
    se,
    batch_col  = "batch",
    group_col  = "condition",
    assay_name = "batch_corrected"
) {
  
  # -------------------------------------------------------------------------
  # 1.  Validate that the required colData columns exist
  # -------------------------------------------------------------------------
  missing_cols <- setdiff(c(batch_col, group_col), colnames(SummarizedExperiment::colData(se)))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following columns were not found in colData(se): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 2.  Warn if the assay name already exists (would silently overwrite)
  # -------------------------------------------------------------------------
  if (assay_name %in% SummarizedExperiment::assayNames(se)) {
    warning(
      "Assay '", assay_name, "' already exists and will be overwritten.",
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 3.  Compute batch-corrected matrix
  #     `group` protects biological signal from being removed -- always pass
  #     it even when groups are balanced, as it improves correction accuracy.
  # -------------------------------------------------------------------------
  corrected_matrix <- limma::removeBatchEffect(
    SummarizedExperiment::assay(se),
    batch = se[[batch_col]],
    group = se[[group_col]]
  )
  
  # -------------------------------------------------------------------------
  # 4.  Append as a new assay; original assay is untouched
  # -------------------------------------------------------------------------
  SummarizedExperiment::assays(se)[[assay_name]] <- corrected_matrix
  
  message(
    "Batch effect removed using '", batch_col, "' as batch and '",
    group_col, "' as group.  ",
    "New assay '", assay_name, "' appended.  ",
    "Available assays: ",
    paste(SummarizedExperiment::assayNames(se), collapse = ", "), "."
  )
  
  se
}




#########
#' Subset a Dataframe Based on Conditions, Contrasts, and Additional Columns
#'
#' This function subsets a dataframe by selecting specified condition columns, contrast columns, 
#' and optionally extra columns. It also supports selecting summary intensity columns for conditions, 
#' and ensures that the intensity summary columns are selected in a specific order.
#'
#' @param df A dataframe from which to subset.
#' @param id_col A string specifying the name of the ID column. Default is "name".
#' @param condition_list A character vector of conditions. Columns matching these conditions followed by numbers 
#'   will be selected.
#' @param contrast_list A character vector of contrasts. Columns starting with these contrasts will be selected.
#' @param get_summaries A logical value indicating whether to include summary intensity columns for the specified conditions.
#'   The intensity summary columns include mean, standard deviation, min, max, quartiles, IQR, missing counts, and measured counts.
#'   Default is TRUE.
#' @param extra_cols A character vector of additional column names to include in the subset. The columns must exist in the dataframe.
#'   Default is NULL, meaning no extra columns will be selected.
#'
#' @return A dataframe containing the selected columns in the following order: ID column, extra columns (if specified), 
#'   condition columns (including intensity summary columns if `get_summaries` is TRUE in the specified order), and contrast columns.
#'
#' @examples
#' # Example dataframe
#' df <- data.frame(
#'   name = c("A", "B", "C"),
#'   cond1_1 = c(1, 2, 3),
#'   cond1_mean_intensity = c(10, 20, 30),
#'   contrastA_1 = c(100, 200, 300),
#'   col2 = c("X", "Y", "Z")
#' )
#'
#' # Subset the dataframe with conditions and contrasts
#' subset_results(df, id_col = "name", condition_list = c("cond1"), contrast_list = c("contrastA"), extra_cols = c("col2"))
#'
#' @export
subset_results <- function(df, id_col = "name", condition_list, contrast_list, get_summaries = TRUE, extra_cols = NULL) {
  
  # Define intensity suffixes to look for in condition columns if get_summaries is TRUE
  intensity_suffixes <- c("mean_intensity", "sd_intensity", "min_intensity", 
                          "max_intensity", "p25_intensity", "p50_intensity", 
                          "p75_intensity", "iqr_intensity", "missing_count", 
                          "measured_count")
  
  # Create a pattern for condition columns: match the condition followed by any number at the end
  condition_cols <- unlist(lapply(condition_list, function(cond) {
    grep(paste0("^", cond, "_[0-9]+$"), names(df), value = TRUE)
  }))
  
  # Optionally add any columns that start with the condition and end with one of the intensity suffixes
  if (get_summaries) {
    intensity_cols <- unlist(lapply(condition_list, function(cond) {
      unlist(lapply(intensity_suffixes, function(suffix) {
        grep(paste0("^", cond, "_", suffix, "$"), names(df), value = TRUE)
      }))
    }))
    all_condition_cols <- c(condition_cols, intensity_cols)
  } else {
    all_condition_cols <- condition_cols
  }
  
  # Create a pattern for contrast columns: match columns starting with contrast strings
  contrast_cols <- unlist(lapply(contrast_list, function(contrast) {
    grep(paste0("^", contrast), names(df), value = TRUE)
  }))
  
  # Check if extra_cols is provided and exists in the dataframe
  if (!is.null(extra_cols)) {
    extra_cols <- extra_cols[extra_cols %in% names(df)]  # Only keep valid column names
  } else {
    extra_cols <- character(0)  # Empty vector if extra_cols is NULL
  }
  
  # Select columns in the desired order: first 'id_col', then extra_cols, then condition columns (including intensity if selected in the correct order), and finally contrast columns
  df_subset <- df %>%
    dplyr::select(all_of(id_col), all_of(extra_cols), all_of(all_condition_cols), all_of(contrast_cols))
  
  # Return the subsetted dataframe
  return(df_subset)
}

