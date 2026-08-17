# Build / annotate SummarizedExperiment objects
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' merge_se
#'
#' @param se A list of SummarizedExperiment objects
#' @param keep_all Logical indicating whether to keep all rows or only those present in all objects. Default is FALSE.
#'
#' @return A SummarizedExperiment object with the merged data.
#' @export
#'
#' @importFrom DEP2 make_se
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom tibble column_to_rownames
#' @import dplyr
merge_se <- function(se = list(), keep_all = FALSE){
  
  old_names <- names(se)
  new_names <- paste0("nof_stripped_sequences_identified_experiment_wide_", names(se))
  mw_names <- paste0("moleculaweight_", names(se))
  
  se <- lapply(seq_along(se), function(x){
    se_rd <- rowData(se[[x]]) %>% as.data.frame() 
    
    colnames(se_rd) <- gsub("nof_stripped_sequences_identified_experiment_wide", new_names[[x]], colnames(se_rd))
    colnames(se_rd) <- gsub("moleculaweight", mw_names[[x]], colnames(se_rd))
    
    se_cd <- colData(se[[x]]) %>% as.data.frame() %>% mutate(set = names(se)[[x]])
    
    se_as <- assay(se[[x]]) %>% as.data.frame()
    
    se_sa <- gsub("(.*)_(quantity|log2quantity)", "\\1", colnames(dplyr::select(se_rd, ends_with("_quantity"), ends_with("_log2quantity"))))
    
    se_rd <- se_rd %>% dplyr::select(name, matches(paste0(se_sa, collapse = "|")), contains("_vs_"),
                                     contains("nof_stripped_sequences_identified_experiment_wide"),
                                     contains("moleculaweight"))
    
    list(cd = se_cd, rd = se_rd, as = se_as, sa = se_sa)
  })
  
  names(se) <- old_names
  
  #Check for duplicated IDs and rename if necessary
  temp <- c()
  for(se_ in se){
    temp <- c(se_$cd$ID, temp)
  }
  
  if(any(duplicated(temp))){
    warning("Duplicated IDs! Renamed to avoid conflicts")
    
    for(i in seq_along(se)){
      se[[i]]$cd$ID <- paste0(names(se)[[i]], "_",se[[i]]$cd$ID)
      se[[i]]$cd$condition <- paste0(names(se)[[i]], "_",se[[i]]$cd$condition)
      colnames(se[[i]]$as) <- paste0(names(se)[[i]], "_", colnames(se[[i]]$as))
      
      
      colnames(se[[i]]$rd) <- gsub("(.*)_vs_(.*)", paste0(names(se)[[i]], "_", "\\1_vs_", names(se)[[i]], "_","\\2"), colnames(se[[i]]$rd))
      
    }
  }
  
  # Merge
  for(i in seq_along(se)){
    if(i == 1){
      rd <- se[[i]]$rd %>% dplyr::select(-name)
      cd <- se[[i]]$cd
      as <- se[[i]]$as
    }else{
      rd <- merge(rd, dplyr::select(se[[i]]$rd, -name), by = "row.names", all = keep_all, sort = FALSE) %>% tibble::column_to_rownames("Row.names") %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      cd <- rbind(cd, se[[i]]$cd) %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      as <- merge(as, se[[i]]$as, by = "row.names", all = TRUE, sort = FALSE, ) %>% tibble::column_to_rownames("Row.names") %>%
        mutate(across(everything(), ~ ifelse(.=="NaN", NA, .)))
      
      colnames(as) <- cd$ID  
      rownames(cd) <- cd$ID
    }
  }
  
  
  se <- SummarizedExperiment::SummarizedExperiment(assay = as.matrix(as), rowData = rd, colData = cd)
  rowData(se)$name <- rownames(se)
  
  return(se)
}



#' Explore type of missingness
#' 
#' Add column to rowData that specifies for each row if values are missing at random or
#' not (MAR vs MNAR). The method is very rudimentary:
#' 1. Missing values are replaced by zeros, measured values by 1.
#' 2. ANOVA between conditions is performed to determine significant differences
#' of the mean.
#' 3. Rows with p-values < 0.05 are termed missing not at random.
#'
#' @param se A SummarizedExperiment object.
#'
#' @return se with added rowData column
#' @importFrom HybridMTest row.oneway.anova
#' @importFrom fdrtool gcmlcm
#' @export
add_randna <- function(se){
  
  # replace na with 0, values with 1
  dummy <- as.matrix(ifelse(is.na(assay(se)), 0,1))
  
  #ANOVA
  oa <- HybridMTest::row.oneway.anova(dummy, colData(se)$condition)
  
  #Add information regarding mode of missingness (MAR = TRUE)
  rowData(se)$randna<-ifelse(is.na(oa$pval), FALSE, ifelse(oa$pval<0.05, FALSE, TRUE))
  
  return(se)
}


#' Add Significance Information to a SummarizedExperiment Object
#'
#' This function takes a SummarizedExperiment object, calculates significance information,
#' and adds it to the rowData.
#'
#' @param se A SummarizedExperiment object.
#' @param p_thr A numeric value for the p-value threshold (default: 0.05).
#' @param diff_thr A numeric value for the difference threshold (default: 1).
#' @return A SummarizedExperiment object with added significance information.
#' @importFrom tidyr pivot_longer pivot_wider
#' @import dplyr
#' @examples
#' \dontrun{
#' # Assuming 'se' is a SummarizedExperiment object
#' se_with_sign <- add_significance(se)
#' }
#' @export
add_significance <- function(se, p_thr = 0.05, diff_thr = 1){
  se_long <- get_rowdata(se) %>% 
  dplyr::select(c(site_id_mult, 
                  ends_with("_diff"), 
                  ends_with("p.adj"), 
                  ends_with("p.val"), 
                  ends_with("CI.L"), 
                  ends_with("CI.R"))) %>% 
    tidyr::pivot_longer(cols = -site_id_mult, 
                        names_to = c("label", "type"), 
                        values_to = "value", 
                        names_pattern = "(.*)_(.*)") %>% 
    filter(!label == "score") %>%
    tidyr::pivot_wider(names_from = "type", 
                       values_from = "value") %>% 
    mutate(significant = ifelse(p.adj < p_thr & (diff > diff_thr | diff < -diff_thr), TRUE, FALSE)) %>% 
    tidyr::pivot_wider(names_from = "label", 
                       values_from = c("p.val","p.adj", "significant", "diff", "CI.L", "CI.R"), 
                       names_glue = "{label}_{.value}") %>% 
    mutate(significant = ifelse(if_any(ends_with("significant"), ~ .x == TRUE), TRUE, FALSE)) %>%
    dplyr::select(-site_id_mult)
  
  rowData(se)[, colnames(se_long)] <- NULL
  
  rowData(se) <- cbind(rowData(se), se_long)
  
  return(se)

}



############################################################################
############################################################################
############################################################################
#' Add statistical summaries to SummarizedExperiment
#'
#' This function calculates various statistical summaries for protein intensities
#' in a SummarizedExperiment object. The user can specify which statistics to compute
#' via the `type` argument, with the default being "all" to calculate all available statistics.
#' The computed statistics are merged back into the `rowData` of the SummarizedExperiment object.
#'
#' @param se A `SummarizedExperiment` object containing the protein data to be summarized.
#' @param type A character vector specifying which statistics to calculate. Options include:
#'   \describe{
#'     \item{"all"}{(default) calculates all available statistics (mean, sd, min, max, percentiles, IQR, counts).}
#'     \item{Custom set, e.g. \code{c("mean", "sd", "min")}}{selects only the specified statistics. Available options are:
#'         \itemize{
#'           \item "mean": Mean intensity
#'           \item "sd": Standard deviation
#'           \item "min": Minimum intensity
#'           \item "max": Maximum intensity
#'           \item "p25": 25th percentile intensity
#'           \item "p50": 50th percentile (median) intensity
#'           \item "p75": 75th percentile intensity
#'           \item "iqr": Interquartile range
#'           \item "missing_count": Count of missing (NA) values
#'           \item "measured_count": Count of non-NA values
#'         }
#'     }
#'   }
#'
#' @return The updated `SummarizedExperiment` object with the selected statistics added to the `rowData`.
#'
#' @examples
#' \dontrun{
#' # Example usage with default (all statistics):
#' updated_se <- add_stats(se)
#'
#' # Example usage to return only mean, sd, and min statistics:
#' updated_se_selected <- add_stats(se, type = c("mean", "sd", "min"))
#'
#' # Example usage to calculate mean and missing_count only:
#' updated_se_mean_missing <- add_stats(se, type = c("mean", "missing_count"))
#'
#' }
#' @export
add_stats <- function(se, type = "all") {
  # Extract data in long format from the SummarizedExperiment object
  df_long <- DEP2::get_df_long(se)
  
  # Convert all non-finite values (NaN, Inf, -Inf) to NA
  df_long <- df_long %>%
    mutate(intensity = ifelse(is.finite(intensity), intensity, NA))
  
  # Define the list of possible statistics and their respective calculation methods with the is.na check
  stats_to_calculate <- list(
    mean_intensity = ~ ifelse(all(is.na(.x)), NA, mean(.x, na.rm = TRUE)),
    sd_intensity = ~ ifelse(all(is.na(.x)), NA, sd(.x, na.rm = TRUE)),
    min_intensity = ~ ifelse(all(is.na(.x)), NA, min(.x, na.rm = TRUE)),
    max_intensity = ~ ifelse(all(is.na(.x)), NA, max(.x, na.rm = TRUE)),
    p25_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.25, na.rm = TRUE)),
    p50_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.50, na.rm = TRUE)),  # Median
    p75_intensity = ~ ifelse(all(is.na(.x)), NA, quantile(.x, probs = 0.75, na.rm = TRUE)),
    iqr_intensity = ~ ifelse(all(is.na(.x)), NA, IQR(.x, na.rm = TRUE)),  # Interquartile range
    missing_count = ~ sum(is.na(.x)),  # Count of NA values
    measured_count = ~ sum(!is.na(.x))  # Count of non-NA values
  )
  
  # If type is "all", calculate all stats; otherwise, filter the stats to calculate
  if (type == "all") {
    selected_stats <- names(stats_to_calculate)
  } else {
    selected_stats <- intersect(names(stats_to_calculate), paste0(type, "_intensity"))
    if ("missing_count" %in% type) selected_stats <- c(selected_stats, "missing_count")
    if ("measured_count" %in% type) selected_stats <- c(selected_stats, "measured_count")
  }
  
  # Summarize the selected statistics
  df_summaries <- df_long %>%
    group_by(name, condition) %>%
    summarize(across(
      .cols = intensity,
      .fns = stats_to_calculate[selected_stats],
      .names = "{.fn}"
    )) %>%
    ungroup()
  
  # Pivot the data to a wide format, appending suffixes to columns
  df_summaries_wide <- df_summaries %>%
    tidyr::pivot_wider(
      names_from = condition, 
      values_from = selected_stats,
      names_glue = "{condition}_{.value}"  # Automatically create column names with suffixes
    )
  
  # Arrange columns alphabetically (with 'name' as the first column)
  df_summaries_wide <- df_summaries_wide %>%
    dplyr::select(name, sort(colnames(df_summaries_wide)[-1]))  # Keep 'name' first
  
  # Extract the existing rowData from the SummarizedExperiment object
  row_data <- rowData(se)
  
  # Merge the new summary data into rowData
  row_data <- as.data.frame(row_data) %>%
    left_join(df_summaries_wide, by = "name")  # Merge by 'name' which is the protein identifier
  
  # Update the rowData of the SummarizedExperiment object
  rowData(se) <- as(row_data, "DataFrame")  # Convert it back to a DataFrame
  
  # Return the updated SummarizedExperiment object
  return(se)
}



######
#' Convert Log2 Differences to Fold Change in SummarizedExperiment Object
#'
#' This function takes a SummarizedExperiment object and adds new columns containing fold-changes. 
#' The function looks for columns ending with "_diff" in the `rowData` slot and uses dplyr::mutate to 'unlog' them from log2 scale. 
#' It preserves `NA` values in the input and returns a modified SummarizedExperiment object.
#' The function also checks to ensure that the rowData of the summarizedExperiment object remains
#' identical before and after the transformation, except for the newly created fold change columns.
#' Upon success, it prints a confirmation message and lists the names of the newly added columns.
#'
#' @param se A SummarizedExperiment object. This object should contain 
#'   rowData with columns that end with "_diff", representing log2 differences.
#' @return A SummarizedExperiment object with the same structure, but with 
#'   additional columns where log2 values are converted to fold change.
#' @examples
#' # Assuming `se` is a SummarizedExperiment object:
#' # se <- log2fc_to_fc(se)
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr mutate across
#' @importFrom digest digest
#' @export
log2fc_to_fc <- function(se) {
  # Check if input is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("The input must be a SummarizedExperiment object.")
  }
  
  # Extract rowData
  rowdata_df <- as.data.frame(rowData(se))
  
  # Check if rowData is empty
  if (nrow(rowdata_df) == 0) {
    stop("rowData is empty. The SummarizedExperiment object must contain data in rowData.")
  }
  
  # Check if any columns end with '_diff'
  diff_cols <- grep('_diff$', colnames(rowdata_df))
  if (length(diff_cols) == 0) {
    stop("No columns ending with '_diff' found in rowData.")
  }
  
  # Store a digest (hash) of the original rowData before modifications
  original_rowdata_hash <- digest::digest(rowdata_df, algo = "sha256")
  
  # Perform mutation: convert log2 to fold change, while handling NA values
  rowdata_df <- rowdata_df %>%
    mutate(across(ends_with('_diff'), 
                  ~ ifelse(is.na(.x), NA, 2^.x),  # Convert log2 to fold change and keep NA
                  .names = "{col}_foldchange"))
  
  # Capture the names of the newly created columns
  new_foldchange_cols <- grep('_foldchange$', colnames(rowdata_df), value = TRUE)
  
  # Remove '_diff' from the newly created columns' names
  colnames(rowdata_df) <- sub('_diff_foldchange$', '_foldchange', colnames(rowdata_df))
  
  # Reassign the mutated rowData back to the SummarizedExperiment object
  rowData(se) <- rowdata_df
  
  # Check that the object is identical apart from the newly created columns
  modified_rowdata_hash <- digest::digest(rowdata_df[ , !grepl('_foldchange$', colnames(rowdata_df))], algo = "sha256")
  
  if (original_rowdata_hash != modified_rowdata_hash) {
    stop("The SummarizedExperiment object has been modified in a way other than adding fold change columns.")
  }
  
  # If successful, print confirmation and list newly added columns
  message("The RowData of the SummarizedExperiment object is identical apart from the newly added columns.")
  message("Newly added columns: ", paste(new_foldchange_cols, collapse = ", "))
  
  # Return the updated SummarizedExperiment object
  return(se)
}

