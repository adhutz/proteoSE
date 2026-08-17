# SummarizedExperiment accessors and modifiers
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' Get column data from a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @return A data.frame containing the column data.
#' @examples
#' \dontrun{
#' # Assuming 'se' is a SummarizedExperiment object
#' get_coldata(se)
#' }
#' @export
get_coldata <- function(se) {
  return(as.data.frame(colData(se)))
}


#' Get row data from a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @return A data.frame containing the row data.
#' @examples
#' \dontrun{
#' # Assuming 'se' is a SummarizedExperiment object
#' get_rowdata(se)
#' }
#' @export
get_rowdata <- function(se) {
  return(as.data.frame(rowData(se)))
}


#' Convert a SummarizedExperiment object to a long format data frame
#'
#' This function takes a SummarizedExperiment object and converts it to a long format data frame.
#' It allows the user to specify the assays to include in the output data frame.
#'
#' @param se A SummarizedExperiment object.
#' @param assays_ A character vector of assay names to include in the output data frame. Default is an empty character vector, which selects all assays.
#' @return A data frame in long format with assay values, rowData, and colData combined.
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @examples
#' \dontrun{
#' # Assuming 'se' is a SummarizedExperiment object
#' result <- se_to_long(se, assays_ = c("assay1", "assay2"))
#'
#' # To select all assays
#' result_all <- se_to_long(se)
#' }
#' @export
se_to_long <- function(se, assays_ = c("")){
  
  if(length(assays_) == 0){
    assays_ <- names(assays(se)) 
  }
  
  df <- lapply(assays(se)[assays_], "as.data.frame")
  
  df <- do.call("cbind", df) %>% cbind(as.data.frame(rowData(se)))
  
  assays <- grep(paste(assays_, "\\.", sep = ""), colnames(df), value = TRUE) %>% gsub("\\..*","",.) %>% unique()
  
  df_ <- df %>% dplyr::select(-ID) %>% tidyr::pivot_longer(grep(paste(assays_,"\\.", sep = ""), colnames(df), value = TRUE), names_to = c("assays", "ID"), names_pattern = "(.*)\\.(.*)")
  
  df_ <- merge(df_, as.data.frame(colData(se)), by.x = "ID")
  #%>% tidyr::pivot_longer(paste0("assay_",names(assays(se))))
  
  return(df_)
}


#' Add assay to se
#' 
#' Convinience function to add a new assay directly as the main assay. The old assay is retained. 
#' @param se Summarized experiment
#' @param new_assay Assay to add as the main assay (first in the assays() list)
#' @param name Name of the added assay
#' @param withDimnames Logical; passed to \code{SummarizedExperiment::assays<-}.
#'   If \code{TRUE} (default), the dimnames of the new assay are checked against
#'   those of \code{se}.
#'
#' @return se with added assay. Old main assay is retained.
#' @export
#'
add_assay <- function(se, new_assay, name, withDimnames = TRUE){
  
    if(name %in% names(assays(se))){
      names(assays(se))[names(assays(se)) == name] <- paste0(name, "_old")
    }
    temp <- assay(se)
    temp_n <- names(assays(se))[1]
    assay(se, withDimnames = withDimnames) <- new_assay
    assays(se, withDimnames = withDimnames)[[1]] <- new_assay
    names(assays(se))[1] <- name
    
    # add old assay
    assays(se)[[temp_n]] <- temp
    
    return(se)
  }


#####
#' Modify rowData or colData of a SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @param datatype A character string specifying whether to modify "rowData" or "colData". Must be one of "rowData" or "colData".
#' @param modify_function A function to apply to the selected data (e.g., dplyr functions like mutate or left_join).
#' @param ... Additional arguments passed to the modify_function.
#'
#' @return A modified SummarizedExperiment object with updated rowData or colData.
#'
#' @examples
#' \dontrun{
#' se <- SummarizedExperiment(assays = list(counts = matrix(1:12, nrow = 3)))
#' colData(se) <- DataFrame(sample = c("S1", "S2", "S3"))
#' rowData(se) <- DataFrame(name = c("gene1", "gene2", "gene3"))
#'
#' # Modify colData by adding a new column
#' se <- se_modify(se, "colData", mutate, new_col = c("A", "B", "C"))
#'
#' # Modify rowData by adding a new column
#' se <- se_modify(se, "rowData", mutate, new_col = c(10, 20, 30))
#'
#' # Example that would trigger an error by scrambling 'sample' column in colData
#' try(se <- se_modify(se, "colData", mutate, sample = sample(sample)))
#'
#' # Example that would trigger an error by scrambling 'name' column in rowData
#' try(se <- se_modify(se, "rowData", mutate, name = sample(name)))
#' }
#'
#' @import SummarizedExperiment
#' @import dplyr
#' @export
se_modify <- function(se, datatype = c("rowData", "colData"), modify_function, ...) {
  # Ensure datatype is either 'rowData' or 'colData'
  datatype <- match.arg(datatype)
  
  # Extract the relevant data (either rowData or colData)
  data <- if (datatype == "rowData") {
    as.data.frame(rowData(se))
  } else {
    as.data.frame(colData(se))
  }
  
  # Store original rownames for safety checks
  original_rownames <- rownames(data)
  
  # Store additional column check depending on datatype
  if (datatype == "rowData") {
    original_name_column <- data$name
  } else {
    original_sample_column <- data$sample
  }
  
  # Apply the specified function to the data
  modified_data <- modify_function(data, ...)
  
  # Safety checks for colData
  if (datatype == "colData") {
    # Check if rownames are still in the same order
    if (!identical(original_rownames, rownames(modified_data))) {
      stop("Error: The rownames of colData have changed. Ensure that the operation preserves rownames.")
    }
    
    # Check if rownames match the 'sample' column
    if (!all(original_rownames == modified_data$sample)) {
      stop("Error: The rownames of colData do not match the 'sample' column after modification.")
    }
  }
  
  # Safety checks for rowData
  if (datatype == "rowData") {
    # Check if rownames are still in the same order
    if (!identical(original_rownames, rownames(modified_data))) {
      stop("Error: The rownames of rowData have changed. Ensure that the operation preserves rownames.")
    }
    
    # Check if rownames match the 'name' column
    if (!all(original_rownames == modified_data$name)) {
      stop("Error: The rownames of rowData do not match the 'name' column after modification.")
    }
  }
  
  # Assign the modified data back to the appropriate slot in SummarizedExperiment
  if (datatype == "rowData") {
    rowData(se) <- S4Vectors::DataFrame(modified_data)
  } else {
    colData(se) <- S4Vectors::DataFrame(modified_data)
  }
  
  return(se)
}



################################################
################################################
################################################
#' Get Data in Wide Format with Specified/Available Assays Information from summarized experiment
#'
#' This function converts a SummarizedExperiment object into a wide-format data frame.
#' It appends all or selected assays data to the regular DEP2 table created with DEP2::get_df_wide()
#'
#' @param se A SummarizedExperiment object containing the assay(s).
#' @param assays Optional character vector specifying the assay(s) to include.
#'        If NULL or not provided, all available assays from the SummarizedExperiment
#'        are used.
#'
#' @return A data frame in wide format with selected/all assays data joined via dplyr::full_join().
#' 
#' @examples
#' \dontrun{
#' # Assuming se_diff is an existing SummarizedExperiment object
#' df_out <- get_df_wide_append_assay(se_diff, assays = c("lfq_raw", "imputed_DEP2"))
#' df_out_all <- get_df_wide_append_assay(se_diff)  # Uses all available assays
#' }
#'
#' @export

get_df_wide_append_assay <- function(se, assays = NULL) {
  # If assays is NULL, retrieve all assay names from SummarizedExperiment
  if (is.null(assays)) {
    assays <- names(assays(se))
    print("Using all available assays.")
  } else {
    print("Using specified assays.")
  }
  
  # Print the assays being processed
  message("Assay names: ", paste(assays, collapse = ", "))
  
  if (length(assays) < 1) {
    stop("No valid assays provided or available in the SummarizedExperiment object.")
  }
  
  # Proceed with the original dataframe and merging process using selected assays
  # Initial df using DEP2::get_df_wide(), and some sorting
  df <- get_df_wide(se) %>% 
    dplyr::select(name,gene_names,protein_ids, protein_descriptions = description, ID, orig_prot_ids = protein, contains('_vs_'),1:ncol(.))
  
  # Initialize an empty tibble for the eventual joins
  result <- tibble(name = df$name)
  
  # Loop through each specified assay, process data, and join
  for (assay_name in assays) {
    if (!assay_name %in% names(assays(se))) {
      warning(paste("Assay", assay_name, "not found in SummarizedExperiment; it will be skipped."))
      next
    }
    
    assay_data <- as.data.frame(assay(se, assay_name)) %>%
      rownames_to_column("name") %>%
      as_tibble() %>%
      rename_with(~paste0(.x, "_", assay_name), -name)
    
    # Join with the main result DataFrame
    result <- result %>%
      full_join(assay_data, by = "name")
  }
  
  # Join all together with initial df
  result2 <- df %>%
    full_join(result, by = "name")
  
  return(result2)
}

