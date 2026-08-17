# Import protein-level search-engine output into a SummarizedExperiment
#
# Part of the 'sev' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.


#' Create se from MaxQuant
#' 
#' Reads in MaxQuant or other output and creates a summarized experiment.
#'
#' @param file path to proteinGroups or similar file
#' @param gene_column name of gene_name column after janitor
#' @param protein_column name of protein column after janitor
#' @param sep character describing the separator between sample name and replicate number (e.g. "_rep_", "_r")
#' @param filt character vector of column names used for filtering. Genes with a "+" will be dropped.
#' @param keep_all_proteins if FALSE, only first protein per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param keep_all_genes if FALSE, only first gene per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param experimental_design dataframe with information regarding samples. If not specified,
#' sample names and groups are read automatically (start with "lfq_intensity"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @export
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all if_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names
#' @details
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#'
se_read_in <- function(file, gene_column = "gene_names", protein_column = "protein_ids", sep="_rep_",
                       filt = c("reverse", "only_identified_by_site", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){

  #ProteinGroups
  data <- read.delim(file, sep="\t")

  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()

  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_column,
           orig_gene_names = gene_column) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes) %>%
    dplyr::rename("perseus_intensity" = "intensity")

  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ .x == ""))
  
  #Make gene_names unique
  data_unique <- make_unique(data, gene_column, protein_column, delim=";")

  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){

    LFQ_labels <- colnames(data)[grep("^lfq_intensity_", colnames(data_unique))]

    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("lfq_intensity_", "", LFQ_labels),
                                    condition=gsub(paste0("lfq_intensity_|",sep, "[0-9].*"), "", LFQ_labels),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE))
  }else{
    LFQ_labels<-experimental_design$label
  }

  data_se <- make_se(data_unique, which(colnames(data_unique) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}



#' Spectronaut to se
#' 
#' Reads in spectronaut or other output and creates a summarized experiment.
#'
#' @param file path to proteinGroups or similar file
#' @param gene_column name of gene_name column after janitor
#' @param protein_column name of protein column after janitor
#' @param sep character describing the separator between sample name and replicate number (e.g. "_rep_", "_r")
#' @param keep_all_proteins if FALSE, only first protein per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param keep_all_genes if FALSE, only first gene per protein-group is kept. If TRUE,
#' entries are split in multiple rows
#' @param experimental_design dataframe with information regarding samples. If not specified,
#' sample names and groups are read automatically (start with "log2quantity_"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP2 make_unique
#' @details
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#' @export

spectronaut_read_in <- function(file, gene_column = "genes", protein_column = "uni_prot_ids", sep="_rep_", keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  #####read in data (spectronaut output from Fatih Demir)
  df <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
    rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
    janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
    dplyr::rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
    dplyr::rename("condition" = "r_condition", "file_name" = "r_file_name", "replicate" = "r_replicate") #remove the strange r_ from those three column names
  
  #turn dataframe (df) to wide format, which is required for summarised experiment (or perseus)
  data <- df %>% 
    mutate(cond_rep = paste0(condition, "_rep", replicate )) %>% # make yourself a new column with the condition and replicate pasted together
    DEP2::make_unique("genes", "uni_prot_ids", delim = ";") %>%#this is a function of DEP2 to make unique names of genes
    select(cond_rep, ID,genes,uni_prot_ids,  protein_descriptions,log2quantity, ibaq) %>% #I just select these columns because otherwise the 'pivot_wider' function is not working because each row contains non-unique info so everything is matched to everything
    tidyr::pivot_wider(names_from = cond_rep, values_from = c(log2quantity, ibaq)) %>% #make wide
    select(-ID) %>% #now I remove the ID col that was added when I ran make_unique, because I need to make unique again, because later I will need both ID and Name column
    DEP2::make_unique("genes", "uni_prot_ids", delim = ";") %>% #re-run
    janitor::clean_names() %>% #clean_up col-names
    dplyr::rename(ID = id) %>% #clean up col names further
    dplyr::rename_all(list(~as.character(gsub("0_2", "0p2", .)))) %>% #clean up further, as this underscore separation is annoying later on 
    mutate_if(is.numeric, function(x) ifelse(is.infinite(x), NA, x)) %>% # turns Infinite to NA
    mutate_if(is.numeric, function(x) ifelse(is.na(x), NA, x)) %>% #turns NaN to NA
    mutate_if(is.numeric, function(x) 2^x) %>% #inverse log2 of all numerics (careful if you change to include anything else, like the ibaqs that have not been inverse log2)
    mutate_if(is.numeric, function(x) ifelse(is.na(x), NA, x))  %>% #turns NaN to NA
    dplyr::rename("gene_names" = gene_column,
                  "protein_ids" = protein_column)
  
    #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = protein_ids,
           orig_gene_names = gene_names) %>%
    split_genes(colname = "protein_ids", keep_all = keep_all_proteins) %>%
    split_genes(colname = "gene_names", keep_all = keep_all_genes)
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    filename <- df %>% 
      mutate(cond_rep = paste0(tolower(condition), "_rep", replicate)) %>%
      select(file_name, cond_rep) %>%
      unique() %>%
      mutate(cond_rep = janitor::make_clean_names(cond_rep))
    
    LFQ_labels <- colnames(data)[grep("^log2quantity_", colnames(data))]
    
    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("^log2quantity_", "", LFQ_labels),
                                    condition=gsub(paste0("^log2quantity_|",sep, "[0-9].*"), "", LFQ_labels),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE)) %>%
      merge(filename, by.x = "sample", by.y = "cond_rep", all = T)
  }else{
    LFQ_labels<-experimental_design$label
  }
  
  data_se <- make_se(data, which(colnames(data) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}





#' fragpipe to se
#' 
#' Reads in fragpipe or other output and creates a summarized experiment.
#'
#' @param file path to proteinGroups or similar file
#' @param gene_column name of gene_name column after janitor
#' @param protein_column name of protein column after janitor
#' @param sep character describing the separator between sample name and replicate number (e.g. "_rep_", "_r")
#' @param experimental_design dataframe with information regarding samples. If not specified,
#' sample names and groups are read automatically (start with "log2quantity_"), end with "_r" followed by
#' replicate number.
#' @return summarized experiment
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP2 make_unique
#' @details
#' Mandatory columns: sample (sample name), condition (treatment), replicate
#' @export

fragpipe_read_in <- function(file, gene_column = "gene", protein_column = "protein_id", sep="_rep_", experimental_design = ""){
  
  #####read in data (spectronaut output from Fatih Demir)
  df <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
    rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
    janitor::clean_names() #make column names clean and unique (makes later coding easier)
    
  #turn dataframe (df) to wide format, which is required for summarised experiment (or perseus)
  data <- df %>% 
    dplyr::rename("gene_names" = gene_column,
                  "protein_ids" = protein_column) %>%
    make_unique("gene_names", "protein_ids")
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    LFQ_labels <- colnames(data)[grep(".*_max_lfq_intensity", colnames(data))]
    
    experimental_design<-data.frame(label=LFQ_labels,
                                    sample=gsub("_max_lfq_intensity", "", LFQ_labels),
                                    condition=gsub(paste0("_max_lfq_intensity|",sep, "[0-9]*"), "", LFQ_labels),
                                    replicate=gsub(paste0("_max_lfq_intensity|^.*",sep,"(?=[0-9])"), "", LFQ_labels, perl = TRUE)) 
    }else{
    LFQ_labels<-experimental_design$label
  }
  
  data_se <- make_se(data, which(colnames(data) %in% LFQ_labels), experimental_design)
  rownames(data_se) <- data$name
  names(assays(data_se)) <- "lfq_raw"
  return(data_se)
}



#' Read in Spectronaut Output
#'
#' This function reads in data from Spectronaut output, cleans and formats the data and finally returns the data as a summarized experiment.
#'
#' @param candidates Character string representing the path to the candidates file. Alternatively, a data frame can be passed.
#' @param report Character string representing the path to the report file. Alternatively, a data frame can be passed.
#' @param contrasts Character vector representing the contrasts.
#' @param conditionSetup Character string representing the path to the condition setup file. Alternatively, a data frame can be passed.
#' @param quant_col Character string representing the column name for the quantification data (default: "log2quantity"). If other than "log2quantiy", data is log2 transformed in the process!
#' 
#' @return A data frame with the processed and merged data.
#'
#' @export
#'
#' @importFrom vroom vroom
#' @importFrom janitor clean_names
#' @importFrom tidyr pivot_wider
#' @importFrom DEP2 make_se
#' @import dplyr

spectronaut_to_se <- function(candidates = NULL, report = NULL, contrasts = NULL, conditionSetup = NULL, quant_col = "log2quantity"){
    #####read in data (spectronaut output from Fatih Demir)
    if(typeof(candidates) == "character"){
      df_candidates <- vroom::vroom(candidates, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("r_", "", .)))
    }else{
      df_candidates <- candidates
    }
    #####read in data (spectronaut output from Fatih Demir)  
    if(typeof(report) == "character"){
      df_wide_report <- vroom(report, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        rename_all(.funs = list(~gsub("pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("r_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
        rename_all(.funs = list(~gsub("x[0-9]*_", "", .)))
    }else{
      df_wide_report <- report
    }
    #####read in colData
    if(typeof(conditionSetup) == "character"){
      coldata <- vroom(conditionSetup, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
        rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
        janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
        mutate(label = paste0(condition, "_", replicate)) %>% 
        dplyr::select(run_label, condition, replicate, file_name) %>%
        dplyr::rename(label = run_label) %>% 
        mutate(label = paste0(janitor::make_clean_names(label),"_", quant_col))
    }else{
      coldata <- conditionSetup
    }
    
    ##### get contrasts in correct "orientation"
    if(!is.null(contrasts)){
      df_candidates_new <- data.frame()
      for(con in contrasts){
        condition_num <- strsplit(con, "_vs_")[[1]][1]
        condition_den <- strsplit(con, "_vs_")[[1]][2]
        temp <- df_candidates %>% filter(condition_numerator == condition_num & condition_denominator == condition_den) %>%
          mutate(diff = avg_log2_ratio,
                 contrast = con)
        temp2 <- df_candidates %>% filter(condition_numerator == condition_den & condition_denominator == condition_num) %>%
          mutate(diff = -1 * avg_log2_ratio,
                 contrast = con)
        df_candidates_new <- rbind(df_candidates_new, temp, temp2)
      }
    }else{
      df_candidates_new <- df_candidates %>% 
        mutate(diff = avg_log2_ratio, 
               contrast = paste0(condition_numerator, "_vs_", condition_denominator))
    }
    ##### rename columns
    df_candidates_new <- df_candidates_new %>% dplyr::rename(p.val = pvalue, p.adj = qvalue) %>% 
      dplyr::select(diff, contrast, p.val, p.adj, protein_groups)
    ##### pivot
    df_candidates_new <- df_candidates_new %>% tidyr::pivot_wider(names_from = "contrast", values_from = c("diff", "p.val", "p.adj"), names_glue = "{contrast}_{.value}")
    ##### merge with report
    combined <- merge(df_candidates_new, df_wide_report, by = "protein_groups", all.y = T) %>%
      make_unique("genes", "protein_groups", delim = ";") %>%
      mutate(gene_names = genes,
             protein_ids = protein_groups)
    
    
    return(DEP2::make_se(combined, grep(paste0(".*_", quant_col), colnames(combined)), coldata, log2transform = ifelse(quant_col == "log2quantity", F, T))) 
  }


############################################################################
############################################################################
############################################################################
#' Create SummarizedExperiment object from Spectronaut reports
#'
#' This function processes Spectronaut reports and returns a SummarizedExperiment object. 
#' It allows for optional transformations like collapsing extra underscores in the 
#' condition column and converting condition values to lowercase.
#'
#' @param candidates Character or dataframe. Path to the candidates file or a dataframe.
#' @param report Character or dataframe. Path to the Spectronaut report file or a dataframe.
#' @param contrasts Character vector. Optional. Specify contrasts for differential analysis. 
#' If NULL, default diff calculations will be used.
#' @param conditionSetup Character or dataframe. Path to the condition setup file or a dataframe.
#' @param quant_col Character. Column name indicating the quantity column. Default is "log2quantity".
#' @param suggest_contrasts Logical. If TRUE, the function will suggest possible contrasts based on conditions. Default is TRUE.
#' @param collapse_conditions Logical. If TRUE, extra underscores in the condition column are collapsed to the minimum number of underscores. Default is FALSE.
#' @param conditions_tolower Logical. If TRUE, the condition column is converted to lowercase. Default is FALSE.
#'
#' @return A SummarizedExperiment object containing the processed data.
#'
#' @details 
#' The function processes Spectronaut reports, condition setups, and contrasts (if provided) 
#' to create a SummarizedExperiment object for downstream analysis. It allows for condition-specific 
#' transformations such as collapsing underscores and converting conditions to lowercase.
#' 
#' If `collapse_conditions` is TRUE, any extra underscores in the condition column will be collapsed 
#' to the minimum number of underscores found in the column.
#' 
#' If `cond_tolower` is TRUE, the values in the condition column will be converted to lowercase.
#'
#' @examples
#' \dontrun{
#' # Example with custom arguments
#' se <- spectronaut_to_se(candidates = "path_to_candidates.txt", 
#'                         report = "path_to_report.txt", 
#'                         conditionSetup = "path_to_conditions.txt", 
#'                         quant_col = "log2quantity", 
#'                         collapse_conditions = TRUE, 
#'                         cond_tolower = TRUE)
#'                         
#' # Example with defaults
#' se <- spectronaut_to_se()
#'
#' }
#' @export
optimized_spectronaut_to_se <- function(candidates = NULL, report = NULL, contrasts = NULL, 
                                        conditionSetup = NULL, quant_col = "log2quantity", 
                                        suggest_contrasts = TRUE, collapse_conditions = FALSE, 
                                        conditions_tolower = FALSE) {
  
  # Check if candidates data is provided and load if necessary
  if (!is.null(candidates)) {
    if (typeof(candidates) == "character") {
      df_candidates <- vroom::vroom(candidates, delim = "\t", col_names = TRUE, 
                                    guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
        dplyr::rename_all(tolower) %>% 
        janitor::clean_names() %>% 
        dplyr::rename_all(~ gsub("pg_", "", .)) %>% 
        dplyr::rename_all(~ gsub("^r_", "", .)) %>%  # Remove 'r_' only at the start of string
        dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA)))  # Replace NaN, -Inf, Inf with NA
    } else {
      df_candidates <- candidates
    }
  } else {
    df_candidates <- NULL
  }
  
  # Load and clean report data if it's a file path
  if (typeof(report) == "character") {
    df_wide_report <- vroom::vroom(report, delim = "\t", col_names = TRUE, 
                                   guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
      dplyr::rename_all(tolower) %>% 
      janitor::clean_names() %>% 
      dplyr::rename_all(~ gsub("pg_", "", .)) %>% 
      dplyr::rename_all(~ gsub("^r_", "", .)) %>%  # Remove 'r_' only at the start of string
      dplyr::rename_all(~ gsub("^x[0-9]*_", "", .)) %>%  # Remove "x" at beginning of colname, followed by numbers
      dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA)))  # Replace NaN, -Inf, Inf with NA
  } else {
    df_wide_report <- report
  }
  
  # Load and clean condition setup data if it's a file path
  if (typeof(conditionSetup) == "character") {
    coldata <- vroom::vroom(conditionSetup, delim = "\t", col_names = TRUE, 
                            guess_max = 30000, .name_repair = janitor::make_clean_names) %>% 
      dplyr::rename_all(tolower) %>% 
      janitor::clean_names() %>% 
      #dplyr::select(run_label, condition, replicate, file_name) %>%
      dplyr::select(-label) %>% #remove in case you have a label column already as we are making a new label col.
      mutate(condition = gsub("[^[:alnum:]]", "_", condition)) %>% 
      dplyr::rename(label = run_label) %>% 
      dplyr::mutate(label = paste0(janitor::make_clean_names(label), "_", quant_col)) %>% 
      dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA))) %>%  # Replace NaN, -Inf, Inf with NA
      dplyr::mutate(renamingcol = sub("^x", "", label)) %>% 
      dplyr::mutate(renamingcol = gsub("_raw_log2quantity$", "", renamingcol)) %>% 
      dplyr::mutate(condition_rep = paste0(condition, "_rep", replicate))
  } else { #load and clean condition setup data if it is a tibble or data frame (same code as above.)
    coldata <- conditionSetup
    coldata <- coldata %>% 
      janitor::clean_names() %>%
      dplyr::rename_all(tolower) %>% 
      janitor::clean_names() %>% 
      #dplyr::select(run_label, condition, replicate, file_name) %>%
      dplyr::select(-label) %>% #remove in case you have a label column already as we are making a new label col.
      mutate(condition = gsub("[^[:alnum:]]", "_", condition)) %>% 
      dplyr::rename(label = run_label) %>% 
      dplyr::mutate(label = paste0(janitor::make_clean_names(label), "_", quant_col)) %>% 
      dplyr::mutate(across(where(is.numeric), ~ replace(.x, !is.finite(.x), NA))) %>%  # Replace NaN, -Inf, Inf with NA
      dplyr::mutate(renamingcol = sub("^x", "", label)) %>% 
      dplyr::mutate(renamingcol = gsub("_raw_log2quantity$", "", renamingcol)) %>% 
      dplyr::mutate(condition_rep = paste0(condition, "_rep", replicate))
  }
  
  # If conditions_tolower is TRUE, convert the condition column to lowercase
  if (conditions_tolower) {
    coldata <- coldata %>%
      dplyr::mutate(condition = tolower(condition))
  }
  
  # If collapse_conditions is TRUE, collapse underscores in the condition column
  if (collapse_conditions) {
    # Find the minimum number of underscores in the condition column
    min_underscores <- min(str_count(coldata$condition, "_"))
    
    # Collapse extra underscores to match the minimum number
    coldata <- coldata %>%
      dplyr::mutate(
        condition = ifelse(
          stringr::str_count(condition, "_") > min_underscores,
          stringr::str_replace(condition, "(.*)_(.*)", "\\1\\2"), 
          condition
        )
      )
  }
  
  # Suggest possible contrasts using expand.grid based on unique conditions
  if (suggest_contrasts) {
    unique_conditions <- coldata %>% dplyr::distinct(condition) %>% dplyr::pull(condition)
    possible_contrasts <- expand.grid(condition_num = unique_conditions, 
                                      condition_den = unique_conditions) %>%
      dplyr::filter(condition_num != condition_den) %>%
      dplyr::mutate(contrast = paste0(condition_num, "_vs_", condition_den)) %>%
      dplyr::arrange(contrast) %>% 
      dplyr::pull(contrast)
    
    # Print possible contrasts to the console
    message("Possible contrasts based on conditions in the condition setup:")
    message(paste(possible_contrasts, collapse = "\n"))
  }
  
  # If candidates is not NULL, process contrasts or use default diff calculations
  if (!is.null(df_candidates)) {
    if (!is.null(contrasts)) {
      df_candidates_new <- data.frame()
      for (con in contrasts) {
        condition_num <- strsplit(con, "_vs_")[[1]][1]
        condition_den <- strsplit(con, "_vs_")[[1]][2]
        temp <- df_candidates %>% 
          dplyr::filter(condition_numerator == condition_num & condition_denominator == condition_den) %>% 
          dplyr::mutate(diff = avg_log2_ratio, contrast = con)
        temp2 <- df_candidates %>% 
          dplyr::filter(condition_numerator == condition_den & condition_denominator == condition_num) %>% 
          dplyr::mutate(diff = -1 * avg_log2_ratio, contrast = con)
        df_candidates_new <- dplyr::bind_rows(df_candidates_new, temp, temp2)
      }
    } else {
      df_candidates_new <- df_candidates %>% 
        mutate(diff = avg_log2_ratio, contrast = paste0(condition_numerator, "_vs_", condition_denominator))
    }
    
    # Prepare the candidates data
    df_candidates_new <- df_candidates_new %>% 
      dplyr::rename(p.val = pvalue, p.adj = qvalue) %>% 
      dplyr::select(diff, contrast, p.val, p.adj, protein_groups) %>% 
      tidyr::pivot_wider(names_from = "contrast", values_from = c("diff", "p.val", "p.adj"), names_glue = "{contrast}_{.value}")
  } else {
    # # If no candidates are provided, create an empty placeholder
    # df_candidates_new <- data.frame(protein_groups = unique(df_wide_report$protein_groups))
    
    if ("protein_groups" %in% names(df_wide_report)) {
      protein_id_col <- "protein_groups"
    } else if ("uni_prot_ids" %in% names(df_wide_report)) {
      protein_id_col <- "uni_prot_ids"
    } else {
      stop("Neither 'protein_groups' nor 'uni_prot_ids' found in report.")
    }
    
    df_candidates_new <- data.frame(protein_groups = unique(df_wide_report[[protein_id_col]]))
  }
  
  # Merge candidates with the report data
  # combined <- df_candidates_new %>%
  #   dplyr::right_join(df_wide_report, by = "protein_groups") %>% 
  #   DEP2::make_unique("genes", "protein_groups", delim = ";") %>% 
  #   dplyr::mutate(gene_names = genes, protein_ids = protein_groups)
  combined <- df_candidates_new %>%
    dplyr::rename(!!protein_id_col := protein_groups) %>%
    dplyr::right_join(df_wide_report, by = protein_id_col) %>% 
    DEP2::make_unique("genes", protein_id_col, delim = ";") %>% 
    dplyr::mutate(gene_names = genes, protein_ids = .data[[protein_id_col]])
  
  
  # Update column names based on condition setup
  for (i in 1:nrow(coldata)) {
    current_prefix <- coldata$renamingcol[i]  # The prefix you want to replace
    new_prefix <- coldata$condition_rep[i]    # The new prefix
    
    # Replace the exact prefix only in the columns that match it fully
    colnames(combined) <- stringr::str_replace(colnames(combined), 
                                               paste0("^", current_prefix, "(?=_|$)"), 
                                               new_prefix)
  }
  
  # Adjust labels in coldata
  coldata <- coldata %>% 
    dplyr::mutate(old_label = label) %>% 
    dplyr::mutate(label = paste0(condition_rep, "_raw_", quant_col))
  
  # Create the SummarizedExperiment object
  se <- DEP2::make_se(combined, grep(paste0(".*_", quant_col), colnames(combined)), coldata, log2transform = ifelse(quant_col == "log2quantity", FALSE, TRUE))
  
  # Add 'sample' column for compatibility with SEV filtering functions
  colData(se)$sample <- colData(se)$ID
  
  # Return the SE object
  return(se)
}




#' Handling of pig derived data
#' 
#' Pig proteins are poorly annotated and have a different "style" when it comes to fasta files. Thus, the output of MaxQuant 
#' is not directly usable. This function takes in a MaxQuant file from an experiment using pig derived proteins and a fasta
#' file and returns MaxQuant like files with correct protein names and gene names. In addition, columns are added that specify
#' the organism in case proteins from multiple organisms are expected. In this case, the fasta input must contain entries for 
#' the other organisms as well in the "pig-style". 
#'
#' @param proteingroups String specifying path to ProteinGroups file
#' @param peptides String specifying path to ProteinGroups file
#' @param fasta String specifying path to ProteinGroups file
#' @param mult_org Logical specifying if fasta and experiment contain multiple organisms.
#' @param obj Logical specifying if results should be returned as one named list. If FALSE, results are saved as individual files. 
#' @importFrom dplyr mutate select
#' @return Either a named list containing all dataframes with correct protein names and gene names or locations where results have been stored.
#' @export
fix_maxquant_pig_annotation <- function(proteingroups, peptides, fasta, mult_org = FALSE, obj = FALSE){
  
  # Read fasta file and create additional file containing only headers
  system2(command = "grep", args = c("\"^>\"", fasta), stdout = paste0(gsub("\\..*", "", fasta), "_fasta_headers.txt"))
  
  # Read newly created header file
  fasta_headers_ <- read.delim(paste0(gsub("\\..*", "", fasta), "_fasta_headers.txt"), header = FALSE, col.names = "fasta")
  
  # Create annotation from fasta headers
  fasta_headers <- fasta_headers_ %>% 
    mutate(
      uniprot = sub("^>.*?\\|(.*?)\\|.*", "\\1", fasta, perl = TRUE), # extract UNIPROT ID (Between "|")
      
      uniprot_name = ifelse(grepl("\\|.*\\|[A-Z0-9]*_[A-Z0-9]* ", fasta), # Extract uniprot name if present
                            sub("^>.*?\\|.*\\|(.*?) .*", "\\1", fasta, perl = TRUE), 
                            NA),
      protein_name = ifelse(grepl("\\|.*\\|[A-Z0-9]*_[A-Z0-9]*", fasta), # Extract protein name
                            sub("^>.*?\\|.*?\\|.*? (.*?) OS.*", "\\1", fasta, perl = TRUE), # If uniprot_name is present 
                            sub("^>.*?\\|.*?\\|(.*?)\\ OS.*", "\\1", fasta, perl = TRUE)), #Else
      organism = ifelse(grepl("OS=", fasta), 
                        sub(".*OS=(.*? [a-z]*).*", "\\1", fasta, perl = TRUE), #Extract Organism ID if present. If no OS is specified, protein_names won't be returned
                        NA),
      organism_id = ifelse(grepl("OX=", fasta), 
                           sub(".* OX=([0-9]*).*", "\\1", fasta, perl = TRUE), #Extract Organism ID
                           NA),
      gene_name = ifelse(grepl("GN=.*? ", fasta), 
                         sub(".*GN=(.*?) .*", "\\1", fasta, perl = TRUE), 
                         ifelse(grepl("GN=.*?", fasta), 
                                sub(".*GN=.*?", "\\1", fasta, perl = TRUE), 
                                NA)
      )
    )
  
  # For peptides
  peptides_first<-read.delim(file = peptides) %>%
    split_genes(., "Proteins", FALSE) %>%
    merge(fasta_headers, by.x = "Proteins", by.y = "uniprot", all.x = T) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "organism_id") 
      else c("gene_name", "protein_name", "uniprot_name")
    )
  
  peptides_all<-read.delim(file = peptides) %>%
    split_genes(., "Proteins", TRUE) %>%
    merge(fasta_headers, by.x="Proteins", by.y = "uniprot", all.x = T) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "organism_id") 
      else c("gene_name", "protein_name", "uniprot_name")
    )
  
  # For protein groups
  protein_groups <- read.delim(file = proteingroups)
  
  protein_groups_first <- split_genes(protein_groups, "Protein.IDs", FALSE) %>%
    merge(fasta_headers, by.x = "Protein.IDs", by.y = "uniprot", all.x = TRUE) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name,
           Fasta.headers = fasta) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "fasta", "organism_id") 
      else c("gene_name", "protein_name", "fasta", "uniprot_name")
    )
  
  protein_groups_all <- split_genes(protein_groups, "Protein.IDs", TRUE) %>%
    merge(fasta_headers, by.x = "Protein.IDs", by.y = "uniprot", all.x = TRUE) %>%
    mutate(Protein.names = protein_name, 
           Gene.names = gene_name,
           Fasta.headers = fasta) %>%
    dplyr::select( 
      - if(!mult_org) c("gene_name", "protein_name", "uniprot_name", "organism", "fasta", "organism_id") 
      else c("gene_name", "protein_name", "fasta", "uniprot_name")
    )
  
  if(obj){
    return(list("fasta_headers" = fasta_headers, 
                "peptides_first" = peptides_first, 
                "peptides_all" = peptides_all, 
                "proteins_groups_first" = protein_groups_first,
                "protein_groups_all" = protein_groups_all, 
                "fasta_file" = fasta)
    )
  }
  
  else{
    write.table(fasta_headers, format(Sys.time(), "%Y%m%d_fasta_headers.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(peptides_first, format(Sys.time(), "%Y%m%d_peptides_first.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(peptides_all, format(Sys.time(), "%Y%m%d_peptides_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(protein_groups_first, format(Sys.time(), "%Y%m%d_protein_groups_first.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(protein_groups_all, format(Sys.time(), "%Y%m%d_protein_groups_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    return(
      paste0("Results are accessible in the folder:\n\"", getwd(), "\".\n\nFile names:\n", 
             format(Sys.time(), "%Y%m%d_fasta_headers.txt\n"),
             format(Sys.time(), "%Y%m%d_peptides_first.txt\n"),
             format(Sys.time(), "%Y%m%d_peptides_all.txt\n"),
             format(Sys.time(), "%Y%m%d_protein_groups_first.txt\n"),
             format(Sys.time(), "%Y%m%d_protein_groups_all.txt\n"),
             "\n\nInput fasta file was:\n", 
             fasta
      ) %>% 
        cat(.)
    ) 
  }
}
