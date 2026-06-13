# Experimental design: blocking and randomisation
#
# Part of the 'sev' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.





############################################################################
############################################################################
############################################################################
#' Evenly spread/randomize samples into manageable blocks for batch-processing
#'
#' The function considers batches, conditions and bioreplicates, and evenly spreads/randomizes samples into manageable experimental blocks.
#' For example, you have an experiment with 200 samples (4 conditions,n=50 per condition).
#' But you cannot process all of these in a single block, because you can only process a maximum of 30 samples at one time.
#' Therefore you need to split your lab work into blocks of a certain size..
#' This function assists in the blocking, so that conditions (and any other prior batching or cohorts) are spread EVENLY across your blocks.
#' This ensures that 
#' i) the effects you are measuring between your conditions are not clouded/overshadowed by blocking-effects and 
#' ii) the effects of the blocking are not accidentally interpreted by you as 'real' effects caused by your conditions
#'
#' @param df Input Data frame that contains your blocked.
#' @param batch_col Column name representing batches (in case the samples come from different batches, e.g. cohorts or prior batching during processing), defaults to "batch".
#' @param condition_col Column name representing your test conditions, defaults to a column called "condition" in your dataframe.
#' @param bio_replicate_col Column name representing biological replicates, defaults to "bio_replicate".
#' @param scramble Logical flag (TRUE or FALSE) to determine whether rows should be shuffled, defaults to TRUE.
#' @param seed Optional numeric seed for reproducibility, defaults to NULL (no seed). But the recommendation is to set the seed to make reproducible.
#' @return A data frame with rows rearranged into blocks. 
#' The main outcome is column called "block". 
#' The "block" column contains a number for your experimental blocks. Importantly, the "block_size" column contains the minimum of samples that you have to run together in one block.
#' If this number exceeds what you can process in one block, then you need to strike some compromises. Compromise options include i) leaving out batches all togehter or,
#' ii) in-silico reducing of any prior batching/cohorting, e.g. if you have 4 prior batches, turn them into 2 evenly mixed in-silico batches).
#' 
#' @examples
#' \dontrun{
#' set.seed(42)
#' 
#' #' Define the conditions, batches, and maximum bioreplicates
#' conditions <- c("KO_ctrl", "WT_ctrl", "KO_treated", "WT_treated")
#' num_batches <- 3
#' max_bioreplicates <- 5
#' 
#' #' Initialize lists to store data
#' conditions_list <- c()
#' batches_list <- c()
#' bioreplicates_list <- c()
#' 
#' #' Generate random data for the sample data frame
#' for (condition in conditions) {
#'   for (batch in 1:num_batches) {
#'     num_bioreplicates <- sample(1:max_bioreplicates, 1)
#' 
#'     conditions_list <- c(conditions_list, rep(condition, num_bioreplicates))
#'     batches_list <- c(batches_list, rep(paste0("Batch", batch), num_bioreplicates))
#'     bioreplicates_list <- c(bioreplicates_list, seq(1, num_bioreplicates))
#'   }
#' }
#' 
#' #' Create a data frame
#' sample_df <- data.frame(
#'   condition = factor(conditions_list),
#'   batch = factor(batches_list),
#'   bio_replicate = bioreplicates_list
#' )
#' 
#' #' Display the sample data frame
#' print(sample_df)
#' 
#' 
#' 
#' result <- block_randomize(df = sample_df,
#'                                 batch_col = "batch",
#'                                 condition_col = "condition",
#'                                 bio_replicate_col = "bio_replicate",
#'                                 scramble = T,
#'                                 seed = 1000)
#'                                 
#' }
#' @export

block_randomize <- function(df, batch_col = "batch", condition_col = "condition", bio_replicate_col = "bio_replicate", scramble = TRUE, seed = NULL) {
  # If batch_col is left empty, create a new 'batch' column with all 1s as values to indicate it is a single batch.
  # If batch_col is left empty, any pre-existing column called 'batch' is renamed to 'user_provided_batch'. (This is to ensure it is not overwritten or discarded)
  if (missing(batch_col) || batch_col == "") {
    batch_col <- "batch"
    if (exists(batch_col, where = df)) {
      df <- df %>% rename(NOTUSED_user_defined_batch_NOTUSED = !!sym(batch_col))
      batch_col <- "batch"
    }
    df <- df %>% mutate(!!batch_col := 1)
  }
  
  
  # Set seed for reproducibility if scramble is TRUE
  if (scramble) {
    if (!is.null(seed)) {
      set.seed(as.numeric(seed))
    } else {
      set.seed(12345)
    }
  }
  
  # Clean and the dataframe up and then scramble if desired
  tryCatch({
    df_cleaned <- df %>%
      dplyr::rename(batch = {{batch_col}},
                    condition = {{condition_col}},
                    bio_replicate = {{bio_replicate_col}}) %>%
      dplyr::filter(if_any(everything(), ~ !is.na(.))) %>%
      dplyr::filter(if_any(condition, ~ !is.na(.))) %>%
      dplyr::arrange(if (scramble) sample(nrow(.)) else NULL) %>% #scramble if desired
      dplyr::arrange(batch, condition) #order rows by batch, condition
    
    # Calculate parameters that are important for the anticlust package
    min_set <- (length(unique(df_cleaned$batch)) * length(unique(df_cleaned$condition))) #this is the number of minimum samples to process simultaneously. The code multiplies the number of unique batches by the number of unique conditions. This operation calculates the total number of unique combinations of batches and conditions
    number_sets <- ceiling(nrow(df_cleaned) / min_set) #this number indicates how many min sets you have. performs a calculation to determine how many minimum sets (min_set) are needed to cover the total number of rows in a data frame (df_cleaned)
    
    # Create a new dataframe with categorical sampling using the anticlust package
    df_new <- df_cleaned %>%
      mutate(block = anticlust::categorical_sampling(., number_sets)) %>% #if number of sets =1 it doesn't work. not sure how it handles not even numbers.
      dplyr::arrange(block) %>%
      mutate(new_sample_number = row_number()) %>%
      dplyr::group_by(block) %>% 
      mutate(block_size = n()) %>% 
      ungroup() %>% 
      dplyr::group_by(block, condition) %>% 
      mutate(block_condition_size = n()) %>% 
      ungroup() %>% 
      mutate(bio_reps_scrambled = scramble) %>% 
      dplyr::select(new_sample_number, block, block_size,dplyr::contains("batch"), condition, block_condition_size, bio_replicate, bio_reps_scrambled, 1:ncol(.) )
    
    # Add seed_used column only if scramble is TRUE
    if (scramble) {
      df_new <- df_new %>%
        mutate(seed_used = ifelse(!is.null(seed), as.character(seed), "12345"))
    }
    
    return(df_new)
  }, error = function(e) {
    cat("\n")
    cat("The following Error occurred in the block_randomiz function: \n")
    cat("\n")
    cat(conditionMessage(e), "\n")
    cat("\n")
    cat("Most likely, you have mis-spelled a column name in the arguments. \n")
    cat("\n")
    cat("In case of error around argument K needing to be bigger than 1, you will have to process all samples in one go. \n")
    cat("If you don't want to do that, leave out the 'batch_col'argument.\n")
    cat("However, doing so jeopordizes your entire future pipeline, because any prior batching of samples will not be accounted for. \n")
    
    return(NULL)
  })
}



#' Plan Experimental Design with Flexible Batching
#'
#' This function generates a full factorial experimental design from specified conditions
#' (e.g., cell type, treatment, timepoint, replicates) and supports optional batching.
#' Replicates are always kept intact across batches. Batch summaries and exports to CSV/XLSX
#' are available as built-in helpers.
#'
#' @param type A character vector of experimental types (e.g. `"cell_line1"`, `"kidney"`, `"organoid"`).
#' @param treatment A character vector of treatments (e.g. `"control"`, `"tnfa"`, `"lps"`).
#' @param timepoint A character vector of timepoints (e.g. `"0h"`, `"24h"`, `"48h"`).
#' @param replicates Integer. Number of replicates (e.g. `5`).
#' @param batched Logical. Whether to assign rows to batches. e.g. batch1 = reps 1:3, batch2 = reps 4:6 or similar (`TRUE`/`FALSE`).
#' @param batch_list Named list of replicate indices for manual batching, e.g. `list(batch1 = 1:3, batch2 = 4:6)`.
#' @param batch_size Approximate maximum number of rows per batch for automatic batching (mutually exclusive with `batch_list`).
#' @param summarize_batches Logical. Print a summary of batches to console. Default is `TRUE`.
#' @param export_csv Logical. Export design to CSV. Default is `TRUE`.
#' @param export_xlsx Logical. Export design to XLSX. Default is `TRUE`.
#' @param csv_name File name for CSV export. Default is `"experiment_design.csv"`.
#' @param xlsx_name File name for XLSX export. Default is `"experiment_design.xlsx"`.
#'
#' @return A `data.frame` containing the full design with columns: `type`, `treatment`, `timepoint`, `rep`, `group`, and optionally `batch`.
#'
#' @examples
#' \dontrun{
#' # Basic full factorial design without batching
#' plan_experiment(
#'   type = c("cell_line1", "kidney"),
#'   treatment = c("control", "tnfa"),
#'   timepoint = c("0h", "24h"),
#'   replicates = 3,
#'   batched = FALSE
#' )
#'
#' # With automatic batching (replicates are grouped together)
#' plan_experiment(
#'   type = c("cell_line1", "kidney", "organoid"),
#'   treatment = c("control", "tnfa", "lps"),
#'   timepoint = c("0h", "24h", "48h"),
#'   replicates = 6,
#'   batched = TRUE,
#'   batch_size = 50
#' )
#'
#' # With manual batching
#' plan_experiment(
#'   type = c("cell_line1", "kidney"),
#'   treatment = c("control", "tnfa"),
#'   timepoint = c("0h", "24h"),
#'   replicates = 6,
#'   batched = TRUE,
#'   batch_list = list(batch1 = 1:3, batch2 = 4:6)
#' )
#'
#' }
#' @importFrom dplyr %>% bind_rows group_by summarise n
#' @importFrom utils write.csv
#' @export
plan_experiment <- function(type, treatment, timepoint, replicates,
                            batched = FALSE,
                            batch_list = NULL,
                            batch_size = NULL,
                            summarize_batches = TRUE,
                            export_csv = TRUE,
                            export_xlsx = TRUE,
                            csv_name = "experiment_design.csv",
                            xlsx_name = "experiment_design.xlsx") {
  
  # dplyr is a package Import; openxlsx (optional, for xlsx export) is checked
  # at the point of use below. No runtime install/library() is performed.

  # --- Build full factorial design ----------------------------------------
  rep_names <- paste0("rep", seq_len(replicates))
  df <- expand.grid(
    type = type,
    treatment = treatment,
    timepoint = timepoint,
    rep = rep_names,
    stringsAsFactors = FALSE
  )
  df$group <- with(df, paste(type, treatment, timepoint, rep, sep = "_"))
  
  # --- Batching logic -----------------------------------------------------
  if (batched) {
    if (!is.null(batch_list) && !is.null(batch_size)) {
      stop("Please specify either 'batch_list' or 'batch_size', not both.")
    }
    
    if (!is.null(batch_list)) {
      # Manual batching
      batch_lookup <- unlist(lapply(names(batch_list), function(batch) {
        reps <- paste0("rep", batch_list[[batch]])
        setNames(rep(batch, length(reps)), reps)
      }))
      df$batch <- batch_lookup[df$rep]
      if (any(is.na(df$batch))) {
        warning("Some replicates were not assigned to any batch.")
      }
    } else if (!is.null(batch_size)) {
      # Auto batching (keep reps together)
      reps <- split(df, df$rep)
      rep_sizes <- sapply(reps, nrow)
      
      batches <- list()
      current_batch <- list()
      current_size <- 0
      batch_index <- 1
      
      for (rep_name in names(rep_sizes)) {
        rep_df <- reps[[rep_name]]
        rep_size <- nrow(rep_df)
        
        if (current_size + rep_size > batch_size && current_size > 0) {
          batches[[paste0("batch", batch_index)]] <- bind_rows(current_batch)
          batch_index <- batch_index + 1
          current_batch <- list()
          current_size <- 0
        }
        
        current_batch[[rep_name]] <- rep_df
        current_size <- current_size + rep_size
      }
      
      if (length(current_batch) > 0) {
        batches[[paste0("batch", batch_index)]] <- bind_rows(current_batch)
      }
      
      df <- bind_rows(lapply(names(batches), function(batch_name) {
        b <- batches[[batch_name]]
        b$batch <- batch_name
        b
      }))
    } else {
      stop("When 'batched = TRUE', provide either 'batch_list' or 'batch_size'.")
    }
  }
  
  # --- Summarize batches --------------------------------------------------
  if (summarize_batches && "batch" %in% names(df)) {
    summary <- df %>%
      group_by(batch, rep) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(batch) %>%
      summarise(
        total_rows = sum(n),
        replicates = paste(rep, collapse = ", ")
      )
    message("\nBatch Summary:")
    print(summary)
  }

  # --- Export files -------------------------------------------------------
  if (export_csv) {
    write.csv(df, file = csv_name, row.names = FALSE)
    message("\nCSV saved to: ", csv_name)
  }
  if (export_xlsx) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for xlsx export; install it or set export_xlsx = FALSE.")
    }
    openxlsx::write.xlsx(df, file = xlsx_name)
    message("XLSX saved to: ", xlsx_name)
  }
  
  return(df)
}

