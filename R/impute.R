# Missing-value imputation
#
# Part of the 'sev' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#'Impute as with perseus
#'
#'Imputes missing values in a summarized experiment similar to Perseus. Missing values are
#'replaced by drawing from a left shifted gaussian distribution that is calculated based on either
#'individual samples or all measurements.
#' @param se summarized experiment
#' @param width width of gaussian distribution
#' @param downshift shift of gaussian distribution
#' @param per_col if distribution should be calculated per column (sample)
#' or accross all samples
#' 
#' @return summarized experiment without missing values
#' @importFrom dplyr group_by summarize filter ungroup select mutate_all mutate
#' @importFrom BiocGenerics unique
#' @importFrom tidyr pivot_wider
#' @export
#'
impute_perseus = function(se, width = 0.3, downshift = 1.8, per_col=T) {


  # 1. transform to long and set lfq_imputed = TRUE for all missing values
  df <- se %>%
    get_df_long() %>%
    select(sample, condition, intensity, name) %>%
    ungroup()%>%
    mutate(lfq_imputed = ifelse(is.na(intensity), TRUE, FALSE))

  # 2. impute with the selected method
  if(per_col){
    cols<-split(df, df$sample)

    for (cols_ in names(cols)){
      set.seed(1)
      temp = cols[[cols_]]$intensity
      temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
      temp.mean = mean(temp, na.rm = TRUE) -
        downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
      n.missing = sum(is.na(temp))
      cols[[cols_]]$intensity[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    }

    df <- do.call(rbind, cols)

  }

  else{
    set.seed(1)
    temp = df$intensity
    temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
    temp.mean = mean(temp, na.rm = TRUE) -
      downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
    n.missing = sum(is.na(temp))
    df$lfq_intensity_norm[is.na(temp)] = stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  }

  # Calculate table with all values and another assay describing which values were imputed
  imp_val <- df  %>%
    select(sample, intensity, name) %>%
    pivot_wider(names_from = "sample", values_from = "intensity") %>%
    as.data.frame() %>%
    select(-name)

  imp_yn <- df  %>%
    select(sample, lfq_imputed, name) %>%
    pivot_wider(names_from = "sample", values_from = "lfq_imputed") %>%
    as.data.frame() %>%
    select(-name)

  colnames(imp_yn) <- paste0(se$sample,"_imputed")

  # Add as assay
  se <- add_assay(se, imp_val %>% as.matrix(), name = "imputed_perseus", withDimnames = FALSE)
  assays(se, withDimnames = FALSE)$imputed <- imp_yn %>% dplyr::mutate_all(~ ifelse(.x, 1,0)) %>% as.matrix()

  # Add rowData
  rowData(se) <- cbind(rowData(se), imp_yn)
  
  return(se)

}



#' Imputation via DEP2 while keeping original data
#' 
#' Utilizes the DEP2::impute() function to replace missing values. In addition, original raw data
#' is retained in an additional assay.
#'
#' @param se with missing data
#' @param fun Imputation function/method passed to \code{DEP2::impute()}; one of
#'   the methods in the default vector (e.g. \code{"bpca"}, \code{"knn"},
#'   \code{"QRILC"}, \code{"MinProb"}). The first element is used by default.
#' @param ... parameters for DEP2::impute()
#'
#' @return se with imputed data in the main assay and raw data in another assay
#' @export
#'
impute_DEP2 <- function(se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
                               "man", "min", "zero", "mixed", "nbavg"), ...){
  
  assays(se, withDimnames = FALSE)$imputed <- assay(se) %>% as.data.frame() %>% dplyr::mutate_all(~ ifelse(is.na(.x), 1,0)) %>% as.matrix()
  
  if(fun == "mixed"){
    
    temp <- DEP2::impute(se, fun = fun, randna = rowData(se)$randna, ...)
    
  }else{
    
    temp <- DEP2::impute(se, fun = fun, ...)
    
  }
  
  se <- add_assay(se, assay(temp), "imputed_DEP2")
  
  
  return(se)
  
}



##############
#' Condition-wise Imputation for Proteomics Data
#'
#' \code{condition_impute} imputes missing values in a proteomics dataset on a condition-wise basis.
#' It supports various imputation methods and allows filtering rows based on the percentage
#' of measured or missing values, with customizable logical operators.
#'
#' @param se SummarizedExperiment
#'   Proteomics data in a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param fun Character
#'   Imputation method to use. Available options include:
#'   \code{"QRILC"}, \code{"bpca"}, \code{"knn"}, \code{"MLE"}, \code{"MinDet"}, \code{"MinProb"},
#'   \code{"man"}, \code{"min"}, \code{"zero"}, \code{"mixed"}, \code{"nbavg"}, \code{"RF"}, \code{"GSimp"}.
#'   See \code{\link[MSnbase]{impute,MSnSet-method}}, \code{\link{manual_impute}}, and
#'   \code{\link[missForest]{missForest}} for details.
#' @param percent Numeric
#'   A value between 0 and 1 specifying the percentage threshold for filtering rows based
#'   on measured or missing values.
#' @param percent_type Character
#'   Specifies whether to filter rows based on the percentage of measured (\code{"measured"}) or
#'   missing (\code{"NA"}) values.
#' @param logical_operator Character
#'   Logical operator to apply for filtering rows. Options are: \code{">"}, \code{"<"},
#'   \code{">="}, \code{"<="}, \code{"=="}, \code{"!="}.
#' @param ... Additional arguments for specific imputation methods. For example:
#'   \itemize{
#'     \item \code{man}: \code{\link{manual_impute}} supports \code{shift} and \code{scale}.
#'     \item \code{RF}: \code{\link[missForest]{missForest}} supports \code{ntree} and \code{parallelize}.
#'     \item Other methods: See \code{\link[MSnbase]{impute,MSnSet-method}}.
#'   }
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with imputed values.
#' @examples
#' \dontrun{
#' # Example usage with SummarizedExperiment data
#' # Impute NA in rows where more than 30% of values are measured
#' imputed_se <- condition_impute(se, 
#'                                fun = "man", 
#'                                percent = 0.3, 
#'                                percent_type = "measured", 
#'                                logical_operator = ">",
#'                                shift = 1.8, scale = 0.3)
#'
#' # Impute NA in rows where more than 70% of values are missing
#' imputed_se <- condition_impute(se, 
#'                                fun = "knn", 
#'                                percent = 0.7, 
#'                                percent_type = "NA", 
#'                                logical_operator = ">",
#'                                k = 5, rowmax = 0.8)
#'
#' # Impute NA in rows where missing values are <=50%
#' imputed_se <- condition_impute(se, 
#'                                fun = "QRILC", 
#'                                percent = 0.5, 
#'                                percent_type = "NA", 
#'                                logical_operator = "<=")
#'
#' }
#' @import SummarizedExperiment
#' @import MSnbase
#' @import assertthat
#' @export
condition_impute <- function(se, 
                             fun = c("QRILC", "bpca", "knn", "MLE", "MinDet", 
                                     "MinProb", "man", "min", "zero", "mixed", "nbavg", "RF", "GSimp"), 
                             percent = 0.3, 
                             percent_type = c("measured", "NA"), 
                             logical_operator = ">", 
                             ...) {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  fun <- match.arg(fun)
  percent_type <- match.arg(percent_type)
  
  if (!any(is.na(assay(se)))) {
    warning("No missing values in the dataset. Returning the unchanged object.")
    return(se)
  }
  
  # Add NA tracking metadata
  SummarizedExperiment::rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  SummarizedExperiment::rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  
  # Group data by condition
  conditions <- SummarizedExperiment::colData(se)$condition
  unique_conditions <- unique(conditions)
  
  # Perform condition-wise imputation
  for (cond in unique_conditions) {
    message(paste("Processing condition:", cond))
    
    # Subset data for the current condition
    condition_indices <- which(conditions == cond)
    condition_data <- SummarizedExperiment::assay(se)[, condition_indices, drop = FALSE]
    
    # Calculate row percentages
    if (percent_type == "measured") {
      row_percent <- rowSums(!is.na(condition_data)) / ncol(condition_data)
    } else if (percent_type == "NA") {
      row_percent <- rowSums(is.na(condition_data)) / ncol(condition_data)
    }
    
    # Apply the logical operator to determine rows to impute
    rows_to_impute <- switch(
      logical_operator,
      ">" = which(row_percent > percent),
      "<" = which(row_percent < percent),
      ">=" = which(row_percent >= percent),
      "<=" = which(row_percent <= percent),
      "==" = which(row_percent == percent),
      "!=" = which(row_percent != percent),
      stop("Invalid logical_operator. Use one of: '>', '<', '>=', '<=', '==', '!='.")
    )
    
    if (length(rows_to_impute) > 0) {
      message(paste("Imputing", length(rows_to_impute), "rows in condition", cond))
      
      # Impute missing values for the selected rows
      condition_data[rows_to_impute, ] <- switch(
        fun,
        "man" = manual_impute_subset(condition_data[rows_to_impute, , drop = FALSE], ...),
        "RF" = {
          args <- list(...)
          if (!"ntree" %in% names(args)) args[["ntree"]] <- 60
          if (!"parallelize" %in% names(args)) args[["parallelize"]] <- "variables"
          args$xmis <- condition_data[rows_to_impute, ]
          doParallel::registerDoParallel(cores = min(6, parallel::detectCores() - 2))
          imputed_data <- do.call(missForest::missForest, args)$ximp
          doParallel::stopImplicitCluster()
          imputed_data
        },
        {
          # Use MSnbase imputation for other methods
          MSnSet_data <- MSnbase::MSnSet(exprs = condition_data[rows_to_impute, ])
          MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
          MSnbase::exprs(MSnSet_imputed)
        }
      )
      
      # Replace the imputed data back into the assay
      SummarizedExperiment::assay(se)[rows_to_impute, condition_indices] <- condition_data[rows_to_impute, ]
    } else {
      message(paste("No rows to impute in condition", cond, "based on the criteria."))
    }
  }
  
  return(se)
}


manual_impute_subset <- function(data, scale = 0.3, shift = 1.8) {
  # Perform manual imputation on a subset of the data
  assertthat::assert_that(is.matrix(data) || is.data.frame(data))
  
  for (col in seq_len(ncol(data))) {
    col_values <- data[, col]
    if (any(is.na(col_values))) {
      # Calculate column statistics
      col_mean <- mean(col_values, na.rm = TRUE)
      col_sd <- sd(col_values, na.rm = TRUE)
      col_median <- median(col_values, na.rm = TRUE)
      
      # Impute missing values
      na_indices <- which(is.na(col_values))
      data[na_indices, col] <- rnorm(
        length(na_indices),
        mean = col_median - shift * col_sd,
        sd = col_sd * scale
      )
    }
  }
  
  return(data)
}



#' Impute missing values with a reproducible, content-derived seed
#'
#' Generates a deterministic seed by hashing the input object, ensuring
#' identical objects always produce identical imputation results regardless
#' of processing order. The seed combines a user-supplied base seed with an
#' xxhash64 hash of the object, keeping results distinct across different
#' datasets while remaining fully reproducible.
#'
#' @param x A \code{SummarizedExperiment} object (or any R object) containing
#'   missing values to impute.
#' @param impute_fun Function. The imputation function to apply to \code{x}.
#'   Must accept \code{x} as its first argument.
#' @param base_seed Integer. Base value added to the object-derived hash seed.
#'   Default \code{42}. Change this to shift all seeds systematically (e.g.
#'   for sensitivity analysis).
#' @param ... Additional arguments passed to \code{impute_fun}.
#'
#' @return The result of \code{impute_fun(x, ...)}, with missing values imputed.
#'
#' @details
#' The seed is computed as:
#' \enumerate{
#'   \item Hash \code{x} using \code{\link[digest]{digest}} with
#'     \code{algo = "xxhash64"}.
#'   \item Convert the first 7 hex characters to an integer via
#'     \code{strtoi(..., base = 16)}.
#'   \item Add \code{base_seed} to ensure the result stays well within
#'     \code{.Machine$integer.max}.
#' }
#' \code{set.seed()} is called immediately before \code{impute_fun}, so any
#' random state set elsewhere in the session does not affect the result.
#'
#' @seealso \code{\link[digest]{digest}}
#'
#' @importFrom digest digest
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' library(imputeLCMD)
#'
#' # Using QRILC imputation
#' se_imputed <- impute_seeded(se, impute_fun = impute.QRILC)
#'
#' # With a different base seed for sensitivity analysis
#' se_imputed_alt <- impute_seeded(se, impute_fun = impute.QRILC, base_seed = 100)
#'
#' # Passing extra arguments to the imputation function
#' se_imputed <- impute_seeded(se, impute_fun = impute.QRILC, tune.sigma = 0.5)
#' }
#'
#' @export
impute_seeded <- function(x,
                          impute_fun,
                          base_seed = 42,
                          ...) {
  
  obj_hash <- digest::digest(x, algo = "xxhash64")
  obj_seed <- base_seed + strtoi(substr(obj_hash, 1, 7), base = 16L)
  
  set.seed(obj_seed)
  
  impute_fun(x, ...)
}


