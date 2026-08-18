# General-purpose helpers
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.

# Guard for a package in Suggests. `purpose` completes "... is required to <purpose>".
.require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required to ", purpose, ". Install it with ",
         "install.packages('", pkg, "') or BiocManager::install('", pkg, "').",
         call. = FALSE)
  }
}


#' Extract the centered substring of a given length from an input string
#'
#' This function extracts a substring of length `n` from the center of the input string. If `n` is greater than the length of the input string,
#' the entire input string is returned.
#'
#' @param input_string A character string from which the centered substring will be extracted.
#' @param n An integer specifying the length of the substring to be extracted.
#' @return A character string representing the centered substring of the input string.
#' @export
#' @examples
#' center_substring("hello", 3) # returns "ell"
#' center_substring("world", 2) # returns "or"
#' center_substring("example", 7) # returns "example"
center_substring <- function(input_string, n) {
  string_length <- nchar(input_string)
  
  start_pos <- ifelse(n > string_length, 1, ceiling((string_length - n) / 2) + 1)
  end_pos <- start_pos + n - 1
  
  return(substr(input_string, start_pos, end_pos))
}


#' Handling of columns with multiple entries
#'
#' Split columns with multiple entries (need to be separated by a semicolon) into multiple rows with one entry each.
#' @param table dataframe
#' @param colname column that needs to be split (e.g., gene_names, protein_ids)
#' @param keep_all if TRUE, all entries are retained in additional rows. If FALSE,
#' only the first name/id is kept.
#' @param sep separator
#' @importFrom splitstackshape cSplit
#' @return table with split entries
split_genes <- function (table, colname = "gene_names", keep_all = FALSE, sep = ";")
{
  if (keep_all) {
    return(splitstackshape::cSplit(table, colname, direction = "long", sep = sep))
  }
  else {
    x <- as.character(table[[colname]])
    first <- vapply(
      strsplit(x, sep, fixed = TRUE),
      function(z) if (length(z)) z[[1]] else NA_character_,
      character(1)
    )
    first[is.na(x)] <- NA_character_
    table[[colname]] <- first
    return(table)
  }
}


#############################################################################
#' Clean String Values Column-Wise
#'
#' This function cleans string values by:
#' - Converting to lowercase.
#' - Replacing all whitespace with underscores.
#' - Replacing non-alphanumeric characters (except underscores) with underscores.
#' - Collapsing multiple underscores into a single underscore.
#' - Removing leading and trailing underscores.
#' - Replacing the micro sign (Unicode U+00B5) with \code{u}, and notifying the user if replacements occur.
#'
#' @param x A character vector to be cleaned.
#'
#' @return A character vector with cleaned values.
#'
#' @importFrom stringr str_detect str_replace_all str_remove_all str_trim
#' @importFrom dplyr %>%
#'
#' @examples
#' library(dplyr)
#' library(stringr)
#'
#' data <- data.frame(
#'   raw_column = c("10\u00b5g", "\u00b5mol/L", "Hello \u00b5World")
#' )
#'
#' data_clean <- data %>%
#'   mutate(cleaned_column = clean_values_colwise(raw_column))
#'
#' @export
clean_values_colwise <- function(x) {
  # Check if \u00b5 (micro) exists in the input
  if (any(stringr::str_detect(x, "\u00b5"))) {
    message("Note: Replaced '\u00b5' with 'u' in the input.")
  }
  
  x %>%
    stringr::str_replace_all("\u00b5", "u") %>%  # Replace 'u' with 'u'
    stringr::str_replace_all("\\s+", "_") %>%   # Replace all whitespace with "_"
    stringr::str_replace_all("[^a-z0-9_]", "_") %>% # Replace non-alphanumeric characters with "_"
    stringr::str_replace_all("_+", "_") %>%     # Replace multiple "_" with a single "_"
    stringr::str_remove_all("^_|_$")             # Remove leading and trailing "_"
}


#####################
#' Compare two data frames and report differences column by column
#'
#' Produces a structured comparison report between two data frames, including
#' dimension checks, column name alignment, per-column match statistics, and a
#' detailed diff via \code{waldo::compare()}. Only columns present in both
#' data frames are compared.
#'
#' @param df1 A data frame. Used as the reference (left-hand side) frame.
#' @param df2 A data frame. Compared against \code{df1} (right-hand side).
#' @param max_diffs Integer or \code{Inf}. Maximum number of differences shown
#'   by \code{waldo::compare()}. Defaults to \code{Inf} (show all).
#'
#' @return A data frame (invisibly) with one row per common column and the
#'   following columns:
#'   \describe{
#'     \item{column}{Column name.}
#'     \item{all_match}{Logical. \code{TRUE} if \code{identical(df1[[col]], df2[[col]])},
#'       meaning all values \emph{and} NA positions are identical.}
#'     \item{matching_rows}{Integer. Number of rows where values are equal,
#'       counting NA-vs-NA as a match.}
#'     \item{na_in_df1}{Integer. Number of \code{NA} values in \code{df1[[col]]}.}
#'     \item{na_in_df2}{Integer. Number of \code{NA} values in \code{df2[[col]]}.}
#'     \item{na_both}{Integer. Number of rows where both \code{df1[[col]]} and
#'       \code{df2[[col]]} are \code{NA} (these are counted as matches).}
#'     \item{total_rows}{Integer. Total number of rows (from \code{df1}).}
#'     \item{match_pct}{Numeric. Percentage of matching rows (including NA-vs-NA)
#'       out of \code{total_rows}, rounded to one decimal. Meaningful only when
#'       \code{nrow(df1) == nrow(df2)}.}
#'   }
#'
#' @details
#' The function prints a human-readable report to the console via \code{cat()}
#' and \code{print()}, covering:
#' \itemize{
#'   \item Row and column dimensions of each data frame.
#'   \item Whether column names are identical.
#'   \item A warning if the two frames have different columns, listing
#'     columns exclusive to each.
#'   \item A per-column summary table (see \strong{Value}).
#'   \item Lists of fully matched and mismatched columns.
#'   \item A detailed diff from \code{\link[waldo]{compare}}.
#' }
#'
#' Note: \code{matching_rows} and \code{match_pct} treat NA-vs-NA as a match,
#' consistent with \code{\link{identical}}. Use \code{na_both} to see how many
#' of those matches are NA pairs rather than equal values.
#'
#' @seealso \code{\link[waldo]{compare}} for the underlying diff engine.
#'
#' @examples
#' df1 <- data.frame(id = 1:3, val = c("a", NA, "c"))
#' df2 <- data.frame(id = 1:3, val = c("a", NA, "c"))
#' result <- compare_dataframes_full(df1, df2)
#'
#' @export
compare_dataframes_full <- function(df1, df2, max_diffs = Inf) {
  .require_pkg("waldo", "show the detailed diff")

  cat("=== DATAFRAME COMPARISON REPORT ===\n\n")
  
  # Dimensions
  cat("Dimensions:\n")
  cat("  df1:", nrow(df1), "rows,", ncol(df1), "columns\n")
  cat("  df2:", nrow(df2), "rows,", ncol(df2), "columns\n\n")
  
  # Column names
  cat("Column names match:", identical(names(df1), names(df2)), "\n\n")
  
  # Find common columns only
  common_cols <- intersect(names(df1), names(df2))
  
  if (length(common_cols) < length(names(df1))) {
    cat("  Only comparing", length(common_cols), "common columns\n")
    cat("   Columns only in df1:", setdiff(names(df1), common_cols), "\n")
    cat("   Columns only in df2:", setdiff(names(df2), common_cols), "\n\n")
  }
  
  # Per-column summary
  cat("=== COLUMN-BY-COLUMN ===\n\n")
  
  summary_df <- data.frame(
    column = common_cols,
    all_match = sapply(common_cols, function(col) {
      identical(df1[[col]], df2[[col]])
    }),
    matching_rows = sapply(common_cols, function(col) {
      x <- df1[[col]]; y <- df2[[col]]
      sum((x == y & !is.na(x) & !is.na(y)) | (is.na(x) & is.na(y)))
    }),
    na_in_df1 = sapply(common_cols, function(col) {
      sum(is.na(df1[[col]]))
    }),
    na_in_df2 = sapply(common_cols, function(col) {
      sum(is.na(df2[[col]]))
    }),
    na_both = sapply(common_cols, function(col) {
      sum(is.na(df1[[col]]) & is.na(df2[[col]]))
    }),
    total_rows = nrow(df1),
    match_pct = sapply(common_cols, function(col) {
      x <- df1[[col]]; y <- df2[[col]]
      n_match <- sum((x == y & !is.na(x) & !is.na(y)) | (is.na(x) & is.na(y)))
      round(100 * n_match / nrow(df1), 1)
    }),
    stringsAsFactors = FALSE
  )
  
  print(summary_df)
  
  # Highlight fully matching columns
  cat("\n=== FULLY MATCHED COLUMNS (100%) ===\n")
  matched_cols <- summary_df$column[summary_df$all_match]
  if (length(matched_cols) > 0) {
    cat("", paste(matched_cols, collapse = ", "), "\n")
  } else {
    cat("None\n")
  }
  
  # Highlight mismatched columns
  cat("\n=== COLUMNS WITH DIFFERENCES ===\n")
  mismatched_cols <- summary_df$column[!summary_df$all_match]
  if (length(mismatched_cols) > 0) {
    cat(" ", paste(mismatched_cols, collapse = ", "), "\n")
  } else {
    cat("None - all common columns are identical!\n")
  }
  
  # Show waldo diff for details
  cat("\n=== DETAILED DIFFERENCES (waldo) ===\n\n")
  print(waldo::compare(df1[common_cols], df2[common_cols], max_diffs = max_diffs))
  
  invisible(summary_df)
}



#' Save R Session Info and Package Usage Report
#'
#' @param root Character. Root directory to scan for package dependencies.
#'   Defaults to current working directory.
#' @param output_subdir Character. Subdirectory under root for output files.
#'   Default "r_version_updates".
#' @param overwrite Logical. Overwrite existing files? Default TRUE.
#'
#' @return Invisibly returns a list of output file paths.
#' @export
#' @examples
#' \dontrun{
#' save_session_report()
#' save_session_report(root = "C:/projects/my_project/")
#' }

save_session_report <- function(
    root         = getwd(),
    output_subdir = "r_session_info",
    overwrite    = TRUE
) {
  
  for (pkg in c("openxlsx", "renv")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required by save_session_report(); please install it.")
    }
  }

  # -- Normalise root --------------------------------------------------------
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  
  # -- Version strings -------------------------------------------------------
  r_version       <- R.version.string
  r_ver <- paste0("R_V", R.version$major, "_", gsub("\\.", "_", R.version$minor), "_")
  rstudio_version <- tryCatch(
    as.character(rstudioapi::versionInfo()$version),
    error = function(e) "RStudio not available"
  )
  library_paths <- paste(.libPaths(), collapse = ", ")
  
  # -- Output paths ----------------------------------------------------------
  out_dir <- file.path(root, output_subdir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  stem                  <- file.path(out_dir, paste0(format(Sys.time(), "%Y%m%d_"), r_ver))
  filename_session_xlsx <- paste0(stem, "session_infos.xlsx")
  filename_all_packs    <- paste0(stem, "all_packages.rda")
  filename_used_packs   <- paste0(stem, "used_packages.rda")
  
  # -- All installed packages ------------------------------------------------
  all_installed_packages <- as.data.frame(installed.packages())
  save(all_installed_packages, file = filename_all_packs)
  message("Saved: ", filename_all_packs)
  
  # -- Packages called in root -----------------------------------------------
  deps <- renv::dependencies(root, progress = FALSE)
  used_packages <- data.frame(
    Package = sort(unique(deps$Package))
  )
  # flag which are actually installed
  used_packages$Installed <- used_packages$Package %in% rownames(installed.packages())
  
  save(used_packages, file = filename_used_packs)
  message("Saved: ", filename_used_packs)
  
  # -- Settings record -------------------------------------------------------
  settings_used <- data.frame(
    Parameter = c("root", "output_subdir", "Run timestamp"),
    Value     = c(root, output_subdir, format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  
  # -- Excel workbook --------------------------------------------------------
  wb <- openxlsx::createWorkbook()
  
  ## Sheet 1 - session info + all installed packages
  openxlsx::addWorksheet(wb, "Session_info_all_packs")
  session_info <- data.frame(
    Field = c("R version", "RStudio version", "Library paths"),
    Value = c(r_version, rstudio_version, library_paths)
  )
  openxlsx::writeData(wb, "Session_info_all_packs", session_info, startRow = 1, startCol = 1)
  
  start_row <- nrow(session_info) + 4
  openxlsx::writeData(wb, "Session_info_all_packs", all_installed_packages,
            startRow = start_row, startCol = 1, withFilter = TRUE)
  openxlsx::freezePane(wb, "Session_info_all_packs", firstActiveRow = start_row + 1)
  openxlsx::setColWidths(wb, "Session_info_all_packs",
               cols = 1:ncol(all_installed_packages), widths = 20)
  
  ## Sheet 2 - packages used in root
  openxlsx::addWorksheet(wb, "Used_packages")
  openxlsx::writeData(wb, "Used_packages", used_packages,
            startRow = 1, startCol = 1, withFilter = TRUE)
  openxlsx::setColWidths(wb, "Used_packages", cols = 1, widths = 30)
  openxlsx::setColWidths(wb, "Used_packages", cols = 2, widths = 15)
  openxlsx::freezePane(wb, "Used_packages", firstActiveRow = 2)
  
  ## Sheet 3 - settings used to generate this report
  openxlsx::addWorksheet(wb, "Report_settings")
  openxlsx::writeData(wb, "Report_settings", settings_used, startRow = 1, startCol = 1)
  openxlsx::setColWidths(wb, "Report_settings", cols = 1, widths = 25)
  openxlsx::setColWidths(wb, "Report_settings", cols = 2, widths = 60)
  
  openxlsx::saveWorkbook(wb, filename_session_xlsx, overwrite = overwrite)
  message("Saved: ", filename_session_xlsx)
  
  invisible(list(
    xlsx          = filename_session_xlsx,
    all_packages  = filename_all_packs,
    used_packages = filename_used_packs
  ))
}



#' Extract values by protein ID position within semicolon-separated columns
#'
#' For each row, finds the position of \code{protein_id} within the
#' semicolon-separated \code{protein_accessions} column and extracts the
#' value at that same position from every other character column that shares
#' the same split-length structure. This converts group-level (multi-entry)
#' columns into single-entry columns aligned to the focal protein.
#'
#' @param df A data frame containing at least the columns specified by
#'   \code{protein_col} and \code{match_col}.
#' @param protein_col Character scalar. Name of the column holding the single
#'   protein identifier to locate (default: \code{"protein_id"}).
#' @param match_col Character scalar. Name of the column whose
#'   semicolon-separated entries are searched for \code{protein_col}
#'   (default: \code{"protein_accessions"}).
#' @param sep Character scalar. Separator used to split multi-value columns
#'   (default: \code{";"}).
#'
#' @return A data frame with the same dimensions as \code{df}. All character
#'   columns whose entry counts (after splitting by \code{sep}) match those of
#'   \code{match_col} are updated so that each cell contains only the value
#'   at the position of the focal protein. Rows where the protein ID cannot
#'   be found in \code{match_col} receive \code{NA} in all updated columns;
#'   a message reports the number of such unmatched rows.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Splits \code{match_col} on \code{sep} to obtain per-row entry
#'         vectors and records their lengths.
#'   \item Uses \code{match()} to find the 1-based position of each
#'         \code{protein_col} value within the corresponding entry vector.
#'   \item Identifies candidate columns (all columns except \code{protein_col})
#'         that are (a) of type character and (b) have the same per-row
#'         split-length as \code{match_col}.
#'   \item Replaces each such column's cell with the value at the matched
#'         position (or \code{NA_character_} when no match or position is
#'         out of bounds).
#' }
#'
#' @importFrom stringr str_split fixed
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   protein_id          = c("P001", "P003"),
#'   protein_accessions  = c("P001;P002", "P003;P004"),
#'   genes               = c("GeneA;GeneB", "GeneC;GeneD"),
#'   stringsAsFactors    = FALSE
#' )
#' extract_by_protein_id(df)
#' }
#'
#' @export
extract_by_protein_id <- function(df,
                                  protein_col = "protein_id",
                                  match_col = "protein_accessions",
                                  sep = ";") {
  
  # split the matching column (e.g. protein_accessions)
  match_entries <- stringr::str_split(df[[match_col]], stringr::fixed(sep))
  match_lengths <- lengths(match_entries)
  
  # find the position of protein_id within protein_accessions
  positions <- mapply(function(entries, pid) {
    match(pid, entries)
  }, match_entries, df[[protein_col]], USE.NAMES = FALSE)
  
  message("NA positions (unmatched rows): ", sum(is.na(positions)), " / ", nrow(df))
  
  # candidate columns: all except the single protein_id column
  candidate_cols <- names(df)[!names(df) %in% protein_col]
  
  # only update character columns with the same split structure as match_col
  cols_to_update <- Filter(function(col) {
    x <- df[[col]]
    
    if (!is.character(x)) return(FALSE)
    
    col_entries <- stringr::str_split(x, stringr::fixed(sep))
    col_lengths <- lengths(col_entries)
    
    all(col_lengths == match_lengths)
  }, candidate_cols)
  
  message("Columns to update: ", paste(cols_to_update, collapse = ", "))
  
  result <- df
  
  for (col in cols_to_update) {
    col_entries <- stringr::str_split(df[[col]], stringr::fixed(sep))
    
    result[[col]] <- mapply(function(entries, pos) {
      if (is.na(pos) || pos > length(entries)) NA_character_ else entries[pos]
    }, col_entries, positions, USE.NAMES = FALSE)
  }
  
  return(result)
}

