# Export results to files / external tools
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.


#' Write phosphoproteomics data to a file
#'
#' This function takes a SummarizedExperiment object, extracts the relevant
#' data, and writes it to a file.
#'
#' @param se A SummarizedExperiment object containing phosphoproteomics data.
#' @param file The filename where the resulting data should be saved (default: current directory).
#' @return relevant data
#' @importFrom DEP2 get_df_wide
#' @import dplyr
#' @export
write_phos <- function(se, file = ""){
  .assert_se(se)

  exp <- se %>% 
    get_df_wide() %>%
    select(name,
           gene_names,
           protein,
           protein_names,
           amino_acid, 
           position,
           multiplicity,
           sequence_window,
           ends_with(c("_diff","_p.val","p.adj","significant")),
           2:(ncol(se))) 
  
  colnames(exp) <- gsub("_diff", "_log2FC", colnames(exp))
  
  if (!file == "") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required to write .xlsx files; install it or leave 'file' empty.")
    }
    openxlsx::write.xlsx(exp, file = file)
  }
  
  return(exp)
  
}


#' Write proteomics data to a file
#'
#' This function takes a SummarizedExperiment object, extracts the relevant
#' data, and writes it to a file.
#'
#' @param se A SummarizedExperiment object containing proteomics data.
#' @param file The filename where the resulting data should be saved (default: current directory).
#' @return relevant data
#' @importFrom DEP2 get_df_wide
#' @import dplyr
#' @export
write_prot <- function(se, file = ""){
  .assert_se(se)

  exp <- se %>% 
    get_df_wide() %>%
    select(name,
           gene_names,
           protein_ids,
           protein_descriptions,
           orig_prot_ids,
           orig_gene_names,
           ends_with(c("_diff","_p.val","p.adj","significant")),
           2:(ncol(se))) 
  
  colnames(exp) <- gsub("_diff", "_log2FC", colnames(exp))
  
  if (!file == "") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required to write .xlsx files; install it or leave 'file' empty.")
    }
    openxlsx::write.xlsx(exp, file = file)
  }
  
  return(exp)
  
}


########################################################################
########################################################################
########################################################################
#' Prepare iPATH3 annotation tables from proteomics results
#'
#' Converts a `SummarizedExperiment` or wide result data frame into one or more
#' iPATH3-ready annotation tables, with line color representing fold change and
#' line width representing statistical significance.
#'
#' This function is designed for DEP2-style proteomics output where each
#' contrast has columns such as `"<contrast>_diff"` together with either
#' `"<contrast>_p.val"` or `"<contrast>_p.adj"`.
#'
#' @param x A `SummarizedExperiment` or a `data.frame`. If a
#'   `SummarizedExperiment` is supplied, it is converted to a wide data frame via
#'   `DEP2::get_df_wide()`.
#' @param contrasts Character vector of contrast names to use, or `"all"` to
#'   use all available contrasts detected from the fold-change columns.
#' @param pval_or_qval Which significance column type to use for filtering and
#'   width scaling. One of `"pval"` or `"qval"`.
#' @param cut_off Numeric scalar. Significance threshold applied to the selected
#'   p-value or adjusted p-value column.
#' @param color_mode Character scalar. One of `"binary"` or `"continuous"`.
#'   In binary mode, significant positive and negative fold changes are assigned
#'   fixed high/low colors and non-significant entries are grey. In continuous
#'   mode, only significant entries are colored by fold change and
#'   non-significant entries are grey.
#' @param range_mode Character scalar. One of `"local"` or `"total"`. In local
#'   mode, the color scale is computed separately for each contrast. In total
#'   mode, one global symmetric fold-change range is used across all selected
#'   contrasts.
#' @param width_range Numeric vector of length 2 giving the minimum and maximum
#'   iPATH3 line widths.
#' @param clip_fc Optional numeric scalar. If supplied, fold changes are clipped
#'   to `[-clip_fc, clip_fc]` for color scaling. Values beyond this range are
#'   shown with the most extreme colors.
#' @param color_pal Character vector of colors used for the fold-change
#'   gradient. Typically low, midpoint, and high.
#' @param color_res Integer scalar giving the number of colors used to build the
#'   gradient.
#' @param id_col Character scalar. Column containing UniProt identifiers.
#'   Default is `"uni_prot_ids"`.
#' @param id_mode How to handle multiple UniProt IDs within one row. One of:
#'   \itemize{
#'     \item `"first"`: keep only the first UniProt ID
#'     \item `"all"`: split all semicolon-separated UniProt IDs into separate
#'       iPATH3 rows
#'   }
#' @param collapse_duplicates How to handle duplicated `ipath3_id` values after
#'   ID expansion. One of:
#'   \itemize{
#'     \item `"none"`: keep all rows
#'     \item `"first"`: keep the first row per `ipath3_id`
#'     \item `"max_width"`: keep the row with the largest line width per
#'       `ipath3_id`
#'   }
#' @param diff_suffix Character scalar giving the suffix for fold-change
#'   columns. Default is `"_diff"`.
#' @param pval_suffix Character scalar giving the suffix for p-value columns.
#'   Default is `"_p.val"`.
#' @param qval_suffix Character scalar giving the suffix for adjusted p-value
#'   columns. Default is `"_p.adj"`.
#' @param open_in_browser Logical; if `TRUE`, open iPATH3 in the browser.
#' @param ipath_url Character scalar. URL to open when `open_in_browser = TRUE`.
#' @param plot_scalebars Logical; if `TRUE`, print one fold-change scalebar per
#'   contrast.
#' @param neutral_color Character scalar used for non-significant entries.
#' @param copy_to_clipboard Logical; if `TRUE`, copy iPATH3 output to the system
#'   clipboard using `clipr::write_clip()`.
#' @param clipboard_mode How to copy output when `copy_to_clipboard = TRUE`.
#'   One of:
#'   \itemize{
#'     \item `"combined"`: combine all selected contrasts into one data frame and
#'       add a `contrast` column
#'     \item `"first"`: copy only the first contrast
#'   }
#'
#' @details
#' For each selected contrast, the function:
#' \enumerate{
#'   \item identifies the fold-change and significance columns
#'   \item assigns colors according to `color_mode`
#'   \item assigns line widths based on significance
#'   \item extracts UniProt IDs from `id_col`
#'   \item optionally collapses duplicated `ipath3_id` values
#'   \item creates an iPATH3-ready table with columns `ipath3_id`, `color`,
#'     and `width`
#' }
#'
#' In `id_mode = "first"`, only the first UniProt ID per row is used.
#' In `id_mode = "all"`, all semicolon-separated UniProt IDs are expanded into
#' separate rows.
#'
#' If `copy_to_clipboard = TRUE` and more than one contrast is selected,
#' `clipboard_mode = "combined"` copies a single combined table with a
#' `contrast` column, which is usually the most practical clipboard format.
#'
#' The returned object is a named list with one entry per contrast. Each entry
#' contains:
#' \itemize{
#'   \item `df`: iPATH3-ready annotation data frame
#'   \item `scalebar`: data frame describing the fold-change color scale
#' }
#'
#' @return A named list, one element per contrast.
#'
#' @examples
#' \dontrun{
#' res_ipath <- prep_ipath3(
#'   x = se,
#'   contrasts = "all",
#'   pval_or_qval = "pval",
#'   cut_off = 0.01,
#'   color_mode = "continuous",
#'   range_mode = "local",
#'   width_range = c(10, 20),
#'   clip_fc = 0.58,
#'   id_mode = "all",
#'   collapse_duplicates = "max_width"
#' )
#'
#' prep_ipath3(
#'   x = se,
#'   contrasts = c("tnfa50_48h_vs_ctrl_48h", "il1b_vs_ctrl"),
#'   copy_to_clipboard = TRUE,
#'   clipboard_mode = "combined"
#' )
#' }
#'
#' @export
prep_ipath3 <- function(
    x,
    contrasts = "all",
    pval_or_qval = c("pval", "qval"),
    cut_off = 0.05,
    color_mode = c("binary", "continuous"),
    range_mode = c("local", "total"),
    width_range = c(10, 20),
    clip_fc = NULL,
    color_pal = c("#3B82F6", "#eeeeee", "#EF4444"),
    color_res = 256,
    id_col = "uni_prot_ids",
    id_mode = c("first", "all"),
    collapse_duplicates = c("none", "first", "max_width"),
    diff_suffix = "_diff",
    pval_suffix = "_p.val",
    qval_suffix = "_p.adj",
    open_in_browser = TRUE,
    ipath_url = "https://pathways.embl.de/ipath3.cgi?map=metabolic",
    plot_scalebars = TRUE,
    neutral_color = "#eeeeee",
    copy_to_clipboard = FALSE,
    clipboard_mode = c("combined", "first")
) {
  pval_or_qval <- match.arg(pval_or_qval)
  color_mode <- match.arg(color_mode)
  range_mode <- match.arg(range_mode)
  id_mode <- match.arg(id_mode)
  collapse_duplicates <- match.arg(collapse_duplicates)
  clipboard_mode <- match.arg(clipboard_mode)
  
  if (!is.numeric(cut_off) || length(cut_off) != 1 || is.na(cut_off)) {
    stop("`cut_off` must be a single numeric value.")
  }
  
  if (!is.numeric(width_range) || length(width_range) != 2 || anyNA(width_range)) {
    stop("`width_range` must be a numeric vector of length 2.")
  }
  
  if (!is.null(clip_fc)) {
    if (!is.numeric(clip_fc) || length(clip_fc) != 1 || is.na(clip_fc) || clip_fc <= 0) {
      stop("`clip_fc` must be NULL or a single positive numeric value.")
    }
  }
  
  if (!is.numeric(color_res) || length(color_res) != 1 || color_res < 2) {
    stop("`color_res` must be a single integer-like value >= 2.")
  }
  
  if (!is.character(color_pal) || length(color_pal) < 2) {
    stop("`color_pal` must contain at least two colors.")
  }
  
  if (inherits(x, "SummarizedExperiment")) {
    df <- DEP2::get_df_wide(x)
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("`x` must be either a SummarizedExperiment or a data.frame.")
  }
  
  if (!is.data.frame(df)) {
    stop("Input could not be converted to a data.frame.")
  }
  
  if (!id_col %in% names(df)) {
    stop("Identifier column not found: ", id_col)
  }
  
  diff_cols <- names(df)[endsWith(names(df), diff_suffix)]
  if (length(diff_cols) == 0) {
    stop("No fold-change columns found ending in '", diff_suffix, "'.")
  }
  
  available_contrasts <- sub(
    paste0(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", diff_suffix), "$"),
    "",
    diff_cols,
    perl = TRUE
  )
  
  if (length(contrasts) == 1 && identical(contrasts, "all")) {
    contrasts <- available_contrasts
  } else {
    if (!is.character(contrasts) || length(contrasts) < 1) {
      stop("`contrasts` must be 'all' or a character vector of contrast names.")
    }
    
    missing_contrasts <- setdiff(contrasts, available_contrasts)
    if (length(missing_contrasts) > 0) {
      stop(
        "The following contrasts were not found: ",
        paste(missing_contrasts, collapse = ", ")
      )
    }
  }
  
  sig_suffix <- if (identical(pval_or_qval, "pval")) pval_suffix else qval_suffix
  
  required_sig_cols <- paste0(contrasts, sig_suffix)
  missing_sig_cols <- setdiff(required_sig_cols, names(df))
  if (length(missing_sig_cols) > 0) {
    stop(
      "The following significance columns were not found: ",
      paste(missing_sig_cols, collapse = ", ")
    )
  }
  
  low_color <- color_pal[1]
  high_color <- color_pal[length(color_pal)]
  
  if (is.null(clip_fc)) {
    message("Scalebar values cover the full observed fold-change range (no clipping).")
  } else {
    message(
      "Scalebar values are clipped at +/-", clip_fc, ". ",
      "This improves visibility of small changes, but values beyond this range ",
      "share the same extreme colors."
    )
  }
  
  clip_values <- function(x, clip_fc = NULL) {
    if (is.null(clip_fc)) {
      return(x)
    }
    pmax(pmin(x, clip_fc), -clip_fc)
  }
  
  make_color_fun <- function(range_use, color_pal, color_res) {
    pal <- grDevices::colorRampPalette(color_pal)(color_res)
    scale_vals <- seq(range_use[1], range_use[2], length.out = color_res)
    
    function(x) {
      x_clipped <- pmax(pmin(x, max(scale_vals)), min(scale_vals))
      pal[findInterval(x_clipped, scale_vals, all.inside = TRUE)]
    }
  }
  
  collapse_ipath_duplicates <- function(df_ipath, method) {
    if (identical(method, "none")) {
      return(df_ipath)
    }
    
    if (identical(method, "first")) {
      return(df_ipath[!duplicated(df_ipath$ipath3_id), , drop = FALSE])
    }
    
    if (identical(method, "max_width")) {
      width_num <- suppressWarnings(as.numeric(sub("^W", "", df_ipath$width)))
      df_ipath$.width_num <- width_num
      df_ipath <- df_ipath[order(df_ipath$ipath3_id, -df_ipath$.width_num), , drop = FALSE]
      df_ipath <- df_ipath[!duplicated(df_ipath$ipath3_id), , drop = FALSE]
      df_ipath$.width_num <- NULL
      return(df_ipath)
    }
    
    df_ipath
  }
  
  build_ipath_df <- function(df_comp, id_col, id_mode, collapse_duplicates) {
    base_df <- df_comp %>%
      dplyr::filter(!is.na(.data[[id_col]]), .data[[id_col]] != "")
    
    if (identical(id_mode, "first")) {
      df_ipath <- base_df %>%
        dplyr::mutate(
          uniprot_id = sub(";.*", "", .data[[id_col]]),
          ipath3_id = paste0("UNIPROT:", .data$uniprot_id),
          width = dplyr::if_else(
            is.na(.data$line_width),
            NA_character_,
            paste0("W", round(.data$line_width))
          )
        ) %>%
        dplyr::select(.data$ipath3_id, .data$color, .data$width)
    } else {
      df_ipath <- base_df %>%
        dplyr::mutate(uniprot_id = stringr::str_split(.data[[id_col]], ";")) %>%
        tidyr::unnest(.data$uniprot_id) %>%
        dplyr::mutate(
          uniprot_id = trimws(.data$uniprot_id),
          ipath3_id = paste0("UNIPROT:", .data$uniprot_id),
          width = dplyr::if_else(
            is.na(.data$line_width),
            NA_character_,
            paste0("W", round(.data$line_width))
          )
        ) %>%
        dplyr::filter(!is.na(.data$uniprot_id), .data$uniprot_id != "") %>%
        dplyr::select(.data$ipath3_id, .data$color, .data$width)
    }
    
    collapse_ipath_duplicates(df_ipath, method = collapse_duplicates)
  }
  
  global_range <- NULL
  if (identical(range_mode, "total")) {
    all_diffs <- unlist(lapply(contrasts, function(con) {
      df[[paste0(con, diff_suffix)]]
    }), use.names = FALSE)
    
    all_diffs <- all_diffs[!is.na(all_diffs)]
    if (length(all_diffs) == 0) {
      stop("No non-missing fold-change values found across selected contrasts.")
    }
    
    absmax <- max(abs(all_diffs), na.rm = TRUE)
    if (!is.null(clip_fc)) {
      absmax <- min(absmax, clip_fc)
    }
    if (absmax == 0) {
      absmax <- 1
    }
    
    global_range <- c(-absmax, absmax)
  }
  
  results_list <- vector("list", length(contrasts))
  names(results_list) <- contrasts
  
  for (contrast in contrasts) {
    diff_col <- paste0(contrast, diff_suffix)
    sig_col <- paste0(contrast, sig_suffix)
    
    if (identical(range_mode, "local")) {
      diffs <- df[[diff_col]]
      diffs <- diffs[!is.na(diffs)]
      
      if (length(diffs) == 0) {
        absmax <- 1
      } else {
        absmax <- max(abs(diffs), na.rm = TRUE)
        if (!is.null(clip_fc)) {
          absmax <- min(absmax, clip_fc)
        }
        if (absmax == 0) {
          absmax <- 1
        }
      }
      
      range_use <- c(-absmax, absmax)
    } else {
      range_use <- global_range
    }
    
    diff_to_color <- make_color_fun(
      range_use = range_use,
      color_pal = color_pal,
      color_res = color_res
    )
    
    scale_vals <- seq(range_use[1], range_use[2], length.out = color_res)
    scalebar <- data.frame(
      value = scale_vals,
      color = diff_to_color(scale_vals),
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(plot_scalebars)) {
      plot_scalebar <- function(scalebar, title = "log2FC") {
        ggplot2::ggplot(scalebar, ggplot2::aes(x = .data$value, y = 1, fill = .data$color)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_identity() +
          ggplot2::theme_void() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(hjust = 0.5)
          ) +
          ggplot2::scale_x_continuous(
            breaks = pretty(scalebar$value, n = 5),
            labels = pretty(scalebar$value, n = 5)
          ) +
          ggplot2::labs(x = title, y = NULL, title = title)
      }
      
      print(plot_scalebar(scalebar, title = paste0(contrast, " log2FC")))
    }
    
    df_comp <- df %>%
      dplyr::mutate(
        .diff = .data[[diff_col]],
        .sig = .data[[sig_col]],
        .diff_clipped = clip_values(.data$.diff, clip_fc = clip_fc),
        color = dplyr::case_when(
          color_mode == "binary" & !is.na(.data$.sig) & !is.na(.data$.diff) &
            .data$.sig <= cut_off & .data$.diff > 0 ~ high_color,
          color_mode == "binary" & !is.na(.data$.sig) & !is.na(.data$.diff) &
            .data$.sig <= cut_off & .data$.diff < 0 ~ low_color,
          color_mode == "binary" ~ neutral_color,
          color_mode == "continuous" & (is.na(.data$.sig) | .data$.sig > cut_off) ~ neutral_color,
          color_mode == "continuous" & .data$.sig <= cut_off & !is.na(.data$.diff_clipped) ~
            diff_to_color(.data$.diff_clipped)
        )
      )
    
    sig_vals <- -log10(df_comp[[sig_col]])
    sig_vals[is.infinite(sig_vals)] <- NA_real_
    
    sig_range <- range(sig_vals, na.rm = TRUE)
    if (any(!is.finite(sig_range)) || diff(sig_range) == 0) {
      sig_range <- c(0, 1)
    }
    
    df_comp <- df_comp %>%
      dplyr::mutate(
        line_width = dplyr::case_when(
          is.na(.data$.sig) ~ width_range[1],
          .data$.sig > cut_off ~ width_range[1],
          TRUE ~ scales::rescale(
            -log10(.data$.sig),
            to = width_range,
            from = sig_range
          )
        )
      )
    
    df_sub <- build_ipath_df(
      df_comp = df_comp,
      id_col = id_col,
      id_mode = id_mode,
      collapse_duplicates = collapse_duplicates
    )
    
    results_list[[contrast]] <- list(
      df = df_sub,
      scalebar = scalebar
    )
  }
  
  if (isTRUE(copy_to_clipboard)) {
    if (!requireNamespace("clipr", quietly = TRUE)) {
      warning("`copy_to_clipboard = TRUE`, but the `clipr` package is not installed.")
    } else {
      if (length(contrasts) == 1 || identical(clipboard_mode, "first")) {
        clipr::write_clip(results_list[[contrasts[1]]]$df)
        message("Copied iPATH3 table for contrast '", contrasts[1], "' to clipboard.")
      } else {
        combined_df <- dplyr::bind_rows(
          lapply(names(results_list), function(con) {
            dplyr::mutate(results_list[[con]]$df, contrast = con, .before = 1)
          })
        )
        clipr::write_clip(combined_df)
        message("Copied combined iPATH3 table for ", length(contrasts), " contrasts to clipboard.")
      }
    }
  }
  
  if (isTRUE(open_in_browser)) {
    utils::browseURL(ipath_url)
    message("iPATH3 opened in browser.")
  }
  
  results_list
}

