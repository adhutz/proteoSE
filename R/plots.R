# General plotting (volcano, scatter, heatmap, correlation)
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.

# Smallest value in `x` closest to `a`. Replaces min(DescTools::Closest(x, a)),
# which cost a 13-package dependency for this one line.
.closest <- function(x, a) {
  d <- abs(x - a)
  min(x[d == min(d, na.rm = TRUE)], na.rm = TRUE)
}

# cbind data frames of unequal length, padding short ones with NA.
# Replaces gdata::cbindX; data-frame row indexing past the end yields NA rows.
.cbind_pad <- function(dfs) {
  n <- max(vapply(dfs, nrow, integer(1)))
  do.call(cbind, lapply(dfs, function(d) d[seq_len(n), , drop = FALSE]))
}


#' Volcano Plot for Proteomics and Phosphoproteomics
#'
#' Creates a volcano plot from a \code{SummarizedExperiment} object, visualising
#' differential abundance results for a given contrast. Supports both standard
#' proteomics (one row per protein) and phosphoproteomics (multiple rows per
#' gene, one per phosphosite) via the \code{unique_id_col} argument.
#'
#' @param se A \code{SummarizedExperiment} object whose \code{rowData} contains
#'   differential abundance results. Expected columns are derived from
#'   \code{contrast}: \code{<contrast>_p.val}, \code{<contrast>_p.adj},
#'   \code{<contrast>_diff}, and \code{<contrast>_significant}.
#' @param contrast Character scalar. Prefix used to identify the relevant
#'   result columns in \code{rowData(se)}, e.g. \code{"condA_vs_condB"}.
#' @param id_col Character scalar. Name of the \code{rowData} column containing
#'   gene or protein identifiers used to match \code{target_names}.
#'   Default: \code{"gene_names"}.
#' @param unique_id_col Character scalar or \code{NULL}. Name of a \code{rowData}
#'   column that uniquely identifies each row (e.g. a phosphosite ID such as
#'   \code{"phosphosite_id"}). When \code{NULL} (default), \code{id_col} is
#'   used for deduplication, which is appropriate for standard proteomics.
#'   For phosphoproteomics, set this to a site-level identifier so that all
#'   sites of the same gene are labelled independently.
#' @param label_col Character scalar or \code{NULL}. Name of the \code{rowData}
#'   column whose values are used as point labels in the plot. When \code{NULL}
#'   (default), the \code{"name"} column is used. For phosphoproteomics, set
#'   this to a column that includes site information (e.g.
#'   \code{"phosphosite_name"}) so labels read as \code{"EGFR_S1234"} rather
#'   than just \code{"EGFR"}.
#' @param target_names Character vector. Gene or protein identifiers (matched
#'   against \code{id_col}) to highlight as targets. All rows whose
#'   \code{id_col} value is in this vector are coloured and labelled as
#'   targets, including multiple phosphosites of the same gene.
#'   Default: \code{c("")} (no targets).
#' @param label_sign Logical. Whether to label significant non-target points.
#'   When \code{FALSE}, no significant points are labelled regardless of
#'   \code{top_n_sign}. When \code{TRUE} (default), the top \code{top_n_sign}
#'   significant points by score are labelled.
#' @param label_targets Logical. Whether to label target points. When
#'   \code{FALSE}, target points are still coloured but not labelled.
#'   Default: \code{TRUE}.
#' @param top_n_sign Integer. Number of top-scoring significant non-target
#'   points to label. The budget is split evenly between left (down-regulated)
#'   and right (up-regulated) when \code{direction = "both"}; ceiling goes to
#'   the left side. Default: \code{20}.
#' @param top_n_nonsign Integer. Number of top-scoring non-significant
#'   non-target points to label. Follows the same directional split logic as
#'   \code{top_n_sign}. Default: \code{0}.
#' @param direction Character scalar. Which side of the volcano to draw labels
#'   from. One of \code{"both"} (default), \code{"left"} (down-regulated
#'   only), or \code{"right"} (up-regulated only).
#' @param max.overlaps Integer. Passed to \code{\link[ggrepel]{geom_text_repel}}
#'   to control the maximum number of overlapping labels allowed before some
#'   are hidden. Default: \code{15}.
#' @param labelsize Numeric. Font size for point labels, in \code{mm} as
#'   interpreted by \code{ggrepel}. Default: \code{2}.
#' @param pvalorqval Character scalar. Whether to apply the significance cutoff
#'   to the raw p-value (\code{"pval"}) or the adjusted p-value / q-value
#'   (\code{"qval"}, default). Controls which column is compared against
#'   \code{pqcutoff}.
#' @param pqcutoff Numeric. Significance threshold applied to the column
#'   selected by \code{pvalorqval}. Points with a value \eqn{\le} this
#'   threshold AND an absolute fold change \eqn{>} \code{diffcutoff} are
#'   considered significant. Default: \code{0.05}.
#' @param diffcutoff Numeric. Absolute log2 fold-change threshold for
#'   significance. Default: \code{0.58} (equivalent to 1.5-fold change).
#' @param sig_from_column Logical. If \code{TRUE}, significance is read directly
#'   from a precomputed \code{<contrast>_significant} column in \code{rowData}
#'   (the behaviour of the former \code{se_volcano()}/\code{lars_volcano()}),
#'   and \code{pvalorqval}, \code{pqcutoff}, and \code{diffcutoff} are ignored.
#'   If \code{FALSE} (default), significance is recomputed from those cutoffs.
#' @param fill_colors Named character vector of length 3. Fill colours for the
#'   three point categories. Names must be \code{"Significant"},
#'   \code{"Target"}, and \code{"None"}.
#'   Default: \code{c(Significant = "#E69F00", Target = "#0072B2", None = "grey90")}.
#' @param shape_values Named numeric vector of length 2. Point shapes for
#'   non-significant (\code{FALSE}) and significant (\code{TRUE}) points.
#'   Values should be \code{pch} codes that support a fill aesthetic (21-25).
#'   Default: \code{c(`FALSE` = 21, `TRUE` = 24)}.
#' @param alpha_values Named numeric vector of length 2. Transparency values
#'   for background points (\code{"None"}) and points of interest (\code{"goi"}).
#'   Default: \code{c(None = 0.2, goi = 0.7)}.
#' @param size_values Named numeric vector of length 2. Point sizes for
#'   non-significant (\code{FALSE}) and significant (\code{TRUE}) points.
#'   Default: \code{c(`FALSE` = 1.5, `TRUE` = 2.5)}.
#' @param target_shape_override Integer or \code{NULL}. If not \code{NULL}, target
#'   points are re-drawn on top of all other layers using this \code{pch} shape
#'   (e.g. \code{23} for a diamond, \code{17} for a solid triangle), making
#'   them visually distinct regardless of their significance status.
#'   Default: \code{NULL} (no override; targets use the same shape scale as
#'   other points).
#' @param target_color Character scalar. Fill colour used for target points
#'   when \code{target_shape_override} is set. Default: \code{"#0072B2"}.
#' @param label_colors Named character vector of length 2. Colours for the two
#'   label types. Names must be \code{"Target"} and \code{"Top"}.
#'   Default: \code{c(Target = "#0072B2", Top = "black")}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{plot_out}}{A \code{ggplot2} object containing the volcano plot.}
#'   \item{\code{df_out}}{A data frame of all labelled points (targets and
#'     top-scoring non-targets), with all \code{rowData} columns plus the
#'     derived columns \code{target_gene}, \code{sign}, \code{target_sign},
#'     \code{alpha}, \code{score}, and \code{label_type}.}
#' }
#'
#' @details
#' ## Scoring and label selection
#' Each point is assigned a score of \eqn{-\log_{10}(p) \times |\Delta|}, where
#' \eqn{p} is the raw p-value and \eqn{\Delta} is the log2 fold change. This
#' balances statistical evidence with effect size. The top-\eqn{n} significant
#' and non-significant non-target points by score are selected for labelling
#' according to \code{top_n_sign}, \code{top_n_nonsign}, and \code{direction}.
#'
#' ## Phosphoproteomics support
#' When \code{unique_id_col} is supplied, deduplication in the label-combining
#' step uses that column instead of \code{id_col}. This ensures that multiple
#' phosphosites of the same gene are all retained as separate labelled points.
#' Pair with \code{label_col} to show site-aware labels (e.g.
#' \code{"EGFR_S1234"}).
#'
#' ## Rendering order
#' Points are drawn in the order: background (non-targets) first, then targets
#' on top. If \code{target_shape_override} is set, an additional layer re-draws
#' target points above everything else.
#'
#' ## Label rendering
#' All labels (targets and top-scoring non-targets) are passed to a single
#' \code{\link[ggrepel]{geom_text_repel}} call so the repulsion algorithm
#' treats them simultaneously and avoids inter-group overlaps.
#'
#' @seealso \code{\link[ggrepel]{geom_text_repel}},
#'   \code{\link[SummarizedExperiment]{rowData}}
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select all_of starts_with mutate case_when arrange filter slice_max bind_rows distinct sym
#' @importFrom rlang quo
#' @importFrom ggplot2 ggplot aes aes_string geom_point scale_fill_manual scale_shape_manual scale_alpha_manual scale_size_manual guides guide_legend labs theme_bw scale_color_manual
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' \dontrun{
#' # --- Standard proteomics ---
#' result <- plot_volcano(
#'   se       = my_se,
#'   contrast = "treated_vs_control",
#'   id_col   = "gene_names"
#' )
#' result$plot_out
#'
#' # --- Highlight specific targets ---
#' result <- plot_volcano(
#'   se           = my_se,
#'   contrast     = "treated_vs_control",
#'   target_names = c("EGFR", "TP53", "MYC"),
#'   top_n_sign   = 15
#' )
#' result$plot_out
#'
#' # --- Phosphoproteomics: all sites of target genes labelled ---
#' result <- plot_volcano(
#'   se             = my_phospho_se,
#'   contrast       = "treated_vs_control",
#'   id_col         = "gene_names",
#'   unique_id_col  = "phosphosite_id",
#'   label_col      = "phosphosite_name",
#'   target_names   = c("EGFR", "AKT1"),
#'   top_n_sign     = 10
#' )
#' result$plot_out
#' }
#'
#' @export
plot_volcano <- function(
    se,
    contrast,
    id_col                = "gene_names",
    unique_id_col         = NULL,
    label_col             = NULL,
    target_names          = c(""),
    label_sign            = TRUE,
    label_targets         = TRUE,
    top_n_sign            = 20,
    top_n_nonsign         = 0,
    direction             = "both",
    max.overlaps          = 15,
    labelsize             = 2,
    pvalorqval            = "qval",
    pqcutoff              = 0.05,
    diffcutoff            = 0.58,
    sig_from_column       = FALSE,
    fill_colors           = c(Significant = "#E69F00", Target = "#0072B2", None = "grey90"),
    shape_values          = c(`FALSE` = 21, `TRUE` = 24),
    alpha_values          = c(None = 0.2, goi = 0.7),
    size_values           = c(`FALSE` = 1.5, `TRUE` = 2.5),
    target_shape_override = NULL,
    target_color          = "#0072B2",
    label_colors          = c(Target = "#0072B2", Top = "black")
) {
  .assert_se(se, require_rowdata = id_col, require_contrast = contrast)


  # -------------------------------------------------------------------------
  # 0.  Resolve helper columns
  #     dedup_col:          column used to deduplicate the combined label data
  #                         frame; site-level for phosphoproteomics, gene-level
  #                         for standard proteomics.
  #     label_col_resolved: column whose values appear as text labels on the
  #                         plot; defaults to "name" if the user did not
  #                         specify label_col.
  # -------------------------------------------------------------------------
  dedup_col          <- if (!is.null(unique_id_col)) unique_id_col else id_col
  label_col_resolved <- if (!is.null(label_col))     label_col     else "name"
  
  # -------------------------------------------------------------------------
  # 1.  Derive column names from the contrast prefix
  #     All result columns are expected to follow the pattern
  #     <contrast>_<suffix> as produced by the upstream DE workflow.
  # -------------------------------------------------------------------------
  pval       <- paste0(contrast, "_p.val")
  padj       <- paste0(contrast, "_p.adj")
  diff       <- paste0(contrast, "_diff")
  sign       <- paste0(contrast, "_significant")
  # Select the column to threshold for significance
  pvalue_col <- ifelse(pvalorqval == "pval", pval, padj)
  
  # -------------------------------------------------------------------------
  # 2.  Build the plotting data frame from rowData
  #
  #     All columns needed downstream are collected upfront in extra_cols so
  #     that a single select() call suffices -- avoiding repeated rowData()
  #     calls and pipe-in-brace workarounds that are incompatible with |>.
  #
  #     Derived columns added here:
  #       target_gene  - "Target" / "Non-Target" flag for rendering order
  #       sign         - logical significance flag (p-value AND fold-change)
  #       target_sign  - three-level factor for fill colour
  #       alpha        - transparency group for background vs. foreground points
  #       score        - ranking score: -log10(p) * |diff|
  # -------------------------------------------------------------------------
  
  # Validate unique_id_col exists in rowData before attempting to select it
  if (!is.null(unique_id_col)) {
    rd_names <- names(SummarizedExperiment::rowData(se))
    if (!unique_id_col %in% rd_names)
      stop(paste0("unique_id_col '", unique_id_col, "' not found in rowData."))
  }

  # Significance can be taken from a precomputed <contrast>_significant column
  # (sig_from_column = TRUE, the behaviour of the older se_volcano()/lars_volcano())
  # or recomputed from the p-value and fold-change cutoffs (default). The
  # expression is built here and unquoted inside the mutate() below.
  if (isTRUE(sig_from_column) && !sign %in% names(SummarizedExperiment::rowData(se)))
    stop(paste0("sig_from_column = TRUE but significance column '", sign,
                "' was not found in rowData."))
  sig_expr <- if (isTRUE(sig_from_column)) {
    rlang::quo(as.logical(!!dplyr::sym(sign)))
  } else {
    rlang::quo(!!dplyr::sym(pvalue_col) <= pqcutoff & abs(!!dplyr::sym(diff)) > diffcutoff)
  }

  # Build the complete set of non-contrast columns to retain
  extra_cols <- unique(c(
    id_col,
    "name",
    if (!is.null(unique_id_col)) unique_id_col,
    if (!is.null(label_col))     label_col
  ))
  
  plot_data <- SummarizedExperiment::rowData(se) |>
    as.data.frame() |>
    dplyr::select(dplyr::all_of(extra_cols), dplyr::starts_with(contrast)) |>
    dplyr::mutate(
      # Flag rows belonging to user-specified targets (matched by id_col)
      target_gene = ifelse(!!dplyr::sym(id_col) %in% target_names, "Target", "Non-Target"),
      
      # Significance flag: from the precomputed column or recomputed (see sig_expr)
      sign = !!sig_expr,
      
      # Three-level category drives the fill colour scale
      target_sign = dplyr::case_when(
        !!dplyr::sym(id_col) %in% target_names ~ "Target",
        sign                                    ~ "Significant",
        TRUE                                    ~ "None"
      ),
      
      # Two-level alpha group: faint background vs. highlighted foreground
      alpha = ifelse(target_sign == "None", "None", "goi"),
      
      # Ranking score used to pick top-n labels; raw p-value is intentional
      # so that the score reflects unadjusted statistical strength
      score = -log10(!!dplyr::sym(pval)) * abs(!!dplyr::sym(diff))
    ) |>
    # Arrange so non-targets are drawn first and targets appear on top
    dplyr::arrange(target_gene == "Target")
  
  # -------------------------------------------------------------------------
  # 3.  Internal helper: select top-n rows by score within a direction
  #
  #     Targets are excluded from this selection; they are added back in
  #     step 4 so they are always labelled regardless of their score rank.
  #
  #     When direction = "both", the label budget (n) is split as:
  #       left  (down-regulated, diff < 0): ceiling(n / 2)
  #       right (up-regulated,   diff > 0): floor(n / 2)
  # -------------------------------------------------------------------------
  .top_by_direction <- function(data, sig_flag, n) {
    if (direction == "left") {
      data |>
        dplyr::filter(!!dplyr::sym(diff) < 0, sign == sig_flag) |>
        dplyr::slice_max(order_by = score, n = n)
      
    } else if (direction == "right") {
      data |>
        dplyr::filter(!!dplyr::sym(diff) > 0, sign == sig_flag) |>
        dplyr::slice_max(order_by = score, n = n)
      
    } else {
      # "both": select independently from each side then combine
      left  <- data |>
        dplyr::filter(!!dplyr::sym(diff) < 0, sign == sig_flag) |>
        dplyr::slice_max(order_by = score, n = ceiling(n / 2))
      right <- data |>
        dplyr::filter(!!dplyr::sym(diff) > 0, sign == sig_flag) |>
        dplyr::slice_max(order_by = score, n = floor(n / 2))
      dplyr::bind_rows(left, right)
    }
  }
  
  # Exclude targets before selecting top-n so target labels are not consumed
  # by the top-n budget and are always shown regardless of their score rank
  non_targets <- plot_data |>
    dplyr::filter(!(!!dplyr::sym(id_col) %in% target_names))
  
  # Top significant non-target points (suppressed entirely when label_sign = FALSE)
  label_data_sign <- non_targets |>
    .top_by_direction(sig_flag = TRUE,  n = if (isTRUE(label_sign)) top_n_sign else 0L) |>
    dplyr::mutate(label_type = "Top")
  
  # Top non-significant non-target points (typically 0; exposed for edge cases)
  label_data_nonsign <- non_targets |>
    .top_by_direction(sig_flag = FALSE, n = top_n_nonsign) |>
    dplyr::mutate(label_type = "Top")
  
  label_data <- dplyr::bind_rows(label_data_sign, label_data_nonsign)
  
  # -------------------------------------------------------------------------
  # 4.  Target labels
  #
  #     All rows whose id_col value is in target_names are retained -- in
  #     phosphoproteomics this means every phosphosite of a target gene gets
  #     its own label, not just the first one encountered.
  # -------------------------------------------------------------------------
  # Target labels (suppressed entirely when label_targets = FALSE)
  label_targets_df <- plot_data |>
    dplyr::filter(if (isTRUE(label_targets)) !!dplyr::sym(id_col) %in% target_names else FALSE) |>
    dplyr::mutate(label_type = "Target")
  
  # -------------------------------------------------------------------------
  # 5.  Combine label sets and deduplicate at the site level
  #
  #     Targets are bound after top-n rows so that distinct() keeps the
  #     "Target" label_type when a site appears in both sets (targets win
  #     because bind_rows puts label_data first and distinct() keeps the
  #     first occurrence -- but targets come second, so we bind targets last
  #     and rely on targets being absent from non_targets above).
  #     Deduplication uses dedup_col (= unique_id_col if supplied, else
  #     id_col) so multiple phosphosites of the same gene are preserved.
  # -------------------------------------------------------------------------
  label_combined <- dplyr::bind_rows(label_data, label_targets_df) |>
    dplyr::distinct(!!dplyr::sym(dedup_col), .keep_all = TRUE)
  
  # -------------------------------------------------------------------------
  # 6.  Base volcano plot
  #
  #     Aesthetic mappings:
  #       fill  -> target_sign  (Significant / Target / None)
  #       shape -> sign         (FALSE = circle, TRUE = triangle)
  #       size  -> sign         (significant points are drawn larger)
  #       alpha -> alpha        (background points are faint)
  # -------------------------------------------------------------------------
  p1 <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes_string(x = diff, y = paste0("-log10(", pval, ")"))
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill  = target_sign,
        alpha = alpha,
        shape = sign,
        size  = sign
      ),
      stroke      = 0.1,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(values = fill_colors,   drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_values, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = alpha_values, drop = FALSE) +
    ggplot2::scale_size_manual(values = size_values,   drop = FALSE) +
    ggplot2::guides(
      alpha = "none",   # alpha is a rendering detail, not informative for the reader
      fill  = ggplot2::guide_legend(
        title        = "Category",
        override.aes = list(shape = 21, size = 2)
      ),
      shape = ggplot2::guide_legend(
        title        = "Significant",
        override.aes = list(size = 2)
      ),
      size = "none"     # size is redundant with shape for significance
    ) +
    ggplot2::labs(
      title = contrast,
      x     = "Log2(Fold Change)",
      y     = "-Log10(P-value)"
    ) +
    ggplot2::theme_bw()
  
  # -------------------------------------------------------------------------
  # 7.  Optional target shape override
  #
  #     When target_shape_override is set, target points are re-drawn as an
  #     additional layer on top of all existing layers using a fixed pch shape
  #     and fill colour. This guarantees targets are always visually distinct
  #     from significant non-target points, which share the same shape scale.
  # -------------------------------------------------------------------------
  if (!is.null(target_shape_override)) {
    p1 <- p1 +
      ggplot2::geom_point(
        data        = subset(plot_data, target_gene == "Target"),
        shape       = target_shape_override,
        fill        = target_color,
        color       = "black",
        size        = max(size_values),
        stroke      = 0.2,
        show.legend = FALSE   # already represented in the fill legend
      )
  }
  
  # -------------------------------------------------------------------------
  # 8.  Unified repelled labels
  #
  #     A single geom_text_repel call handles all label types simultaneously
  #     so the repulsion algorithm treats targets and top-n points together,
  #     avoiding inter-group overlaps that would occur with separate calls.
  #
  #     Label styling:
  #       color    -> label_type  (Target = blue, Top = black)
  #       fontface -> bold for targets, plain for top-n
  # -------------------------------------------------------------------------
  if (nrow(label_combined) > 0) {
    p1 <- p1 +
      ggrepel::geom_text_repel(
        data = label_combined,
        ggplot2::aes(
          label    = !!dplyr::sym(label_col_resolved),
          color    = label_type,
          fontface = ifelse(label_type == "Target", "bold", "plain")
        ),
        size               = labelsize,
        min.segment.length = 0,        # always draw a segment, even for close labels
        max.overlaps       = max.overlaps,
        segment.size       = 0.2,
        segment.color      = "black",
        box.padding        = 0.3,
        point.padding      = 0.3,
        bg.color           = "white",  # white halo improves legibility on dense plots
        bg.r               = 0.15
      ) +
      ggplot2::scale_color_manual(values = label_colors, drop = FALSE) +
      ggplot2::guides(color = "none")  # label colour is conveyed by the fill legend
  }
  
  # -------------------------------------------------------------------------
  # 9.  Return plot and label data frame
  # -------------------------------------------------------------------------
  list(
    plot_out = p1,   # ggplot2 object; render with print() or ggsave()
    df_out   = label_combined   # data frame of all labelled points for inspection
  )
}




#' Create a 2D plot_scatter
#' 
#' Plots two columns against each other and marks entries that differ more than standard_dev SD. This calculation is performed as a moving window.
#' @param df data for plotting
#' @param col_x string column name of x-axis data 
#' @param col_y string column name of y-axis data
#' @param col_label string column name of labels
#' @param show_labels Boolean show labels on plot
#' @param title string plot title
#' @param standard_dev int difference to sd that specifies highlighted entries 
#' @param window rolling window for which sd is calculated
#' @return list containing the scatterplot and dfs for highlighted entries
#' @importFrom dplyr group_by summarize filter ungroup select mutate mutate_all rename rename_all
#' @importFrom stats cor.test
#' @importFrom ggrepel geom_text_repel
#' @importFrom stringr str_wrap
#' @export

plot_scatter <- function(df, col_x, col_y, col_label = "gene_names", show_labels = TRUE, title = "", standard_dev=2, window=1) {
  x <- df[,col_x]
  y <- df[,col_y]
  col_label <- df[,col_label]
  mod = lm(unlist(y)~unlist(x))
  slope=mod$coefficients["x"]
  intercept=if (is.na(mod$coefficients["(Intercept)"])) 0 else mod$coefficients["(Intercept)"]
  angle = atan(slope)
  
  rotate = function (x,y, theta, intercept){
    c(cos(theta) * x + sin(theta) * (y-intercept),
      -sin(theta) * x + cos(theta) * (y-intercept))
  }
  
  ## Rotate so that slope is 0
  rot = mapply(rotate, x, y, angle, intercept)
  xp = rot[1,]
  yp = rot[2,]
  
  ## Compute standard deviation along the x axis at position pos,
  ## for data in a given window (a lower window will be less smooth)
  sdWindow = function (pos, window, x, y) {
    sel = which (x >= pos-window/2 & x <= pos+window/2) 
    sd(y[sel], na.rm = TRUE)
  }
  
  xsd = seq(min(xp, na.rm = TRUE),max(xp, na.rm = TRUE), length=100)
  ysd = sapply(xsd, sdWindow, window, xp, yp)
  
  ## Transform back the data, and two curves for the + and - standard dev
  xpp = mapply(rotate, xp, yp, -angle, -intercept)
  sdp = mapply(rotate, xsd, standard_dev*ysd, -angle, -intercept)
  sdm = mapply(rotate, xsd, -standard_dev*ysd, -angle, -intercept)
  
  sdpdf = data.frame(x=sdp[1,], y=sdp[2,])
  sdmdf = data.frame(x=sdm[1,], y=sdm[2,])
  
  cbPalette <- c("#dddddd","#0000ff", "#00ff00")
  df = data.frame(x=x, y=y)
  
  up<-data.frame("x"=lapply(df$x, function(i) .closest(sdpdf$x, i))%>%unlist(),"x_old"=df$x,"y"=df$y, "col_label"=col_label)
  up_col_labels<-merge(up,sdpdf, by="x")%>%filter(y.x>=y.y)

  down<-data.frame("x"=lapply(df$x, function(i) .closest(sdmdf$x, i))%>%unlist(),"x_old"=df$x,"y"=df$y, "col_label"=col_label)
  down_col_labels<-merge(down,sdmdf, by="x")%>%filter(y.x<=y.y)
  
  pearson<-cor.test(x,y)
  
  p = ggplot(df, aes(x=x, y=y))
  
  if(show_labels == TRUE){
    p <- p + geom_point (shape=20, color="black", size=3, alpha = 0.5) +
      geom_abline (intercept=intercept, slope=mod$coefficients["x"], color="darkgrey", linetype=2) +
      geom_line (data=sdpdf, aes(x=x, y=y), color="red", alpha=0.5,linetype=2, size=1) +
      geom_line (data=sdmdf, aes(x=x, y=y), color="blue", alpha=0.5, linetype=2, size=1) +
      
      #scale_colour_manual(values=cbPalette) +
      scale_size_manual(values=c(2,3)) +
      geom_point(data=up_col_labels, aes(x=x_old, y=y.x), colour="red")+
      geom_text_repel(data=up_col_labels, aes(x=x_old, y=y.x,label=col_label), nudge_y = 0.1,fill = alpha(c("white"),0.3), label.padding = 0.5, size=3, max.overlaps = Inf)+
      geom_point(data=down_col_labels, aes(x=x_old, y=y.x), colour="blue")+
      geom_text_repel(data=down_col_labels, aes(x=x_old, y=y.x,label=str_wrap(col_label,30)),nudge_y=-0.2,fill = alpha(c("white"),0.3), label.padding = 0.5, size=3, max.overlaps = Inf)+
      geom_label(aes(x=ifelse(any(up_col_labels$x_old<= -0.5), 0.7,-1),y=0.8),label=paste("r: ",as.character(signif(pearson$estimate,digits=3)) , "\np-value: ",as.character(signif(pearson$p.value, digits=3))), hjust=0, fontface="bold", max.overlaps = Inf)+
      theme_bw () + 
      theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5), legend.position = "none", 
            axis.title = element_text(face="bold", lineheight=0.6))+
      labs(x=col_x, y=col_y, title=title)
  }else{
    p <- p + geom_point (shape=20, color="black", size=3, alpha = 0.5) +
      geom_abline (intercept=intercept, slope=mod$coefficients["x"], color="darkgrey", linetype=2) +
      geom_line (data=sdpdf, aes(x=x, y=y), color="red", alpha=0.5,linetype=2, size=1) +
      geom_line (data=sdmdf, aes(x=x, y=y), color="blue", alpha=0.5, linetype=2, size=1) +
      
      #scale_colour_manual(values=cbPalette) +
      scale_size_manual(values=c(2,3)) +
      geom_point(data=up_col_labels, aes(x=x_old, y=y.x), colour="red")+
      geom_point(data=down_col_labels, aes(x=x_old, y=y.x), colour="blue")+
      geom_label(aes(x=ifelse(any(up_col_labels$x_old<= -0.5), 0.7,-1),y=0.8),label=paste("r: ",as.character(signif(pearson$estimate,digits=3)) , "\np-value: ",as.character(signif(pearson$p.value, digits=3))), hjust=0, fontface="bold", max.overlaps = Inf)+
      theme_bw () + 
      theme(plot.title = element_text(lineheight=.8, face="bold", hjust=0.5), legend.position = "none", 
            axis.title = element_text(face="bold", lineheight=0.6))+
      labs(x=col_x, y=col_y, title=title)    
  }
  return(list(plot = p, up_col_labels = up_col_labels, down_col_labels = down_col_labels))
}


#' Create a correlation plot with different conditions and labels
#'
#' This function takes a dataframe and creates a scatter plot with correlation
#' analysis, highlighting points based on certain criteria such as being in a
#' gene list or having a rank among the top N differences.
#'
#' @param df A data frame with columns for x, y, and label names.
#' @param x_column Character string specifying the column name for the x-axis values.
#' @param y_column Character string specifying the column name for the y-axis values.
#' @param label_names Character string specifying the column name for the labels.
#' @param gene_list A character vector of gene names to be highlighted. Default is an empty vector.
#' @param top_genes Integer specifying the number of top differing proteins to consider. Default is  20.
#' @param max.overlaps Integer specifying the maximum number of overlapping labels. Default is Inf.
#' 
#' @return A list containing the plot and subsets of the data used for each part of the plot.
#' @export
#'
#' @importFrom patchwork plot_layout
#' @importFrom ggpubr stat_cor
#' @importFrom ggrepel geom_label_repel
#' @import ggplot2
#' @import dplyr

plot_correlation <- function(df, x_column, y_column, label_names, gene_list = c(""), top_genes = 20, max.overlaps = Inf) {
  x_char <- df %>% dplyr::select({ {x_column} }) %>% colnames()
  y_char <- df %>% dplyr::select({ {y_column} }) %>% colnames()
  
  df <- df %>% dplyr::select({ {x_column} }, { {y_column} }, { {label_names} }) %>%
    filter(!is.na({ {x_column} }) | !is.na({ {y_column} })) %>%
    mutate(abs_diff = ifelse((is.na({{x_column}}) | is.na({{y_column}})), NA, abs({{x_column}} - {{y_column}})),
           in_gene_list = {{label_names}} %in% gene_list)
  
  
  sub_1 <- df %>% filter(!is.na({ {x_column} }) & !is.na({ {y_column} })) %>% 
    mutate(top_genes_common = rank(-abs_diff, ties.method ="first") < top_genes,
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_common, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_2 <- sub_1 %>% filter(top_genes_common & !in_gene_list)
  sub_3 <- sub_1 %>% filter(in_gene_list)
  
  x_min = min(dplyr::select(sub_1, c({ {x_column} }, { {y_column} })))
  x_max = max(dplyr::select(sub_1, c({ {x_column} }, { {y_column} })))
  
  p1 <- ggplot(sub_1, aes(x = { {x_column} }, y = { {y_column} }, label = { {label_names} })) +
    geom_point(aes(fill = target, shape = factor(top_genes_common, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    scale_shape_manual(values = c("FALSE"= 21, "TRUE" = 23), drop = FALSE) +
    geom_label_repel(data = filter(sub_3, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = max.overlaps, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_2, top_genes_common), color = "firebrick", show.legend = FALSE, max.overlaps = max.overlaps, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    geom_abline(intercept = 0, linetype = "solid", color = "black") +
    ggpubr::stat_cor(data = sub_1, aes(label = ..r.label..)) +
    labs(title = "Common proteins", fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(x_min, x_max), y = c(x_min, x_max)) +
    theme_bw() 
  
  sub_4 <- df %>% filter(!is.na({{x_column}}) & is.na({{y_column}})) %>% 
    mutate(`:=`({{y_column}}, "NA")) %>% 
    mutate(top_genes_x = rank(-{{x_column}}, ties.method ="first") < top_genes,
           jittered = runif(nrow(.), min = -0.4, max = 0.4),
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_x, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_5 <- sub_4 %>% filter(!in_gene_list)
  sub_6 <- sub_4 %>% filter(in_gene_list)
  
  p2 <- ggplot(sub_4, aes(x = jittered, y = {{x_column}}, label = {{label_names}})) +
    geom_point(aes(fill = target, shape = factor(top_genes_x, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    scale_shape_manual(values = c(`FALSE`= 21, `TRUE` = 23), drop = FALSE) +
    
    geom_label_repel(data = filter(sub_6, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = Inf,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_5, top_genes_x), color = "firebrick", show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    labs(title = paste0("Only in ", x_char), fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(-0.4, 0.4)) +
    theme_bw()+
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()    
    )
  
  sub_7 <- df %>% 
    filter(is.na({{x_column}}) & !is.na({{y_column}})) %>% 
    mutate(`:=`({{x_column}}, "NA")) %>% 
    mutate(top_genes_y = rank(-{{y_column}}, ties.method ="first") < top_genes,
           jittered = runif(nrow(.), min = -0.4, max = 0.4),
           target = factor(ifelse(in_gene_list, "in_gene_list", ifelse(top_genes_y, "top_gene", "gene")), levels = c("in_gene_list", "top_gene", "gene")))
  
  sub_8 <- sub_7 %>% filter(!in_gene_list)
  sub_9 <- sub_7 %>% filter(in_gene_list)
  
  p3 <- ggplot(sub_7, aes(x = jittered, y = {{y_column}}, label = {{label_names}})) +
    geom_point(aes(fill = target, shape = factor(top_genes_y, levels = c("TRUE", "FALSE"))), size = 2, alpha = 0.4, show.legend = TRUE) +
    scale_shape_manual(values = c(`FALSE`= 21, `TRUE` = 23), drop = FALSE) +
    scale_fill_manual(values = c("gene" = "grey90", "top_gene" = "firebrick", "in_gene_list" = "darkorange"), drop = FALSE) +
    geom_label_repel(data = filter(sub_9, in_gene_list), color = "darkorange", show.legend = FALSE, max.overlaps = Inf,
                     fill = alpha(c("white"),0.7)) +
    geom_label_repel(data = filter(sub_8, top_genes_y), color = "firebrick", show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                     fill = alpha(c("white"),0.7)) +
    
    labs(title = paste0("Only in ", y_char), fill = "Targets", shape = "Top Genes") +
    guides(fill = guide_legend( 
      override.aes=list(shape = 21))) +
    lims(x = c(-0.4, 0.4)) +
    theme_bw()+
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  panel_list <- list()
  panel_list[["common"]] <- p1
  panel_list[[x_char]] <- p2
  panel_list[[y_char]] <- p3
  
  l <- list()
  l[["plot"]] <- (p1 + p2 + p3) + patchwork::plot_layout(guides = "collect", widths = c(0.5, 0.3, 0.3))
  l[["common"]] <- sub_1
  l[[x_char]] <- sub_4
  l[[y_char]] <- sub_7
  l[["top_genes"]] <- sub_2
  l[["panels"]] <- panel_list
  
  return(l)
}


#' Create a clustered heatmap from a given SummarizedExperiment object
#'
#' @param se A SummarizedExperiment object.
#' @param indicate A character string indicating the variable to be used for annotation.
#' @param type A character string indicating the type of heatmap.
#' @param k An integer specifying the number of clusters to be used in k-means clustering.
#' @param ... Additional arguments to be passed to the underlying ComplexHeatmap::Heatmap function.
#' 
#' @return A list containing the heatmap plot, a data frame with the heatmap data, and a data frame with the protein clusters.
#'
#' @importFrom DEP2 plot_heatmap
#' @import dplyr
#' @export
#' @examples
#' \dontrun{
#' # Assuming se_diff is a SummarizedExperiment object with appropriate data
#' result <- plot_heatmap_clustered(se_diff, indicate = "condition", type = "centered", k = 5)
#' result$plot
#' result$df
#' result$clusters
#' }
plot_heatmap_clustered <- function(se, indicate = "condition", type = "centered", k = 3, ...){
  .assert_se(se, require_coldata = indicate)

  p_heatmap <- se %>% DEP2::plot_heatmap(indicate = indicate, type = type, kmeans = TRUE, k = k)  

  p_heatmap_data <- se %>% DEP2::plot_heatmap(indicate = indicate, type = type, kmeans = TRUE, k = k, plot = FALSE)  
  
  split <- p_heatmap_data %>% select(k, protein) %>% split.data.frame(.$k)
  split <- lapply(split, function(x) select(x,protein))
  
  split_df <- .cbind_pad(split)
  colnames(split_df) <- paste0("cluster_", 1:length(split))
  rownames(split_df) <- NULL
  
  return(list("plot" = p_heatmap, "df" = p_heatmap_data, "clusters" = split_df))
}

########################################################################
########################################################################
########################################################################
#' 2-D Scatter Plot Comparing Two Differential Expression Results
#'
#' Joins two differential expression data sets on a shared key, classifies
#' each feature into a directional group (e.g. \code{"both-up"},
#' \code{"x-up"}), and draws a scatter plot of log2 fold-changes with
#' optional correlation annotation, threshold lines, and repelled labels.
#'
#' @param x,y Data inputs.  Each may be a \code{data.frame}, a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}, or any object
#'   coercible by \code{as_df}.
#' @param x_key,y_key Character scalar.  Column name used as the join key in
#'   \code{x} and \code{y} respectively (e.g. protein or gene identifier).
#' @param x_col,y_col Character scalar.  Column containing the log2
#'   fold-change in \code{x} and \code{y}.
#' @param x_anno,y_anno Character scalar.  Column used as the point label /
#'   display name in \code{x} and \code{y}.  When both sides carry identical
#'   values the label is shown once; otherwise the two are concatenated with
#'   \code{" / "}.
#' @param x_padj,y_padj Character scalar.  Column containing the adjusted
#'   p-value in \code{x} and \code{y}.
#' @param x_lfc_diff,y_lfc_diff Numeric scalar.  Absolute log2 fold-change
#'   threshold for significance in \code{x} and \code{y}.
#'   Default: \code{0.58} (~1.5-fold).
#' @param x_alpha_diff,y_alpha_diff Numeric scalar.  Adjusted p-value
#'   threshold for significance in \code{x} and \code{y}.
#'   Default: \code{0.05}.
#' @param require_sig Logical.  If \code{TRUE} (default), directional group
#'   assignment requires both the fold-change and p-value thresholds to be
#'   met.  If \code{FALSE}, only the fold-change threshold is used.
#' @param as_df Function or \code{NULL}.  Coercion function applied to
#'   \code{x} and/or \code{y} when they are not already data frames or
#'   \code{SummarizedExperiment} objects (e.g. \code{DEP2::get_df_wide}).
#'   Also applied to \code{SummarizedExperiment} inputs when supplied,
#'   overriding the default \code{DEP2::get_df_wide} fallback.
#'   Default: \code{NULL}.
#' @param join Character scalar.  Join type used to merge \code{x} and
#'   \code{y}.  One of \code{"inner"} (default), \code{"left"},
#'   \code{"right"}, or \code{"full"}.
#' @param label_top_n Integer scalar.  Number of top-scoring features to
#'   label across all groups listed in \code{label_groups}.  Ranking is
#'   performed globally across all eligible groups (including \code{"both-up"}
#'   / \code{"both-down"} when their constituent groups are eligible), so
#'   the top-n are always the genuinely highest-scoring points regardless of
#'   which specific group they fell into.  \code{0} disables top-N labelling.
#'   Default: \code{0}.
#' @param label_groups Character vector or \code{NULL}.  Directional groups
#'   eligible for top-N labelling.  Valid values: \code{"x-up"},
#'   \code{"x-down"}, \code{"y-up"}, \code{"y-down"}, \code{"both-up"},
#'   \code{"both-down"}.  \code{"both-up"} is automatically included when
#'   either \code{"x-up"} or \code{"y-up"} is present; similarly for
#'   \code{"both-down"}.  Default: \code{NULL} (no top-N labels).
#' @param label_keys Character vector or \code{NULL}.  Specific feature keys
#'   (matched against the join key column) to label regardless of score.
#'   Default: \code{NULL}.
#' @param label_keys_mode Character scalar.  \code{"always"} (default) labels
#'   \code{label_keys} regardless of their directional group;
#'   \code{"if_visible"} only labels them when their group is listed in
#'   \code{highlight_groups}.
#' @param highlight_groups Character vector or \code{NULL}.  Directional
#'   groups to render in colour; all others are collapsed to \code{"neither"}
#'   (grey).  \code{NULL} (default) colours all groups.
#' @param label_score Character scalar.  Scoring method used to rank features
#'   for top-N labelling and to scale point sizes.  One of:
#'   \describe{
#'     \item{\code{"diff"}}{Absolute difference between x and y LFC.}
#'     \item{\code{"sum_abs"}}{Sum of absolute LFCs.}
#'     \item{\code{"product_abs"}}{Product of absolute LFCs.}
#'     \item{\code{"max_abs"}}{Maximum absolute LFC across x and y.}
#'   }
#'   Default: \code{"diff"}.
#' @param title,subtitle,caption Character scalar or \code{NULL}.  Plot
#'   annotations.  \code{subtitle} is auto-generated from matched feature
#'   counts when \code{NULL}.
#' @param x_label,y_label Character scalar or \code{NULL}.  Axis labels.
#'   When \code{NULL} (default), labels are derived from the name of the
#'   input object and the column name (e.g. \code{"imp_se_std\ntnfa_vs_pbs_diff"}).
#'   Supply a string to override completely, which is useful when both
#'   \code{x_col} and \code{y_col} share the same name.
#' @param x_lim,y_lim Numeric vector of length 2 or \code{NULL}.  Manual
#'   axis limits (e.g. \code{c(-3, 3)}).  Partial override is supported:
#'   supply only one to override a single axis while the other is computed
#'   automatically.  When both are \code{NULL} (default), limits are derived
#'   from the data, optionally with \code{symmetric_axes = TRUE}.
#' @param group_colors Named character vector.  Point colours for each
#'   directional group.  Must include \code{"neither"} for background points.
#'   Any groups not present in the named vector receive auto-generated colours.
#' @param point_alpha Numeric scalar in \code{(0, 1]}.  Point transparency.
#'   Default: \code{0.8}.
#' @param point_size_range Numeric vector of length 2.  Minimum and maximum
#'   point sizes mapped to \code{label_score}.  Default: \code{c(1.2, 4.5)}.
#' @param add_diagonal Logical.  Draw a \eqn{y = x} reference line.
#'   Default: \code{TRUE}.
#' @param add_threshold_lines Logical.  Draw dashed lines at the LFC
#'   thresholds.  Default: \code{TRUE}.
#' @param symmetric_axes Logical.  Force both axes to the same absolute
#'   range, centred on zero.  Ignored for any axis with a manual limit
#'   supplied via \code{x_lim}/\code{y_lim}.  Default: \code{TRUE}.
#' @param expand Numeric scalar.  Fractional padding added to each
#'   automatically computed axis limit.  Default: \code{0.06}.
#' @param cor_method Character scalar.  Correlation method passed to
#'   \code{\link[stats]{cor.test}}.  One of \code{"pearson"} (default),
#'   \code{"spearman"}, or \code{"kendall"}.
#' @param show_cor Logical.  Annotate the plot with the correlation
#'   coefficient and p-value.  Default: \code{TRUE}.
#' @param show_counts Logical.  Include per-group feature counts in the
#'   auto-generated subtitle.  Default: \code{TRUE}.
#' @param highlight_label_keys Logical.  Re-draw \code{label_keys} points on
#'   top of all other layers with a black outline for visibility.
#'   Default: \code{TRUE}.
#' @param label_key_point_size Numeric scalar.  Size of the re-drawn key
#'   points.  Default: \code{3.0}.
#' @param labelsize Numeric scalar.  Text size for repelled labels (in
#'   \pkg{ggplot2} units).  Default: \code{3.1}.
#' @param max.overlaps Integer.  Passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.  Default: \code{Inf}.
#' @param label_colors Named character vector with elements \code{"Target"}
#'   and \code{"Top"}.  Text colours for key-label and top-N labels
#'   respectively.
#' @param return_data Logical.  If \code{TRUE}, returns a named list with
#'   elements \code{plot}, \code{data} (the full joined data frame), and
#'   \code{label_data} (the subset used for labels) instead of just the
#'   plot object.  Default: \code{FALSE}.
#'
#' @details
#' ## Directional groups
#' Each feature is assigned to exactly one mutually exclusive group in the
#' following priority order:
#' \enumerate{
#'   \item \code{"both-up"}   -- up in both x and y
#'   \item \code{"both-down"} -- down in both x and y
#'   \item \code{"x-up"}      -- up in x only
#'   \item \code{"x-down"}    -- down in x only
#'   \item \code{"y-up"}      -- up in y only
#'   \item \code{"y-down"}    -- down in y only
#'   \item \code{"neither"}   -- does not meet any threshold
#' }
#' When \code{require_sig = TRUE}, "up/down" means the feature passes both
#' the p-value and fold-change thresholds.  When \code{FALSE}, only the
#' fold-change threshold is applied.
#'
#' ## Top-N label ranking
#' Ranking for \code{label_top_n} is performed \emph{globally} across all
#' eligible groups rather than per-group.  \code{"both-up"} and
#' \code{"both-down"} are automatically included in the eligible pool when
#' their constituent directional groups appear in \code{label_groups}, so
#' they are never inadvertently excluded from a small within-group pool.
#'
#' ## Colour scale
#' \code{scale_color_manual} is called once with the combined point and label
#' text palette so that the two do not silently override each other.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or when
#'   \code{return_data = TRUE} a named list with elements \code{plot},
#'   \code{data}, and \code{label_data}.
#'
#' @examples
#' \dontrun{
#' plot_2D_scatter(
#'   x      = diff_final[["imp_se_std"]],
#'   y      = diff_final[["nonimp_se_unmatched_phos"]],
#'   x_key  = "name",              y_key  = "final_id",
#'   x_col  = "tnfa_vs_pbs_diff",  y_col  = "tnfa_vs_pbs_diff",
#'   x_anno = "name",              y_anno = "name",
#'   x_padj = "tnfa_vs_pbs_p.adj", y_padj = "tnfa_vs_pbs_p.adj",
#'   as_df          = DEP2::get_df_wide,
#'   label_top_n    = 10,
#'   label_groups   = c("y-up", "y-down"),
#'   highlight_groups = c("y-up", "y-down"),
#'   label_keys     = c("NPHS1", "THSD7A", "C3", "PIEZO1"),
#'   x_label        = "Imputed proteome\nTNFa vs ctrl",
#'   y_label        = "Non-imputed phospho\nTNFa vs ctrl",
#'   x_lim          = c(-4, 4),
#'   title          = "TNFa 48h vs ctrl 48h",
#'   require_sig    = FALSE
#' )
#' }
#'
#' @importFrom dplyr transmute mutate case_when if_else coalesce group_by ungroup arrange distinct filter inner_join left_join right_join full_join min_rank desc
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_abline geom_point annotate coord_cartesian scale_color_manual scale_fill_manual scale_shape_manual scale_size_continuous guides guide_legend labs theme_classic theme element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales hue_pal pretty_breaks
#' @importFrom stats cor.test
#'
#' @export
plot_2D_scatter <- function(
    x, y,
    x_key, y_key,
    x_col, y_col,
    x_anno, y_anno,
    x_padj, y_padj,
    x_lfc_diff           = 0.58,
    y_lfc_diff           = 0.58,
    x_alpha_diff         = 0.05,
    y_alpha_diff         = 0.05,
    require_sig          = TRUE,
    as_df                = NULL,
    join                 = c("inner", "left", "right", "full"),
    label_top_n          = 0,
    label_groups         = NULL,
    label_keys           = NULL,
    label_keys_mode      = c("always", "if_visible"),
    highlight_groups     = NULL,
    label_score          = c("diff", "sum_abs", "product_abs", "max_abs"),
    title                = NULL,
    subtitle             = NULL,
    caption              = NULL,
    x_label              = NULL,
    y_label              = NULL,
    x_lim                = NULL,
    y_lim                = NULL,
    group_colors         = c(
      "x-up"      = "#2166AC",
      "x-down"    = "#92C5DE",
      "y-up"      = "#B2182B",
      "y-down"    = "#F4A582",
      "both-up"   = "#762A83",
      "both-down" = "#C51B7D",
      "neither"   = "grey80"
    ),
    point_alpha          = 0.8,
    point_size_range     = c(1.2, 4.5),
    add_diagonal         = TRUE,
    add_threshold_lines  = TRUE,
    symmetric_axes       = TRUE,
    expand               = 0.06,
    cor_method           = c("pearson", "spearman", "kendall"),
    show_cor             = TRUE,
    show_counts          = TRUE,
    highlight_label_keys = TRUE,
    label_key_point_size = 3.0,
    labelsize            = 3.1,
    max.overlaps         = Inf,
    label_colors         = c(Target = "#0072B2", Top = "black"),
    return_data          = FALSE
) {
  
  # ---- Capture axis label names FIRST ----------------------------------------
  # deparse(substitute()) must be called before x/y are touched in any way,
  # otherwise it captures the internal variable name rather than the user's
  # expression.  width.cutoff keeps long list-indexing expressions on one line.
  x_deparsed <- deparse(substitute(x), width.cutoff = 60L)
  y_deparsed <- deparse(substitute(y), width.cutoff = 60L)
  
  # Clean up common noisy list-indexing patterns:
  #   diff_final[["imp_se_std"]]  ->  imp_se_std
  #   diff_final[['imp_se_std']]  ->  imp_se_std
  .clean_deparse <- function(s) {
    # If the expression contains [[ ]], extract just the quoted content inside
    if (grepl("\\[\\[", s)) {
      # Remove everything up to and including [[, and everything from ]] onward
      # Also strip any quote characters and backslashes left behind
      inner <- gsub(".*\\[\\[|\\]\\].*", "", s)          # extract between [[ and ]]
      inner <- gsub('["\'\\ ]', "", inner)                 # strip quotes, backslashes, spaces
      if (nzchar(inner)) return(inner)
    }
    # Fallback: truncate if too long
    if (nchar(s) > 40L) paste0(substr(s, 1L, 37L), "...") else s
  }
  
  # ---- Match arguments -------------------------------------------------------
  join            <- match.arg(join)
  label_score     <- match.arg(label_score)
  cor_method      <- match.arg(cor_method)
  label_keys_mode <- match.arg(label_keys_mode)
  
  # ---- Expand NULL sentinels -------------------------------------------------
  all_groups <- c("x-up", "x-down", "y-up", "y-down", "both-up", "both-down")
  if (is.null(highlight_groups)) highlight_groups <- all_groups
  if (is.null(label_groups))     label_groups     <- character(0)
  
  # ---- Input validation ------------------------------------------------------
  .check_scalar_pos <- function(val, name) {
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop("`", name, "` must be a single positive number.", call. = FALSE)
  }
  .check_scalar_nn <- function(val, name) {
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 0)
      stop("`", name, "` must be a single non-negative number.", call. = FALSE)
  }
  .check_lim <- function(val, name) {
    if (!is.null(val) && (!is.numeric(val) || length(val) != 2L || anyNA(val)))
      stop("`", name, "` must be a numeric vector of length 2 or NULL.", call. = FALSE)
  }
  
  .check_scalar_nn(label_top_n,  "label_top_n")
  .check_scalar_nn(expand,       "expand")
  .check_scalar_pos(label_key_point_size, "label_key_point_size")
  .check_scalar_pos(labelsize,            "labelsize")
  .check_lim(x_lim, "x_lim")
  .check_lim(y_lim, "y_lim")
  
  if (!is.numeric(point_size_range) || length(point_size_range) != 2L || anyNA(point_size_range))
    stop("`point_size_range` must be a numeric vector of length 2.", call. = FALSE)
  if (!is.null(label_keys) && !is.character(label_keys))
    stop("`label_keys` must be NULL or a character vector.", call. = FALSE)
  
  # -------------------------------------------------------------------------
  # 1.  Coerce inputs to data frames
  #     Priority: is.data.frame > SummarizedExperiment > as_df fallback.
  #     as_df overrides the SE default when explicitly supplied.
  # -------------------------------------------------------------------------
  .coerce <- function(obj, req_cols, side) {
    out <- if (is.data.frame(obj)) {
      obj
    } else if (inherits(obj, "SummarizedExperiment")) {
      if (!is.null(as_df)) as_df(obj) else DEP2::get_df_wide(obj)
    } else if (!is.null(as_df)) {
      as_df(obj)
    } else {
      stop("`", side, "` must be a data.frame or SummarizedExperiment, ",
           "or supply `as_df`.", call. = FALSE)
    }
    if (!is.data.frame(out))
      stop("Coercion of `", side, "` did not return a data.frame.", call. = FALSE)
    miss <- setdiff(req_cols, names(out))
    if (length(miss) > 0)
      stop("Missing columns in `", side, "`: ", paste(miss, collapse = ", "),
           call. = FALSE)
    out
  }
  
  x_df <- .coerce(x, c(x_key, x_col, x_anno, x_padj), "x")
  y_df <- .coerce(y, c(y_key, y_col, y_anno, y_padj), "y")
  
  # -------------------------------------------------------------------------
  # 2.  Standardise to a common schema and join
  # -------------------------------------------------------------------------
  x_std <- dplyr::transmute(x_df,
                            key    = .data[[x_key]],
                            x_lfc  = .data[[x_col]],
                            x_padj = .data[[x_padj]],
                            x_anno = .data[[x_anno]]
  )
  y_std <- dplyr::transmute(y_df,
                            key    = .data[[y_key]],
                            y_lfc  = .data[[y_col]],
                            y_padj = .data[[y_padj]],
                            y_anno = .data[[y_anno]]
  )
  
  join_fn <- list(
    inner = dplyr::inner_join,
    left  = dplyr::left_join,
    right = dplyr::right_join,
    full  = dplyr::full_join
  )[[join]]
  
  df <- join_fn(x_std, y_std, by = "key")
  
  if (nrow(df) == 0)
    stop(
      "No rows remained after the join.  ",
      "Check that `x_key` and `y_key` share common values.",
      call. = FALSE
    )
  
  # -------------------------------------------------------------------------
  # 3.  Significance flags and directional group assignment
  #
  #     .dir_flag() encapsulates the per-axis up/down test so the
  #     require_sig branch is not repeated four times.
  #
  #     Group priority (mutually exclusive, first match wins):
  #       both-up > both-down > x-up > x-down > y-up > y-down > neither
  # -------------------------------------------------------------------------
  .dir_flag <- function(sig, lfc, threshold, direction) {
    # direction: +1 = up, -1 = down
    lfc_ok <- if (direction == 1L) lfc >= threshold else lfc <= -threshold
    if (require_sig) sig & lfc_ok else lfc_ok
  }
  
  df <- dplyr::mutate(df,
                      sig_x  = !is.na(.data$x_padj) & .data$x_padj < x_alpha_diff & abs(.data$x_lfc) >= x_lfc_diff,
                      sig_y  = !is.na(.data$y_padj) & .data$y_padj < y_alpha_diff & abs(.data$y_lfc) >= y_lfc_diff,
                      x_up   = .dir_flag(.data$sig_x, .data$x_lfc, x_lfc_diff,  1L),
                      x_down = .dir_flag(.data$sig_x, .data$x_lfc, x_lfc_diff, -1L),
                      y_up   = .dir_flag(.data$sig_y, .data$y_lfc, y_lfc_diff,  1L),
                      y_down = .dir_flag(.data$sig_y, .data$y_lfc, y_lfc_diff, -1L),
                      group = dplyr::case_when(
                        .data$x_up   & .data$y_up                       ~ "both-up",
                        .data$x_down & .data$y_down                     ~ "both-down",
                        .data$x_up   & !.data$y_up & !.data$y_down      ~ "x-up",
                        .data$x_down & !.data$y_up & !.data$y_down      ~ "x-down",
                        .data$y_up   & !.data$x_up & !.data$x_down      ~ "y-up",
                        .data$y_down & !.data$x_up & !.data$x_down      ~ "y-down",
                        TRUE                                             ~ "neither"
                      ),
                      sig_cat = dplyr::case_when(
                        .data$sig_x & .data$sig_y ~ "Both",
                        .data$sig_x               ~ "X only",
                        .data$sig_y               ~ "Y only",
                        TRUE                      ~ "Neither"
                      ),
                      # Points not in highlight_groups are demoted to "neither" for colour only;
                      # their true group is preserved in `group` for all other logic.
                      group_plot = dplyr::if_else(
                        .data$group %in% highlight_groups, .data$group, "neither"
                      )
  )
  
  # -------------------------------------------------------------------------
  # 4.  Label score -- drives point sizing and top-N ranking.
  #     Named function list avoids a repeated case_when and is trivially
  #     extensible (add new scoring modes here without touching anything else).
  # -------------------------------------------------------------------------
  .score_fns <- list(
    diff        = function(x, y) abs(x - y),
    sum_abs     = function(x, y) abs(x) + abs(y),
    product_abs = function(x, y) abs(x) * abs(y),
    max_abs     = function(x, y) pmax(abs(x), abs(y))
  )
  
  df <- dplyr::mutate(df,
                      .score = dplyr::coalesce(
                        .score_fns[[label_score]](.data$x_lfc, .data$y_lfc),
                        -Inf
                      )
  )
  
  # -------------------------------------------------------------------------
  # 5.  Annotation label text
  #     When x_anno and y_anno carry identical values, show only one copy;
  #     otherwise concatenate with " / " to disambiguate.
  # -------------------------------------------------------------------------
  anno_ok   <- !is.na(df$x_anno) & !is.na(df$y_anno)
  same_anno <- any(anno_ok) && all(df$x_anno[anno_ok] == df$y_anno[anno_ok])
  
  # -------------------------------------------------------------------------
  # 6.  Top-N and keyword label selection
  #
  #     Key design decision for top-N ranking:
  #       Ranking is performed GLOBALLY across all eligible groups rather
  #       than per-group.  "both-up" / "both-down" are automatically added
  #       to the eligible pool when their constituent groups appear in
  #       label_groups.  This prevents "both-up" points from being ranked
  #       against a tiny pool of their own kind rather than competing
  #       fairly with the larger "x-up" / "y-up" populations.
  #
  #     label_type priority: "Target" (label_keys) beats "Top" (top-N).
  #     One row per key: duplicates are resolved by priority then score.
  # -------------------------------------------------------------------------
  
  # Expand eligible groups: both-up/down inherit eligibility from constituents
  eligible_groups <- c(
    label_groups,
    if (any(c("x-up",   "y-up")   %in% label_groups)) "both-up",
    if (any(c("x-down", "y-down") %in% label_groups)) "both-down"
  )
  
  # Step 1: add .eligible and .score are already on df at this point,
  # so compute .rank outside mutate where df$.eligible is real
  df <- dplyr::mutate(df,
                      .eligible = .data$group %in% eligible_groups
  )
  
  # Step 2: compute global rank vector outside mutate, then bind it in
  rank_vec <- rep(NA_integer_, nrow(df))
  elig_idx <- which(df$.eligible)
  if (length(elig_idx) > 0) {
    rank_vec[elig_idx] <- dplyr::min_rank(dplyr::desc(df$.score[elig_idx]))
  }
  
  df <- dplyr::mutate(df,
                      .rank = rank_vec,
                      .label_top = .data$.eligible &
                        label_top_n > 0 &
                        !is.na(.data$.rank) &
                        .data$.rank <= label_top_n,
                      .label_key = {
                        in_keys <- if (!is.null(label_keys)) .data$key %in% label_keys else FALSE
                        if (identical(label_keys_mode, "if_visible")) {
                          in_keys & .data$group %in% highlight_groups
                        } else {
                          in_keys
                        }
                      },
                      label_type = dplyr::case_when(
                        .data$.label_key ~ "Target",
                        .data$.label_top ~ "Top",
                        TRUE             ~ NA_character_
                      ),
                      .label_final = !is.na(.data$label_type),
                      label = dplyr::if_else(
                        .data$.label_final,
                        if (same_anno) as.character(.data$x_anno)
                        else           paste0(.data$x_anno, " / ", .data$y_anno),
                        NA_character_
                      )
  )
  
  # One label row per key; Target beats Top, then highest score breaks ties
  label_df <- df |>
    dplyr::filter(.data$.label_final) |>
    dplyr::mutate(.priority = dplyr::if_else(.data$label_type == "Target", 1L, 2L)) |>
    dplyr::arrange(.data$.priority, dplyr::desc(.data$.score)) |>
    dplyr::distinct(.data$key, .keep_all = TRUE) |>
    dplyr::select(-.data$.priority)
  
  # -------------------------------------------------------------------------
  # 7.  Colour palette
  #     Auto-generate colours for any unexpected groups not in group_colors.
  #     The combined palette (point colours + label text colours) is built
  #     once here so that scale_color_manual is called only once in the plot,
  #     avoiding the silent override bug from two competing scale calls.
  # -------------------------------------------------------------------------
  present_groups <- unique(df$group_plot)
  missing_cols   <- setdiff(present_groups, names(group_colors))
  if (length(missing_cols) > 0) {
    extra <- scales::hue_pal()(length(missing_cols))
    names(extra) <- missing_cols
    group_colors <- c(group_colors, extra)
  }
  active_colors   <- group_colors[intersect(names(group_colors), present_groups)]
  combined_colors <- c(active_colors, label_colors)
  
  # -------------------------------------------------------------------------
  # 8.  Axis limits
  #     Manual x_lim / y_lim override computed limits on a per-axis basis.
  #     symmetric_axes applies only to axes without a manual override.
  # -------------------------------------------------------------------------
  .pad_range <- function(r, frac = expand) {
    if (!all(is.finite(r))) return(c(-1, 1))
    if (diff(r) == 0)        return(r + c(-1, 1))
    r + c(-1, 1) * diff(r) * frac
  }
  
  if (symmetric_axes && is.null(x_lim) && is.null(y_lim)) {
    # Both axes computed symmetrically from the joint range
    max_abs <- max(abs(c(df$x_lfc, df$y_lfc)), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    x_lim_final <- y_lim_final <- .pad_range(c(-max_abs, max_abs))
  } else {
    x_lim_final <- .pad_range(range(df$x_lfc, na.rm = TRUE))
    y_lim_final <- .pad_range(range(df$y_lfc, na.rm = TRUE))
  }
  # Manual overrides always win, per-axis
  if (!is.null(x_lim)) x_lim_final <- x_lim
  if (!is.null(y_lim)) y_lim_final <- y_lim
  
  # -------------------------------------------------------------------------
  # 9.  Correlation annotation
  # -------------------------------------------------------------------------
  ok <- !is.na(df$x_lfc) & !is.na(df$y_lfc)
  cor_label <- NULL
  
  if (show_cor && sum(ok) >= 3L) {
    ct       <- stats::cor.test(df$x_lfc[ok], df$y_lfc[ok], method = cor_method)
    est_name <- c(pearson = "italic(r)", spearman = "italic(rho)",
                  kendall = "italic(tau)")[[cor_method]]
    p_str    <- if (ct$p.value < 0.001) " < 0.001" else sprintf(" == %.3f", ct$p.value)
    cor_label <- sprintf(
      "%s == %.3f*','~~italic(p)%s", est_name, unname(ct$estimate), p_str
    )
  }
  
  # -------------------------------------------------------------------------
  # 10. Axis labels
  #     Defaults use the cleaned deparse of the user's input expression +
  #     the column name on a new line, giving e.g.:
  #       "imp_se_std\ntnfa_vs_pbs_diff"
  #     This retains dataset identity even when x_col == y_col (a common
  #     case when comparing the same contrast across two experiments).
  #     Supply x_label / y_label to override completely.
  # -------------------------------------------------------------------------
  if (is.null(x_label)) x_label <- paste0(.clean_deparse(x_deparsed), "\n", x_col)
  if (is.null(y_label)) y_label <- paste0(.clean_deparse(y_deparsed), "\n", y_col)
  
  # -------------------------------------------------------------------------
  # 11. Subtitle
  # -------------------------------------------------------------------------
  if (is.null(subtitle)) {
    subtitle <- paste0("n = ", sum(ok), " matched features")
    if (show_counts) {
      counts   <- table(df$group)
      subtitle <- paste(
        subtitle,
        paste(paste0(names(counts), ": ", as.integer(counts)), collapse = " | "),
        sep = " | "
      )
    }
  }
  
  # -------------------------------------------------------------------------
  # 12. Assemble the plot
  #
  #     Layer order:
  #       reference lines -> threshold lines -> diagonal ->
  #       background (neither) points -> highlighted points ->
  #       target key redraw -> correlation text -> repelled labels
  #
  #     scale_color_manual is called ONCE with combined_colors so that
  #     directional group colours and label text colours share a single
  #     scale and cannot silently override each other.
  # -------------------------------------------------------------------------
  
  # Sort so background ("neither") points are drawn first, highlighted on top
  df_plot <- dplyr::arrange(df, .data$group_plot == "neither")
  
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$x_lfc, y = .data$y_lfc, color = .data$group_plot)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55")
  
  if (add_diagonal)
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linewidth = 0.35, linetype = "dotted", color = "grey40"
    )
  
  if (add_threshold_lines)
    p <- p +
    ggplot2::geom_hline(
      yintercept = c(-y_lfc_diff, y_lfc_diff),
      linetype = "dashed", linewidth = 0.4, color = "grey65"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-x_lfc_diff, x_lfc_diff),
      linetype = "dashed", linewidth = 0.4, color = "grey65"
    )
  
  # Base scatter points
  p <- p + ggplot2::geom_point(
    ggplot2::aes(shape = .data$sig_cat, size = .data$.score),
    alpha = point_alpha, stroke = 0.2
  )
  
  # Re-draw label_keys on top with a filled shape + black outline
  if (isTRUE(highlight_label_keys) && nrow(label_df) > 0) {
    target_df <- dplyr::filter(label_df, .data$label_type == "Target")
    if (nrow(target_df) > 0)
      p <- p + ggplot2::geom_point(
        data         = target_df,
        ggplot2::aes(
          x    = .data$x_lfc,
          y    = .data$y_lfc,
          fill = .data$group_plot
        ),
        inherit.aes  = FALSE,
        shape        = 21,
        size         = label_key_point_size,
        stroke       = 0.35,
        color        = "black",
        alpha        = 1,
        show.legend  = FALSE
      )
  }
  
  # Repelled labels -- single geom for both Target (bold) and Top (plain)
  if (nrow(label_df) > 0)
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        x        = .data$x_lfc,
        y        = .data$y_lfc,
        label    = .data$label,
        color    = .data$label_type,
        fontface = ifelse(.data$label_type == "Target", "bold", "plain")
      ),
      inherit.aes        = FALSE,
      size               = labelsize,
      min.segment.length = 0,
      max.overlaps       = max.overlaps,
      segment.size       = 0.25,
      segment.color      = "black",
      box.padding        = 0.35,
      point.padding      = 0.3,
      bg.color           = "white",
      bg.r               = 0.12,
      show.legend        = FALSE,
      na.rm              = TRUE
    )
  
  p <- p +
    ggplot2::coord_cartesian(xlim = x_lim_final, ylim = y_lim_final) +
    ggplot2::scale_color_manual(
      name   = "Directional group",
      values = combined_colors,
      breaks = setdiff(names(active_colors), "neither"),
      drop   = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = active_colors,
      guide  = "none",
      drop   = TRUE
    ) +
    ggplot2::scale_shape_manual(
      name   = "Significance",
      values = c(Both = 18, "X only" = 15, "Y only" = 17, Neither = 20),
      drop   = TRUE
    ) +
    ggplot2::scale_size_continuous(
      name   = paste0("Score (", label_score, ")"),
      range  = point_size_range,
      breaks = scales::pretty_breaks(n = 4)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(size = 3.5, alpha = 1)),
      shape = ggplot2::guide_legend(order = 2),
      size  = ggplot2::guide_legend(order = 3, override.aes = list(shape = 16, alpha = 1))
    ) +
    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = title,
      subtitle = subtitle,
      caption  = caption
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle   = ggplot2::element_text(size = 10, color = "grey30"),
      plot.caption    = ggplot2::element_text(size = 9,  color = "grey40"),
      axis.title      = ggplot2::element_text(face = "bold"),
      legend.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  if (!is.null(cor_label))
    p <- p + ggplot2::annotate(
      "text",
      x     = x_lim_final[1],
      y     = y_lim_final[2],
      label = cor_label,
      hjust = 0, vjust = 1,
      size  = 3.4,
      parse = TRUE,
      color = "grey30"
    )
  
  if (return_data) list(plot = p, data = df, label_data = label_df) else p
}


#### added 23 April 2026 (inspired by Arvids function)


#' Custom theme for ggplot2
#'
#' This function applies a custom theme to a plot by combining the DEP2::theme_DEP1()
#' function with additional theme elements.
#'
#' @return A ggplot2 theme object.
#' @importFrom DEP2 theme_DEP1
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
#' p + theme_proteoSE()
#' @export
theme_proteoSE <- function() {
  DEP2::theme_DEP1() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}



#' Plot to pdf
#' 
#' Convenience function to print plots as pdf document. 
#'
#' @param .data plot object
#' @param filename String specifying filename. Can also specify path. 
#' @param w width in inches
#' @param h height in inches
#'
#' @return plot
#' @export
#'
to_pdf<-function(.data, filename, w=7,h=7){
  
  pdf(  sub("(.*)(?<=/)(.*)", paste0("\\1",gsub(":","-",format(Sys.time(), "%Y%m%d_")), "\\2", ".pdf"), filename, perl=TRUE ), # File name
        width = w, height = h, # Width and height in inches
        bg = "white",          # Background color
        colormodel = "srgb",    # Color model (cmyk is required for most publications)
        onefile=T)          # Paper size
  
  # Creating a plot
  plot(.data)
  # Closing the graphical device
  dev.off() 
}

