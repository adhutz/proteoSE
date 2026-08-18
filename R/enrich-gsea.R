# GSEA / ssGSEA / KSEA preparation and GO-term reduction
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.


#' Convert a GCT object to a long format data frame
#'
#' This function takes a GCT object (or a file containing GCT data) and converts it to a long format data frame.
#'
#' @param gct A GCT object. If the file parameter is provided, this parameter will be ignored.
#' @param file A file path to a GCT file. If provided, the function will parse the GCT data from the file. Default is an empty string.
#' @return A long format data frame with columns for id.y, id.x, and metadata.
#' 
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom plyr colwise
#' @importFrom dplyr filter
#' 
#' @examples
#' \dontrun{
#' # Use an example GCT object from cmapR package
#' example_gct <- cmapR::example_gct
#' long_format_data <- gct_to_long(example_gct)
#' head(long_format_data)
#' }
#' @export
gct_to_long <- function(gct, file = ""){
  .require_pkg("cmapR", "read GCT files")

  if(!file == ""){
    gct <- cmapR::parse_gctx(file)
  }
  
  m <- cmapR::melt_gct(gct)
  
  m <- m %>% 
    tidyr::pivot_longer(grep(paste(unique(.$id.y), 
                                   collapse = "|"), 
                             colnames(.)), 
                        names_to = c("type", "set"), 
                        values_to = "val", 
                        names_pattern = paste0("^(.*)\\.(.*_[0-9]*)"), 
                        values_transform = list(val = as.character)) %>% 
    filter(set == id.y) %>% 
    tidyr::pivot_wider(names_from = "type", 
                       values_from = "val")
  
  m <- plyr::colwise(type.convert)(m)
  
}



#' Prepare Data for Single Sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' This function prepares the data from a SummarizedExperiment object for ssGSEA analysis. It removes duplicated rows based on the rowData sequence_window, creates a GCT object, and optionally writes the GCT object to a file.
#'
#' @param se A SummarizedExperiment object.
#' @param file A character string indicating the file name to save the GCT object. Defaults to an empty string, which means the GCT object will not be saved to a file.
#' @return A list with the following elements:
#' \describe{
#'   \item{"data"}{A GCT object with the unique rows from the input SummarizedExperiment object.}
#'   \item{"path"}{A character string indicating the path of the saved GCT file, if applicable.}
#'   \item{"duplicated"}{A SummarizedExperiment object containing the duplicated rows removed from the input object.}
#' }
#' @examples
#' \dontrun{
#' # Create a SummarizedExperiment object (se) here
#' # ...
#' result <- prep_ssgsea2(se)
#' result <- prep_ssgsea2(se, file = "output")
#' }
#' @export
prep_ssgsea2 <- function(se, file = ""){
  # cmapR defines the GCT S4 class used below; requireNamespace() loads it.
  .require_pkg("cmapR", "build GCT objects")

  rid <- gsub("^.{8}(.*).{8}$", 
              "\\1-p", 
              rowData(se)$sequence_window)
  
  dupl <- se[duplicated(rid),]
  uni <- se[!duplicated(rid),]
  
  temp <- new("GCT", 
              mat=assay(uni), 
              rdesc=as.data.frame(rowData(uni)), 
              cdesc=as.data.frame(colData(uni)), 
              rid = gsub("^.{8}(.*).{8}$", 
                         "\\1-p", 
                         rowData(uni)$sequence_window))
  
  if(!file == ""){
    cmapR::write_gct(temp, file)
    return(list("data" = temp, "path" = paste0(file, "_n", ncol(se), "x", nrow(uni), ".gct"), "duplicated" = dupl))
  } else{
    return(list("data" = temp, "path" = paste0(file, "_n", ncol(se), "x", nrow(uni), ".gct"), "duplicated" = dupl))
  }
  
}


#' Prepare KSEA input data from SummarizedExperiment object
#'
#' This function processes a SummarizedExperiment object and extracts relevant information
#' for the given contrast. The output is a data frame containing Protein, Gene, Peptide, Residue.Both,
#' p-value, and fold change (FC) columns, filtered to exclude rows with missing values.
#'
#' @param se A SummarizedExperiment object containing the input data.
#' @param contrast A character string specifying the contrast of interest.
#' @return A data frame with columns: Protein, Gene, Peptide, Residue.Both, p, and FC.
#' @importFrom rlang !! 
#' @importFrom ggplot2 sym
#' @import dplyr
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'se' is a SummarizedExperiment object and 'contrast' is a character string
#' # representing the contrast of interest:
#' prep_ksea_data <- prep_ksea(se, contrast)
#' }
prep_ksea <- function(se, contrast){
  diff_ <- paste0(contrast, "_diff")
  p_ <- paste0(contrast, "_p.val")
  
  rowData(se) %>% 
    as.data.frame() %>%
    select(proteins, 
           gene_names, 
           amino_acid, 
           position,
           contains(!! p_),
           contains(!! diff_)) %>%
    mutate(Residue.Both = paste0(amino_acid, position), 
           Peptide = "NULL",
           FC = 2^!!sym(diff_)) %>% 
    dplyr::rename("Protein" = "proteins", 
           "Gene" = "gene_names", 
           "p" = p_) %>% 
    select(Protein, 
           Gene, 
           Peptide, 
           Residue.Both, 
           p, 
           FC)  %>% 
    filter(!is.na(Gene)) %>%
    return()
}

#########################################################
#########################################################
#########################################################

#' Prepare ranked gene lists from a SummarizedExperiment for GSEA
#'
#' Extracts ranking statistics from a wide data frame derived from a
#' `SummarizedExperiment` and returns a named list of sorted numeric vectors,
#' suitable for GSEA workflows such as `clusterProfiler::compareCluster()`.
#'
#' @param se A `SummarizedExperiment` object.
#' @param contrasts Character vector of contrast names, or `"all"` to use all
#'   available contrasts.
#' @param stat_suffix Suffix identifying statistic columns in the wide data
#'   frame. Default is `"_diff"`.
#' @param gene_col Name of the gene identifier column. Default is `"name"`.
#' @param deduplicate_method Method to handle duplicated gene identifiers.
#'   One of:
#'   \itemize{
#'     \item `"first"`: keep first occurrence
#'     \item `"max_abs"`: keep value with largest absolute statistic
#'     \item `"max"`: keep largest statistic
#'     \item `"mean"`: average statistics per gene
#'   }
#'   Default is `"first"`.
#'
#' @return A named list of ranked numeric vectors.
#' @export
prep_gse <- function(se,
                     contrasts = "all",
                     stat_suffix = "_diff",
                     gene_col = "name",
                     deduplicate_method = c("first", "max_abs", "max", "mean")) {
  
  deduplicate_method <- match.arg(deduplicate_method)
  
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }
  
  df_wide <- get_df_wide(se)
  
  if (!gene_col %in% names(df_wide)) {
    stop("`gene_col` not found in wide data frame: ", gene_col)
  }
  
  stat_cols <- names(df_wide)[grepl(stat_suffix, names(df_wide), fixed = TRUE)]
  
  if (length(stat_cols) == 0) {
    stop("No columns found with suffix '", stat_suffix, "'.")
  }
  
  if (length(contrasts) == 1 && identical(contrasts, "all")) {
    selected_cols <- stat_cols
  } else {
    selected_cols <- paste0(contrasts, stat_suffix)
    
    missing <- selected_cols[!selected_cols %in% stat_cols]
    if (length(missing) > 0) {
      stop("The following contrasts were not found: ",
           paste(missing, collapse = ", "))
    }
  }
  
  result_list <- lapply(selected_cols, function(col) {
    
    # Remove NA rows
    valid_rows <- !is.na(df_wide[[gene_col]]) & !is.na(df_wide[[col]])
    genes <- df_wide[[gene_col]][valid_rows]
    stats <- df_wide[[col]][valid_rows]
    
    # Handle duplicates
    if (anyDuplicated(genes)) {
      
      if (deduplicate_method == "first") {
        
        keep <- !duplicated(genes)
        genes <- genes[keep]
        stats <- stats[keep]
        
      } else if (deduplicate_method == "max_abs") {
        
        ord <- order(abs(stats), decreasing = TRUE)
        genes <- genes[ord]
        stats <- stats[ord]
        
        keep <- !duplicated(genes)
        genes <- genes[keep]
        stats <- stats[keep]
        
      } else if (deduplicate_method == "max") {
        
        ord <- order(stats, decreasing = TRUE)
        genes <- genes[ord]
        stats <- stats[ord]
        
        keep <- !duplicated(genes)
        genes <- genes[keep]
        stats <- stats[keep]
        
      } else if (deduplicate_method == "mean") {
        
        df_tmp <- data.frame(gene = genes, stat = stats)
        
        df_agg <- aggregate(stat ~ gene, data = df_tmp, FUN = mean)
        
        genes <- df_agg$gene
        stats <- df_agg$stat
      }
    }
    
    ranked_vec <- stats
    names(ranked_vec) <- genes
    
    sort(ranked_vec, decreasing = TRUE)
  })
  
  # Clean names (remove suffix only at end)
  names(result_list) <- sub(paste0(stat_suffix, "$"), "", selected_cols)
  
  result_list
}


########################################################################
########################################################################
########################################################################
#' Pre-filter GSEA Results Across Clusters
#'
#' Applies a set of numerical and categorical filters to a
#' \pkg{clusterProfiler} GSEA result object and returns a data frame that
#' retains any GO term that passes the filter in \emph{at least one}
#' cluster.  Terms that pass in one cluster are therefore kept for
#' \emph{all} clusters in which they were identified, enabling
#' consistent cross-cluster visualisation.
#'
#' @param gsea_res A \pkg{clusterProfiler} result object.  Accepted classes
#'   are \code{compareClusterResult} (multi-cluster), objects with a
#'   \code{result} slot (single-comparison GSEA/ORA), or anything directly
#'   coercible to a data frame.  The underlying data frame must contain at
#'   least the columns \code{Cluster}, \code{ID}, \code{ONTOLOGY},
#'   \code{p.adjust}, \code{NES}, \code{core_enrichment}, and
#'   \code{setSize}.
#' @param ont Character vector. GO ontology branches to retain.  Valid
#'   values are \code{"BP"}, \code{"CC"}, \code{"MF"}, or \code{"all"}
#'   (default), which expands to all three.
#' @param comparisons Character vector. Cluster names to include.
#'   \code{"all"} (default) includes every cluster present in the data.
#' @param p.adj Numeric scalar. Adjusted p-value threshold.  Only terms
#'   with \code{p.adjust <= p.adj} in at least one cluster are retained.
#'   Default: \code{0.05}.
#' @param min_nes Numeric scalar. Minimum absolute Normalised Enrichment
#'   Score.  Default: \code{0} (no NES filter).
#' @param n_genes_cut Integer scalar. Minimum number of leading-edge genes
#'   (\code{core_enrichment}) required.  Default: \code{0} (no filter).
#' @param gene_ratio_cut Numeric scalar. Minimum ratio of leading-edge
#'   genes to gene-set size (\code{n_genes / setSize}).  Default: \code{0}
#'   (no filter).
#'
#' @details
#' ## Computed columns added to the output
#' \describe{
#'   \item{\code{n_genes}}{Number of leading-edge genes in
#'     \code{core_enrichment}, derived by counting \code{"/"} separators
#'     and adding 1.  Assumes no leading or trailing slash.}
#'   \item{\code{gene_ratio}}{\code{n_genes / setSize}.}
#'   \item{\code{gse_score}}{\eqn{-\log_{10}(\text{p.adjust}) \times
#'     |\text{NES}|}: a composite score balancing statistical significance
#'     with effect size, useful for ranking terms.}
#'   \item{\code{is_significant}}{Logical; \code{TRUE} if
#'     \code{p.adjust <= p.adj}.}
#'   \item{\code{survived_filter}}{Logical; \code{TRUE} if the term passed
#'     \emph{all} filters in that specific cluster, \code{FALSE} otherwise.}
#'   \item{\code{n_survived}}{Integer; number of clusters in which this GO
#'     term passed all filters.}
#' }
#'
#' ## Filtering logic
#' Filters are applied per cluster x GO term row.  A term is then retained
#' across \emph{all} clusters in which it appears if it passes the filter
#' in \emph{any} one cluster.  This means a term significant in cluster A
#' but not cluster B will still appear for cluster B (with
#' \code{survived_filter = FALSE}), preserving the ability to compare
#' enrichment patterns across the full set of clusters.
#'
#' @return A \code{\link[tibble]{tibble}} / data frame containing all rows
#'   from the input that belong to GO terms which passed the filter in at
#'   least one cluster.  All original columns are retained, plus the
#'   computed columns described above.
#'
#' @examples
#' \dontrun{
#' # Basic usage -- default filters
#' res <- gse_filt(gsea_res = my_gsea)
#'
#' # Stricter filtering, BP only, subset of clusters
#' res <- gse_filt(
#'   gsea_res       = my_gsea,
#'   ont            = "BP",
#'   comparisons    = c("Cluster1", "Cluster2"),
#'   p.adj          = 0.01,
#'   min_nes        = 1.5,
#'   n_genes_cut    = 5,
#'   gene_ratio_cut = 0.1
#' )
#'
#' # Inspect how many terms survived per cluster
#' res |>
#'   dplyr::filter(survived_filter) |>
#'   dplyr::count(Cluster)
#' }
#'
#' @importFrom dplyr mutate filter select pull left_join group_by ungroup if_else
#' @importFrom stringr str_count
#' @importFrom tidyr replace_na
#'
#' @export
gse_filt <- function(
    gsea_res,
    ont            = "all",
    comparisons    = "all",
    p.adj          = 0.05,
    min_nes        = 0,
    n_genes_cut    = 0,
    gene_ratio_cut = 0
) {
  
  # -------------------------------------------------------------------------
  # 1.  Coerce clusterProfiler result object to a plain data frame.
  #     compareClusterResult (multi-cluster) takes priority over the generic
  #     result slot; direct coercion is the fallback for other classes.
  # -------------------------------------------------------------------------
  if ("compareClusterResult" %in% methods::slotNames(gsea_res)) {
    df <- as.data.frame(gsea_res@compareClusterResult)
  } else if ("result" %in% methods::slotNames(gsea_res)) {
    df <- as.data.frame(gsea_res@result)
  } else {
    df <- as.data.frame(gsea_res)
  }
  
  # -------------------------------------------------------------------------
  # 2.  Expand sentinel values to their full sets
  # -------------------------------------------------------------------------
  
  # "all" ontologies -> all three branches
  if (identical(ont, "all")) ont <- c("BP", "CC", "MF")
  
  # Derive computed columns before resolving the comparisons sentinel so
  # Cluster is already a character when we unique() it below
  df <- df |>
    dplyr::mutate(
      Cluster = as.character(Cluster),
      
      # Count "/" separators in core_enrichment and add 1 to get gene count.
      # NOTE: assumes no leading/trailing slash in the core_enrichment string.
      n_genes = stringr::str_count(core_enrichment, "/") + 1,
      
      # Fraction of the gene set represented by leading-edge genes
      gene_ratio = n_genes / setSize,
      
      # Composite ranking score: penalises weak effect sizes and high p-values
      gse_score = -log10(p.adjust) * abs(NES),
      
      # Simple significance flag at the chosen threshold
      is_significant = dplyr::if_else(p.adjust <= p.adj, TRUE, FALSE)
    )
  
  # "all" comparisons -> every cluster present in the data
  if (identical(comparisons, "all")) {
    comparisons <- unique(df$Cluster)
  }
  
  # -------------------------------------------------------------------------
  # 3.  Build a per-cluster filter table: one row per (Cluster x ID) that
  #     passes every threshold.  This is joined back to the full data later.
  # -------------------------------------------------------------------------
  go_filtered <- df |>
    dplyr::filter(ONTOLOGY  %in% ont)          |>
    dplyr::filter(Cluster   %in% comparisons)  |>
    dplyr::filter(p.adjust  <= p.adj)          |>
    dplyr::filter(abs(NES)  >= min_nes)        |>
    dplyr::filter(n_genes   >= n_genes_cut)    |>
    dplyr::filter(gene_ratio >= gene_ratio_cut)|>
    dplyr::mutate(survived_filter = TRUE)      |>
    dplyr::select(Cluster, ID, survived_filter)
  
  # -------------------------------------------------------------------------
  # 4.  Report how many unique GO terms survived filtering
  # -------------------------------------------------------------------------
  n_input    <- dplyr::n_distinct(df$ID)
  n_survived <- dplyr::n_distinct(go_filtered$ID)
  
  message(
    "Input contained ", n_input, " unique GO terms.  ",
    "After filtering, ", n_survived, " GO terms remain.\n",
    "Any term significant in at least one cluster is retained for all clusters."
  )
  
  # -------------------------------------------------------------------------
  # 5.  Expand filtered IDs back to all clusters.
  #     Logic:
  #       - left_join brings in survived_filter = TRUE only for rows that
  #         passed; all other rows get NA.
  #       - group_by(ID) then any(survived_filter) keeps a GO term for every
  #         cluster if it passed in *any* cluster.
  #       - n_survived counts how many clusters each term passed in.
  #       - replace_na converts remaining NA -> FALSE for clean downstream use.
  # -------------------------------------------------------------------------
  res <- df |>
    dplyr::left_join(go_filtered, by = c("Cluster", "ID")) |>
    dplyr::group_by(ID) |>
    dplyr::filter(any(survived_filter == TRUE, na.rm = TRUE)) |>
    dplyr::mutate(
      n_survived     = sum(survived_filter == TRUE, na.rm = TRUE),
      survived_filter = tidyr::replace_na(survived_filter, FALSE)
    ) |>
    dplyr::ungroup()
  
  res
}

########################################################################
########################################################################
########################################################################
#' Reduce GO Term Redundancy by Removing Ancestral Terms
#'
#' Removes redundant (broad/ancestral) GO terms from a
#' \pkg{clusterProfiler} GSEA result, retaining the most specific
#' (leaf-like) terms.  Specificity is determined using a pre-built
#' parent-offspring table: if a GO term acts as an ancestor of other
#' terms present in the data across \emph{all} its occurrences, it is
#' removed.  The NES score plays no role in this reduction.
#'
#' @param gsea_res A \pkg{clusterProfiler} result object or data frame.
#'   Must contain at least columns \code{ID} and \code{is_significant}.
#'   If a \code{Cluster} column is present (i.e. a multi-cluster
#'   \code{compareClusterResult}), parent detection is performed within
#'   each cluster independently before aggregating across the full data
#'   set.
#' @param ont Character vector.  GO ontology branches to retain before
#'   reduction.  One or more of \code{"BP"}, \code{"CC"}, \code{"MF"},
#'   or \code{"all"} (default), which expands to all three.  Note: if
#'   \code{gsea_res} does not contain an \code{ONTOLOGY} column this
#'   filter is silently skipped.
#' @param terms A data frame with at least columns \code{ID},
#'   \code{parents}, and \code{n_children}, as produced by
#'   \code{\link{build_go_terms}}.  \code{parents} should be a
#'   pipe-separated (\code{"|"}) string of ancestor GO IDs, or
#'   \code{NA} for root terms.
#'
#' @details
#' ## Reduction algorithm
#' \enumerate{
#'   \item The GSEA result is joined to \code{terms} to bring in the
#'     \code{parents} column.
#'   \item Optionally filtered to the requested \code{ont} branches.
#'   \item Within each group (cluster x dataset, or dataset-wide if no
#'     \code{Cluster} column), the set of all parent IDs referenced by
#'     any term in the group is collected.
#'   \item Each term is flagged \code{is_parent_in_cluster = TRUE} if its
#'     own \code{ID} appears in that parent set -- i.e. it is an ancestor
#'     of at least one other term in the same group.
#'   \item Across the full data set (grouped by \code{ID}), the fraction
#'     of occurrences where the term was flagged as a parent is computed
#'     (\code{perc_is_parent_in_df}).
#'   \item Terms where \code{perc_is_parent_in_df == 1} -- i.e. the term
#'     was \emph{always} a parent, never a leaf -- are removed.
#' }
#'
#' ## What is and is not removed
#' A term is only removed if it is a parent in \emph{every} cluster it
#' appears in.  A term that is a parent in some clusters but a leaf in
#' others is retained, preserving cross-cluster comparability.
#'
#' ## NES not considered
#' The reduction is purely topology-based.  Two terms with the same
#' \code{core_enrichment} genes but different specificity will always
#' resolve to the more specific (smaller \code{n_children}) term,
#' regardless of NES magnitude.
#'
#' @return The input data frame with redundant ancestral GO terms removed
#'   and the following columns added:
#' \describe{
#'   \item{\code{parents_in_group}}{Character. Pipe-separated set of all
#'     parent IDs observed in the same group (cluster or full data set).}
#'   \item{\code{is_parent_in_cluster}}{Logical. \code{TRUE} if this
#'     term's ID appears as a parent of another term in its group.}
#'   \item{\code{n_term_in_df}}{Integer. Number of times this GO term
#'     appears across the full data set.}
#'   \item{\code{n_is_parent_in_df}}{Integer. Number of times this term
#'     was flagged as a parent.}
#'   \item{\code{perc_is_parent_in_df}}{Numeric. Fraction of occurrences
#'     where the term was a parent (\code{n_is_parent_in_df /
#'     n_term_in_df}).}
#'   \item{\code{n_is_significant_in_df}}{Integer. Number of clusters in
#'     which this term was significant.}
#' }
#'
#' @examples
#' \dontrun{
#' terms <- build_go_terms()
#'
#' reduced <- reduce_terms(
#'   gsea_res = my_gsea_filtered,
#'   ont      = "BP",
#'   terms    = terms
#' )
#'
#' # How many terms were removed?
#' nrow(my_gsea_filtered) - nrow(reduced)
#' }
#'
#' @importFrom dplyr mutate filter group_by ungroup pull left_join n
#' @importFrom stringr str_split
#'
#' @export
reduce_terms <- function(gsea_res, ont = "all", terms) {
  
  # -------------------------------------------------------------------------
  # 1.  Expand sentinel and coerce input to data frame
  # -------------------------------------------------------------------------
  if (identical(ont, "all")) ont <- c("BP", "CC", "MF")
  
  df <- as.data.frame(gsea_res) |>
    dplyr::left_join(terms, by = "ID")
  
  # -------------------------------------------------------------------------
  # 2.  Optionally filter to requested ontology branches
  #     (silently skipped if the ONTOLOGY column is absent)
  # -------------------------------------------------------------------------
  if ("ONTOLOGY" %in% colnames(df)) {
    df <- dplyr::filter(df, ONTOLOGY %in% ont)
  }
  
  # -------------------------------------------------------------------------
  # 3.  Report input size
  # -------------------------------------------------------------------------
  n_rows_in         <- nrow(df)
  n_unique_terms_in <- dplyr::n_distinct(df$ID)
  message(n_rows_in, " input rows with ", n_unique_terms_in, " unique input terms.")
  
  # -------------------------------------------------------------------------
  # 4.  Core reduction logic
  #
  #     The Cluster / no-Cluster branches only differ in the first grouping
  #     level, so we handle that with a single conditional group_by rather
  #     than duplicating the entire pipeline.
  #
  #     Step A -- within-group parent detection:
  #       For each group, collect every parent ID referenced by any term
  #       (by splitting the pipe-separated `parents` strings), then flag
  #       each row whose own ID appears in that set.
  #
  #     Step B -- dataset-level aggregation:
  #       Group by ID and compute how often each term was flagged as a
  #       parent across its total occurrences.
  #
  #     Step C -- filter:
  #       Remove terms that were *always* a parent (ratio == 1), i.e. never
  #       a leaf node in any group.
  # -------------------------------------------------------------------------
  has_cluster <- "Cluster" %in% colnames(df)
  
  # Step A: within-group parent flagging
  df <- {
    if (has_cluster) dplyr::group_by(df, Cluster) else df
  } |>
    dplyr::mutate(
      # Collect all parent IDs cited by any term in this group into one vector
      parents_in_group = paste(
        unique(unlist(stringr::str_split(na.omit(parents), "\\|"))),
        collapse = "|"
      ),
      # Flag: is this term's ID an ancestor of another term in this group?
      is_parent_in_cluster = ID %in% unlist(stringr::str_split(parents_in_group, "\\|"))
    ) |>
    dplyr::ungroup()
  
  # Step B: dataset-level parent frequency per GO term
  df <- df |>
    dplyr::group_by(ID) |>
    dplyr::mutate(
      n_term_in_df           = dplyr::n(),
      n_is_parent_in_df      = sum(is_parent_in_cluster),
      perc_is_parent_in_df   = n_is_parent_in_df / n_term_in_df,
      n_is_significant_in_df = sum(is_significant)
    ) |>
    dplyr::ungroup()
  
  # Step C: remove terms that were a parent in every one of their occurrences
  df <- dplyr::filter(df, perc_is_parent_in_df != 1)
  
  # -------------------------------------------------------------------------
  # 5.  Report output size
  # -------------------------------------------------------------------------
  n_rows_out         <- nrow(df)
  n_unique_terms_out <- dplyr::n_distinct(df$ID)
  
  message(
    "Reduced to ", n_rows_out, " output rows with ",
    n_unique_terms_out, " unique output terms.\n",
    "Reduction kept the most specific GO term (leaf over branch); ",
    "NES magnitude was not considered."
  )
  
  df
}

########################################################################
########################################################################
########################################################################
#' Select a Diverse and Balanced Set of GO Terms
#'
#' Chooses a representative, non-redundant subset of GO terms from a
#' filtered GSEA result by stratifying terms into \code{n_quantiles}
#' bins of a chosen diversity variable, picking the top
#' \code{n_per_stratum} terms per bin, and then optionally upsampling
#' under-represented ontologies so that BP/CC/MF contribute equally.
#' Specific terms can be pinned to the output via \code{keep_keywords}.
#'
#' @param gsea_res A \pkg{clusterProfiler} result object or data frame.
#'   Must contain columns \code{ONTOLOGY}, \code{ID}, \code{Description},
#'   and the column named by \code{which_diverse}.  A \code{Cluster}
#'   column is optional; its presence activates multi-cluster mode.
#' @param ont Character vector.  Ontology branches to retain.  One or
#'   more of \code{"BP"}, \code{"CC"}, \code{"MF"}, or \code{"all"}
#'   (default), which expands to all three.
#' @param terms A data frame as produced by \code{\link{build_go_terms}},
#'   containing at least \code{ID} and \code{n_children}.  Required only
#'   if \code{n_children} is not already a column in \code{gsea_res}.
#' @param which_diverse Character scalar.  Column to use for diversity
#'   stratification.  One of \code{"NES"} (default), \code{"rank"},
#'   \code{"gene_ratio"}, or \code{"n_children"}.  Terms are split into
#'   \code{n_quantiles} bins by the absolute value of this column so
#'   that the selection spans the full distribution rather than being
#'   dominated by only the top-scoring terms.
#' @param n_per_stratum Integer scalar.  Number of terms to select per
#'   stratum, where one stratum = one ontology x one cluster (if present)
#'   x one quantile bin.  The approximate total number of unique IDs
#'   selected before balancing is therefore:
#'   \deqn{n\_ontologies \times n\_clusters \times n\_quantiles \times
#'   n\_per\_stratum}
#'   Default: \code{3}.
#' @param n_quantiles Integer scalar.  Number of quantile bins used to
#'   stratify \code{which_diverse}.  Higher values create finer strata
#'   (more diverse output) but also increase the total number of terms
#'   selected.  Default: \code{3} (tertiles).
#' @param keep_keywords Character vector or \code{NULL}.  Terms whose
#'   \code{Description} matches any of these strings (case-insensitive
#'   regex) are pinned to the top of the selection priority within their
#'   stratum, ensuring they consume selection slots before non-keyword
#'   terms.  \code{NULL} (default) disables keyword pinning.
#' @param keywords_only Logical.  If \code{TRUE}, only keyword-matched
#'   terms are considered for selection.  Requires \code{keep_keywords}
#'   to be non-\code{NULL}; a warning is issued and the flag is ignored
#'   otherwise.  Default: \code{FALSE}.
#' @param balance_ont Logical.  If \code{TRUE} (default), ontologies with
#'   fewer unique IDs than the target are upsampled (guided by
#'   \code{which_diverse}, not randomly) until all present ontologies
#'   reach the same count.  If \code{FALSE}, the raw per-stratum
#'   selection is returned without balancing.
#'
#' @details
#' ## Selection algorithm
#' \enumerate{
#'   \item Terms are filtered to the requested \code{ont} branches.
#'   \item Each term's \code{Description} is matched against
#'     \code{keep_keywords}; matched terms are flagged \code{force_keep
#'     = TRUE} (best row per ID x keyword only).
#'   \item Terms are split into \code{n_quantiles} bins by
#'     \code{abs(which_diverse)} within each ontology x cluster group.
#'   \item Within each stratum, \code{force_keep} terms sort to the top,
#'     then by \code{abs(which_diverse)}.  Up to \code{n_per_stratum}
#'     unique IDs are taken per stratum.
#'   \item If \code{balance_ont = TRUE}, under-represented ontologies are
#'     topped up from the remaining (not-yet-selected) terms in priority
#'     order until all ontologies reach the target count (see below).
#' }
#'
#' ## Upsampling target
#' The per-ontology target unique-ID count is:
#' \deqn{n\_ont\_present \times n\_clusters \times n\_quantiles \times
#' n\_per\_stratum}
#' divided equally across the ontologies present.  This means passing
#' \code{ont = "BP"} targets the same total as passing all three, but
#' assigns it entirely to BP -- so you will see more BP terms than when
#' all three are active.  The message printed by the function shows the
#' exact arithmetic so the output size is never a surprise.
#'
#' ## Keyword pinning
#' Within each stratum, \code{force_keep} rows are sorted before all
#' others, so keyword-matched terms consume selection slots first.  Only
#' the single strongest row per \code{ID x keyword} combination
#' (ordered by \code{keyword_match} index then \code{abs(which_diverse)})
#' is marked \code{force_keep = TRUE}, preventing one ID from occupying
#' multiple slots.
#'
#' @return A data frame of selected rows from \code{gsea_res}, with the
#'   following additional columns:
#' \describe{
#'   \item{\code{quantile}}{Integer 1-\code{n_quantiles}. Bin of
#'     \code{abs(which_diverse)} within each ontology x cluster group.}
#'   \item{\code{keyword_match}}{Integer or \code{NA}. Index of the first
#'     matching element of \code{keep_keywords}, or \code{NA} if none.}
#'   \item{\code{matched_keyword}}{Character or \code{NA}. The matched
#'     keyword string, or \code{NA}.}
#'   \item{\code{best_row_for_keyword_match}}{Logical. \code{TRUE} for
#'     the single strongest row per ID x keyword group.}
#'   \item{\code{force_keep}}{Logical. \code{TRUE} if this row matched a
#'     keyword and is the best representative for that ID.}
#' }
#'
#' @examples
#' \dontrun{
#' # Diverse NES-stratified terms, balanced across all three ontologies
#' selected <- select_diverse_terms(
#'   gsea_res = my_gsea_reduced,
#'   terms    = go_terms
#' )
#'
#' # Finer stratification: 5 bins, 2 terms per stratum, BP only
#' selected <- select_diverse_terms(
#'   gsea_res      = my_gsea_reduced,
#'   terms         = go_terms,
#'   ont           = "BP",
#'   n_quantiles   = 5,
#'   n_per_stratum = 2,
#'   balance_ont   = FALSE
#' )
#'
#' # Pin metabolic terms to the top across all ontologies
#' selected <- select_diverse_terms(
#'   gsea_res      = my_gsea_reduced,
#'   terms         = go_terms,
#'   keep_keywords = c("metabol", "biosynthes"),
#'   balance_ont   = TRUE
#' )
#'
#' # Return only keyword-matched terms
#' selected <- select_diverse_terms(
#'   gsea_res      = my_gsea_reduced,
#'   terms         = go_terms,
#'   keep_keywords = c("lysosome", "autophagy"),
#'   keywords_only = TRUE
#' )
#' }
#'
#' @importFrom dplyr mutate filter group_by ungroup arrange distinct slice_head pull semi_join anti_join bind_rows inner_join group_modify n_distinct count across all_of if_else row_number ntile
#' @importFrom tibble tibble
#'
#' @export
select_diverse_terms <- function(
    gsea_res,
    ont           = "all",
    terms         = NULL,
    which_diverse = c("NES", "rank", "gene_ratio", "n_children"),
    n_per_stratum = 3,
    n_quantiles   = 3,
    keep_keywords = NULL,
    keywords_only = FALSE,
    balance_ont   = TRUE
) {
  
  which_diverse <- match.arg(which_diverse)
  
  # -------------------------------------------------------------------------
  # 1.  Coerce to data frame and optionally join GO term metadata.
  #     The column check for which_diverse is intentionally deferred to
  #     step 2, after the join, so that columns arriving via `terms`
  #     (e.g. n_children) are present before validation.
  # -------------------------------------------------------------------------
  df <- as.data.frame(gsea_res)
  
  if (!"n_children" %in% colnames(df)) {
    if (is.null(terms)) {
      stop(
        "Column 'n_children' not found and 'terms' not provided. ",
        "Please supply 'terms' or run reduce_terms() first.",
        call. = FALSE
      )
    }
    df <- dplyr::left_join(df, terms, by = "ID")
    message("Joined with terms data frame.")
  } else {
    message("Terms already joined; skipping merge.")
  }
  
  # -------------------------------------------------------------------------
  # 2.  Validate which_diverse (after potential join)
  # -------------------------------------------------------------------------
  if (!which_diverse %in% colnames(df)) {
    stop("Column '", which_diverse, "' not found in gsea_res.", call. = FALSE)
  }
  
  # -------------------------------------------------------------------------
  # 3.  Expand sentinel and filter to requested ontology branches
  # -------------------------------------------------------------------------
  if (identical(ont, "all")) ont <- c("BP", "CC", "MF")
  
  df <- dplyr::filter(df, ONTOLOGY %in% ont)
  
  if (nrow(df) == 0) {
    message("No terms passed the ontology filter.")
    return(df)
  }
  
  # -------------------------------------------------------------------------
  # 4.  Keyword pinning
  #
  #     For each Description, find the index of the first matching keyword
  #     (preserving keep_keywords order, so earlier keywords take priority).
  #     Only the single strongest row per ID x keyword combination is
  #     flagged force_keep = TRUE, preventing one ID from occupying
  #     multiple selection slots within the same stratum.
  # -------------------------------------------------------------------------
  if (!is.null(keep_keywords) && length(keep_keywords) > 0) {
    
    df <- df |>
      dplyr::mutate(
        # Index of first matching keyword; NA if no match
        keyword_match = sapply(Description, function(desc) {
          idx <- which(vapply(
            keep_keywords,
            function(k) grepl(k, desc, ignore.case = TRUE),
            logical(1)
          ))
          if (length(idx) == 0L) NA_integer_ else min(idx)
        }),
        matched_keyword = dplyr::if_else(
          !is.na(keyword_match),
          keep_keywords[keyword_match],
          NA_character_
        )
      ) |>
      # Mark the single best row per ID x keyword (keyword order first,
      # then strongest which_diverse) to prevent duplicate slot consumption
      dplyr::group_by(ONTOLOGY, ID) |>
      dplyr::arrange(
        keyword_match,
        dplyr::desc(abs(.data[[which_diverse]])),
        .by_group = TRUE
      ) |>
      dplyr::mutate(best_row_for_keyword_match = dplyr::row_number() == 1L) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        force_keep = !is.na(keyword_match) & best_row_for_keyword_match
      )
    
  } else {
    # No keywords: fill columns with neutral defaults for a consistent schema
    df <- dplyr::mutate(
      df,
      keyword_match              = NA_integer_,
      matched_keyword            = NA_character_,
      best_row_for_keyword_match = FALSE,
      force_keep                 = FALSE
    )
  }
  
  # -------------------------------------------------------------------------
  # 5.  Optionally restrict to keyword-matched terms only
  # -------------------------------------------------------------------------
  if (keywords_only) {
    if (is.null(keep_keywords) || length(keep_keywords) == 0L) {
      warning(
        "keywords_only = TRUE but no keep_keywords supplied; ignoring.",
        call. = FALSE
      )
    } else {
      message("keywords_only = TRUE: restricting to keyword-matched terms.")
      df <- dplyr::filter(df, !is.na(keyword_match))
      if (nrow(df) == 0L) {
        message("No keyword-matched terms found.")
        return(df)
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # 6.  Diversity stratification and per-stratum selection
  #
  #     Terms are binned into n_quantiles bins by abs(which_diverse).
  #     Up to n_per_stratum unique IDs are taken per stratum, with
  #     force_keep terms prioritised within each slice.
  #
  #     Stratum definition:
  #       - with Cluster column:    ontology x cluster x quantile bin
  #       - without Cluster column: ontology x quantile bin
  #
  #     Both cases are handled by a single pipeline using dynamic grouping
  #     via dplyr::across(all_of(...)), avoiding duplicated code blocks.
  # -------------------------------------------------------------------------
  has_cluster <- "Cluster" %in% colnames(df)
  
  group_vars_bin       <- if (has_cluster) c("ONTOLOGY", "Cluster")           else "ONTOLOGY"
  group_vars_selection <- if (has_cluster) c("ONTOLOGY", "Cluster", "quantile") else c("ONTOLOGY", "quantile")
  
  # Assign quantile bins within each ontology x [cluster] group
  df <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_bin))) |>
    dplyr::mutate(quantile = dplyr::ntile(abs(.data[[which_diverse]]), n_quantiles)) |>
    dplyr::ungroup()
  
  # Select top n_per_stratum unique IDs per stratum
  selected_ids <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_selection))) |>
    dplyr::arrange(
      dplyr::desc(force_keep),
      dplyr::desc(abs(.data[[which_diverse]])),
      .by_group = TRUE
    ) |>
    dplyr::distinct(ID, .keep_all = TRUE) |>   # one slot per unique ID per stratum
    dplyr::slice_head(n = n_per_stratum) |>
    dplyr::ungroup() |>
    dplyr::pull(ID) |>
    unique()
  
  df_selected <- dplyr::filter(df, ID %in% selected_ids)
  
  # -------------------------------------------------------------------------
  # 7.  Guided upsampling to balance ontology representation
  #     (skipped entirely when balance_ont = FALSE)
  #
  #     Target unique IDs per ontology:
  #       n_ont_present x n_clusters x n_quantiles x n_per_stratum
  #       / n_ont_present
  #     = n_clusters x n_quantiles x n_per_stratum
  #
  #     Under-represented ontologies are topped up in priority order
  #     (force_keep first, then abs(which_diverse)) from the remaining
  #     not-yet-selected terms -- never randomly.
  # -------------------------------------------------------------------------
  if (balance_ont) {
    
    n_clusters    <- if (has_cluster) dplyr::n_distinct(df$Cluster) else 1L
    n_ont_present <- dplyr::n_distinct(df$ONTOLOGY)
    target_n      <- n_clusters * n_quantiles * n_per_stratum
    
    message(
      "Upsampling target: ", n_ont_present, " ontolog",
      if (n_ont_present == 1L) "y" else "ies", " \u00d7 ",
      n_clusters,  " cluster",  if (n_clusters  == 1L) "" else "s", " \u00d7 ",
      n_quantiles, " quantile", if (n_quantiles  == 1L) "" else "s", " \u00d7 ",
      n_per_stratum, " per stratum = ",
      target_n, " unique IDs per ontology."
    )
    
    sel_unique <- dplyr::distinct(df_selected, ONTOLOGY, ID)
    
    ont_counts <- dplyr::count(sel_unique, ONTOLOGY, name = "n_selected_unique")
    
    need_tbl <- dplyr::mutate(
      ont_counts,
      n_needed = pmax(0L, target_n - n_selected_unique)
    )
    
    # Candidate pool: IDs not yet selected, one row per ID, priority-ordered
    remaining_unique <- df |>
      dplyr::anti_join(sel_unique, by = c("ONTOLOGY", "ID")) |>
      dplyr::arrange(
        ONTOLOGY,
        dplyr::desc(force_keep),
        dplyr::desc(abs(.data[[which_diverse]]))
      ) |>
      dplyr::distinct(ONTOLOGY, ID, .keep_all = TRUE)
    
    # Take exactly n_needed additional IDs per ontology
    guided_add_ids <- remaining_unique |>
      dplyr::inner_join(need_tbl, by = "ONTOLOGY") |>
      dplyr::group_by(ONTOLOGY) |>
      dplyr::group_modify(~ {
        n_to_take <- unique(.x$n_needed)
        if (is.na(n_to_take) || n_to_take <= 0L)
          return(tibble::tibble(ID = character()))
        dplyr::slice_head(.x, n = min(n_to_take, nrow(.x))) |>
          dplyr::select(ID)
      }) |>
      dplyr::ungroup()
    
    # Combine originally selected + upsampled IDs, then filter df to match
    final_ids <- dplyr::bind_rows(sel_unique, guided_add_ids) |>
      dplyr::distinct(ONTOLOGY, ID)
    
    df_out <- dplyr::semi_join(df, final_ids, by = c("ONTOLOGY", "ID"))
    
  } else {
    df_out <- df_selected
  }
  
  # -------------------------------------------------------------------------
  # 8.  Report
  # -------------------------------------------------------------------------
  report_type <- if (balance_ont) "balanced" else "unbalanced"
  n_forced    <- sum(df_out$force_keep, na.rm = TRUE)
  
  id_counts <- df_out |>
    dplyr::distinct(ONTOLOGY, ID) |>
    dplyr::count(ONTOLOGY, name = "n_unique_ids")
  
  message(
    "Selected ", nrow(df_out), " rows total (", report_type, ").  ",
    n_forced, " force-kept rows.  ",
    "Unique IDs by ontology: ",
    paste(
      sprintf("%s=%d", id_counts$ONTOLOGY, id_counts$n_unique_ids),
      collapse = ", "
    ), "."
  )
  
  if (keywords_only)
    message("Note: only keyword-matched terms were considered for selection.")
  
  df_out
}

