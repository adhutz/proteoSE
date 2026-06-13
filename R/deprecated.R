# Deprecated function aliases
#
# Part of the 'sev' package. When a function is renamed or several near-duplicate
# functions are consolidated, the old exported names are kept here as thin
# wrappers that still work but emit .Deprecated(). Each wrapper forwards to the
# canonical function with arguments chosen to reproduce the old behaviour, so
# existing user scripts keep running for one release cycle.


# ---------------------------------------------------------------------------
# Volcano plots -> plot_volcano()
#
# plot_volcano() is the consolidated superset of the three former volcano
# functions. The aliases below map each old function's arguments onto it:
#   - se_volcano() / lars_volcano() read significance from the precomputed
#     <contrast>_significant column (sig_from_column = TRUE) and labelled ALL
#     significant points (top_n_sign = Inf); they returned the ggplot object
#     directly, so the aliases return plot_out rather than the full list.
#   - volcano() was simply renamed; it forwards verbatim and returns the list.
# ---------------------------------------------------------------------------

#' Deprecated volcano plot functions
#'
#' These functions are deprecated. Use \code{\link{plot_volcano}} instead, which
#' consolidates \code{se_volcano}, \code{lars_volcano}, and \code{volcano} into a
#' single function. The wrappers below forward to \code{plot_volcano} with
#' arguments that reproduce the original output.
#'
#' @inheritParams plot_volcano
#' @param ... Passed on to \code{\link{plot_volcano}}.
#' @return For \code{se_volcano} and \code{lars_volcano}, a \code{ggplot} object
#'   (the \code{plot_out} element). For \code{volcano}, the full list returned by
#'   \code{\link{plot_volcano}}.
#' @name sev-deprecated-volcano
#' @keywords internal
NULL

#' @rdname sev-deprecated-volcano
#' @export
se_volcano <- function(se, contrast, id_col = "gene_names", target_names = c(""),
                       label_sign = TRUE, label_targets = TRUE, max.overlaps = Inf) {
  .Deprecated("plot_volcano")
  plot_volcano(
    se, contrast,
    id_col          = id_col,
    target_names    = target_names,
    label_sign      = label_sign,
    label_targets   = label_targets,
    max.overlaps    = max.overlaps,
    sig_from_column = TRUE,
    top_n_sign      = Inf,
    fill_colors     = c(Significant = "firebrick", Target = "darkorange", None = "grey90"),
    shape_values    = c(`FALSE` = 21, `TRUE` = 23),
    label_colors    = c(Target = "darkorange", Top = "firebrick")
  )$plot_out
}

#' @rdname sev-deprecated-volcano
#' @export
lars_volcano <- function(se, contrast, id_col = "gene_names", target_names = c(""),
                         label_sign = TRUE, label_targets = TRUE, max.overlaps = 15,
                         labelsize = 2) {
  .Deprecated("plot_volcano")
  plot_volcano(
    se, contrast,
    id_col          = id_col,
    target_names    = target_names,
    label_sign      = label_sign,
    label_targets   = label_targets,
    max.overlaps    = max.overlaps,
    labelsize       = labelsize,
    sig_from_column = TRUE,
    top_n_sign      = Inf,
    fill_colors     = c(Significant = "#F8756D", Target = "#7570b3", None = "grey90"),
    shape_values    = c(`FALSE` = 21, `TRUE` = 23),
    label_colors    = c(Target = "#7570b3", Top = "black")
  )$plot_out
}

#' @rdname sev-deprecated-volcano
#' @export
volcano <- function(...) {
  .Deprecated("plot_volcano")
  plot_volcano(...)
}


# ---------------------------------------------------------------------------
# PTM report -> PhosphoExperiment -> read_ptm_to_ppe()
#
# read_ptm_to_ppe(source = ...) is the consolidated entry point.
#   - maxq_to_ppe() and ptm_to_ppe() forward to it with no behaviour change
#     (their former bodies are now the internal read_ptm_maxquant() /
#     read_ptm_spectronaut() engines).
#   - sn_ptm_to_ppe() is the older Spectronaut converter; its body is preserved
#     verbatim below because it uses a different intermediate schema, so it
#     cannot be reproduced exactly by the current engine.
# ---------------------------------------------------------------------------

#' Deprecated PTM-to-PhosphoExperiment converters
#'
#' These functions are deprecated. Use \code{\link{read_ptm_to_ppe}} instead.
#'
#' @param ... Passed on to \code{\link{read_ptm_to_ppe}}.
#' @return A \code{\link[PhosR]{PhosphoExperiment}} object (see
#'   \code{\link{read_ptm_to_ppe}}).
#' @name sev-deprecated-ptm
#' @keywords internal
NULL

#' @rdname sev-deprecated-ptm
#' @export
maxq_to_ppe <- function(...) {
  .Deprecated("read_ptm_to_ppe(source = \"maxquant\")")
  read_ptm_to_ppe(source = "maxquant", ...)
}

#' @rdname sev-deprecated-ptm
#' @export
ptm_to_ppe <- function(...) {
  .Deprecated("read_ptm_to_ppe(source = \"spectronaut\")")
  read_ptm_to_ppe(source = "spectronaut", ...)
}

#' Convert Spectronaut PTM Report to PhosphoExperiment (deprecated)
#'
#' @description
#' \strong{Deprecated.} Use \code{\link{read_ptm_to_ppe}} with
#' \code{source = "spectronaut"} instead. This older Spectronaut converter is
#' retained unchanged because it builds a different intermediate schema (notably
#' the \code{use_id = "collapse_key_gene"} identifier) than the current engine.
#'
#' This function reads in a long-format Spectronaut PTM report (Moritz style),
#' processes it, and outputs a PhosphoExperiment object, a long dataframe, a wide dataframe,
#' or a list of all three depending on `output_type`.
#'
#' @param file Path to Spectronaut export file (TSV) or a data.frame
#' @param sep Separator between condition and replicate in sample labels. Default: "_rep_"
#' @param filt Character vector of columns used to filter out contaminants. Default: c("reverse", "potential_contaminant")
#' @param experimental_design Optional data.frame with columns: label, sample, condition, replicate
#' @param output_type Output type: one of "long_df", "wide_df", "ppe", or "all". Default: "ppe"
#' @param assay_name Name of the assay in the final PhosphoExperiment object. Default: "orig_assay"
#' @param mod_type Character vector of modification types to retain (e.g., "Phospho (STY)"). Default: "Phospho (STY)"
#' @param site_prob numeric vector for modifications probabilities to retain (e.g., 0.75). Default: 0.75
#' @param use_id Character vector to choose for your "ID" in all subsequent analysis and plots default: "collapse_key_gene")
#'
#' @return Depending on `output_type`, returns:
#'   - `"ppe"`: a PhosphoExperiment object
#'   - `"long_df"`: a long-format cleaned tibble
#'   - `"wide_df"`: a wide-format tibble after pivot
#'   - `"all"`: a named list with all three components
#'
#' @importFrom vroom vroom
#' @importFrom dplyr %>% mutate group_by ungroup filter rename_all rename_with select if_all
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom janitor make_clean_names clean_names
#' @importFrom DEP2 make_unique
#' @importFrom PhosR PhosphoExperiment
#' @importFrom purrr when
#' @keywords internal
#' @export
sn_ptm_to_ppe <- function(file, sep = "_rep_", filt = c("reverse", "potential_contaminant"),
                          experimental_design = NA, output_type = "ppe",
                          assay_name = "orig_assay", mod_type = "Phospho (STY)",
                          site_prob = 0.75, use_id = "collapse_key_gene") {
  .Deprecated("read_ptm_to_ppe(source = \"spectronaut\")")

  # Load data from file or use supplied data frame
  if (is.character(file)) {
    data <- vroom::vroom(file, delim = "\t", col_names = TRUE, guess_max = 30000,
                         .name_repair = janitor::make_clean_names)
  } else {
    data <- file %>% janitor::clean_names()
  }
  
  #####read in data (spectronaut output from Moritz (long format) - this is in contrast to wide format report style from Fatih Demir)
  #long_df <- vroom::vroom(file, delim = "\t", col_names = TRUE, guess_max = 30000, .name_repair = janitor::make_clean_names) %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
  long_df <- data %>% #https://www.rdocumentation.org/packages/readxl/versions/1.3.1/topics/cell-specification
    
    dplyr::rename_all(tolower) %>% #turn all column names to lower case (makes it easier for later code writing)
    janitor::clean_names() %>% #make column names clean and unique (makes later coding easier)
    dplyr::rename_all(.funs = list(~gsub("^pg_", "", .))) %>%  # pg probably stands for protein group, just remove it from column names
    dplyr::rename_all(.funs = list(~gsub("^r_", "", .))) %>% 
    dplyr::rename_all(.funs = list(~gsub("^ptm_", "", .))) %>% 
    dplyr::mutate(sample = paste0(condition, "_rep_", replicate)) %>%  #paste together a sample column
    dplyr::group_by(collapse_key) %>% 
    dplyr::mutate(max_prob_across_exp = max(site_probability, na.rm = TRUE)) %>% #calculate maximum site probability by phosphosite - I think this is what maxquant does 
    dplyr::ungroup() %>% 
    dplyr::filter(site_probability >= site_prob) %>% #no need really, because spectrnoaut does this by default
    #dplyr::filter(modification_title %in% mod_type) %>% #kepp only phosphosite (gets rid of  carbamidomethylation, acetylation, oxidation)
    purrr::when(!is.null(mod_type) && !all(is.na(mod_type)) ~ dplyr::filter(., modification_title %in% mod_type),~ .) %>% 
    dplyr::mutate(collapse_key_simplified = stringr::str_replace(collapse_key, "_M\\d+$", "")) %>%  #this creates a column with the site sans the _M1, _M2, _M3, indicating Multiplicity
    { if (!"reverse" %in% names(.)) mutate(., reverse = NA) else . } %>% #check if "reverse" col exists, if not, create it
    { if (!"potential_contaminant" %in% names(.)) mutate(., potential_contaminant = NA) else . } #check if "potential_contaminant" col exists, if not, create it
  
  
  #next bit of code is to make a new column (collapse_key_gene) where we basically 
  #replace uniprot ids with gene ids in the collapsed_key column (if there is no gene id, then the uniprotid will be used)
  # Step 1: Extract and split once
  protein_ids <- stringr::str_extract(long_df$collapse_key, "^[^_]+")
  protein_lists <- strsplit(long_df$protein_groups, ";")
  gene_lists <- strsplit(long_df$genes, ";")
  
  # Step 2: Match and extract gene_id (NA if no match)
  gene_ids <- vapply(
    seq_along(protein_ids),
    function(i) {
      pid <- protein_ids[i]
      plist <- protein_lists[[i]]
      glist <- gene_lists[[i]]
      idx <- match(pid, plist)
      if (!is.na(idx) && idx <= length(glist)) glist[idx] else NA_character_
    },
    character(1)
  )
  
  # Step 3: Final ID (gene if available, else protein)
  final_ids <- ifelse(!is.na(gene_ids), gene_ids, protein_ids)
  
  # Step 4: Replace in collapse_key
  collapse_key_genes <- stringr::str_replace(long_df$collapse_key, protein_ids, final_ids)
  
  # Step 5: Add everything into the final dataframe
  df_final <- long_df %>%
    dplyr::mutate(
      protein_id = protein_ids,
      gene_id = gene_ids,
      final_id = final_ids,
      collapse_key_gene = collapse_key_genes
    )
  
  #pivot wider and keep some ID cols - maybe revisit this in case we need the #commented lines 
  wide_df <- df_final %>% 
    tidyr::pivot_wider(
      id_cols = c(genes,
                  protein_id,
                  gene_id,
                  final_id,
                  protein_groups,
                  protein_accessions,
                  uni_prot_ids,
                  protein_names,
                  protein_descriptions,
                  organisms,
                  molecular_weight,
                  nr_of_stripped_sequences_measured,
                  nr_of_stripped_sequences_identified_experiment_wide,
                  # nr_of_collapsed_peptides,# add these to the pivot_wider command because it is sample specific. if you keep it here it creates extra rows compared to fatih's export (in Fatih's wide report, each sample gets a column with values)
                  collapse_key,
                  collapse_key_simplified,
                  collapse_key_gene,
                  flanking_region,
                  #group # add these to the pivot_wider command because it is sample specific. if you keep it here it creates extra rows compared to fatih's export (in Fatih's wide report, each sample gets a column with values)
                  site_aa,
                  site_location,
                  #site_probability, # add these to the pivot_wider command because it is sample specific. if you keep it here it creates extra rows compared to fatih's export (in Fatih's wide report, each sample gets a column with values)
                  max_prob_across_exp,#you created this column above, this col is analogous to Maxquant probability
                  modification_title,
                  multiplicity,
                  reverse,
                  potential_contaminant
      ),
      names_from = sample,#the new wide columns will be assigned these names (you created the "sample" column above)
      values_from = c(quantity, site_probability, group, nr_of_collapsed_peptides), #these are sample specific columns that I want to add as wide columns
      names_glue = "{.value}_{sample}" # ensures a pre-fix in front ot the columns, not quite sure how pivoting works this out.
    )
  
  #rename some cols
  data <- wide_df %>%
    dplyr::rename_with(~ sub("^quantity_", "intensity_", .x), starts_with("quantity_")) %>% #rename the "quantity" columns to "intensity" to match Arvid's SEV pacakge
    dplyr::rename_with(~ sub("^group", "ptm_group_", .x), starts_with("group"))
  
  # Standardize column names
  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()
  
  # Ensure key columns are character and missing values are "" (to match read.delim behavior)
  data <- data %>%
    dplyr::mutate(dplyr::across(
      c(genes, site_aa, site_location, reverse, potential_contaminant),
      ~ ifelse(is.na(.), "", as.character(.))
    ))
  
  # Clean up ID columns
  data <- data %>%
    dplyr::mutate(
      #orig_prot_ids = proteins,
      #orig_gene_names = gene_names,
      #orig_positions_within_proteins = positions_within_proteins,
      #orig_sequence_window = sequence_window,
      uniprot_first = gsub("(.*?);.*", "\\1", uni_prot_ids),
      gene_first = gsub("(.*?);.*", "\\1", genes),
      site_location_first = gsub("(.*?);.*", "\\1", site_location),
      flanking_region_first = gsub("(.*?);.*", "\\1", flanking_region)
    )
  
  # Filter out reverse and contaminant hits
  data <- data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(filt), ~.x == ""))
  
  # Extract intensity matrices
  int <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]*$"), colnames(data))]) %>% log2()
  # mult_1 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_1"), colnames(data))]) %>% log2()
  # mult_2 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_2"), colnames(data))]) %>% log2()
  # mult_3 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_3"), colnames(data))]) %>% log2()
  
  # Replace -Inf with NA
  int[is.infinite(int)] <- NA
  # mult_1[is.infinite(mult_1)] <- NA
  # mult_2[is.infinite(mult_2)] <- NA
  # mult_3[is.infinite(mult_3)] <- NA
  
  # Format column names and IDs
  label <- colnames(int)
  ID <- gsub(paste0("intensity_(.*)", sep, "(.*)"), "\\1_\\2", label)
  colnames(int) <- ID
  # colnames(mult_1) <- ID
  # colnames(mult_2) <- ID
  # colnames(mult_3) <- ID
  
  # Experimental design fallback
  if (!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))) {
    experimental_design <- data.frame(
      label = label,
      sample = gsub("intensity_", "", label),
      ID = ID,
      condition = gsub(paste0("intensity_|", sep, "[0-9].*"), "", label),
      replicate = gsub(paste0("^.*", sep, "(?=[0-9])"), "", label, perl = TRUE)
    )
  }
  
  # Construct rowData with proper types
  rowdata <- data %>% dplyr::select(-grep("intensity", colnames(data))) #this will be rowdata in your summarized experiment
  rowdata <- DEP2::make_unique(rowdata, names = use_id, ids = use_id, delim = "/") #make unique ensures unique IDs, probably not necessary, because the collapse_key is unique in spectronaut (but i do it anyway.
  
  #add columns required later on 
  rowdata <- rowdata %>% 
    mutate(gene_names = genes, protein_ids = protein_groups, ID = final_id)
  #dplyr::select(name,genes, gene_names,uni_prot_ids, protein_ids, ID, protein_descriptions, contains('_vs_'),1:ncol(.)))
  
  
  # Assemble PhosphoExperiment object
  ppe <- PhosR::PhosphoExperiment(
    assays = setNames(list(int), assay_name),
    rowData = rowdata,
    colData = experimental_design,
    Site = as.numeric(data$site_location_first),
    GeneSymbol = data$gene_first,
    Residue = data$site_aa,
    Sequence = data$flanking_region_first,
    UniprotID = data$uniprot_first,
    Localisation = as.numeric(data$max_prob_across_exp)
  )
  
  # Row names = unique site identifier
  # rownames(ppe) <- paste(
  #   rowdata$uniprot_first,
  #   rowdata$gene_first,
  #   paste0(rowdata$site_aa, rowdata$site_location_first),
  #   rowdata$flanking_region_first,
  #   sep = ";"
  # )
  rownames(ppe) <-rowdata$name
  
  # Return based on output type
  if (output_type == "long_df") return(DEP2::get_df_long(ppe))
  if (output_type == "wide_df") return(DEP2::get_df_wide(ppe))
  if (output_type == "all") return(list(long_df = DEP2::get_df_long(ppe), wide_df = DEP2::get_df_wide(ppe), ppe = ppe))
  return(ppe)
}


# ---------------------------------------------------------------------------
# Renamed functions (Phase 3) -> snake_case verb_noun names
#
# Each old name forwards to its renamed counterpart. See AUDIT.md section 6.
# ---------------------------------------------------------------------------

#' Deprecated function names (renamed in Phase 3)
#'
#' These functions have been renamed for consistency (snake_case, predictable
#' \code{verb_noun} prefixes). The old names still work but emit a deprecation
#' warning and forward to the new function.
#'
#' \describe{
#'   \item{\code{scatterPlot}}{use \code{\link{plot_scatter}}}
#'   \item{\code{corr_plot}}{use \code{\link{plot_correlation}}}
#'   \item{\code{clustered_heatmap}}{use \code{\link{plot_heatmap_clustered}}}
#'   \item{\code{my_theme}}{use \code{\link{theme_sev}}}
#'   \item{\code{colab_subset_df}}{use \code{\link{subset_results}}}
#'   \item{\code{plot_indiviuals}}{use \code{\link{plot_individuals}}}
#'   \item{\code{fix_maxq_pig}}{use \code{\link{fix_maxquant_pig_annotation}}}
#'   \item{\code{se_GOE}}{use \code{\link{enrich_go_se}}}
#'   \item{\code{to_long}}{use \code{\link{se_to_long}}}
#'   \item{\code{advanced_test}}{use \code{\link{test_diff_limma}}}
#'   \item{\code{long_test}}{use \code{\link{test_diff_long}}}
#'   \item{\code{add_sign}}{use \code{\link{add_significance}}}
#' }
#'
#' @param ... Passed on to the renamed function.
#' @return The value returned by the renamed function.
#' @name sev-deprecated-renamed
#' @keywords internal
NULL

#' @rdname sev-deprecated-renamed
#' @export
scatterPlot <- function(...) { .Deprecated("plot_scatter"); plot_scatter(...) }

#' @rdname sev-deprecated-renamed
#' @export
corr_plot <- function(...) { .Deprecated("plot_correlation"); plot_correlation(...) }

#' @rdname sev-deprecated-renamed
#' @export
clustered_heatmap <- function(...) { .Deprecated("plot_heatmap_clustered"); plot_heatmap_clustered(...) }

#' @rdname sev-deprecated-renamed
#' @export
my_theme <- function(...) { .Deprecated("theme_sev"); theme_sev(...) }

#' @rdname sev-deprecated-renamed
#' @export
colab_subset_df <- function(...) { .Deprecated("subset_results"); subset_results(...) }

#' @rdname sev-deprecated-renamed
#' @export
plot_indiviuals <- function(...) { .Deprecated("plot_individuals"); plot_individuals(...) }

#' @rdname sev-deprecated-renamed
#' @export
fix_maxq_pig <- function(...) { .Deprecated("fix_maxquant_pig_annotation"); fix_maxquant_pig_annotation(...) }

#' @rdname sev-deprecated-renamed
#' @export
se_GOE <- function(...) { .Deprecated("enrich_go_se"); enrich_go_se(...) }

#' @rdname sev-deprecated-renamed
#' @export
to_long <- function(...) { .Deprecated("se_to_long"); se_to_long(...) }

#' @rdname sev-deprecated-renamed
#' @export
advanced_test <- function(...) { .Deprecated("test_diff_limma"); test_diff_limma(...) }

#' @rdname sev-deprecated-renamed
#' @export
long_test <- function(...) { .Deprecated("test_diff_long"); test_diff_long(...) }

#' @rdname sev-deprecated-renamed
#' @export
add_sign <- function(...) { .Deprecated("add_significance"); add_significance(...) }
