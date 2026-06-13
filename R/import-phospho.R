# Import phospho / PTM output into (Phospho)Experiment objects
#
# Part of the 'sev' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' Read in Phosphoproteomics Data and Perform Initial Processing
#'
#' This function reads in phosphoproteomics data from the "Phospho (STY)Sites.txt" MaxQuant output, processes it by splitting protein groups,
#' filtering low-quality hits, and expanding the site table. It also creates site IDs and makes gene names unique. In the end, it returns a 
#' SummarizedExperiment.
#' Experimental design can be provided or generated automatically if not provided.
#'
#' @param file A string representing the file path of the tab-delimited phosphoproteomics data.
#' @param gene_column A string representing the name of the gene column in the input data. Default is "gene_names".
#' @param protein_column A string representing the name of the protein column in the input data. Default is "proteins".
#' @param sep A string used as the separator for intensity columns in the input data. Default is "_rep_".
#' @param filt A character vector containing the filters to apply for false and low-quality hits. Default is c("reverse", "potential_contaminant").
#' @param keep_all_proteins A boolean indicating whether to keep all proteins when splitting protein groups. Default is FALSE.
#' @param keep_all_genes A boolean indicating whether to keep all genes when splitting protein groups. Default is FALSE.
#' @param experimental_design A data.frame containing the experimental design information. Default is an empty string.
#' 
#' @importFrom vroom vroom
#' @importFrom janitor make_clean_names
#' @importFrom dplyr select filter
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom DEP2 make_se make_unique
#' 
#' @return A summarizedExperiment containing the phosphoproteomics data.
#' @examples
#' \dontrun{
#' # Assuming "Phospho (STY)Sites.txt" is a valid MaxQuant phosphoproteomics data file
#' processed_data <- phos_read_in_int("Phospho (STY)Sites.txt")
#' }
#' @export
phos_read_in_int <- function(file, gene_column = "gene_names", protein_column = "proteins", sep="_rep_", filt = c("reverse", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  data <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names)
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = .[[protein_column]],
           orig_gene_names =  .[[gene_column]]) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes) %>%
    dplyr::select(-contains("score"))
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ is.na(.x)))
  
  # intensity columns
  cols <- grep(paste0(".*",sep,"[0-9]*_[0-9]*$"), colnames(data), value = TRUE)
  
  #expand site_table
  data <-  data %>% tidyr::pivot_longer(cols = cols, names_to = c("set", "multiplicity"), values_to = "int", names_pattern = "(.*)_([0-9]*)$") %>% mutate(set = gsub("intensity", "value", set)) %>% tidyr::pivot_wider(values_from = "int", names_from = "set")
  
  #create site id
  data <- data %>% mutate(site_id = paste0(gene_names, "_", amino_acid, position), 
                          site_id_mult = paste0(gene_names, "_", amino_acid, position, "_", multiplicity))
  
  #Make gene_names unique
  data_unique <- make_unique(data, "site_id_mult", protein_column, delim=";")
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    
    labels <- grep(paste0("intensity_.*", sep, "[0-9]*$"), colnames(data), value = TRUE)
    
    experimental_design <- data.frame(label=labels, 
                                      sample=gsub("intensity_", "", labels),
                                      ID = gsub("intensity_", "", labels),
                                      condition=gsub(paste0("intensity_|",sep, "[0-9]*"), "", labels),
                                      replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", labels, perl = TRUE))
    
  }else{
    experimental_design <- data.frame(experimental_design)
  }
  
  data_se <- DEP2::make_se(data_unique, which(colnames(data_unique) %in% gsub("intensity", "value",labels)), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "intensity_raw"
  return(data_se)
}


#' Read in phosphoproteomics occupancy data and process it
#'
#' This function reads in phosphoproteomics occupancy data from the "Phospho (STY)Sites.txt" MaxQuant output,
#' processes the data, filters false and low-quality hits, and returns a SummarizedExperiment object.
#'
#' @param file A character string specifying the input file path.
#' @param gene_column A character string specifying the column name for gene names in the input file (default: "gene_names").
#' @param protein_column A character string specifying the column name for protein names in the input file (default: "proteins").
#' @param sep A character string specifying the separator for the replicate information in the column names (default: "_rep_").
#' @param filt A character vector specifying the columns to filter for NA values (default: c("reverse", "potential_contaminant")).
#' @param keep_all_proteins A logical value indicating whether to keep all proteins when splitting protein groups (default: FALSE).
#' @param keep_all_genes A logical value indicating whether to keep all genes when splitting gene names (default: FALSE).
#' @param experimental_design A data frame specifying the experimental design. If not provided, it will be generated based on the input data.
#'
#' @importFrom vroom vroom
#' @importFrom janitor make_clean_names
#' @importFrom dplyr select filter
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom DEP2 make_se make_unique
#' 
#' @return A SummarizedExperiment object containing the processed phosphoproteomics occupancy data.
#' @examples
#' \dontrun{
#' # Example input file path (replace with your own file path)
#' input_file <- "path/to/your/input_file.tsv"
#'
#' # Read and process the data
#' data_se <- phos_read_in_occ(input_file)
#' }
#'@export
phos_read_in_occ <- function(file, gene_column = "gene_names", protein_column = "proteins", sep="_rep_", filt = c("reverse", "potential_contaminant"), keep_all_proteins = F, keep_all_genes = F, experimental_design = ""){
  
  data <- vroom::vroom(file, delim = "\t", col_names = T,guess_max = 30000,  .name_repair = janitor::make_clean_names)
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = .[[protein_column]],
           orig_gene_names =  .[[gene_column]]) %>%
    split_genes(colname = gene_column, keep_all = keep_all_proteins) %>%
    split_genes(colname = protein_column, keep_all = keep_all_genes)
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ is.na(.x)))
  
  # intensity columns
  labels <- colnames(data)[grepl("occupancy_", colnames(data)) &
                             !grepl("occupancy_(ratio|error)", colnames(data))]
  
  #create site id
  data <- data %>% mutate(site_id = paste0(gene_names, "_", amino_acid, position))
  
  #Make gene_names unique
  data_unique <- make_unique(data, "site_id", protein_column, delim=";")
  
  #2^x because it is later log2 transformed
  data_unique <- data_unique %>% mutate(across(all_of(labels), ~2^.x))
  
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    
    experimental_design <- data.frame(label=labels, 
                                      sample=gsub("occupancy_", "", labels),
                                      ID = gsub("occupancy_", "", labels),
                                      condition=gsub(paste0("occupancy_|",sep, "[0-9]*"), "", labels),
                                      replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", labels, perl = TRUE))
    
  }else{
    experimental_design <- data.frame(experimental_design)
  }
  
  data_se <- make_se(data_unique, which(colnames(data_unique) %in% labels), experimental_design)
  rownames(data_se) <- data_unique$name
  names(assays(data_se)) <- "occupancy"
  
  assay(data_se)[assay(data_se) == "NaN"] <- NA
  return(data_se)
}



#' Read a PTM / phospho report into a PhosphoExperiment
#'
#' Unified entry point for converting a post-translational-modification (PTM)
#' search-engine report into a \code{\link[PhosR]{PhosphoExperiment}} object.
#' Dispatches to a source-specific reader based on \code{source}. This
#' consolidates the former \code{maxq_to_ppe()} (MaxQuant) and
#' \code{ptm_to_ppe()} / \code{sn_ptm_to_ppe()} (Spectronaut) functions, which
#' are retained as deprecated aliases.
#'
#' @param file Path to the report file (tab-delimited) or a pre-loaded
#'   data frame.
#' @param source Character scalar selecting the input format: \code{"maxquant"}
#'   for a MaxQuant \code{Phospho (STY)Sites.txt} table, or \code{"spectronaut"}
#'   for a Spectronaut PTM report. Default: \code{"maxquant"}.
#' @param ... Further arguments passed to the source-specific reader.
#'   \describe{
#'     \item{\code{source = "maxquant"}}{\code{sep}, \code{filt},
#'       \code{experimental_design}.}
#'     \item{\code{source = "spectronaut"}}{\code{sep}, \code{filt},
#'       \code{experimental_design}, \code{output_type}, \code{assay_name},
#'       \code{mod_type}, \code{site_prob}, \code{aggregate_multiplicity},
#'       \code{use_id} (the \code{rowData} column used as the unique site
#'       identifier, default \code{"site"}).}
#'   }
#'
#' @return A \code{\link[PhosR]{PhosphoExperiment}} object. For the Spectronaut
#'   reader with \code{output_type} other than \code{"ppe"}, a data frame or
#'   list as documented for that reader.
#'
#' @seealso \code{\link{aggregate_phospho}}, \code{\link[PhosR]{PhosphoExperiment}}
#'
#' @importFrom vroom vroom
#' @importFrom janitor clean_names make_clean_names
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom purrr when
#' @importFrom PhosR PhosphoExperiment
#' @importFrom DEP2 make_unique get_df_long get_df_wide
#' @importFrom BiocGenerics colnames rownames
#' @importFrom rlang .data
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # MaxQuant
#' ppe <- read_ptm_to_ppe("Phospho (STY)Sites.txt", source = "maxquant")
#'
#' # Spectronaut, aggregating peptide multiplicities by sum
#' ppe <- read_ptm_to_ppe(
#'   "Spectronaut_PTM_Report.tsv",
#'   source                 = "spectronaut",
#'   aggregate_multiplicity = "sum"
#' )
#' }
#'
#' @export
read_ptm_to_ppe <- function(file, source = c("maxquant", "spectronaut"), ...) {
  source <- match.arg(source)
  switch(source,
    maxquant    = read_ptm_maxquant(file, ...),
    spectronaut = read_ptm_spectronaut(file, ...)
  )
}


#' Convert MaxQuant output to PhosphoExperiment object
#'
#' This function reads a MaxQuant output file, filters the data, and converts it into a PhosphoExperiment object.
#'
#' @param file The input file in tab-delimited format.
#' @param sep The separator used in column names for replicates (default: "_rep_").
#' @param filt A vector of filters to remove unwanted rows (default: c("reverse", "potential_contaminant")).
#' @param experimental_design An optional data frame to provide custom experimental design information.
#' @return A PhosphoExperiment object containing the processed data.
#' @importFrom janitor make_clean_names
#' @importFrom DEP2 make_unique
#' @importFrom BiocGenerics colnames rownames
#' @importFrom PhosR PhosphoExperiment
#' @import dplyr
#' @keywords internal
#' @noRd
read_ptm_maxquant <- function(file, sep="_rep_",
                        filt = c("reverse", "potential_contaminant"), experimental_design = NA){
  
  data <- read.delim(file, sep="\t")
  
  colnames(data) <- colnames(data) %>%
    tolower() %>%
    janitor::make_clean_names()
  
  #Split protein groups to single proteins, keep all
  data <- data %>%
    mutate(orig_prot_ids = proteins,
           orig_gene_names = gene_names, 
           orig_positions_within_proteins = positions_within_proteins,
           orig_sequence_window = sequence_window,
           protein_ids = gsub("(.*?);.*", "\\1", proteins),
           gene_names = gsub("(.*?);.*", "\\1", gene_names),
           positions_within_proteins = gsub("(.*?);.*", "\\1", positions_within_proteins),
           sequence_window = gsub("(.*?);.*", "\\1", sequence_window))
  
  
  #Filter false and low quality hits
  data <- data %>% filter(if_all(filt, ~ .x == ""))
  

  # Create assays for ppe
  int <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]*$"), colnames(data))]) %>% log2() 
  mult_1 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_1"), colnames(data))]) %>% log2()
  mult_2 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_2"), colnames(data))]) %>% log2()
  mult_3 <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]_3"), colnames(data))]) %>% log2()
  
  # replace missing values with NA
  int[is.infinite(int)] <- NA 
  mult_1[is.infinite(mult_1)] <- NA
  mult_2[is.infinite(mult_2)] <- NA
  mult_3[is.infinite(mult_3)] <- NA
  
  # Change colnames
  label<- colnames(int)
  
  ID <- gsub(paste0("intensity_(.*)", sep, "(.*)"), "\\1_\\2", label)
  colnames(int) <- ID
  colnames(mult_1) <- ID
  colnames(mult_2) <- ID
  colnames(mult_3) <- ID
  
  # Create colData if not provided
  if(!all(c("label", "sample", "condition", "replicate") %in% colnames(experimental_design))){
    experimental_design<-data.frame(label=label,
                                    sample=gsub("intensity_", "", label),
                                    ID = ID,
                                    condition=gsub(paste0("intensity_|",sep, "[0-9].*"), "", label),
                                    replicate=gsub(paste0("^.*",sep,"(?=[0-9])"), "", label, perl = TRUE))
  }
  
  #Create rowData
  rowdata <- data %>% select(-grep("intensity", colnames(data)), -starts_with("score"), -contains("_diff")) %>% mutate(gene_names = ifelse(gene_names == "", "NA", gene_names)) %>% mutate(psite = paste0(gene_names, ";", amino_acid, ";", positions_within_proteins)) %>% make_unique(names = "psite", ids = "psite", delim = "/")
  
  # Create PPE
  ppe <- PhosphoExperiment(assays = list(Quantification = int, 
                                         multiplicity_1 = mult_1,
                                         multiplicity_2 = mult_2,
                                         multiplicity_3 = mult_3),
                           rowData = rowdata,
                           colData = experimental_design,
                           Site = as.numeric(data$positions_within_proteins), 
                           GeneSymbol = data$gene_names, 
                           Residue = data$amino_acid, 
                           Sequence = data$sequence_window, 
                           UniprotID = data$protein_ids,
                           Localisation = as.numeric(data$localization_prob))
  
  rownames(ppe) <- paste(rowdata$protein_ids, rowdata$gene_names, paste0(rowdata$amino_acid, rowdata$positions_within_proteins), rowdata$sequence_window, sep = ";")
  return(ppe)
}



#' Convert a Spectronaut PTM report to a PhosphoExperiment object
#'
#' Reads a Spectronaut phosphoproteomics report (tab-delimited file or
#' pre-loaded data frame), filters and cleans the data, resolves
#' protein-group columns to single-protein entries, builds unique site
#' identifiers, optionally aggregates peptide multiplicity, pivots to wide
#' format, log2-transforms intensities, and assembles a
#' \code{\link[PhosR]{PhosphoExperiment}} object ready for downstream analysis
#' with \pkg{PhosR} or \pkg{DEP2}.
#'
#' @param file Either a path (character scalar) to a tab-delimited Spectronaut
#'   PTM report, or a data frame / tibble that has already been loaded into R.
#' @param sep Character scalar. Separator used to split sample names into
#'   condition and replicate parts when constructing the experimental design
#'   automatically (default: \code{"_rep_"}).
#' @param filt Character vector of column names whose non-empty values flag
#'   rows to be removed (e.g. reverse hits, contaminants).
#'   Default: \code{c("reverse", "potential_contaminant")}.
#' @param experimental_design Optional data frame with at least columns
#'   \code{file_name} and \code{sample}. When supplied it overrides the
#'   automatically derived sample labels. Default: \code{NULL}.
#' @param output_type Character scalar controlling what is returned:
#'   \describe{
#'     \item{\code{"ppe"} (default)}{A \code{PhosphoExperiment} object.}
#'     \item{\code{"long_df"}}{Long-format data frame via
#'       \code{DEP2::get_df_long()}.}
#'     \item{\code{"wide_df"}}{Wide-format data frame via
#'       \code{DEP2::get_df_wide()}.}
#'     \item{\code{"all"}}{A named list with elements \code{long_df},
#'       \code{wide_df}, and \code{ppe}.}
#'   }
#' @param assay_name Character scalar. Name given to the intensity assay
#'   stored in the \code{PhosphoExperiment} (default: \code{"orig_assay"}).
#' @param mod_type Character vector of modification titles to retain
#'   (matched against the \code{modification_title} column). Set to
#'   \code{NULL} to skip modification filtering.
#'   Default: \code{"Phospho (STY)"}.
#' @param site_prob Numeric scalar in \[0, 1\]. Minimum site localisation
#'   probability required to keep a row (default: \code{0.75}).
#' @param aggregate_multiplicity Character scalar passed to
#'   \code{\link{aggregate_phospho}}. Controls how peptide entries
#'   with the same site but different multiplicities are collapsed:
#'   \code{"sum"}, \code{"max"}, \code{"mean"}, or \code{"none"}.
#'   When \code{"none"} (or \code{NULL}), multiplicity is preserved as part
#'   of the site identifier (\code{_M<n>} suffix).
#'
#' @return Depends on \code{output_type}:
#' \itemize{
#'   \item \code{"ppe"}: a \code{\link[PhosR]{PhosphoExperiment}} object whose
#'     assay contains log2-transformed intensities and whose \code{rowData}
#'     holds site-level metadata.
#'   \item \code{"long_df"} / \code{"wide_df"}: a data frame.
#'   \item \code{"all"}: a list of length 3.
#' }
#'
#' @details
#' The processing pipeline proceeds in the following order:
#' \enumerate{
#'   \item \strong{Read}: The input is read with \code{vroom::vroom()} (if a
#'     file path) or coerced with \code{janitor::clean_names()} (if a data
#'     frame). All column names are lower-cased and common prefixes
#'     (\code{pg_}, \code{r_}, \code{ptm_}) are stripped.
#'   \item \strong{Filter}: Rows failing the site probability threshold, not
#'     matching \code{mod_type}, or flagged by \code{filt} columns are
#'     removed.
#'   \item \strong{Resolve protein groups}: \code{\link{extract_by_protein_id}}
#'     selects the correct entry from each semicolon-delimited group column.
#'   \item \strong{Site identifier}: A \code{site} column is constructed as
#'     \code{<gene|protein_id>_<aa><location>} (with an optional \code{_M<n>}
#'     suffix when multiplicity is not aggregated).
#'   \item \strong{Aggregate}: \code{\link{aggregate_phospho}} collapses
#'     multiple peptide rows per site using the chosen strategy.
#'   \item \strong{Pivot}: Data are pivoted to wide format; one column per
#'     sample for each of \code{quantity}, \code{site_probability}, and
#'     \code{nr_of_collapsed_peptides}.
#'   \item \strong{Log2 transform}: Intensity columns are log2-transformed;
#'     \code{-Inf} values (from zero intensities) are replaced with \code{NA}.
#'   \item \strong{Assemble PPE}: A \code{PhosphoExperiment} is built with
#'     site, residue, gene symbol, UniProt ID, flanking sequence, and
#'     localisation probability stored in the appropriate slots.
#' }
#'
#' @seealso
#'   \code{\link{extract_by_protein_id}},
#'   \code{\link{aggregate_phospho}},
#'   \code{\link[PhosR]{PhosphoExperiment}},
#'   \code{\link[DEP2]{make_unique}}
#'
#' @importFrom vroom vroom
#' @importFrom janitor clean_names make_clean_names
#' @importFrom dplyr rename_all mutate group_by ungroup filter select left_join rename_with across any_of if_all all_of
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom purrr when
#' @importFrom PhosR PhosphoExperiment
#' @importFrom DEP2 make_unique get_df_long get_df_wide
#'
#' @examples
#' \dontrun{
#' # Minimal usage - read from file, keep default filters and phospho-STY only
#' ppe <- ptm_to_ppe("Spectronaut_PTM_Report.tsv")
#'
#' # Aggregate by summing multiplicities, return wide data frame
#' wide <- ptm_to_ppe(
#'   file                   = "Spectronaut_PTM_Report.tsv",
#'   aggregate_multiplicity = "sum",
#'   output_type            = "wide_df"
#' )
#'
#' # Supply a custom experimental design
#' ed <- data.frame(
#'   file_name = c("run1.raw", "run2.raw"),
#'   sample    = c("ctrl_rep_1", "treat_rep_1")
#' )
#' ppe <- ptm_to_ppe("report.tsv", experimental_design = ed)
#' }
#'
#' @keywords internal
#' @noRd
read_ptm_spectronaut <- function(file, sep = "_rep_", filt = c("reverse", "potential_contaminant"),
                       experimental_design = NULL, output_type = "ppe", assay_name = "orig_assay",
                       mod_type = "Phospho (STY)", site_prob = 0.75, aggregate_multiplicity = c("sum", "max", "mean", "none"),
                       use_id = "site")
{
  
  # read in file
  if (is.character(file)) {
    data <- vroom::vroom(file, delim = "\t", col_names = TRUE, 
                         guess_max = 30000, .name_repair = janitor::make_clean_names)
  }
  else {
    data <- file %>% janitor::clean_names()
  }
  
  # rename and clean up
  long_df <- data %>% dplyr::rename_all(tolower) %>% janitor::clean_names() %>% 
    dplyr::rename_all(.funs = list(~gsub("^pg_", "", .))) %>% 
    dplyr::rename_all(.funs = list(~gsub("^r_", "", .))) %>% 
    dplyr::rename_all(.funs = list(~gsub("^ptm_", "", .))) %>% 
    dplyr::mutate(sample = paste0(condition, "_rep_", replicate)) %>% 
    dplyr::group_by(collapse_key) %>% dplyr::mutate(max_prob_across_exp = max(site_probability, 
                                                                              na.rm = TRUE)) %>% dplyr::ungroup() %>% dplyr::filter(site_probability >= 
                                                                                                                                      site_prob) %>% purrr::when(!is.null(mod_type) && !all(is.na(mod_type)) ~ 
                                                                                                                                                                   dplyr::filter(., modification_title %in% mod_type), ~.) %>% 
    dplyr::mutate(collapse_key_simplified = stringr::str_replace(collapse_key, 
                                                                 "_M\\d+$", ""))
  
  
  # filter
  long_df <- long_df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(intersect(filt, names(data))), 
                                                     ~.x == ""))
  
  # update site identifiers from uniprot to gene based and only keep correct values for each line
  # This is necessary because all orginial protein_group columns contain multiple values separated by ";"
  df_final <- long_df %>% extract_by_protein_id()
  
  # create a beautiful site identifier that does not contain multiplicity as we are aggregating this in a moment
  # df_final <- df_final %>%
  #   mutate(site = paste0(
  #     ifelse(
  #       !is.na(genes), genes, protein_id), "_", site_aa, site_location, ifelse(is.null(aggregate_multiplicity), paste0("_M", multiplicity), "")))
  
  #I (Moritz) added this on 22 April 2026
  df_final <- df_final %>%
    mutate(
      site = paste0(
        dplyr::coalesce(genes, protein_id),
        "_", site_aa, site_location,
        if (is.null(aggregate_multiplicity) || aggregate_multiplicity == "none") {
          paste0("_M", multiplicity)
        } else {
          ""
        }
      )
    )
  
  # aggregate multiplicity
  df_final <- df_final %>% aggregate_phospho(method = aggregate_multiplicity)
  
  
  # add sample column from experimental_design
  if(!is.null(experimental_design)){
    df_final <- df_final %>% select(-sample) %>% 
      dplyr::left_join(dplyr::select(experimental_design, c(file_name, sample)), by = "file_name")
  }
  
  # transform to wide
  wide_df <- df_final %>% tidyr::pivot_wider(id_cols = any_of(c(
    "genes", "protein_id", "site", "fasta_headers", "protein_groups", 
    "protein_accessions", "uni_prot_ids", "protein_names", 
    "protein_descriptions", "organisms", "molecular_weight", 
    "nr_of_stripped_sequences_measured", 
    "nr_of_stripped_sequences_identified_experiment_wide", "flanking_region", 
    "site_aa", "site_location", "max_prob_across_exp", 
    "modification_title", "multiplicity", # <--- any_of handles this
    "reverse", "potential_contaminant"
  )),
  names_from = sample, values_from = c(quantity, site_probability, nr_of_collapsed_peptides), names_glue = "{.value}_{sample}")
  
  # clean up 
  data <- wide_df %>% dplyr::rename_with(~sub("^quantity_", 
                                              "intensity_", .x), starts_with("quantity_")) %>% dplyr::rename_with(~sub("^group", 
                                                                                                                       "ptm_group_", .x), starts_with("group"))
  colnames(data) <- colnames(data) %>% tolower() %>% janitor::make_clean_names()
  
  data <- data %>% dplyr::mutate(dplyr::across(any_of(c("genes", "site_aa", 
                                                        "site_location", "reverse", "potential_contaminant")), ~ifelse(is.na(.), 
                                                                                                                       "", as.character(.))))
  
  # log2 transformation
  int <- as.matrix(data[grep(paste0("intensity_.*", sep, "[1-9]*$"), 
                             colnames(data))]) %>% log2()
  int[is.infinite(int)] <- NA
  
  label <- colnames(int)
  ID <- gsub(paste0("intensity_(.*)", sep, "(.*)"), "\\1_\\2", 
             label)
  colnames(int) <- ID
  if (is.null(experimental_design)){
    experimental_design <- data.frame(label = label, sample = gsub("intensity_", 
                                                                   "", label), ID = ID, condition = gsub(paste0("intensity_|", 
                                                                                                                sep, "[0-9].*"), "", label), replicate = gsub(paste0("^.*", 
                                                                                                                                                                     sep, "(?=[0-9])"), "", label, perl = TRUE))
  }
  rowdata <- data %>% dplyr::select(-grep("intensity", colnames(data)))
  if (!use_id %in% names(rowdata))
    stop(paste0("use_id = '", use_id, "' is not a column of the assembled rowData."))
  rowdata <- DEP2::make_unique(rowdata, names = use_id, ids = use_id,
                               delim = "/")
  rowdata <- rowdata %>% mutate(gene_names = genes, protein_ids = protein_id,
                                ID = .data[[use_id]])
  ppe <- PhosR::PhosphoExperiment(assays = setNames(list(int), 
                                                    assay_name), rowData = rowdata, colData = experimental_design, 
                                  Site = as.numeric(data$site_location), GeneSymbol = data$genes, 
                                  Residue = data$site_aa, Sequence = data$flanking_region, 
                                  UniprotID = data$protein_id, Localisation = as.numeric(data$max_prob_across_exp))
  rownames(ppe) <- rowdata$name
  if (output_type == "long_df") 
    return(DEP2::get_df_long(ppe))
  if (output_type == "wide_df") 
    return(DEP2::get_df_wide(ppe))
  if (output_type == "all") 
    return(list(long_df = DEP2::get_df_long(ppe), wide_df = DEP2::get_df_wide(ppe), 
                ppe = ppe))
  return(ppe)
}



#' Aggregate phosphopeptide rows to unique phosphosites
#'
#' Collapses multiple rows that map to the same phosphosite
#' (identified by \code{site}, \code{file_name}, and \code{sample}) into a
#' single row using a user-chosen aggregation strategy. Metadata columns are
#' carried forward by taking the first observed value per group.
#'
#' @param df A long-format data frame produced upstream (e.g., inside
#'   \code{\link{ptm_to_ppe}}). Must contain at minimum the columns
#'   \code{site}, \code{file_name}, \code{sample}, \code{quantity},
#'   \code{site_probability}, and \code{nr_of_collapsed_peptides}.
#' @param method Character scalar specifying how \code{quantity} values from
#'   multiple peptides mapping to the same site are combined. One of:
#'   \describe{
#'     \item{\code{"none"} (default)}{No aggregation; \code{df} is returned
#'       unchanged.}
#'     \item{\code{"sum"}}{Sum of non-\code{NA} quantities; returns
#'       \code{NA_real_} when all values are \code{NA}.}
#'     \item{\code{"max"}}{Maximum of non-\code{NA} quantities.}
#'     \item{\code{"mean"}}{Mean of non-\code{NA} quantities.}
#'   }
#'   Passing \code{NULL} is equivalent to \code{"none"}.
#'
#' @return A data frame (tibble) with one row per unique
#'   \code{(site, file_name, sample)} combination. Additional columns retained
#'   via \code{dplyr::first()}:
#'   \code{genes}, \code{protein_id}, \code{protein_accessions},
#'   \code{site_aa}, \code{site_location}, \code{flanking_region},
#'   \code{protein_names}, \code{condition}, \code{replicate}.
#'   A column \code{aggregation_method} recording the chosen \code{method} is
#'   appended. Row order matches the order in which sites first appeared in
#'   \code{df}.
#'
#' @details
#' Regardless of the \code{quantity} aggregation method chosen:
#' \itemize{
#'   \item \code{site_probability} is always aggregated by maximum.
#'   \item \code{nr_of_collapsed_peptides} is always aggregated by sum.
#'   \item All-\code{NA} groups produce \code{NA_real_} rather than raising
#'     warnings (safe wrappers are used internally).
#' }
#' The original row order is preserved by recording row numbers before
#' grouping and re-joining them after summarisation.
#'
#' @importFrom dplyr group_by summarize mutate left_join arrange select distinct across any_of first
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Aggregate by summing quantities of co-eluting multiply phosphorylated
#' # peptides mapping to the same site
#' df_agg <- aggregate_phospho(long_df, method = "sum")
#' }
#'
#' @export
aggregate_phospho <- function(df, method = "none") {
  
  # Validate method upfront
  valid_methods <- c("none", "sum", "max", "mean")
  if (!is.null(method) && !method %in% valid_methods) {
    stop(sprintf("method must be one of: %s (or NULL)", 
                 paste(valid_methods, collapse = ", ")))
  }
  
  # 0. NULL / no-op option -> return df unchanged
  if (is.null(method) || method == "none") {
    return(df)
  }
  
  # 1. Capture original site order before any grouping/summarizing
  site_order <- df %>%
    distinct(site, file_name, sample) %>%
    mutate(.row_order = row_number())
  
  # 2. Define safe aggregation helpers (all return NA_real_ for all-NA input)
  safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x,  na.rm = TRUE)
  safe_sum  <- function(x) if (all(is.na(x))) NA_real_ else sum(x,  na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  
  # 3. Assign aggregation function
  agg_func <- switch(method,
                     sum  = safe_sum,
                     max  = safe_max,
                     mean = safe_mean
  )
  
  # 4. Perform the aggregation
  df_agg <- df %>%
    group_by(site, file_name, sample) %>%
    summarize(
      quantity                 = agg_func(quantity),
      site_probability         = safe_max(site_probability),
      nr_of_collapsed_peptides = safe_sum(nr_of_collapsed_peptides),
      nr_of_stripped_sequences_identified_experiment_wide = safe_max(nr_of_stripped_sequences_identified_experiment_wide),#Moritz added
      across(
        any_of(c(
          "genes", "protein_id", "protein_accessions",
          "site_aa", "site_location", "flanking_region",
          "protein_names", "condition", "replicate"
        )),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) %>%
    mutate(aggregation_method = method)
  
  # 5. Restore original row order and drop the helper column
  df_agg <- df_agg %>%
    left_join(site_order, by = c("site", "file_name", "sample")) %>%
    arrange(.row_order) %>%
    select(-.row_order)
  
  return(df_agg)
}


#' Filter Proteins Based on Minimum Peptide Count
#'
#' Sets protein quantifications to \code{NA} when fewer than the minimum number
#' of stripped peptide sequences were used for quantification in a given sample.
#' Optionally removes proteins that become all-\code{NA} across all filtered
#' assays after filtering.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param min_pep Integer. Minimum number of peptides required for a
#'   quantification to be retained. Values below this threshold are set to
#'   \code{NA}. Default \code{2}.
#' @param peptide_count_pattern Character. Regular expression used to identify
#'   columns in the wide-format data that contain peptide counts. Default matches
#'   Spectronaut output
#'   (\code{"Nr_Of_Stripped_Sequences_Used_For_Quantification$"}).
#' @param remove_all_na Logical. If \code{TRUE}, removes proteins that are
#'   all-\code{NA} across \emph{all} filtered assays after the masking step.
#'   Proteins with data surviving in at least one filtered assay are retained.
#'   Default \code{TRUE}.
#' @param assays_to_filter Character vector of assay names to apply filtering
#'   to. If \code{NULL}, all assays in \code{se} are filtered. Default
#'   \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints filtering statistics at each
#'   step including the number of new \code{NA}s introduced per assay and
#'   proteins removed. Default \code{TRUE}.
#'
#' @return A \code{SummarizedExperiment} object with assay values below the
#'   peptide threshold replaced by \code{NA} in place. No new assays are
#'   created. If \code{remove_all_na = TRUE}, rows that are entirely
#'   \code{NA} across all filtered assays are dropped.
#'
#' @details
#' Sample-to-peptide-column matching tolerates a \code{_rep} infix in column
#' names (e.g. \code{ko_ang2_rep1_raw_...} is normalised to
#' \code{ko_ang2_1_raw_...} before matching against sample names). If multiple
#' peptide count columns match a sample, the first match is used with a
#' warning. \code{NA} peptide counts are treated conservatively as failing the
#' threshold.
#'
#' Requires the helper function \code{get_df_wide} to be available in the
#' session; the function checks for this and stops with an informative message
#' if it is not found.
#'
#' @seealso \code{\link[SummarizedExperiment]{assay}},
#'   \code{\link[SummarizedExperiment]{assayNames}}
#'
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importFrom dplyr select all_of
#'
#' @examples
#' \dontrun{
#' # Filter using default settings (min 2 peptides, all assays)
#' se_filtered <- censor_by_peptide_count(se)
#'
#' # Require at least 3 peptides, only applied to the log2 assay
#' se_filtered <- censor_by_peptide_count(
#'   se,
#'   min_pep          = 3,
#'   assays_to_filter = "log2_intensity",
#'   remove_all_na    = FALSE
#' )
#' }
#'
#' @export
censor_by_peptide_count <- function(
    se,
    min_pep = 2,
    peptide_count_pattern = "Nr_Of_Stripped_Sequences_Used_For_Quantification$",
    remove_all_na = TRUE,
    assays_to_filter = NULL,
    verbose = TRUE
) {
  
  # --- Input validation ---
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  if (!is.numeric(min_pep) || length(min_pep) != 1 || min_pep < 1) {
    stop("min_pep must be a single positive number")
  }
  if (!exists("get_df_wide")) {
    stop("Function 'get_df_wide' not found. Please ensure it is loaded.")
  }
  
  # --- Resolve assays to filter ---
  if (is.null(assays_to_filter)) {
    assays_to_filter <- assayNames(se)
  } else {
    missing_assays <- setdiff(assays_to_filter, assayNames(se))
    if (length(missing_assays) > 0) {
      stop(paste("Assays not found:", paste(missing_assays, collapse = ", ")))
    }
  }
  
  # --- Build peptide count matrix ---
  df_wide <- get_df_wide(se)
  
  peptide_cols <- grep(
    peptide_count_pattern,
    colnames(df_wide),
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(peptide_cols) == 0) {
    stop(paste("No columns found matching pattern:", peptide_count_pattern))
  }
  
  if (verbose) {
    message(paste("Found", length(peptide_cols), "peptide count columns"))
    message(paste("Minimum peptides required:", min_pep))
  }
  
  # --- Align peptide count columns to assay sample order ---
  # Peptide count columns are expected to encode sample names. Extract and match.
  se_samples <- colnames(se)
  
  # --- Align peptide count columns to assay sample order ---
  # Normalize both sides by removing the "_rep" infix, then match.
  # e.g. "ko_ang2_rep1_raw_..." -> "ko_ang2_1_raw_..."
  # and  "ko_ang2_1"            -> "ko_ang2_1"
  
  strip_rep <- function(x) sub("_rep([0-9]+)", "_\\1", x)
  
  peptide_cols_normalized <- strip_rep(peptide_cols)
  
  sample_to_pep_col <- vapply(se_samples, function(s) {
    # Match: normalized peptide col must start with the sample name
    # (followed by "_" to avoid partial matches, e.g. "ctrl_1" vs "ctrl_10")
    match_idx <- which(startsWith(peptide_cols_normalized, paste0(s, "_")))
    
    if (length(match_idx) == 0) {
      stop(paste(
        "No peptide count column found for sample:", s,
        "\nNormalized peptide cols:\n",
        paste(peptide_cols_normalized, collapse = "\n ")
      ))
    }
    if (length(match_idx) > 1) {
      warning(paste(
        "Multiple peptide count columns match sample:", s,
        "-- using first:", peptide_cols[match_idx[1]]
      ))
    }
    peptide_cols[match_idx[1]]  # return the ORIGINAL (un-normalized) column name
  }, character(1))
  
  # Extract peptide counts as a matrix aligned to se columns (proteins x samples)
  pep_count_matrix <- df_wide %>%
    dplyr::select(dplyr::all_of(unname(sample_to_pep_col))) %>%
    as.matrix()
  
  colnames(pep_count_matrix) <- se_samples  # rename to match se column order
  
  if (nrow(pep_count_matrix) != nrow(se)) {
    stop(paste(
      "Row mismatch: peptide count matrix has", nrow(pep_count_matrix),
      "rows but se has", nrow(se)
    ))
  }
  
  # --- Build NA mask: TRUE where peptide count is below threshold OR is NA ---
  # NA peptide counts are treated conservatively as failing the threshold.
  na_mask <- is.na(pep_count_matrix) | (pep_count_matrix < min_pep)
  
  # --- Apply mask to each assay in place ---
  for (assay_name in assays_to_filter) {
    assay_data <- assay(se, assay_name)
    
    if (verbose) {
      new_nas <- sum(na_mask & !is.na(assay_data))
      total_quantified <- sum(!is.na(assay_data))
      message(sprintf(
        "Assay '%s': setting %d new NAs (%.1f%% of %d quantified values)",
        assay_name, new_nas, 100 * new_nas / max(total_quantified, 1), total_quantified
      ))
    }
    
    assay_data[na_mask] <- NA
    assays(se)[[assay_name]] <- assay_data
  }
  
  # --- Optionally remove proteins that are all-NA across ALL filtered assays ---
  # A row is only removed if it is all-NA in every one of the filtered assays.
  # Rows with data surviving in any filtered assay are retained.
  if (remove_all_na) {
    all_na_per_assay <- vapply(assays_to_filter, function(assay_name) {
      rowSums(!is.na(assay(se, assay_name))) == 0
    }, logical(nrow(se)))
    
    # all_na_per_assay is a (proteins x assays) logical matrix
    rows_all_na_in_all_filtered <- rowSums(!all_na_per_assay) == 0
    rows_to_keep <- !rows_all_na_in_all_filtered
    n_removed <- sum(rows_all_na_in_all_filtered)
    
    se <- se[rows_to_keep, ]
    
    if (verbose) {
      message(paste(
        "Removed", n_removed,
        "proteins that were all-NA across all filtered assays after filtering"
      ))
      message(paste("Retained", sum(rows_to_keep), "proteins"))
    }
  }
  
  return(se)
}

