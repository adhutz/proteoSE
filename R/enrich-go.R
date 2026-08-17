# GO-term enrichment and pathway lookups
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.



#' GO term enrichment
#' 
#' Performs GO-term enrichment via clusterprofiler for all contrasts of a summarized experiment.
#' Results are added to the metadata.
#'
#' @param se with t-test results (e.g. via DEP2 test_diff())
#' @param col_names name of columns on that enrichment is performed. If left empty, all columns ending with "_diff" are selected.
#' @param simplify Should similar terms be merged (sim > 0.7)
#' @param ont ontology
#' @param keyType Type of identifier
#' @param OrgDb Organism
#' @param pvalueCutoff significance threshold
#'
#' @return se with GO-term enrichment in metadata
#' @export
#'
#' @importFrom clusterProfiler simplify gseGO
enrich_go_se <- function(se, col_names =c(), simplify = TRUE, ont="all", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pvalueCutoff = 1){

  tab_ <- rowData(se)%>%as.data.frame()

  if(length(col_names) == 0){
    col_names <- colnames(tab_)[grep("_diff",colnames(tab_))]
  }

  m<-list()
  GO_raw <- list()

  for(comp_ in col_names){

    l <- tab_[,comp_]
    names(l)<-tab_$gene_names
    l<-sort(l, decreasing=TRUE)

    #Calculate GO-term enrichment
    temp <- clusterProfiler::gseGO(l, ont="all", keyType = "SYMBOL", OrgDb = OrgDb, pvalueCutoff = 1)

    #Reduce similar terms

    if(simplify == TRUE){
      tryCatch(
        expr = {

          temp <- temp %>% clusterProfiler::simplify(., 0.7, by="NES", select_fun=max)

        },

        error = function(e){
          print("Not enough terms for reduction ")
        }

      )
    }
    #Split genes to list that can be compared to rownames of se
    gene_sets <- temp$core_enrichment %>% strsplit(.,split="/")

    #Assign GO-term names
    names(gene_sets)<-temp$Description

    #Get unique gene names to allow filtering of se by GO-terms
    for(names_ in names(gene_sets)){
      gene_sets[[names_]]<- tab_ %>% filter(gene_names %in% gene_sets[[names_]])%>%rownames()
    }

    gene_sets <- IRanges::CharacterList(gene_sets)

    mcols(gene_sets) <- temp

    metadata(se)$GO_enrichment[[comp_]] <- gene_sets

  }

  return(se)
}


#' Perform Over-representation Analysis (ORA) on Phosphoproteomics Data
#'
#' This function takes a SummarizedExperiment object and conducts Over-representation Analysis (ORA) on proteins with at least one significant phosphosite, using the clusterProfiler package.
#'
#' @param se A SummarizedExperiment object containing phosphoproteomics data.
#' @param contr A character vector specifying the contrasts to be analyzed. Default is "all".
#' @param OrgDb A character string specifying the organism database to be used. Default is "org.Hs.eg.db".
#' @param pvalueCutoff Numeric value for the p-value cutoff. Default is 0.4.
#' @param qvalueCutoff Numeric value for the q-value cutoff. Default is 0.8.
#' @param ont A character vector specifying the Gene Ontology (GO) categories to be analyzed. Default is c("BP", "MF", "CC").
#' @return A list containing the ora object, a combined dataframe and a ggplot DEP2icting terms with q < 0.2
#' @import dplyr
#' @import ggplot2
#' @importFrom clusterProfiler enrichGO
#' @importFrom DOSE parse_ratio
#' @export
phospho_ora <- function(se, contr = "all", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.4, qvalueCutoff = 0.8, ont = c("BP", "MF", "CC")){
  
  se_diff_test <- se %>% test_diff_long()
  if(!contr == "all"){
    se_diff_test <- se_diff_test %>% filter(contrast %in% contr)
  }
  
  # Get significant genes
  sign_genes <- se_diff_test %>% filter(significant == TRUE) %>% select(gene_names, contrast, diff)  %>% distinct(gene_names, .keep_all = TRUE) %>% group_by(contrast) %>% summarize(n = n(), gene_names = list(gene_names), diff = list(diff))
  
  # Get background
  background <- se_diff_test$gene_names %>% unique()
  
  #Create list for ora
  names(sign_genes$gene_names) <- sign_genes$contrast
  
  # Perform ORA via clusterProfiler for each contrast and ontology
  ora_res <- list()
  
  for(ont_ in ont){
    ora_res[[ont_]] <- lapply(sign_genes$gene_names, 
                              clusterProfiler::enrichGO, 
                              OrgDb = OrgDb, 
                              keyType = "SYMBOL", 
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff,
                              ont = ont_, 
                              universe = background)
    
    for(sets_ in names(ora_res[[ont_]])){
      ora_res[[ont_]][[sets_]] <- ora_res[[ont_]][[sets_]]@result %>% mutate(contrast = sets_, ont = ont_)
    }
  }
  
  # Create one large table storing the results
  ora_res_df <- do.call("rbind", 
                        lapply(ora_res, 
                               function(x) do.call("rbind", x))
  ) %>%
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
  
  # Test if faceting variables have at least one value
  if(nrow(ora_res_df %>% filter(qvalue < 0.2)) > 0){
    #Keep only Terms that are significant at least once
    ora_plot <- ora_res_df %>% filter(qvalue < 0.2) %>%
      ggplot(aes(x = FoldEnrichment, 
                 y = reorder(Description, FoldEnrichment),
                 color = ifelse(pvalue < 0.05 & qvalue < 0.5, "Significant", "Non-significant"),
                 fill = contrast,
                 size = p.adjust))+
      geom_point(shape = 21,
                 alpha = 0.5) +
      scale_color_manual(values = c("grey", "black"))+
      scale_size_continuous(range = c(5,2)) +
      labs(title = "q < 0.5", color = "p < 0.05 & q < 0.5", fill = "contrast")+
      facet_grid(ont ~., scales = "free_y", space = "free_y") +
      # scale_color_manual(values = c("grey", "firebrick"))+
      geom_segment(aes(xend=0, yend = Description), color = "black", size = 0.1) +
      theme_bw()
    } else{
      ora_plot <- ggplot() + theme_void()
      }
  
  return(list("res" = ora_res, "df" = ora_res_df, "plot" = ora_plot))
}

#' Retrieve Gene Data from a Specified KEGG Pathway
#'
#' This function uses the R-package KEGGREST to retrieve data for a specified KEGG pathway and then returns a tidy data frame
#' that includes: "kegg_gene_id", "description", "gene_symbol", "kegg_ontology_number", "enzyme_commission_number", "kegg_path_number" and "kegg_path_name" 
#'
#' @param keggpath_id A character string specifying the KEGG pathway ID.
#' Default is "hsa00190" for oxidative phosphorylation in Homo sapiens.
#'
#' @return A data frame containing Gene IDs, Descriptions, Gene symbols, KEGG ontology numbers,
#' Enzyme Commission numbers, the pathway ID queried, and the pathway name.
#'
#' @examples
#' \dontrun{
#' genes_from_kegg("hsa00190")
#' oxphos <- genes_from_kegg(keggpath_id = "hsa00190")  # default, or specify another pathway ID https://www.kegg.jp/entry/map00190
#' cellcycle <- genes_from_kegg(keggpath_id = "hsa04110") #https://www.genome.jp/pathway/map04110
#' mtor <- genes_from_kegg(keggpath_id = "hsa04150") #https://www.kegg.jp/pathway/map=map04150
#' ampk <- genes_from_kegg(keggpath_id = "hsa04152") #https://www.kegg.jp/pathway/map=map04152
#' print(cellcycle)
#' 
#' 
#' # Or mutliple keggpaths, assuming you have a vector of 4 keggpath_ids
#' keggpath_ids <- c("hsa00190", "hsa04110", "hsa04150", "hsa04152")
#' 
#' # Initialize an empty list to store the data frames
#' list_of_dataframes <- list()
#' 
#' # Loop over each keggpath_id
#' for (keggpath_id in keggpath_ids) {
#'   # Call the genes_from_kegg function for each keggpath_id
#'   df <- genes_from_kegg(keggpath_id)
#' 
#'   # Name the data frame according to the keggpath_id and add it to the list
#'   list_of_dataframes[[keggpath_id]] <- df
#' }
#' 
#' # Now, list_of_dataframes is a list of data frames, each named after a keggpath_id
#' # Use bind_rows to combine all data frames into a single data frame
#' combined_df <- bind_rows(list_of_dataframes, .id = "keggpath_id")
#' }
#' @importFrom dplyr mutate
#' @export
genes_from_kegg <- function(keggpath_id = "hsa00190") {
  # KEGGREST, magrittr and dplyr are package Imports (see DESCRIPTION); no
  # runtime install/library() is needed.

  # Fetching KEGG pathway data only once
  pathway_data <- KEGGREST::keggGet(keggpath_id)[[1]]
  
  # Get human genes
  genes <- pathway_data$GENE
  
  # Separate IDs and Descriptions
  kegg_gene_id <- genes[seq(1, length(genes), by = 2)]
  description <- genes[seq(2, length(genes), by = 2)]
  
  # Combine into a data frame
  gene_data <- data.frame(kegg_gene_id, description, stringsAsFactors = FALSE)
  
  # Extract additional details
  gene_data <- gene_data %>%
    dplyr::mutate(
      gene_symbol = sapply(description, function(x) {
        strsplit(x, ";")[[1]][1]
      }),
      kegg_ontology_number = sapply(description, function(x) {
        matches <- regmatches(x, regexpr("(?<=\\[KO:)[^\\]]+", x, perl = TRUE))
        if (length(matches) > 0) matches[1] else NA
      }),
      enzyme_commission_number = sapply(description, function(x) {
        matches <- regmatches(x, regexpr("(?<=\\[EC:)[^\\]]+", x, perl = TRUE))
        if (length(matches) > 0) matches[1] else NA
      }),
      kegg_path_number = keggpath_id, # Add the new column 'kegg'
      kegg_path_name = pathway_data$NAME # Directly add pathway name from fetched data
    )
  
  return(gene_data)
}




###########
#' Build a GO Term Parent-Offspring Relationship Table
#'
#' Retrieves the full offspring and ancestor relationships for Gene Ontology
#' (GO) terms across one or more ontologies using internal
#' \pkg{GOSemSim} functions, and returns them as a single tidy
#' \code{\link[tibble]{tibble}} with one row per GO term per ontology.
#'
#' @param onts Character vector. GO ontologies to process. Valid values are
#'   \code{"BP"} (Biological Process), \code{"CC"} (Cellular Component),
#'   and \code{"MF"} (Molecular Function). Default: \code{c("BP", "CC", "MF")}.
#' @param outfile Character scalar or \code{NULL}. If non-\code{NULL}, the
#'   result is saved as an \code{.Rds} file at this path. The string
#'   \code{"\{date\}"} in the path is replaced with today's date in
#'   \code{YYYYMMDD} format, e.g.
#'   \code{"results/\{date\}_go_terms_full_list.Rds"}.
#'   Default: \code{NULL} (no file written).
#'
#' @details
#' ## Data source
#' Offspring and ancestor vectors are retrieved via the internal
#' \pkg{GOSemSim} functions \code{getOffsprings()} and
#' \code{getAncestors()}.  As these are unexported, they may change without
#' notice; the \pkg{GOSemSim} version this function was written against
#' should be recorded in your project's \code{renv.lock} or
#' \code{DESCRIPTION} file.
#'
#' ## Output columns
#' \describe{
#'   \item{\code{ID}}{GO term accession (e.g. \code{"GO:0008150"}).}
#'   \item{\code{ontology}}{Ontology abbreviation (\code{"BP"}, \code{"CC"},
#'     or \code{"MF"}).}
#'   \item{\code{n_children}}{Number of direct + transitive offspring terms.}
#'   \item{\code{children}}{Pipe-separated offspring IDs, or \code{NA} for
#'     leaf nodes.}
#'   \item{\code{n_parents}}{Number of ancestor terms (root node
#'     \code{"all"} excluded).}
#'   \item{\code{parents}}{Pipe-separated ancestor IDs, or \code{NA} for
#'     the root node.}
#' }
#'
#' ## Saving to disk
#' Use the \code{"\{date\}"} token in \code{outfile} to embed today's date
#' automatically:
#' \preformatted{
#' build_go_terms(outfile = "results/{date}_go_terms_full_list.Rds")
#' }
#'
#' @return A \code{\link[tibble]{tibble}} with columns \code{ID},
#'   \code{ontology}, \code{n_children}, \code{children}, \code{n_parents},
#'   and \code{parents}.
#'
#' @examples
#' \dontrun{
#' # In-memory only
#' terms <- build_go_terms()
#'
#' # Save with a datestamped filename
#' terms <- build_go_terms(
#'   outfile = "results/{date}_go_terms_full_list.Rds"
#' )
#' }
#'
#' @importFrom purrr map_dfr
#' @importFrom dplyr full_join
#' @importFrom tibble tibble
#'
#' @export
build_go_terms <- function(
    onts    = c("BP", "CC", "MF"),
    outfile = NULL
) {
  
  # -------------------------------------------------------------------------
  # Helper 1: collapse a GO ID vector to a pipe-separated string or NA.
  # NAs and zero-length vectors both map to NA_character_.
  # -------------------------------------------------------------------------
  .collapse_ids <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_character_ else paste(x, collapse = "|")
  }
  
  # -------------------------------------------------------------------------
  # Helper 2: turn a named GO list (from getOffsprings / getAncestors) into
  # a tibble with count and pipe-collapsed ID columns.
  # `count_col` and `string_col` are injected as column names via tidy eval.
  # -------------------------------------------------------------------------
  .process_go_list <- function(go_list, count_col, string_col, ont) {
    tibble::tibble(
      ID           = names(go_list),
      ontology     = ont,
      !!count_col  := lengths(lapply(go_list, function(x) x[!is.na(x)])),
      !!string_col := vapply(go_list, .collapse_ids, character(1))
    )
  }
  
  # -------------------------------------------------------------------------
  # 1. Offspring (children): one tibble per ontology, bound into one
  # -------------------------------------------------------------------------
  all_offsprings <- purrr::map_dfr(onts, function(ont) {
    .process_go_list(
      go_list    = GOSemSim:::getOffsprings(ont = ont),
      count_col  = "n_children",
      string_col = "children",
      ont        = ont
    )
  })
  
  # -------------------------------------------------------------------------
  # 2. Ancestors (parents): strip the synthetic root node "all" first
  # -------------------------------------------------------------------------
  all_ancestors <- purrr::map_dfr(onts, function(ont) {
    go_list <- GOSemSim:::getAncestors(ont = ont)
    go_list <- lapply(go_list, function(x) x[x != "all"])
    .process_go_list(
      go_list    = go_list,
      count_col  = "n_parents",
      string_col = "parents",
      ont        = ont
    )
  })
  
  # -------------------------------------------------------------------------
  # 3. Join offspring and ancestor tables on GO ID + ontology
  # -------------------------------------------------------------------------
  terms <- dplyr::full_join(
    all_offsprings,
    all_ancestors,
    by = c("ID", "ontology")
  )
  
  # -------------------------------------------------------------------------
  # 4. Optionally save to disk with {date} token substitution
  # -------------------------------------------------------------------------
  if (!is.null(outfile)) {
    outfile <- gsub(
      "\\{date\\}",
      format(Sys.time(), "%Y%m%d"),
      outfile
    )
    saveRDS(terms, file = outfile)
    message("Saved GO term table to: ", outfile)
  }
  
  terms
}


######################################
######################################
######################################

#' Search GO Terms by Keyword
#'
#' Retrieves all Gene Ontology (GO) terms associated with a given organism,
#' then filters them to those whose \code{term} and/or \code{definition}
#' fields match one or more keywords.  Returns both the full gene-level
#' mapping and a deduplicated table of matching GO terms.
#'
#' @param keywords Character vector.  One or more search terms to match
#'   against GO metadata columns.  Terms are combined into a single
#'   \code{|}-separated regex, so partial matches are supported
#'   (e.g. \code{"lyso"} matches \code{"lysosome"} and
#'   \code{"lysophospholipid"}).
#' @param organism_db An \pkg{AnnotationDbi}-based organism annotation
#'   package object, e.g. \code{org.Hs.eg.db} (human, default) or
#'   \code{org.Mm.eg.db} (mouse).  Must support keytypes
#'   \code{"ENTREZID"} and columns \code{"GO"}, \code{"EVIDENCE"},
#'   \code{"ONTOLOGY"}, and \code{"SYMBOL"}.
#' @param ontology_db An \pkg{AnnotationDbi}-based GO database object.
#'   Default: \code{GO.db}.  Must support keytype \code{"GOID"} and
#'   columns \code{"GOID"}, \code{"TERM"}, \code{"ONTOLOGY"}, and
#'   \code{"DEFINITION"}.
#' @param cols_to_search Character vector.  Cleaned column names (i.e. after
#'   \code{\link[janitor]{clean_names}}) to search for keyword matches.
#'   Default: \code{c("term", "definition")}.
#'
#' @details
#' ## Workflow
#' \enumerate{
#'   \item All GO IDs are fetched from \code{ontology_db} with their term,
#'     ontology, and definition metadata.
#'   \item All Entrez Gene IDs are fetched from \code{organism_db} and
#'     mapped to their GO annotations (evidence code, ontology, gene symbol).
#'   \item The two tables are joined on \code{GO}/\code{GOID} and
#'     \code{ONTOLOGY}, then column names are snake_cased via
#'     \code{\link[janitor]{clean_names}}.
#'   \item Rows are retained where any of \code{cols_to_search} matches the
#'     keyword regex (case-insensitive).
#' }
#'
#' ## Evidence codes
#' The gene-to-GO mapping includes all evidence codes as returned by
#' \code{organism_db}.  If you want to restrict to experimentally validated
#' annotations, filter the \code{filtered} element of the result by the
#' \code{evidence} column (e.g. keep \code{"EXP"}, \code{"IDA"},
#' \code{"IPI"}, \code{"IMP"}, \code{"IGI"}, \code{"IEP"}).
#'
#' ## Performance note
#' This function fetches \emph{all} Entrez IDs and GO IDs from the
#' annotation packages on every call, which can be slow for large organisms.
#' Consider caching the result with \code{\link{saveRDS}} if you call it
#' repeatedly within a session.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{filtered}}{A \code{\link[tibble]{tibble}} of all
#'     gene-level GO annotations matching the keyword regex, with columns
#'     \code{entrezid}, \code{go}, \code{evidence}, \code{ontology},
#'     \code{symbol}, \code{term}, and \code{definition}.}
#'   \item{\code{unique_terms}}{A deduplicated \code{\link[tibble]{tibble}}
#'     of the matching GO terms with columns \code{go}, \code{ontology},
#'     \code{term}, and \code{definition}. Suitable for use as an
#'     annotation reference table.}
#' }
#'
#' @examples
#' \dontrun{
#' library(org.Hs.eg.db)
#' library(GO.db)
#'
#' result <- get_go_terms(
#'   keywords    = c("lysosome", "autophagy"),
#'   organism_db = org.Hs.eg.db,
#'   ontology_db = GO.db
#' )
#'
#' # Gene-level mapping
#' result$filtered
#'
#' # Unique GO terms only
#' result$unique_terms
#'
#' # Restrict to experimental evidence
#' experimental <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
#' result$filtered |>
#'   dplyr::filter(evidence %in% experimental)
#' }
#'
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr filter left_join select
#' @importFrom janitor clean_names
#'
#' @export
get_go_terms <- function(
    keywords,
    organism_db    = org.Hs.eg.db,
    ontology_db    = GO.db,
    cols_to_search = c("term", "definition")
) {
  
  # -------------------------------------------------------------------------
  # 1.  Fetch GO term metadata from the ontology database
  #     (GOID, human-readable term, ontology branch, free-text definition)
  # -------------------------------------------------------------------------
  go_ids <- AnnotationDbi::keys(ontology_db, keytype = "GOID")
  
  go_metadata <- AnnotationDbi::select(
    ontology_db,
    keys     = go_ids,
    columns  = c("GOID", "TERM", "ONTOLOGY", "DEFINITION"),
    keytype  = "GOID"
  )
  
  # -------------------------------------------------------------------------
  # 2.  Fetch gene-to-GO mappings from the organism annotation package
  #     Drop unannotated genes (GO == NA) immediately to keep the table lean
  # -------------------------------------------------------------------------
  entrez_ids <- AnnotationDbi::keys(organism_db, keytype = "ENTREZID")
  
  gene_go <- AnnotationDbi::select(
    organism_db,
    keys    = entrez_ids,
    columns = c("GO", "EVIDENCE", "ONTOLOGY", "SYMBOL"),
    keytype = "ENTREZID"
  ) |>
    dplyr::filter(!is.na(GO))
  
  # -------------------------------------------------------------------------
  # 3.  Join gene-GO mapping with GO metadata, then snake_case all columns
  # -------------------------------------------------------------------------
  complete_go_dataset <- dplyr::left_join(
    gene_go,
    go_metadata,
    by = c("GO" = "GOID", "ONTOLOGY" = "ONTOLOGY")
  ) |>
    janitor::clean_names()
  
  # -------------------------------------------------------------------------
  # 4.  Filter rows where any of `cols_to_search` matches the keyword regex
  #     apply() works row-wise across the selected columns; any() means a
  #     match in *either* term or definition (or both) retains the row.
  # -------------------------------------------------------------------------
  pattern <- paste(keywords, collapse = "|")
  
  row_matches <- apply(
    complete_go_dataset[cols_to_search],
    MARGIN = 1,
    FUN    = function(row) any(grepl(pattern, row, ignore.case = TRUE))
  )
  
  df_filtered <- complete_go_dataset[row_matches, ]
  
  # -------------------------------------------------------------------------
  # 5.  Deduplicate to a clean GO-term reference table (no gene-level info)
  # -------------------------------------------------------------------------
  df_unique_terms <- dplyr::select(df_filtered, go, ontology, term, definition) |>
    unique()
  
  list(
    filtered     = df_filtered,
    unique_terms = df_unique_terms
  )
}



###########################################
###########################################
###########################################
#' Shorten biological keywords or GO terms
#'
#' Abbreviates long biological phrases and common words for cleaner display in
#' plots, tables, and labels.
#'
#' Replacements are applied in two stages:
#' \enumerate{
#'   \item phrase-level replacements
#'   \item word-level replacements
#' }
#'
#' Matching is case-insensitive. Word-level replacements use word boundaries to
#' reduce accidental replacement inside larger words.
#'
#' @param terms Character vector of terms to shorten.
#' @param extra_phrases Optional named character vector of additional
#'   phrase-level replacements. Names are patterns to match, values are
#'   replacements.
#' @param extra_words Optional named character vector of additional word-level
#'   replacements. Names are words or short expressions to match, values are
#'   replacements.
#' @param extra_first Logical; if `TRUE`, user-supplied replacements are applied
#'   before the built-in dictionary. If `FALSE` (default), built-in replacements
#'   are applied first.
#'
#' @return A character vector of the same length as `terms`.
#'
#' @details
#' Matching is case-insensitive, but replacement text is inserted exactly as
#' provided in the replacement dictionary.
#'
#' Phrase replacements should generally be more specific than word replacements.
#' For example, `"transmembrane"` should be replaced before `"membrane"`.
#'
#' `NA` values are preserved.
#'
#' @examples
#' x <- c(
#'   "transmembrane transport",
#'   "antigen processing",
#'   "positive regulation of process",
#'   "Positive Regulation of Antigen Processing",
#'   NA
#' )
#'
#' shorten_terms(x)
#'
#' shorten_terms(
#'   x,
#'   extra_words = c("antigen" = "antig."),
#'   extra_first = TRUE
#' )
#'
#' @export
shorten_terms <- function(
    terms,
    extra_phrases = NULL,
    extra_words = NULL,
    extra_first = FALSE
) {
  if (!is.character(terms)) {
    stop("`terms` must be a character vector.")
  }
  
  validate_replacements <- function(x, arg_name) {
    if (is.null(x)) {
      return(x)
    }
    
    if (!is.character(x)) {
      stop(sprintf("`%s` must be a named character vector.", arg_name))
    }
    
    if (is.null(names(x)) || anyNA(names(x)) || any(names(x) == "")) {
      stop(sprintf(
        "`%s` must be a named character vector with non-empty names.",
        arg_name
      ))
    }
    
    x
  }
  
  extra_phrases <- validate_replacements(extra_phrases, "extra_phrases")
  extra_words <- validate_replacements(extra_words, "extra_words")
  
  # Phrase-level replacements:
  # keep the most specific / longest phrases first
  phrase_replacements <- c(
    "positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" =
      "Pos. reg. of adap. imm. resp. via IgSF recep. recomb.",
    "antigen processing and presentation of endogenous peptide antigen via MHC class I" =
      "antigen proc. & pres. of endog. peptide antigen via MHC-I",
    "antigen processing and presentation of exogenous peptide antigen via MHC class II" =
      "antigen proc. & pres. of exog. peptide antigen via MHC-II",
    "endoplasmic reticulum" = "ER",
    "tumor necrosis factor" = "TNF",
    "immunoglobulin superfamily" = "IgSF",
    "amino acid" = "AA"
  )
  
  # Word-level replacements:
  # put more specific terms before more general terms
  word_replacements <- c(
    "NF-kappaB" = "NF-kB",
    "transmembrane" = "transmembr.",
    "membrane" = "membr.",
    "transporter" = "transp.",
    "transport" = "transp.",
    "regulation" = "reg.",
    "positive" = "pos.",
    "negative" = "neg.",
    "processing" = "proc.",
    "process" = "proc.",
    "presentation" = "pres.",
    "endogenous" = "endog.",
    "exogenous" = "exog.",
    "response" = "resp.",
    "activity" = "act.",
    "protein" = "prot.",
    "cellular" = "cell.",
    "organization" = "org.",
    "development" = "dev.",
    "apoptotic" = "apop.",
    "metabolic" = "metab.",
    "pathway" = "path.",
    "superfamily" = "SF",
    "signaling" = "sign.",
    "production" = "prod.",
    "mediated" = "med.",
    "receptor" = "recep.",
    "recombination" = "recomb.",
    "antigen" = "AG",
    "detoxification" = "detox.",
    "compound" = "compd.",
    "intermediate" = "intermed.",
    "compartment" = "compartm.",
    "dehydrogenase" = "DH",
    "synthetic" = "synth.",
    "transduction" = "transd.",
    "expression" = "expr."
  )
  
  if (isTRUE(extra_first)) {
    phrase_replacements <- c(extra_phrases, phrase_replacements)
    word_replacements <- c(extra_words, word_replacements)
  } else {
    phrase_replacements <- c(phrase_replacements, extra_phrases)
    word_replacements <- c(word_replacements, extra_words)
  }
  
  shortened <- terms
  not_na <- !is.na(shortened)
  
  # Escape regex metacharacters in user-supplied or built-in patterns
  escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }
  
  # Apply phrase replacements first
  for (pattern in names(phrase_replacements)) {
    safe_pattern <- escape_regex(pattern)
    
    shortened[not_na] <- gsub(
      pattern = safe_pattern,
      replacement = phrase_replacements[[pattern]],
      x = shortened[not_na],
      ignore.case = TRUE,
      perl = TRUE
    )
  }
  
  # Apply word replacements second, with word boundaries
  for (pattern in names(word_replacements)) {
    safe_pattern <- escape_regex(pattern)
    
    shortened[not_na] <- gsub(
      pattern = paste0("\\b", safe_pattern, "\\b"),
      replacement = word_replacements[[pattern]],
      x = shortened[not_na],
      ignore.case = TRUE,
      perl = TRUE
    )
  }
  
  shortened
}

