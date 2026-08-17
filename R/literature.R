# Literature search query builders (PubMed / Scopus)
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.

#' Construct a PubMed Query String
#'
#' This function constructs a PubMed query string using optional keyword sets,
#' high-impact journal filters, and positive control PMIDs. If no date range is supplied,
#' it defaults to the last \code{years_back} years. All components are optional, but at least one of
#' \code{keywords1} or \code{pmid_ctrl} must be provided.
#'
#' Keywords can be combined using logical operators: \code{"AND"}, \code{"OR"},
#' \code{"COMBINE"} (pairwise combinations of keyword1 + keyword2), or \code{"NOT"}
#' (to exclude keyword2 terms).
#'
#' The resulting query is printed to the console, copied to the clipboard, and optionally
#' opened in your default web browser on PubMed.
#'
#' @param date_range Optional character vector of length 2 for publication date filter.
#' @param years_back Integer; used if \code{date_range} is NULL. Sets default start year as current year minus \code{years_back}. Default is 5.
#' @param authors Optional character vector of author names (e.g., "Rinschen MM").
#' @param author_position Position to match: "any", "first", "last", or "first_last". Default is "any".
#' @param keywords1 Optional character vector of primary keywords.
#' @param keywords2 Optional character vector of secondary keywords.
#' @param operator Optional keywords1 and keywords2 combination operator: "AND", "OR", "COMBINE", "NOT". Default is "AND" if needed.
#' @param high_impact_journals Optional character vector of journal names to restrict the query.
#' @param pmid_ctrl Optional numeric vector of PubMed IDs to include as positive controls.
#' @param open_in_browser Logical; if TRUE, opens the query in PubMed in your default browser. Default is TRUE.
#' @param exclude_reviews Logical; if TRUE, excludes articles with Publication Type "review". Default is TRUE.
#' @return The constructed PubMed query string (also copied to clipboard and optionally opened).
#' @export
#'
#' @examples
#' \dontrun{
#' pubmed_query(keywords1 = c("itaconate", "macrophage"))
#'
#' pubmed_query(
#'   #date_range = c("2020", "3000"),
#'   years_back = 5,
#'   keywords1 = c("itaconate", "itaconylation", "itaconatylation", "itaconation","itaconic acid"),
#'   keywords2 = c("proteomics", "metabolomics", "PTM", "proteome", "metabolite", "metabolome", "modification"),
#'   operator = "AND",
#'   authors = c("Rinschen MM", "Smith AB"),
#'   author_position = "first_last",
#'   high_impact_journals = c("Nature", "Nat Metab", "Nature Immunology", "Nature Medicine", "Nature Communications", "Nat Rev Immunol", 
#'                            "Cell", "Cell Metab", "Cell Reports", "BMJ", "Annu Rev Pathol", "Annu Rev Immunol", "Trends Endocrinol Metab", "Metabolism",
#'                            "Science", "Science Immunology", "Science Translational Medicine",
#'                            "Signal Transduct Target Ther", "Immunity", "J Clin Invest", "EMBO J", "Proc Natl Acad Sci U S A",
#'                            "J Am Soc Nephrol", "Kidney Int", "N Engl J Med", "Circulation", "Nat Rev Nephrol"),
#'   pmid_ctrl = c(40011782, 37271942),
#'   open_in_browser = T
#' )
#' }
pubmed_query <- function(date_range = NULL,
                         years_back = 5,
                         authors = NULL,
                         author_position = "any",
                         keywords1 = NULL,
                         keywords2 = NULL,
                         operator = NULL,
                         high_impact_journals = NULL,
                         pmid_ctrl = NULL,
                         exclude_reviews = FALSE,
                         open_in_browser = TRUE) {
  
  if (!requireNamespace("clipr", quietly = TRUE)) {
    stop("The 'clipr' package is required. Install it with install.packages('clipr').")
  }
  
  if (!is.null(operator) && !operator %in% c("AND", "OR", "COMBINE", "NOT")) {
    stop('operator must be one of "AND", "OR", "COMBINE", or "NOT".')
  }
  

  # Date range
  date_part <- ""
  if (is.null(date_range)) {
    end_year <- as.numeric(format(Sys.Date(), "%Y"))
    start_year <- end_year - years_back
    date_range <- c(as.character(start_year), as.character(end_year))
    message(sprintf("Using default 5-year date range: %s to %s", start_year, end_year))
  }
  
  if (length(date_range) != 2 || !all(sapply(date_range, is.character))) {
    stop("date_range must be a character vector of two elements.")
  }
  
  date_part <- paste0('("', date_range[1], '"[Publication Date] : "', date_range[2], '"[Publication Date]) AND ')
  
  # Keyword logic
  keyword_elements <- character(0)
  if (!is.null(keywords1)) {
    if (is.null(keywords2)) {
      keyword_elements <- paste0('("', keywords1, '")')
    } else {
      op <- ifelse(is.null(operator), "AND", operator)
      if (op == "COMBINE") {
        keyword_elements <- unlist(lapply(keywords1, function(k1) {
          sapply(keywords2, function(k2) paste0('("', k1, ' ', k2, '")'))
        }))
      } else if (op == "AND") {
        k1 <- paste0('("', keywords1, '")', collapse = " OR ")
        k2 <- paste0('("', keywords2, '")', collapse = " OR ")
        keyword_elements <- paste0("(", k1, ") AND (", k2, ")")
      } else if (op == "OR") {
        keyword_elements <- paste0('("', c(keywords1, keywords2), '")')
      } else if (op == "NOT") {
        k1 <- paste0('("', keywords1, '")', collapse = " OR ")
        k2 <- paste0('("', keywords2, '")', collapse = " OR ")
        keyword_elements <- paste0("(", k1, ") NOT (", k2, ")")
      }
    }
  }
  
  search_part <- paste(keyword_elements, collapse = " OR ")
  
  # Journal filter
  journal_part <- ""
  if (!is.null(high_impact_journals)) {
    journal_filters <- paste0('"', high_impact_journals, '"[Journal]')
    journal_part <- paste(journal_filters, collapse = " OR ")
    journal_part <- paste0(" AND (", journal_part, ")")
  }
  
  # Author filtering
  author_part <- ""
  if (!is.null(authors)) {
    if (!author_position %in% c("any", "first", "last", "first_last")) {
      stop('author_position must be one of "any", "first", "last", "first_last"')
    }
    
    author_queries <- switch(author_position,
                             any = paste0('"', authors, '"[Author]'),
                             first = paste0('"', authors, '"[1au]'),
                             last = paste0('"', authors, '"[lastau]'),
                             first_last = paste0('("', authors, '"[1au] OR "', authors, '"[lastau])')
    )
    
    author_part <- paste(author_queries, collapse = " OR ")
    author_part <- paste0(" AND (", author_part, ")")
  }
  
  
  
  
  # Positive control PMIDs
  pmid_part <- ""
  if (!is.null(pmid_ctrl)) {
    pmid_strings <- paste0(pmid_ctrl, '[PMID]')
    pmid_part <- paste(pmid_strings, collapse = " OR ")
    pmid_part <- paste0(" OR (", pmid_part, ")")
  }
  
  # Exclude reviews if requested
  review_exclude <- ""
  if (exclude_reviews) {
    review_exclude <- ' NOT "review"[Publication Type]'
  }
  
  # Final query
  # Build the main blocks
  main_blocks <- list()
  
  # Add date part
  main_blocks <- c(main_blocks, date_part)
  
  # Add keyword search block only if it exists
  if (nzchar(search_part)) {
    main_blocks <- c(main_blocks, paste0("(", search_part, ")"))
  }
  
  # Add optional components
  optional_blocks <- c(journal_part, author_part, pmid_part, review_exclude)
  optional_blocks <- optional_blocks[nzchar(optional_blocks)]
  
  # Combine everything
  query <- paste(c(main_blocks, optional_blocks), collapse = " ")
  
  
  if (any(!nzchar(c(main_blocks, optional_blocks)))) {
    warning("Some query components are empty and have been removed.")
  }
  
  
  # Output
  cat("\nPubMed query:\n", query, "\n\n")
  clipr::write_clip(query)
  message("PubMed query copied to clipboard.")
  
  # Open in browser
  if (open_in_browser) {
    query_url <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", utils::URLencode(query, reserved = TRUE))
    utils::browseURL(query_url)
    message("PubMed opened in browser.")
  }
  
  message("Link to journal abbreviations on PubMed: https://www.ncbi.nlm.nih.gov/nlmcatalog")
  message("Ranked journals here: https://wos-journal.info/")
  message("Or from UKE intranet: https://jcr.clarivate.com/jcr/home?app=jcr&Init=Yes&authCode=null&SrcApp=IC2LS")
  
  
  
  return(invisible(query))
}




############################################################################
############################################################################
############################################################################
#' Construct a Scopus Query String
#'
#' This function constructs a Scopus query string based on the provided
#' publication date range and combinations of two sets of keywords.
#'
#' @param date_range A character vector of length 2 indicating the start and end publication years.
#' @param keywords1 A character vector of the first set of keywords.
#' @param keywords2 A character vector of the second set of keywords.
#' @return A character string representing the Scopus query.
#' @examples
#' date_range <- c("2015", "2025")
#' keywords1 <- c("kidney", "renal")
#' keywords2 <- c("organoid", "spheroid", "tubuloid", "organ-on-chip", "microfluidic")
#' query <- scopus_query(date_range, keywords1, keywords2)
#' cat(query)
#' @export
scopus_query <- function(date_range, keywords1, keywords2) {
  # Validate input
  if (length(date_range) != 2) {
    stop("date_range should be a vector of two elements (start and end year).")
  }
  if (!all(sapply(date_range, is.character))) {
    stop("date_range should be a vector of characters.")
  }
  
  # Construct date part
  date_part <- paste0('PUBYEAR AFT ', date_range[1], ' AND PUBYEAR BEF ', date_range[2])
  
  # Construct keyword combinations
  keyword_combinations <- unlist(lapply(keywords1, function(k1) {
    sapply(keywords2, function(k2) {
      paste0('TITLE-ABS-KEY("', k1, ' ', k2, '")')
    })
  }))
  
  # Combine keyword combinations with OR
  keyword_part <- paste(keyword_combinations, collapse = " OR ")
  
  # Combine date part and keyword part
  query <- paste0(date_part, " AND (", keyword_part, ")")
  
  return(query)
}

