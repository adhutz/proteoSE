# UniProt keyword queries
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.




#' Retrieve UniProt Keywords by Search Term
#'
#' Queries the UniProt Keywords REST API for one or more search terms and
#' returns a tidy data frame of matching keyword entries.  Hits from
#' multiple terms are collapsed so that each unique keyword ID appears only
#' once, with all matching query terms recorded in a single
#' semicolon-separated field.
#'
#' @param terms Character vector.  One or more search terms to look up
#'   against the UniProt Keywords endpoint.
#' @param limit Integer scalar.  Maximum number of results to retrieve
#'   \emph{per term}.  Default: \code{20}.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per unique keyword
#'   ID and the following columns:
#' \describe{
#'   \item{\code{id}}{UniProt keyword accession (e.g. \code{"KW-0597"}).}
#'   \item{\code{label}}{Human-readable keyword name.}
#'   \item{\code{category}}{Keyword category as returned by the API
#'     (e.g. \code{"Biological process"}).}
#'   \item{\code{synonyms}}{Always \code{NA_character_}; reserved for
#'     future enrichment.}
#'   \item{\code{definition}}{Free-text definition, or \code{NA} if the
#'     API does not provide one.}
#'   \item{\code{source}}{Fixed value \code{"UniProtKW"}.}
#'   \item{\code{query}}{Semicolon-separated list of the input
#'     \code{terms} that matched this keyword.}
#' }
#' Returns \code{NULL} (invisibly, via \code{map_dfr} empty-bind) if no
#' hits are found for any term.
#'
#' @examples
#' \dontrun{
#' kws <- find_kws(c("glycolysis", "mitochondrion"), limit = 10)
#' kws
#' }
#'
#' @importFrom httr2 request req_url_query req_perform resp_body_json
#' @importFrom purrr map_dfr
#' @importFrom dplyr group_by summarise
#' @importFrom tibble tibble
#'
#' @export
find_kws <- function(terms, limit = 20) {
  
  kw_search_url <- "https://rest.uniprot.org/keywords/search"
  
  # -------------------------------------------------------------------------
  # 1.  Query the API for each term and flatten into a single data frame
  # -------------------------------------------------------------------------
  results <- purrr::map_dfr(terms, function(term) {
    
    kw_hits <- httr2::request(kw_search_url) |>
      httr2::req_url_query(query = term, format = "json", size = limit) |>
      httr2::req_perform() |>
      httr2::resp_body_json() |>
      _$results   # extract the results list from the parsed JSON
    
    if (length(kw_hits) == 0) return(NULL)   # skip terms with no hits
    
    purrr::map_dfr(kw_hits, function(hit) {
      tibble::tibble(
        query      = term,
        source     = "UniProtKW",
        id         = hit$keyword$id,
        label      = hit$keyword$name,
        category   = hit$category$name,
        synonyms = dplyr::na_if(paste(hit$synonyms %||% character(0), collapse = "; "),""),
        definition = hit$definition %||% NA_character_
      )
    })
  })
  
  # -------------------------------------------------------------------------
  # 2.  Collapse duplicate keyword IDs
  #     A keyword may be returned by more than one query term; merge those
  #     rows so each ID appears once, with all matching terms recorded.
  # -------------------------------------------------------------------------
  results |>
    dplyr::group_by(id, label, category, synonyms, definition, source) |>
    dplyr::summarise(
      query      = paste(unique(query), collapse = "; "),
      .groups    = "drop"
    )
}


################
################
#' Fetch UniProt Accessions for a Set of Keywords
#'
#' For each keyword in a data frame produced by \code{\link{find_kws}},
#' queries UniProtKB via \code{\link[uniprotREST]{uniprot_search}} and
#' returns either a flat gene-keyword mapping or the full per-keyword
#' result list.
#'
#' @param keywords A tibble as returned by \code{find_kws}, containing at
#'   least columns \code{id} and \code{label}.
#' @param proteome Character scalar. UniProt proteome ID to restrict the
#'   search.  Default: \code{"UP000005640"} (Human).
#' @param fields Character vector of UniProt field names to retrieve.
#'   Default covers accession, gene names, organism, GO terms, keywords,
#'   subcellular location, length, and mass.
#' @param flatten Logical. If \code{TRUE} (default), returns a single
#'   two-column tibble of \code{gene}/\code{keyword} pairs suitable for
#'   row annotation.  If \code{FALSE}, returns a named list of per-keyword
#'   tibbles for inspection or further processing.
#' @param shorten Logical. If \code{TRUE} (default) and \code{flatten =
#'   TRUE}, applies \code{shorten_terms()} to the keyword column.
#'
#' @return See \code{flatten}: either a named list of tibbles or a flat
#'   gene-keyword tibble.
#'
#' @importFrom purrr map set_names
#' @importFrom dplyr bind_rows select mutate rename_with
#' @importFrom tibble as_tibble
#' @importFrom janitor clean_names
#' @importFrom uniprotREST uniprot_search
#'
#' @export
fetch_kw_accessions <- function(
    keywords,
    proteome = "UP000005640",
    fields   = c(
      "accession", "id", "gene_names", "gene_primary", "protein_name",
      "organism_id", "organism_name", "reviewed",
      "go_id", "go", "go_p", "go_f", "go_c",
      "keyword", "keywordid",
      "cc_subcellular_location", "length", "mass"
    ),
    flatten  = TRUE,
    shorten  = TRUE
) {
  
  # -------------------------------------------------------------------------
  # 1.  Fetch and clean results for a single keyword row
  # -------------------------------------------------------------------------
  .fetch_one <- function(id_val, label_val) {
    uniprotREST::uniprot_search(
      query    = sprintf(
        "(proteome:%s) AND (keyword:%s) AND (reviewed:true)",
        proteome, id_val
      ),
      database = "uniprotkb",
      format   = "tsv",
      fields   = fields
    ) |>
      tibble::as_tibble() |>
      janitor::clean_names() |>           # lowercases + snake_cases in one step
      dplyr::mutate(
        uniprot_kw_id    = id_val,
        uniprot_kw_label = label_val
      ) |>
      dplyr::select(uniprot_kw_label, uniprot_kw_id, dplyr::everything())
  }
  
  # -------------------------------------------------------------------------
  # 2.  Map over all keywords, naming the list by label
  # -------------------------------------------------------------------------
  results <- purrr::map2(
    keywords$id,
    keywords$label,
    .fetch_one
  ) |>
    purrr::set_names(keywords$label)
  
  # -------------------------------------------------------------------------
  # 3.  Optionally flatten to a gene-keyword annotation tibble
  # -------------------------------------------------------------------------
  if (!flatten) return(results)
  
  flat <- dplyr::bind_rows(results) |>
    dplyr::select(gene = gene_names_primary, keyword = uniprot_kw_label)
  
  if (shorten) flat <- dplyr::mutate(flat, keyword = shorten_terms(keyword))
  
  flat
}

