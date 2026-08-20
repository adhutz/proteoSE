# On-disk cache for the functions that hit an external API.
#
# The memoisation itself happens in .onLoad() (R/zzz.R): memoise() has to run at
# load time, or the cache object would be serialised into the installed package.

# Functions rebound to memoised versions at load. Only real network fetches
# belong here -- pubmed_query() and scopus_query() build query strings offline,
# and pubmed_query() writes the clipboard and opens a browser, so caching them
# would be wrong rather than merely useless.
.cached_fns <- c("get_network", "genes_from_kegg", "find_kws", "fetch_kw_accessions")

# Option override exists so the tests can point at a temp dir.
.cache_dir <- function() {
  getOption("proteoSE.cache_dir", tools::R_user_dir("proteoSE", "cache"))
}

# ponytail: one 30-day TTL for every service. Add a per-service ttl only if a
# stale KEGG or UniProt answer actually bites.
.proteoSE_cache <- function() {
  cachem::cache_disk(
    dir      = .cache_dir(),
    max_age  = 30 * 24 * 3600,
    max_size = 512 * 1024^2
  )
}

.cache_files <- function() list.files(.cache_dir(), full.names = TRUE)

.fmt_bytes <- function(n) {
  format(structure(as.numeric(n), class = "object_size"), units = "auto")
}


#' Inspect the proteoSE API cache
#'
#' Reports where the on-disk cache of external API responses lives and what it
#' currently holds. See the "Caching" section of [proteoSE-package].
#'
#' @return Invisibly, a list with elements `dir`, `n` (number of entries),
#'   `size` (total bytes), `oldest` and `newest` (entry modification times, or
#'   `NA` when the cache is empty). The same information is printed as a
#'   message.
#' @seealso [proteoSE_cache_clear()]
#' @export
#' @examples
#' proteoSE_cache_info()
proteoSE_cache_info <- function() {
  dir   <- .cache_dir()
  files <- .cache_files()
  info  <- file.info(files)

  out <- list(
    dir    = dir,
    n      = length(files),
    size   = sum(info$size, na.rm = TRUE),
    oldest = if (length(files)) min(info$mtime) else as.POSIXct(NA),
    newest = if (length(files)) max(info$mtime) else as.POSIXct(NA)
  )

  message("proteoSE cache: ", out$n, " entries, ", .fmt_bytes(out$size), "\n  ", dir)
  if (out$n) {
    message("  oldest: ", format(out$oldest), "   newest: ", format(out$newest))
  }
  invisible(out)
}


#' Empty the proteoSE API cache
#'
#' Deletes every cached API response. The next call to a network function
#' fetches fresh data.
#'
#' @return Invisibly, the number of bytes freed.
#' @seealso [proteoSE_cache_info()]
#' @export
#' @examples
#' \dontrun{
#' proteoSE_cache_clear()
#' }
proteoSE_cache_clear <- function() {
  files <- .cache_files()
  freed <- sum(file.info(files)$size, na.rm = TRUE)
  unlink(files, recursive = TRUE)
  message("Cleared ", length(files), " cache entries (", .fmt_bytes(freed), ") from ",
          .cache_dir())
  invisible(freed)
}
