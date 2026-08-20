.onLoad <- function(libname, pkgname) {
  # Cache is a convenience: no memoise/cachem, or opted out, means the network
  # functions stay exactly as defined.
  if (!isTRUE(getOption("proteoSE.cache", TRUE))) return(invisible())
  if (!requireNamespace("memoise", quietly = TRUE) ||
      !requireNamespace("cachem", quietly = TRUE)) return(invisible())

  ns    <- topenv(environment())
  cache <- .proteoSE_cache()
  for (fn in .cached_fns) {
    assign(fn, memoise::memoise(get(fn, envir = ns), cache = cache), envir = ns)
  }
  invisible()
}
