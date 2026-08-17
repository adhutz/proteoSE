#' Flatten point layers of a ggplot for lightweight vector export
#'
#' Volcano plots (and other scatter-heavy ggplots) become bloated as PDF/EPS
#' because every point is stored as its own vector path. Illustrator then has
#' to handle thousands of objects, which tanks performance. This function
#' rasterises only the point-type layer(s) (geom_point/geom_jitter by default)
#' at a high DPI, while leaving axes, titles, legends, and text labels
#' (e.g. ggrepel gene names) as normal editable vector objects.
#'
#' At dpi >= 600 the rasterised points are visually indistinguishable from
#' vector points in print, but the exported file collapses from thousands of
#' paths down to one embedded image for that layer.
#'
#' Requires the `ggrastr` package: install.packages("ggrastr")
#'
#' @param plot A ggplot object
#' @param dpi Raster resolution for the flattened layer(s). 600 is a safe
#'   default for print figures; use 900-1200 if you need to zoom in a lot.
#' @param layer_geoms Character vector of geom classes to rasterise.
#'   Defaults to point-type geoms only -- text/labels are left vector on
#'   purpose so gene names etc. stay editable.
#' @return A ggplot object, identical in appearance, with the target layers
#'   replaced by rasterised equivalents.
#' @export
flatten_plot_points <- function(plot,
                                 dpi = 600,
                                 layer_geoms = c("GeomPoint", "GeomJitter")) {
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    stop("Please install ggrastr first: install.packages('ggrastr')")
  }
  if (!inherits(plot, "ggplot")) stop("plot must be a ggplot object")

  is_target <- vapply(plot$layers, function(l) class(l$geom)[1] %in% layer_geoms, logical(1))

  if (!any(is_target)) {
    warning("No layers matching ", paste(layer_geoms, collapse = ", "),
            " found -- returning plot unchanged")
    return(plot)
  }

  plot$layers[is_target] <- lapply(plot$layers[is_target], ggrastr::rasterise, dpi = dpi)
  plot
}

#' Flatten and save a ggplot as an Illustrator-friendly PDF
#'
#' Convenience wrapper: flattens point layers with flatten_plot_points(),
#' then saves via cairo_pdf so the remaining vector elements (text, axes,
#' legend) stay crisp and editable.
#'
#' @param plot A ggplot object
#' @param filename Output path, e.g. "volcano.pdf"
#' @param dpi Raster resolution for point layer(s). Default 600.
#' @param width,height Figure size in inches (matches ggsave defaults)
#' @param layer_geoms Passed to flatten_plot_points()
#' @return Invisibly, the flattened ggplot object.
#' @export
save_flattened_plot <- function(plot, filename, dpi = 600,
                                 width = 7, height = 6,
                                 layer_geoms = c("GeomPoint", "GeomJitter")) {
  flat <- flatten_plot_points(plot, dpi = dpi, layer_geoms = layer_geoms)
  ggplot2::ggsave(filename, flat, width = width, height = height,
                   device = grDevices::cairo_pdf)
  invisible(flat)
}

# --- Example usage on a typical volcano plot ---------------------------------
# library(ggplot2); library(ggrepel)
#
# p <- ggplot(res, aes(log2FC, -log10(pval), color = sig)) +
#   geom_point(size = 1) +
#   geom_text_repel(data = subset(res, label != ""), aes(label = label)) +
#   theme_bw()
#
# # Points get rasterised at 600 dpi; geom_text_repel labels stay vector.
# save_flattened_plot(p, "volcano_flat.pdf", dpi = 600)
#
# # Or just get the modified ggplot object back and inspect/save yourself:
# p_flat <- flatten_plot_points(p, dpi = 600)
# p_flat
