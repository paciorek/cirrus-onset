# maps.R
#
# Map specs, summaries across years, and plotting.  Shared by make_maps and
# the Shiny app.
#
# A map spec has
#   label   title
#   needs   array names required (the map is skipped if any is missing)
#   fun     function(r) returning a lon x lat x year array from the list of
#           arrays r (so derived quantities are computed here, not stored)
#   zlim    list(mean = , sd = ) colour limits; NULL = automatic.  Region
#           config can override these (cfg$zlim[[direction]][[name]]).

suppressPackageStartupMessages({
  library(fields)
  library(maps)
})

map_specs <- function() {
  list(
    chpt = list(
      label = "Changepoint day of year", needs = "chpt",
      fun = function(r) r$chpt, zlim = list(mean = NULL, sd = NULL)),
    dev_expl = list(
      label = "Deviance explained by changepoint model", needs = "dev_expl",
      fun = function(r) r$dev_expl, zlim = list(mean = c(0, 1), sd = NULL)),
    dev_full = list(
      label = "Deviance explained by spline", needs = "dev_full",
      fun = function(r) r$dev_full, zlim = list(mean = c(0, 1), sd = NULL)),
    dev_rel = list(
      label = "Relative deviance: changepoint / spline",
      needs = c("dev_expl", "dev_full"),
      fun = function(r) r$dev_expl / r$dev_full,
      zlim = list(mean = c(0, 1), sd = c(0, 0.4))),
    pre_mean = list(
      label = "Mean daily precip before changepoint (mm/day)", needs = "pre_mean",
      fun = function(r) r$pre_mean, zlim = list(mean = NULL, sd = NULL)),
    post_mean = list(
      label = "Mean daily precip after changepoint (mm/day)", needs = "post_mean",
      fun = function(r) r$post_mean, zlim = list(mean = NULL, sd = NULL)),
    post_frac = list(
      label = "Post-changepoint precip / window total",
      needs = c("pre_total", "post_total"),
      fun = function(r) r$post_total / (r$pre_total + r$post_total),
      zlim = list(mean = c(0, 1), sd = NULL)),
    transition_sd = list(
      label = "Transition sharpness: 7-day running-mean SD",
      needs = "transition_sd",
      fun = function(r) r$transition_sd, zlim = list(mean = NULL, sd = NULL)),
    transition_edf = list(
      label = "Transition complexity: GAM smooth EDF",
      needs = "transition_edf",
      fun = function(r) r$transition_edf, zlim = list(mean = NULL, sd = NULL))
  )
}

# Map specs for gold standard datasets (array names prefixed "gold_").
# They share the changepoint map's colour limits so the two compare directly.
gold_map_specs <- function(cfg, direction, chpt_zlim) {
  gs <- gold_standards_for(cfg, direction)
  if (!length(gs)) return(list())
  specs <- lapply(names(gs), function(nm) {
    key <- paste0("gold_", nm)
    list(label = paste("Gold standard:", gs[[nm]]$label, "(day of year)"),
         needs = key,
         fun = eval(bquote(function(r) r[[.(key)]])),
         zlim = chpt_zlim)
  })
  setNames(specs, paste0("gold_", names(gs)))
}

# All specs for a run, with region overrides applied and specs whose inputs
# are missing dropped.
specs_for_run <- function(r, cfg, direction) {
  specs <- map_specs()
  over  <- cfg$zlim[[direction]]
  for (nm in intersect(names(over), names(specs)))
    specs[[nm]]$zlim <- modifyList(specs[[nm]]$zlim, over[[nm]])
  specs <- c(specs, gold_map_specs(cfg, direction, specs$chpt$zlim))
  Filter(function(s) all(s$needs %in% names(r)), specs)
}

# Summaries over the year dimension.  Add entries to map more summaries.
year_summaries <- list(
  mean = list(label = "mean",    fun = function(x) mean(x, na.rm = TRUE)),
  sd   = list(label = "std dev", fun = function(x) sd(x, na.rm = TRUE))
)

summarize_years <- function(arr, summary = "mean") {
  f <- year_summaries[[summary]]$fun
  z <- suppressWarnings(apply(arr, c(1, 2), f))
  z[!is.finite(z)] <- NA
  z
}

# Read the region's optional shapefile overlay once (NULL if none or if sf is
# unavailable).
load_overlay <- function(cfg) {
  if (is.null(cfg$shapefile)) return(NULL)
  f <- resolve_path(cfg$shapefile)
  if (!file.exists(f)) {
    warning("Shapefile not found, skipping overlay: ", f)
    return(NULL)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    warning("Package 'sf' not installed, skipping shapefile overlay.")
    return(NULL)
  }
  sf::st_geometry(sf::st_read(f, quiet = TRUE))
}

# One map panel.  Values outside zlim are clamped to the end colours.  On
# return the plot's user coordinates are lon/lat (useful for Shiny clicks).
plot_map_panel <- function(z, lon, lat, main = "", zlim = NULL, col_n = 32,
                           map_regions = NULL, overlay = NULL, mark = NULL) {
  z[!is.finite(z)] <- NA
  if (is.null(zlim)) zlim <- if (all(is.na(z))) c(0, 1) else range(z, na.rm = TRUE)
  if (!all(is.finite(zlim)) || zlim[1] == zlim[2]) zlim <- zlim[1] + c(-0.5, 0.5)
  if (!all(is.finite(zlim))) zlim <- c(0, 1)
  z <- pmin(pmax(z, zlim[1]), zlim[2])

  old_cex <- par()$cex  # `cex` being changed to 1 in image.plot somehow.
  on.exit(par(cex = old_cex))
  image.plot(x = lon, y = lat, z = z, col = tim.colors(col_n), zlim = zlim,
             xlab = "", ylab = "", main = main, useRaster = TRUE, asp = 1)
  if (!is.null(map_regions))
    maps::map("world", regions = map_regions, add = TRUE, col = "black", lwd = 1.2)
  if (!is.null(overlay))
    plot(overlay, add = TRUE, border = "gray30", lwd = 0.4)
  box()
  if (!is.null(mark)) {
    points(mark[1], mark[2], pch = 3, cex = 3, lwd = 4, col = "white")
    points(mark[1], mark[2], pch = 3, cex = 3, lwd = 2, col = "magenta")
  }
  invisible(NULL)
}

# A row of panels, one per summary, with an outer title.
plot_map_summaries <- function(arr, lon, lat, title, zlim = list(),
                               summaries = names(year_summaries),
                               map_regions = NULL, overlay = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, length(summaries)), mai = c(0.4, 0.3, 0.5, 0.35),
      mgp = c(1, 0.5, 0), cex = 0.7, cex.main = 1.2, omi = c(0, 0, 0.3, 0))
  for (s in summaries)
    plot_map_panel(summarize_years(arr, s), lon, lat,
                   main = year_summaries[[s]]$label, zlim = zlim[[s]],
                   map_regions = map_regions, overlay = overlay)
  title(main = title, outer = TRUE, line = 0.5)
}
