# compare_gold.R
#
# Side analysis: compare a run's changepoints with a gold standard onset
# (e.g. ROMP for Ethiopia), one point per grid cell and year.
#
# This file is sourced by onset_setup() along with the rest of R/, so it only
# defines functions.  Example (from the base directory):
#
#   source("R/setup.R"); onset_setup(".")
#   run <- load_run("ethiopia", "ethiopia_1998-2024_doy1-240_onset")
#   compare_gold_chpt(run)                                  # to the screen
#   compare_gold_chpt(run, file = "gold_vs_chpt_ethiopia.pdf")
#   compare_gold_chpt(run, bin = 2)                         # finer bins
#   compare_gold_chpt(run, type = "points")                 # scatter plot
#
#   # or step by step
#   compare_gold_chpt(run, dev_range = c(0.3, Inf))         # dev_expl >= 0.3
#   d <- gold_chpt_pairs(run, lat_range = c(8, Inf), lon_range = c(-Inf, 40))
#   plot_gold_vs_chpt(d, main = "lat >= 8, lon <= 40")


# One row per cell-year with both a gold standard date and a changepoint.
#
# run        : from load_run() (must include the gold standard array)
# gold       : name of the gold standard array in run$r
# lon_range, lat_range : c(min, max), inclusive; NULL = no restriction
# dev_range  : c(min, max) for the cell-year's dev_expl (deviance explained
#              by the changepoint model), inclusive; NULL = no restriction.
#              Use -Inf / Inf for a one-sided bound, e.g. c(0.3, Inf).
# years      : years to use (NULL = all years in the run)
#
# Returns a data frame with lon, lat, year, gold, chpt, dev_expl.
gold_chpt_pairs <- function(run, gold = "gold_romp",
                            lon_range = NULL, lat_range = NULL,
                            dev_range = NULL, years = NULL) {
  if (is.null(run$r[[gold]]))
    stop("No array '", gold, "' in the run; available gold standards: ",
         paste(run$gold, collapse = ", "), call. = FALSE)
  d <- data.frame(
    lon      = run$lon[slice.index(run$r$chpt, 1)],
    lat      = run$lat[slice.index(run$r$chpt, 2)],
    year     = run$years[slice.index(run$r$chpt, 3)],
    gold     = as.vector(run$r[[gold]]),
    chpt     = as.vector(run$r$chpt),
    dev_expl = as.vector(run$r$dev_expl)
  )
  keep <- is.finite(d$gold) & is.finite(d$chpt)
  if (!is.null(lon_range)) keep <- keep & d$lon >= lon_range[1] & d$lon <= lon_range[2]
  if (!is.null(lat_range)) keep <- keep & d$lat >= lat_range[1] & d$lat <= lat_range[2]
  if (!is.null(dev_range))
    keep <- keep & is.finite(d$dev_expl) &
      d$dev_expl >= dev_range[1] & d$dev_expl <= dev_range[2]
  if (!is.null(years))     keep <- keep & d$year %in% years
  d <- d[keep, ]
  rownames(d) <- NULL
  d
}

# Changepoint (y) against gold standard (x), with the 1:1 line.  The line
# under the title gives the number of points, the correlation and the median
# and mean absolute difference (changepoint - gold, in days).
#
# d        : from gold_chpt_pairs()
# type     : "density" (default): number of cell-years in square bins of
#            `bin` x `bin` days, colour on a log scale, empty bins blank;
#            "points": one semi-transparent point per cell-year
# bin      : bin width in days (density only).  Values are whole days, and
#            bin edges fall half-way between days, so each day lies in
#            exactly one bin.
# lim      : common axis limits (NULL = range of both variables)
# zlim     : count limits for the colour scale (density only; NULL = 1 to
#            the largest count).  Give the same zlim to compare panels.
# xlab, ylab, main : labels (main goes above the summary line)
# col      : point colour (points only)
# ...      : passed to plot() (points only)
plot_gold_vs_chpt <- function(d, type = c("density", "points"), bin = 4,
                              lim = NULL, zlim = NULL,
                              xlab = "Gold standard onset (day of year)",
                              ylab = "Changepoint (day of year)",
                              main = "", col = "#1F4E7930", ...) {
  type <- match.arg(type)
  if (!nrow(d)) {
    plot.new(); title(main = main); mtext("no cell-years with both values")
    return(invisible(NULL))
  }
  if (is.null(lim)) lim <- range(d$gold, d$chpt)

  if (type == "points") {
    plot(d$gold, d$chpt, pch = 16, cex = 0.4, col = col,
         xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, main = "",
         asp = 1, ...)
  } else {
    b <- gold_chpt_bins(d, bin, lim)
    if (is.null(zlim)) zlim <- c(1, max(b$counts, 1, na.rm = TRUE))
    z <- log10(pmin(b$counts, zlim[2]))
    # Colour bar ticks at 1, 2, 5, 10, 20, 50, ... within zlim
    ticks <- outer(c(1, 2, 5), 10^(0:6))
    ticks <- sort(ticks[ticks >= zlim[1] & ticks <= zlim[2]])
    old_cex <- par()$cex  # image.plot can reset cex
    on.exit(par(cex = old_cex), add = TRUE)
    image.plot(b$mids, b$mids, z, zlim = log10(zlim),
               col = hcl.colors(64, "viridis"),
               xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, main = "",
               asp = 1, useRaster = TRUE,
               legend.lab = sprintf("cell-years per bin (%d days)", bin),
               legend.line = 2.6, legend.cex = 0.8, legend.mar = 7,
               axis.args = list(at = log10(ticks), labels = ticks))
    box()
  }
  title(main = main, line = 1.5, cex.main = 1.05)
  abline(0, 1, col = if (type == "density") "black" else "red", lwd = 1.5,
         lty = if (type == "density") 2 else 1)
  diffs <- d$chpt - d$gold
  mtext(sprintf(paste0("n = %d; r = %.2f; ",
                       "chpt - gold: median %.0f d, mean |diff| %.0f d"),
                nrow(d), cor(d$gold, d$chpt), median(diffs), mean(abs(diffs))),
        side = 3, line = 0.3, cex = 0.7)
  invisible(d)
}

# Counts of cell-years in square bins (used by plot_gold_vs_chpt()).
# Returns list(mids = bin centres, counts = matrix [gold bin x chpt bin],
# NA for empty bins).
gold_chpt_bins <- function(d, bin, lim) {
  lo <- floor(lim[1] / bin) * bin - 0.5
  hi <- ceiling((lim[2] + 1) / bin) * bin - 0.5
  breaks <- seq(lo, hi, by = bin)
  mids <- breaks[-1] - bin / 2
  ix <- findInterval(d$gold, breaks, rightmost.closed = TRUE)
  iy <- findInterval(d$chpt, breaks, rightmost.closed = TRUE)
  counts <- matrix(0, length(mids), length(mids))
  tab <- table(factor(ix, seq_along(mids)), factor(iy, seq_along(mids)))
  counts[] <- tab
  counts[counts == 0] <- NA
  list(mids = mids, counts = counts)
}

# Both plots side by side: all cells, and a sub-region.
# Both panels use the same axis limits and, for
# density plots, the same colour scale.
#
# run      : from load_run()
# gold     : gold standard array name
# sub_lon, sub_lat : sub-region bounds, c(min, max)
# dev_range : keep only cell-years with dev_expl in c(min, max), in both
#            panels (NULL = all).  The axis limits stay those of the
#            unrestricted data, so plots with different dev_range compare
#            directly; the colour scale adapts to the counts kept.
# type, bin : as for plot_gold_vs_chpt() (default: density, 4-day bins)
# file     : PDF to write (NULL = current device)
#
# Returns (invisibly) list(all = , sub = ) of the data frames plotted.
compare_gold_chpt <- function(run, gold = "gold_romp",
                              sub_lon = c(-Inf, Inf), sub_lat = c(-Inf, Inf),
                              dev_range = NULL,
                              type = c("density", "points"), bin = 4,
                              file = NULL) {
  type <- match.arg(type)
  full <- gold_chpt_pairs(run, gold)
  lim  <- range(full$gold, full$chpt)
  all  <- gold_chpt_pairs(run, gold, dev_range = dev_range)
  sub  <- gold_chpt_pairs(run, gold, lon_range = sub_lon, lat_range = sub_lat,
                          dev_range = dev_range)
  # Shared colour scale: the sub-region's bins are a subset of the full
  # region's, so the full region has the largest count
  zlim <- if (type == "density" && nrow(all))
    c(1, max(gold_chpt_bins(all, bin, lim)$counts, na.rm = TRUE)) else NULL

  if (!is.null(file)) {
    pdf(file, width = 11, height = 5.8)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(mfrow = c(1, 2), mai = c(0.8, 0.8, 0.6, 0.2) +
                   if (type == "density") c(0, 0, 0, 1.0) else 0,
                 mgp = c(2, 0.6, 0),
                 oma = c(0, 0, 1.8, if (type == "density") 1.5 else 0))
  on.exit(par(old_par), add = TRUE)

  bounds <- function(nm, r) paste(c(
    if (is.finite(r[1])) sprintf("%s >= %g", nm, r[1]),
    if (is.finite(r[2])) sprintf("%s <= %g", nm, r[2])), collapse = ", ")
  plot_gold_vs_chpt(all, type = type, bin = bin, lim = lim, zlim = zlim,
                    main = "All grid cells")
  plot_gold_vs_chpt(sub, type = type, bin = bin, lim = lim, zlim = zlim,
                    main = paste0(paste(c(bounds("lat", sub_lat),
                                          bounds("lon", sub_lon)),
                                        collapse = ", ")))
  dev_txt <- if (is.null(dev_range)) "" else
    paste0("  |  ", bounds("dev_expl", dev_range))
  title(sprintf("%s: changepoint vs %s%s", run$meta$tag, sub("^gold_", "", gold),
                dev_txt),
        outer = TRUE, line = 0.5, cex.main = 1.1)
  invisible(list(all = all, sub = sub))
}


