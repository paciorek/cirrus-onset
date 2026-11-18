# plot_utils.R
#
# Time series plots for one grid cell.  Map panels use the shared
# plot_map_panel() and summarize_years() in ../R/maps.R.
#
#   spline_fits()           -- GAM spline fit to each year's window
#   plot_location_panels()  -- one panel per year
#   plot_spline_overlay()   -- all years' spline fits in one panel
#   nearest_spline_year()   -- the year whose fit passes nearest a click
#   panel_at()              -- the year of the panel at a click position

suppressPackageStartupMessages(library(mgcv))

# Fit the spline (Gamma GAM with log link, s(t, k = spline_k)) to each
# year's analysis window, on values floored at gamma_floor as in
# metric_dev_full().
#
# rain : [day x year] matrix of daily precip for one cell
#
# Returns a [window day x year] matrix of fitted values (rows = doy_start to
# doy_end), NA for years with missing days in the window or a failed fit.
spline_fits <- function(rain, doy_start, doy_end, spline_k = 20,
                        gamma_floor = 1) {
  n_window <- doy_end - doy_start + 1L
  out <- matrix(NA_real_, n_window, ncol(rain),
                dimnames = list(doy = seq(doy_start, doy_end),
                                year = colnames(rain)))
  t <- seq_len(n_window)
  for (k in seq_len(ncol(rain))) {
    w <- extract_window(as.numeric(rain[, k]), doy_start, doy_end)
    if (!all(is.finite(w))) next
    tmp <- pmax(w, gamma_floor)
    fit <- try(suppressWarnings(
      gam(tmp ~ s(t, k = spline_k), family = Gamma(link = "log"))),
      silent = TRUE)
    if (!inherits(fit, "try-error"))
      out[, k] <- predict(fit, type = "response")
  }
  out
}

# Plot one panel per year for a single grid cell: the full annual series in
# gray, the analysis window shaded, the spline fit (red), the changepoint
# model's segment means (blue), a dashed line at the changepoint and,
# optionally, a green line at the gold standard date.  Panel titles are
# "<year>, dev expl: <deviance explained by the changepoint model>".
#
# rain        : [day x year] matrix of daily precip for the cell
# years       : year labels (length ncol(rain))
# chpt        : changepoints (day of year), one per year
# dev_expl    : deviance explained by the changepoint model, per year
# dev_full    : deviance explained by the spline, per year (not shown in the
#               titles; kept for compatibility)
# doy_start, doy_end : analysis window
# gold        : optional gold standard day of year per year (NULL = none)
# y_max       : y-axis max (mm); NULL = data-driven per panel
# spline_k, gamma_floor : model settings (should match the run)
# spline      : optional precomputed spline_fits() for these years (so the
#               fits can be shared with plot_spline_overlay())
# highlight   : optional year to outline (e.g. the year selected in the
#               spline overlay plot)
# panels_per_row : NULL = square layout
# main        : outer title (NULL = no outer title)
#
# Returns (invisibly) a data frame giving each panel's figure region in
# normalized device coordinates (0-1, origin bottom left): year, x0, x1, y0,
# y1.  Used by panel_at() to find which year's panel was clicked.
# cex         : text/symbol size for the panels
plot_location_panels <- function(rain, years,
                                 chpt, dev_expl, dev_full = NULL,
                                 doy_start   = 1,
                                 doy_end     = 240,
                                 gold        = NULL,
                                 y_max       = 50,
                                 spline_k    = 20,
                                 gamma_floor = 1,
                                 spline      = NULL,
                                 highlight   = NULL,
                                 panels_per_row = NULL,
                                 main        = NULL,
                                 cex         = 0.6) {
  n_panels <- ncol(rain)
  n_window <- doy_end - doy_start + 1L
  doy_window <- seq(doy_start, doy_end)
  if (is.null(spline))
    spline <- spline_fits(rain, doy_start, doy_end, spline_k, gamma_floor)

  n_col <- if (is.null(panels_per_row)) ceiling(sqrt(n_panels)) else
    as.integer(panels_per_row)
  n_row <- ceiling(n_panels / n_col)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(n_row, n_col),
      mai   = c(0.3, 0.3, 0.25, 0.1) * max(1, cex / 0.6),
      mgp   = c(1.5, 0.5, 0),
      omi   = c(0, 0, if (is.null(main)) 0.05 else 0.45, 0),
      cex   = cex)

  # Figure region of the current panel, in normalized device coordinates
  regions <- data.frame(year = as.character(years), x0 = NA_real_,
                        x1 = NA_real_, y0 = NA_real_, y1 = NA_real_)
  record_region <- function(k) {
    regions[k, c("x0", "x1")] <<- grconvertX(c(0, 1), "nfc", "ndc")
    regions[k, c("y0", "y1")] <<- grconvertY(c(0, 1), "nfc", "ndc")
  }

  for (k in seq_len(n_panels)) {
    rain_full   <- as.numeric(rain[, k])
    rain_window <- extract_window(rain_full, doy_start, doy_end)

    if (all(is.na(rain_full))) {
      plot.new()
      record_region(k)
      title(main = years[k], cex.main = 1.0)
      mtext("no data", side = 3, line = -2, cex = 0.7)
      next
    }

    ylim <- if (!is.null(y_max) && is.finite(y_max)) c(0, y_max) else
      c(0, max(rain_full, na.rm = TRUE))

    plot(seq_along(rain_full), rain_full, type = "n",
         xlim = c(1, length(rain_full)), ylim = ylim,
         xlab = paste0("doy ", years[k]), ylab = "mm")
    record_region(k)
    rect(doy_start, 0, doy_end, ylim[2], col = "#3333FF0A", border = NA)
    lines(seq_along(rain_full), rain_full, col = "gray45", lwd = 0.6)

    chpt_abs <- chpt[k]
    chpt_rel <- chpt_abs - doy_start + 1L
    complete <- all(is.finite(rain_window))

    # Spline fit (red), on the floored window as in metric_dev_full()
    if (any(is.finite(spline[, k])))
      lines(doy_window, spline[, k], col = "red", lwd = 0.9)

    if (complete && !is.na(chpt_rel) && chpt_rel >= 1 && chpt_rel < n_window) {
      # Changepoint model: the two segment means of the floored window
      lines(doy_window, changepoint_fitted(rain_window, chpt_rel, gamma_floor),
            col = "blue", lwd = 0.9)
      abline(v = chpt_abs, col = "blue", lty = 2, lwd = 0.8)
      title(main = sprintf("%s, dev expl: %.2f", years[k], dev_expl[k]),
            cex.main = 1.0)
    } else {
      title(main = years[k], cex.main = 1.0)
    }

    # Gold standard date (green vertical line)
    if (!is.null(gold) && is.finite(gold[k]))
      abline(v = gold[k], col = "green3", lty = 1, lwd = 0.8)

    # Outline the highlighted year's panel
    if (!is.null(highlight) && identical(as.character(years[k]),
                                         as.character(highlight)))
      box(which = "figure", col = highlight_col, lwd = 3)
  }

  if (!is.null(main))
    title(main = main, outer = TRUE, line = 0.5, cex.main = 1.3)
  invisible(regions)
}

# The year of the panel containing a point, given in normalized device
# coordinates (0-1, origin bottom left), or NULL if none.
#
# regions : the data frame returned by plot_location_panels()
panel_at <- function(regions, x, y) {
  if (is.null(regions) || is.null(x) || is.null(y)) return(NULL)
  k <- which(x >= regions$x0 & x <= regions$x1 &
             y >= regions$y0 & y <= regions$y1)
  if (length(k)) regions$year[k[1]] else NULL
}


# Colour for the selected year (overlay line and panel outline)
highlight_col <- "#8B0000"

# All years' spline fits for one cell in a single panel, over the analysis
# window only, so patterns across years can be seen at once.
#
# spline    : spline_fits() matrix [window day x year]
# doy_start, doy_end : analysis window (x-axis range)
# years     : year labels (default: the column names of spline)
# highlight : optional year to emphasize: drawn thick and dark on top, the
#             other years faded, and labelled at its peak
# chpt, gold : optional changepoint and gold standard day of year, one per
#             year (aligned with years).  For the highlighted year they are
#             drawn as vertical lines (dashed blue, green), as in the
#             per-year panels, with a legend giving the days.
# gold_label : name of the gold standard for the legend
# col       : line colour when nothing is highlighted
# main      : title (NULL = automatic)
# cex       : text size
#
# Graphics settings are deliberately NOT restored on exit: Shiny reads the
# plot coordinates after drawing to translate clicks into day / value.
plot_spline_overlay <- function(spline, doy_start, doy_end,
                                years = colnames(spline), highlight = NULL,
                                chpt = NULL, gold = NULL,
                                gold_label = "gold standard",
                                col = "#E0000080", main = NULL, cex = 0.85) {
  par(mfrow = c(1, 1), mai = c(0.6, 0.7, 0.4, 0.2), mgp = c(1.8, 0.5, 0),
      cex = cex)

  doy <- seq(doy_start, doy_end)
  ok  <- apply(spline, 2, function(x) any(is.finite(x)))
  if (!any(ok)) {
    plot.new()
    mtext("no spline fits for the displayed years", side = 3, line = -2)
    return(invisible(NULL))
  }
  ylim <- spline_overlay_ylim(spline)
  plot(NA, xlim = c(doy_start, doy_end), ylim = ylim, xaxs = "i",
       xlab = "day of year", ylab = "spline fit (mm/day)", main = "")

  hi <- if (!is.null(highlight))
    match(as.character(highlight), as.character(years)) else NA
  if (!is.na(hi) && !ok[hi]) hi <- NA
  if (is.na(hi)) {
    for (k in which(ok)) lines(doy, spline[, k], col = col, lwd = 1.2)
  } else {
    for (k in setdiff(which(ok), hi))
      lines(doy, spline[, k], col = "#E0909066", lwd = 1)

    # The selected year's changepoint and gold standard date
    in_window <- function(v) !is.null(v) && is.finite(v) &&
      v >= doy_start && v <= doy_end
    cp <- if (!is.null(chpt)) chpt[hi] else NA
    gs <- if (!is.null(gold)) gold[hi] else NA
    leg <- character(0); leg_col <- character(0); leg_lty <- numeric(0)
    if (is.finite(cp)) {
      if (in_window(cp)) abline(v = cp, col = "blue", lty = 2, lwd = 1.5)
      leg <- c(leg, sprintf("changepoint: doy %d", round(cp)))
      leg_col <- c(leg_col, "blue"); leg_lty <- c(leg_lty, 2)
    }
    if (!is.null(gold)) {
      leg <- c(leg, if (!is.finite(gs)) sprintf("%s: none this year", gold_label)
               else sprintf("%s: doy %d%s", gold_label, round(gs),
                            if (in_window(gs)) "" else " (outside window)"))
      if (in_window(gs)) abline(v = gs, col = "green3", lwd = 1.5)
      leg_col <- c(leg_col, "green3"); leg_lty <- c(leg_lty, 1)
    }

    lines(doy, spline[, hi], col = "white", lwd = 5)   # halo
    lines(doy, spline[, hi], col = highlight_col, lwd = 2.8)
    i_max <- which.max(spline[, hi])
    text(doy[i_max], spline[i_max, hi], years[hi], pos = 3, font = 2,
         col = highlight_col, xpd = NA)
    if (length(leg))
      legend("topleft", legend = leg, col = leg_col, lty = leg_lty, lwd = 1.5,
             bg = "#FFFFFFCC", box.col = "gray70", cex = 0.85, inset = 0.01)
  }

  yrs <- years[ok]
  title(main = if (!is.null(main)) main else
    sprintf("Spline fits for %d years (%s-%s), doy %d-%d%s", sum(ok),
            min(yrs), max(yrs), doy_start, doy_end,
            if (is.na(hi)) "" else sprintf("  |  selected: %s", years[hi])),
    cex.main = 1.0)
  invisible(NULL)
}

spline_overlay_ylim <- function(spline) {
  m <- suppressWarnings(max(spline, na.rm = TRUE))
  c(0, if (is.finite(m)) m else 1)
}

# The year whose spline fit passes nearest the clicked point.
#
# spline   : spline_fits() matrix
# doy_start: first day of the window (row 1 of spline)
# x, y     : click position (day of year, value)
# tol      : largest vertical distance accepted, as a fraction of the
#            y-axis range; farther clicks return NULL
#
# Returns the year label (column name), or NULL.
nearest_spline_year <- function(spline, doy_start, x, y, tol = 0.05) {
  if (is.null(x) || is.null(y) || !ncol(spline)) return(NULL)
  i <- round(x) - doy_start + 1L
  if (i < 1 || i > nrow(spline)) return(NULL)
  d <- abs(spline[i, ] - y)
  if (all(is.na(d))) return(NULL)
  k <- which.min(d)
  if (d[k] > tol * diff(spline_overlay_ylim(spline))) return(NULL)
  colnames(spline)[k]
}
