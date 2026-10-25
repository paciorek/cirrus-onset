# plot_utils.R
#
# Per-year time series panels for one grid cell.  Map panels use the shared
# plot_map_panel() and summarize_years() in ../R/maps.R.

suppressPackageStartupMessages(library(mgcv))

# Plot one panel per year for a single grid cell: the full annual series in
# gray, the analysis window shaded, the spline fit (red), the changepoint
# model's segment means (blue), a dashed line at the changepoint and,
# optionally, a green line at the gold standard date.  Panel titles give
# "changepoint dev.expl : spline dev.expl" (the latter only when dev_full has
# been computed).
#
# rain        : [day x year] matrix of daily precip for the cell
# years       : year labels (length ncol(rain))
# chpt        : changepoints (day of year), one per year
# dev_expl    : deviance explained by the changepoint model, per year
# dev_full    : deviance explained by the spline, per year (NULL = not shown)
# doy_start, doy_end : analysis window
# gold        : optional gold standard day of year per year (NULL = none)
# y_max       : y-axis max (mm); NULL = data-driven per panel
# spline_k, gamma_floor : model settings (should match the run)
# panels_per_row : NULL = square layout
# main        : outer title (NULL = no outer title)
# cex         : text/symbol size for the panels
plot_location_panels <- function(rain, years,
                                 chpt, dev_expl, dev_full = NULL,
                                 doy_start   = 1,
                                 doy_end     = 240,
                                 gold        = NULL,
                                 y_max       = 50,
                                 spline_k    = 20,
                                 gamma_floor = 1,
                                 panels_per_row = NULL,
                                 main        = NULL,
                                 cex         = 0.6) {
  n_panels <- ncol(rain)
  n_window <- doy_end - doy_start + 1L
  doy_window <- seq(doy_start, doy_end)

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

  for (k in seq_len(n_panels)) {
    rain_full   <- as.numeric(rain[, k])
    rain_window <- extract_window(rain_full, doy_start, doy_end)

    if (all(is.na(rain_full))) {
      plot.new()
      title(main = years[k], cex.main = 1.0)
      mtext("no data", side = 3, line = -2, cex = 0.7)
      next
    }

    ylim <- if (!is.null(y_max) && is.finite(y_max)) c(0, y_max) else
      c(0, max(rain_full, na.rm = TRUE))

    plot(seq_along(rain_full), rain_full, type = "n",
         xlim = c(1, length(rain_full)), ylim = ylim,
         xlab = paste0("doy ", years[k]), ylab = "mm")
    rect(doy_start, 0, doy_end, ylim[2], col = "#3333FF0A", border = NA)
    lines(seq_along(rain_full), rain_full, col = "gray45", lwd = 0.6)

    chpt_abs <- chpt[k]
    chpt_rel <- chpt_abs - doy_start + 1L
    complete <- all(is.finite(rain_window))

    # Spline fit (red), on the floored window as in metric_dev_full()
    if (complete) {
      tmp <- pmax(rain_window, gamma_floor)
      t   <- seq_len(n_window)
      fullmod <- try(suppressWarnings(
        gam(tmp ~ s(t, k = spline_k), family = Gamma(link = "log"))),
        silent = TRUE)
      if (!inherits(fullmod, "try-error"))
        lines(doy_window, predict(fullmod, type = "response"),
              col = "red", lwd = 0.9)
    }

    if (complete && !is.na(chpt_rel) && chpt_rel >= 1 && chpt_rel < n_window) {
      # Changepoint model: the two segment means of the floored window
      lines(doy_window, changepoint_fitted(rain_window, chpt_rel, gamma_floor),
            col = "blue", lwd = 0.9)
      abline(v = chpt_abs, col = "blue", lty = 2, lwd = 0.8)
      ttl <- paste0(years[k], "   ", round(dev_expl[k], 2))
      if (!is.null(dev_full)) ttl <- paste0(ttl, " : ", round(dev_full[k], 2))
      title(main = ttl, cex.main = 1.0)
    } else {
      title(main = years[k], cex.main = 1.0)
    }

    # Gold standard date (green vertical line)
    if (!is.null(gold) && is.finite(gold[k]))
      abline(v = gold[k], col = "green3", lty = 1, lwd = 0.8)
  }

  if (!is.null(main))
    title(main = main, outer = TRUE, line = 0.5, cex.main = 1.3)
  invisible(NULL)
}
