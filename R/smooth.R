# smooth.R
#
# Spatial smoothing of deviance explained, used to pick out regions with a
# consistently sharp onset (or retreat) by thresholding the smoothed map.
# Smoothing makes region boundaries less ragged and fills small holes.
#
# Steps:
#   1. For each grid cell, average the yearly dev_expl values over the years
#      used, giving one value per location (a lon x lat field):
#        z <- mean_dev_expl(run$r$dev_expl, run$years, years = 2000:2024)
#   2. Smooth that field with a thin-plate spline in (lon, lat), with
#      user-chosen basis size k and GAM family:
#        sm <- smooth_dev_expl(z, run$lon, run$lat, k = 500, family = "betar")
#   3. Threshold the smoothed (or raw) field to get a region mask, and save it:
#        m <- make_region_mask(run, threshold = 0.3, source = "smooth",
#                              years = 2000:2024, smooth = sm)
#        save_region_mask(m)       # output/<region>/<tag>/masks/...
#        head(mask_cells(m))       # lon, lat, value, keep
#
# Fitting time grows quickly with k and is much larger for betar and
# quasibinomial than for gaussian (India, ~5000 cells: gaussian k = 500 about
# 20 s; betar k = 500 several minutes; betar k = 3000 hours).  Fits can be
# saved and reloaded with fit_or_load_smooth(), which the Shiny app also uses,
# and computed offline with scripts/smooth_dev_expl.R.

suppressPackageStartupMessages(library(mgcv))

# GAM families offered for smoothing (name -> family object).  dev_expl lies
# in [0, 1]: betar (beta regression, logit link) and quasibinomial (logit
# link) respect that; gaussian is fastest but can give values outside [0, 1].
smooth_families <- list(
  gaussian      = function() gaussian(),
  betar         = function() betar(link = "logit"),
  quasibinomial = function() quasibinomial()
)

smooth_family <- function(family) {
  if (inherits(family, "family")) return(family)
  if (!family %in% names(smooth_families))
    stop("Unknown smoothing family '", family, "'. Choose from: ",
         paste(names(smooth_families), collapse = ", "), call. = FALSE)
  smooth_families[[family]]()
}

# Per-cell mean over years: averages the yearly dev_expl values at each grid
# cell over the given years (all years if NULL), giving one value per
# location.  Returns a lon x lat matrix, NA where a cell has no data.
mean_dev_expl <- function(dev_expl, all_years, years = NULL) {
  idx <- if (is.null(years)) seq_along(all_years) else match(years, all_years)
  if (anyNA(idx)) stop("Requested years not in the run: ",
                       paste(years[is.na(idx)], collapse = ", "))
  z <- suppressWarnings(apply(dev_expl[, , idx, drop = FALSE], c(1, 2), mean,
                              na.rm = TRUE))
  z[!is.finite(z)] <- NA
  z
}

# Smooth a lon x lat matrix z with a thin-plate spline s(lon, lat, k = k).
#
# z        : lon x lat matrix (e.g. from mean_dev_expl()); NA cells are
#            not used and stay NA in the result
# lon, lat : coordinate vectors
# k        : basis dimension of the spatial smooth
# family   : "gaussian", "betar", "quasibinomial", or a family object
# engine   : "gam" (default) or "bam" (bam(..., discrete = TRUE); faster for
#            gaussian with large k, not for betar)
# method   : smoothing-parameter method passed to gam()/bam(); NULL uses
#            mgcv's default (as in the original scratch/smooth.R).  "REML"
#            gives near-identical maps here but was ~3x slower for gaussian.
#
# Returns a list:
#   fitted   lon x lat matrix of smoothed values (response scale)
#   raw      z
#   model    the fitted gam/bam object
#   k, family (name), engine, edf (total), n (cells used), seconds
smooth_dev_expl <- function(z, lon, lat, k = 500, family = "gaussian",
                            engine = c("gam", "bam"), method = NULL) {
  engine <- match.arg(engine)
  fam <- smooth_family(family)
  fam_name <- if (is.character(family)) family else fam$family

  d <- expand.grid(lon = lon, lat = lat)
  d$y <- c(z)
  ok <- is.finite(d$y)
  if (sum(ok) <= k)
    stop("k = ", k, " must be smaller than the number of cells with data (",
         sum(ok), ").", call. = FALSE)
  dd <- d[ok, ]
  if (fam_name == "betar" && any(dd$y <= 0 | dd$y >= 1))
    message("betar: values at 0 or 1 are moved just inside (0, 1) by mgcv.")

  t0 <- Sys.time()
  args <- list(formula = y ~ s(lon, lat, k = k), data = dd, family = fam)
  if (!is.null(method)) args$method <- method
  model <- if (engine == "gam") do.call(gam, args) else
    do.call(bam, c(args, discrete = TRUE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  fitted <- matrix(NA_real_, length(lon), length(lat), dimnames = dimnames(z))
  fitted[ok] <- predict(model, newdata = dd, type = "response")

  list(fitted = fitted, raw = z, model = model, k = k, family = fam_name,
       engine = engine, edf = sum(model$edf), n = sum(ok), seconds = secs)
}

# Logical lon x lat matrix: TRUE where values >= threshold, FALSE where below
# or missing.
threshold_mask <- function(values, threshold) {
  m <- !is.na(values) & values >= threshold
  dim(m) <- dim(values)
  dimnames(m) <- dimnames(values)
  m
}

# -- Saving and reloading fits -----------------------------------------------

# Fits are saved as output/<region>/<tag>/smooth/dev_expl_<family>_k<k>_<years>.rds
# holding the fitted and raw matrices and the settings (not the model, which
# can be hundreds of MB for large k).
smooth_years_label <- function(years) {
  years <- sort(years)
  lab <- sprintf("%d-%d", min(years), max(years))
  if (length(years) != max(years) - min(years) + 1)
    lab <- paste0(lab, "_n", length(years))
  lab
}

smooth_file <- function(region, tag, family, k, years, engine = "gam") {
  file.path(run_dir(region, tag, create = FALSE), "smooth",
            sprintf("dev_expl_%s_k%d_%s%s.rds", family, as.integer(k),
                    smooth_years_label(years),
                    if (engine == "gam") "" else paste0("_", engine)))
}

# Load a saved fit if there is one; otherwise fit (unless fit = FALSE, in which
# case NULL is returned) and save it.  A loaded fit has model = NULL.
#
# run    : list from load_run()
# years  : years to average over (NULL = all run years)
fit_or_load_smooth <- function(run, k = 500, family = "gaussian", years = NULL,
                               engine = "gam", fit = TRUE, save = TRUE,
                               refit = FALSE) {
  if (is.null(years)) years <- run$years
  f <- smooth_file(run$meta$region, run$meta$tag, family, k, years, engine)
  if (!refit && file.exists(f)) {
    sm <- readRDS(f)
    sm$model <- NULL
    sm$source <- "loaded"
    return(sm)
  }
  if (!fit) return(NULL)
  z  <- mean_dev_expl(run$r$dev_expl, run$years, years)
  sm <- smooth_dev_expl(z, run$lon, run$lat, k = k, family = family,
                        engine = engine)
  sm$years <- sort(years)
  if (save) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    saveRDS(sm[setdiff(names(sm), "model")], f)
  }
  sm$source <- "fitted"
  sm
}


# -- Region masks -------------------------------------------------------------

# A region mask from thresholding the per-cell mean dev_expl, raw or smoothed.
#
# run       : list from load_run()
# threshold : cells with value >= threshold are kept
# source    : "smooth" or "raw"
# years     : years averaged over (NULL = all run years)
# smooth    : an existing fit (from smooth_dev_expl()/fit_or_load_smooth());
#             if NULL and source = "smooth", one is loaded or fitted with k,
#             family and engine
#
# Returns a list of class "region_mask": keep (lon x lat logical), values
# (the thresholded field), lon, lat, and the settings used.
make_region_mask <- function(run, threshold, source = c("smooth", "raw"),
                             years = NULL, smooth = NULL, k = 500,
                             family = "gaussian", engine = "gam") {
  source <- match.arg(source)
  if (is.null(years)) years <- run$years
  years <- sort(years)
  if (source == "raw") {
    values <- mean_dev_expl(run$r$dev_expl, run$years, years)
    settings <- list()
  } else {
    if (is.null(smooth))
      smooth <- fit_or_load_smooth(run, k = k, family = family, years = years,
                                   engine = engine)
    if (!is.null(smooth$years) && !identical(as.integer(smooth$years),
                                             as.integer(years)))
      stop("The smooth was fit for different years than requested.")
    values <- smooth$fitted
    settings <- list(k = smooth$k, family = smooth$family, engine = smooth$engine,
                     edf = smooth$edf)
  }
  structure(c(
    list(keep = threshold_mask(values, threshold), values = values,
         lon = run$lon, lat = run$lat,
         region = run$meta$region, tag = run$meta$tag,
         direction = run$meta$direction,
         doy_start = run$meta$doy_start, doy_end = run$meta$doy_end,
         years = years, source = source, threshold = threshold),
    settings,
    list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
    class = "region_mask")
}

# Default file name for a mask:
# output/<region>/<tag>/masks/mask_<raw|smooth_<family>_k<k>>_thr<thr>_<years>.rds
region_mask_file <- function(m) {
  src <- if (m$source == "raw") "raw" else
    sprintf("smooth_%s_k%d", m$family, as.integer(m$k))
  file.path(run_dir(m$region, m$tag, create = FALSE), "masks",
            sprintf("mask_%s_thr%s_%s.rds", src, format(m$threshold, nsmall = 2),
                    smooth_years_label(m$years)))
}

# Save a mask (as .rds) and return the file path.
save_region_mask <- function(m, file = region_mask_file(m)) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(m, file)
  file
}

load_region_mask <- function(file) readRDS(file)

# Cells with data as a data frame: lon, lat, value (the thresholded field)
# and keep.
mask_cells <- function(m) {
  d <- expand.grid(lon = m$lon, lat = m$lat)
  d$value <- c(m$values)
  d$keep  <- c(m$keep)
  d[!is.na(d$value), , drop = FALSE]
}

print.region_mask <- function(x, ...) {
  cat(sprintf("Region mask: %s, %s mean dev_expl >= %g (years %s)%s\n",
              x$tag, if (x$source == "raw") "raw" else "smoothed", x$threshold,
              smooth_years_label(x$years),
              if (x$source == "smooth")
                sprintf(", %s k = %d", x$family, as.integer(x$k)) else ""))
  cat(sprintf("%d of %d cells with data kept\n", sum(x$keep),
              sum(!is.na(x$values))))
  invisible(x)
}
