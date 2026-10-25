# gold_standard.R
#
# Optional "gold standard" onset or retreat dates (e.g. O2 onset for India).
# Each dataset is a single NetCDF file holding all years on the same grid as
# the precipitation data; see config/regions.R for the fields.
#
# Values are converted to day of year per year: doy = doy(origin date in that
# year) + value - 1.  With origin "04-01", value 1 is April 1, i.e. doy 91 in
# normal years and doy 92 in leap years.

suppressPackageStartupMessages(library(ncdf4))

origin_doy <- function(origin, years) {
  as.integer(format(as.Date(sprintf("%d-%s", years, origin)), "%j"))
}

# entry      : one gold standard entry from the region config
# lon, lat   : the run's grid (increasing)
# years      : the run's years
# tol        : allowed coordinate mismatch (degrees)
#
# Returns a lon x lat x year array of day of year (NA where missing or for
# years not in the file), with dimnames matching the run's arrays.
read_gold_standard <- function(entry, lon, lat, years, tol = 1e-3) {
  f <- resolve_path(entry$file)
  if (!file.exists(f)) stop("Gold standard file not found: ", f, call. = FALSE)

  nc <- nc_open(f)
  on.exit(nc_close(nc))
  vobj <- nc$var[[entry$var_name]]
  if (is.null(vobj)) stop("Variable '", entry$var_name, "' not in ", f)
  id <- nc_dim_index(vobj, entry$dim_names)

  g_lon  <- as.numeric(vobj$dim[[id["lon"]]]$vals)
  g_lat  <- as.numeric(vobj$dim[[id["lat"]]]$vals)
  g_year <- as.integer(round(as.numeric(vobj$dim[[id["year"]]]$vals)))

  x <- ncvar_get(nc, entry$var_name, collapse_degen = FALSE)
  x <- nc_clean_missing(x, nc, entry$var_name)
  x <- aperm(x, c(id["lon"], id["lat"], id["year"]))

  o_lon <- order(g_lon); o_lat <- order(g_lat)
  g_lon <- g_lon[o_lon]; g_lat <- g_lat[o_lat]
  x <- x[o_lon, o_lat, , drop = FALSE]

  if (length(g_lon) != length(lon) || length(g_lat) != length(lat) ||
      max(abs(g_lon - lon)) > tol || max(abs(g_lat - lat)) > tol)
    stop("Gold standard grid in ", basename(f),
         " does not match the precipitation grid.", call. = FALSE)

  out <- array(NA_real_, c(length(lon), length(lat), length(years)),
               dimnames = list(lon  = sprintf("%.4f", lon),
                               lat  = sprintf("%.4f", lat),
                               year = as.character(years)))
  shared <- intersect(years, g_year)
  if (!length(shared))
    warning("No years in common between ", basename(f), " and the run.")
  for (yr in shared) {
    out[, , as.character(yr)] <-
      x[, , match(yr, g_year)] + origin_doy(entry$origin, yr) - 1
  }
  out
}

# The gold standard entries that apply to a run's direction (possibly empty).
gold_standards_for <- function(cfg, direction) {
  gs <- cfg$gold_standard
  if (is.null(gs)) return(list())
  Filter(function(e) identical(e$direction, direction), gs)
}

# Read every gold standard dataset that applies to the run.
# Returns a named list of lon x lat x year arrays (empty list if none).
read_gold_standards <- function(cfg, direction, lon, lat, years) {
  gs <- gold_standards_for(cfg, direction)
  lapply(gs, read_gold_standard, lon = lon, lat = lat, years = years)
}
