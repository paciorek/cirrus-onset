# read_data.R
#
# read_region(): read daily precipitation for one region, a set of years and
# a day-of-year window into a lon x lat x day x year array.
#
# Replaces load_india_array() and load_chirps_array().  Only the window's days
# are read from each yearly file (via start/count), and dimensions are found
# by name, so the order in which they are stored in the file does not matter.
#
# Day slots are day-of-year 1..366 with Feb 29 kept (doy 60 in leap years);
# in non-leap years doy 366 does not exist.  Windows are limited to doy <= 365
# (see check_window()), so every year has data for every day of the window.

suppressPackageStartupMessages(library(ncdf4))

# Locate dimension indices of a NetCDF variable by name (case-insensitive).
nc_dim_index <- function(vobj, dim_names) {
  have <- toupper(vapply(vobj$dim, function(d) d$name, character(1)))
  idx <- vapply(dim_names, function(nm) {
    i <- which(have == toupper(nm))
    if (length(i) != 1)
      stop("Dimension '", nm, "' not found in variable '", vobj$name,
           "' (dimensions: ", paste(have, collapse = ", "), ")")
    i
  }, integer(1))
  setNames(idx, names(dim_names))
}

# Replace _FillValue / missing_value (and NaN) by NA.
nc_clean_missing <- function(x, nc, var_name) {
  for (att in c("_FillValue", "missing_value")) {
    a <- ncatt_get(nc, var_name, att)
    if (isTRUE(a$hasatt) && length(a$value) == 1 && is.finite(a$value))
      x[x == a$value] <- NA_real_
  }
  x[is.nan(x)] <- NA_real_
  x
}

region_files <- function(cfg, years) {
  files <- file.path(resolve_path(cfg$data_dir), sprintf(cfg$file_fmt, years))
  missing <- years[!file.exists(files)]
  if (length(missing))
    stop("Missing data files for years: ", paste(missing, collapse = ", "),
         " (looked in ", resolve_path(cfg$data_dir), ")", call. = FALSE)
  files
}

# All years that have a data file in the region's data folder, found by
# matching file names against cfg$file_fmt (e.g. "data_%d.nc" -> data_1901.nc).
# If doy_end is given, years whose file has fewer than doy_end days (e.g. a
# partial current year) are dropped with a message.
data_years <- function(cfg, doy_end = NULL, verbose = TRUE) {
  data_dir <- resolve_path(cfg$data_dir)
  if (!dir.exists(data_dir)) stop("Data folder not found: ", data_dir, call. = FALSE)
  esc   <- gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", cfg$file_fmt, perl = TRUE)
  regex <- paste0("^", sub("%d", "([0-9]{4})", esc, fixed = TRUE), "$")
  files <- list.files(data_dir, pattern = regex)
  if (!length(files))
    stop("No files matching '", cfg$file_fmt, "' in ", data_dir, call. = FALSE)
  years <- sort(as.integer(sub(regex, "\\1", files)))

  if (!is.null(doy_end)) {
    n_days <- vapply(years, function(yr) {
      nc <- nc_open(file.path(data_dir, sprintf(cfg$file_fmt, yr)))
      on.exit(nc_close(nc))
      vobj <- nc$var[[cfg$var_name]]
      vobj$varsize[nc_dim_index(vobj, cfg$dim_names)["time"]]
    }, numeric(1))
    short <- years[n_days < doy_end]
    if (length(short) && verbose)
      msg("skipping years whose files end before doy", doy_end, ":",
          paste(short, collapse = ", "))
    years <- years[n_days >= doy_end]
    if (!length(years)) stop("No data files cover the window.", call. = FALSE)
  }
  years
}

# cfg       : region config (get_region())
# years     : integer vector of years
# doy_start, doy_end : window (inclusive), 1 <= doy_start < doy_end <= 365
#
# Returns a list:
#   rain   lon x lat x day x year array (dimnames lon, lat, doy, year)
#   lon, lat (increasing), years, doy (window day numbers), files,
#   n_days_by_year (days in each file)
read_region <- function(cfg, years, doy_start, doy_end, verbose = TRUE) {
  check_window(doy_start, doy_end)
  files  <- region_files(cfg, years)
  n_win  <- doy_end - doy_start + 1L
  dn     <- cfg$dim_names

  nc0  <- nc_open(files[1])
  v0   <- nc0$var[[cfg$var_name]]
  if (is.null(v0)) stop("Variable '", cfg$var_name, "' not in ", files[1])
  id0  <- nc_dim_index(v0, dn)
  lon  <- as.numeric(v0$dim[[id0["lon"]]]$vals)
  lat  <- as.numeric(v0$dim[[id0["lat"]]]$vals)
  nc_close(nc0)

  # Arrays are stored with increasing lon and lat (needed by image()).
  o_lon <- order(lon); o_lat <- order(lat)

  rain <- array(NA_real_, dim = c(length(lon), length(lat), n_win, length(years)),
                dimnames = list(lon  = sprintf("%.4f", lon[o_lon]),
                                lat  = sprintf("%.4f", lat[o_lat]),
                                doy  = as.character(doy_start:doy_end),
                                year = as.character(years)))
  n_days_by_year <- setNames(integer(length(years)), years)

  for (k in seq_along(files)) {
    nc   <- nc_open(files[k])
    vobj <- nc$var[[cfg$var_name]]
    id   <- nc_dim_index(vobj, dn)
    lon_k <- as.numeric(vobj$dim[[id["lon"]]]$vals)
    lat_k <- as.numeric(vobj$dim[[id["lat"]]]$vals)
    if (length(lon_k) != length(lon) || length(lat_k) != length(lat) ||
        max(abs(lon_k - lon)) > 1e-6 || max(abs(lat_k - lat)) > 1e-6) {
      nc_close(nc)
      stop("Grid in ", basename(files[k]), " differs from ", basename(files[1]))
    }
    n_days <- vobj$varsize[id["time"]]
    n_days_by_year[k] <- n_days
    if (n_days < doy_end) {
      nc_close(nc)
      stop(basename(files[k]), " has only ", n_days, " days; window ends at doy ",
           doy_end, call. = FALSE)
    }

    start <- rep(1L, length(vobj$dim)); count <- vobj$varsize
    start[id["time"]] <- doy_start
    count[id["time"]] <- n_win
    x <- ncvar_get(nc, cfg$var_name, start = start, count = count,
                   collapse_degen = FALSE)
    x <- nc_clean_missing(x, nc, cfg$var_name)
    nc_close(nc)

    # Keep only the lon/lat/time dimensions (drop any degenerate extras),
    # then permute to lon x lat x time.
    keep <- c(id["lon"], id["lat"], id["time"])
    if (length(dim(x)) > 3) {
      extra <- setdiff(seq_along(dim(x)), keep)
      if (any(dim(x)[extra] != 1))
        stop("Unexpected non-degenerate dimensions in ", basename(files[k]))
    }
    x <- array(aperm(x, c(keep, setdiff(seq_along(dim(x)), keep))),
               dim = dim(x)[keep])
    rain[, , , k] <- x[o_lon, o_lat, , drop = FALSE]

    if (verbose) msg("read", basename(files[k]))
  }

  list(rain = rain, lon = lon[o_lon], lat = lat[o_lat], years = as.integer(years),
       doy = doy_start:doy_end, files = files, n_days_by_year = n_days_by_year)
}
