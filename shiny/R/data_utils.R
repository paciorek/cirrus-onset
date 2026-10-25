# data_utils.R
#
# App-specific data utilities.  Loading a run (changepoints, stats and gold
# standard) uses the shared load_run() and list_runs() in ../R/run_io.R.
#
#   extract_window()     -- window of a full-year series
#   nearest_cell()       -- nearest (optionally non-missing) grid cell
#   read_cell_series()   -- one grid cell's daily series for given years,
#                           read straight from the yearly NetCDF files

suppressPackageStartupMessages(library(ncdf4))

# Extract the analysis window from a full-year daily precipitation vector.
# Days outside the year's actual length are returned as NA.
extract_window <- function(rain_series, doy_start, doy_end) {
  idx  <- seq(doy_start, doy_end)
  out  <- rep(NA_real_, length(idx))
  keep <- idx <= length(rain_series)
  if (any(keep)) out[keep] <- rain_series[idx[keep]]
  out
}

# Nearest grid cell to (lon0, lat0).  If `valid` (a logical [lon x lat]
# matrix) is given, only cells where valid is TRUE are considered, so a
# click just offshore snaps to the nearest land cell.
#
# Returns list(i_lon, i_lat, lon, lat).
nearest_cell <- function(lon0, lat0, lon, lat, valid = NULL) {
  if (is.null(valid)) {
    i_lon <- which.min(abs(lon - lon0))
    i_lat <- which.min(abs(lat - lat0))
  } else {
    idx <- which(valid, arr.ind = TRUE)
    if (nrow(idx) == 0) stop("No valid grid cells.")
    d2  <- (lon[idx[, 1]] - lon0)^2 + (lat[idx[, 2]] - lat0)^2
    k   <- which.min(d2)
    i_lon <- idx[k, 1]
    i_lat <- idx[k, 2]
  }
  list(i_lon = i_lon, i_lat = i_lat, lon = lon[i_lon], lat = lat[i_lat])
}

# Read the daily precipitation series for a single grid cell from the yearly
# NetCDF files.  Only the requested cell is read from each file (via
# start/count), so this takes well under a second for 25 years.
#
# cfg      : region config (data_dir, file_fmt, var_name, dim_names)
# years    : integer vector of years
# lon0, lat0 : grid-cell coordinates (matched to the nearest grid point)
# n_days   : number of rows in the returned matrix
#
# Returns a [n_days x n_years] matrix (columns named by year), NA where a
# file or day is missing.
read_cell_series <- function(cfg, years, lon0, lat0, n_days = 366) {
  out <- matrix(NA_real_, n_days, length(years),
                dimnames = list(day = NULL, year = as.character(years)))
  data_dir <- resolve_path(cfg$data_dir)

  for (k in seq_along(years)) {
    f <- file.path(data_dir, sprintf(cfg$file_fmt, years[k]))
    if (!file.exists(f)) next

    nc   <- nc_open(f)
    vobj <- nc$var[[cfg$var_name]]
    id   <- nc_dim_index(vobj, cfg$dim_names)

    i_lon <- which.min(abs(vobj$dim[[id["lon"]]]$vals - lon0))
    i_lat <- which.min(abs(vobj$dim[[id["lat"]]]$vals - lat0))

    start <- rep(1L, length(vobj$dim)); count <- rep(1L, length(vobj$dim))
    start[id["lon"]] <- i_lon
    start[id["lat"]] <- i_lat
    count[id["time"]] <- vobj$varsize[id["time"]]

    x <- as.numeric(ncvar_get(nc, cfg$var_name, start = start, count = count))
    x <- nc_clean_missing(x, nc, cfg$var_name)
    nc_close(nc)

    n_copy <- min(length(x), n_days)
    out[seq_len(n_copy), k] <- x[seq_len(n_copy)]
  }
  out
}

# Years with a data file present for this region.
available_data_years <- function(cfg, years) {
  data_dir <- resolve_path(cfg$data_dir)
  years[file.exists(file.path(data_dir, sprintf(cfg$file_fmt, years)))]
}
