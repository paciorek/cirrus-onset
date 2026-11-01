# combine_gold_standard.R
#
# Combine one-file-per-year gold standard onset files into a single NetCDF
# with dimensions lon x lat x year, the layout read by read_gold_standard()
# (R/gold_standard.R), like the India O2 files.
#
# Defaults are set for the Ethiopia ROMP onset files
# (data/ethiopia_onset/climatology_onset_doy_by_year/onset_doy_YYYY.nc, one
# lon x lat field of calendar day of year per file, NaN = no onset detected),
# so for those just run
#
#   Rscript scripts/combine_gold_standard.R
#
# Options (paths relative to the base directory unless absolute):
#   --in_dir=<dir>     folder with the yearly files
#   --pattern=<regex>  file-name regex with the year as its first group
#   --var=<name>       variable to combine
#   --out=<file>       output NetCDF (default: <var>_<first>-<last year>.nc in
#                      the parent of in_dir, e.g.
#                      data/ethiopia_onset/onset_doy_2000-2022.nc)
#
# All input files must share the same grid.  Years are taken from the file
# names and need not be consecutive; years with no file are simply absent
# from the output (the app then draws no gold standard line for them).
# Missing values are written as -9999 with _FillValue/missing_value set.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

args <- parse_cli()
in_dir_arg <- if (!is.null(args$in_dir)) args$in_dir else
  "data/ethiopia_onset/climatology_onset_doy_by_year"
in_dir  <- resolve_path(in_dir_arg)
pattern <- if (!is.null(args$pattern)) args$pattern else "^onset_doy_([0-9]{4})\\.nc$"
var     <- if (!is.null(args$var)) args$var else "onset_doy"

files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
if (!length(files)) stop("No files matching ", pattern, " in ", in_dir, call. = FALSE)
years <- as.integer(sub(pattern, "\\1", basename(files)))
o <- order(years); files <- files[o]; years <- years[o]
if (anyDuplicated(years)) stop("Duplicate years in file names.", call. = FALSE)

out <- if (!is.null(args$out)) resolve_path(args$out) else
  file.path(dirname(in_dir), sprintf("%s_%d-%d.nc", var, min(years), max(years)))

# Read the grid and variable attributes from the first file
nc0 <- nc_open(files[1])
v0  <- nc0$var[[var]]
if (is.null(v0)) stop("Variable '", var, "' not in ", files[1], call. = FALSE)
id  <- nc_dim_index(v0, c(lon = "lon", lat = "lat"))
lon <- as.numeric(v0$dim[[id["lon"]]]$vals)
lat <- as.numeric(v0$dim[[id["lat"]]]$vals)
atts <- ncatt_get(nc0, var)
nc_close(nc0)

arr <- array(NA_real_, c(length(lon), length(lat), length(years)))
for (k in seq_along(files)) {
  nc <- nc_open(files[k])
  v  <- nc$var[[var]]
  idk <- nc_dim_index(v, c(lon = "lon", lat = "lat"))
  lon_k <- as.numeric(v$dim[[idk["lon"]]]$vals)
  lat_k <- as.numeric(v$dim[[idk["lat"]]]$vals)
  if (length(lon_k) != length(lon) || length(lat_k) != length(lat) ||
      max(abs(lon_k - lon)) > 1e-6 || max(abs(lat_k - lat)) > 1e-6)
    stop("Grid in ", basename(files[k]), " differs from ", basename(files[1]),
         call. = FALSE)
  x <- ncvar_get(nc, var, collapse_degen = FALSE)
  x <- nc_clean_missing(x, nc, var)
  nc_close(nc)
  arr[, , k] <- aperm(x, c(idk["lon"], idk["lat"]))
  cat(sprintf("%d: %d cells with a value, range %s\n", years[k], sum(!is.na(x)),
              if (all(is.na(x))) "-" else paste(range(x, na.rm = TRUE), collapse = "-")))
}

# Write lon x lat x year
fill  <- -9999
d_lon <- ncdim_def("lon", "degrees_east", lon)
d_lat <- ncdim_def("lat", "degrees_north", lat)
d_yr  <- ncdim_def("year", "year", years)
v_out <- ncvar_def(var, units = if (isTRUE(atts$units != "")) atts$units else "",
                   dim = list(d_lon, d_lat, d_yr), missval = fill, prec = "double")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
if (file.exists(out)) file.remove(out)
nco <- nc_create(out, v_out, force_v4 = TRUE)
ncvar_put(nco, v_out, arr)
ncatt_put(nco, var, "missing_value", fill)
# Carry over descriptive attributes (method settings etc.) from the inputs
skip <- c("_FillValue", "missing_value", "units", "missing")
for (nm in setdiff(names(atts), skip))
  ncatt_put(nco, var, nm, atts[[nm]])
ncatt_put(nco, var, "missing", sprintf("%g = no onset detected", fill))
ncatt_put(nco, 0, "source_files", paste(basename(files), collapse = " "))
ncatt_put(nco, 0, "source_dir", in_dir_arg)
ncatt_put(nco, 0, "history", sprintf("%s: combined by scripts/combine_gold_standard.R",
                                     format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
nc_close(nco)

cat(sprintf("Wrote %s: %d x %d x %d (years %d-%d)\n", out, length(lon),
            length(lat), length(years), min(years), max(years)))
