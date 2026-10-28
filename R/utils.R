# utils.R
#
# Run naming, paths, argument parsing and metadata checks shared by all steps.

base_dir <- function() {
  d <- getOption("onset.base_dir")
  if (is.null(d)) stop("onset_setup() has not been called.")
  d
}

# Resolve a path from the config: absolute paths are returned unchanged,
# relative ones are taken relative to the root directory.
resolve_path <- function(p) {
  if (is.null(p)) return(NULL)
  if (grepl("^(/|~|[A-Za-z]:)", p)) return(path.expand(p))
  file.path(base_dir(), p)
}

get_region <- function(region) {
  if (!region %in% names(regions))
    stop("Unknown region '", region, "'. Known regions: ",
         paste(names(regions), collapse = ", "))
  cfg <- regions[[region]]
  cfg$key <- region
  cfg
}

# -- Command-line arguments ---------------------------------------------------

# Parse --key=value arguments into a named list of strings.
parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (a in args) {
    if (!grepl("^--[A-Za-z_]+=", a))
      stop("Arguments must look like --key=value; got: ", a, call. = FALSE)
    key <- sub("^--([A-Za-z_]+)=.*$", "\\1", a)
    out[[key]] <- sub("^--[A-Za-z_]+=", "", a)
  }
  out
}

require_args <- function(args, keys, usage) {
  missing <- setdiff(keys, names(args))
  if (length(missing))
    stop("Missing argument(s): ", paste0("--", missing, collapse = ", "),
         "\n", usage, call. = FALSE)
}

# "1951:2025", "1990,1995,2000:2010" -> integer vector
parse_years <- function(x) {
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  yrs <- unlist(lapply(parts, function(p) {
    if (grepl(":", p)) {
      ab <- as.integer(strsplit(p, ":", fixed = TRUE)[[1]])
      seq(ab[1], ab[2])
    } else as.integer(p)
  }))
  if (anyNA(yrs) || !length(yrs)) stop("Could not parse years: ", x)
  sort(unique(yrs))
}

# "210:365" -> c(210, 365)
parse_doy <- function(x) {
  ab <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
  if (length(ab) != 2 || anyNA(ab)) stop("--doy must look like 210:365")
  check_window(ab[1], ab[2])
  ab
}

check_window <- function(doy_start, doy_end) {
  if (doy_start < 1 || doy_end > 365 || doy_end <= doy_start)
    stop("Window must satisfy 1 <= doy_start < doy_end <= 365 ",
         "(got ", doy_start, "-", doy_end, ").", call. = FALSE)
  invisible(TRUE)
}

check_direction <- function(direction) {
  if (!direction %in% c("onset", "retreat"))
    stop("direction must be 'onset' or 'retreat'", call. = FALSE)
  direction
}

# -- Run naming ---------------------------------------------------------------

run_tag <- function(region, years, doy_start, doy_end, direction) {
  sprintf("%s_%d-%d_doy%d-%d_%s", region, min(years), max(years),
          doy_start, doy_end, direction)
}

# Inverse of run_tag(); returns NULL if `tag` does not match the pattern.
parse_run_tag <- function(tag) {
  m <- regmatches(tag, regexec(
    "^(.+)_([0-9]{4})-([0-9]{4})_doy([0-9]+)-([0-9]+)_(onset|retreat)$", tag))[[1]]
  if (!length(m)) return(NULL)
  list(region = m[2], year_start = as.integer(m[3]), year_end = as.integer(m[4]),
       doy_start = as.integer(m[5]), doy_end = as.integer(m[6]),
       direction = m[7])
}

output_root <- function() file.path(base_dir(), "output")

run_dir <- function(region, tag, create = TRUE) {
  d <- file.path(output_root(), region, tag)
  if (create) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Common setup for the step scripts: parse the run inputs.  --years is
# optional; without it, all years with a data file covering the window are
# used (see data_years()).
run_spec_from_args <- function(args, usage) {
  require_args(args, c("region", "doy", "direction"), usage)
  doy <- parse_doy(args$doy)
  years <- if (!is.null(args$years)) parse_years(args$years) else NULL
  spec <- make_spec(args$region, years, doy[1], doy[2], args$direction)
  if (is.null(years))
    msg("no --years given; using all", length(spec$years), "years with data files:",
        min(spec$years), "-", max(spec$years))
  spec
}

# -- Metadata -----------------------------------------------------------------

make_meta <- function(spec, ...) {
  c(spec[c("region", "years", "doy_start", "doy_end", "direction", "tag")],
    list(...),
    list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
         R_version = R.version.string))
}

# Stop if the metadata stored in a file does not match the run requested.
check_meta <- function(meta, spec, what = "input file") {
  for (k in c("region", "doy_start", "doy_end", "direction")) {
    if (!identical(as.character(meta[[k]]), as.character(spec[[k]])))
      stop(sprintf("%s has %s = %s but %s was requested.", what, k,
                   paste(meta[[k]], collapse = ","),
                   paste(spec[[k]], collapse = ",")), call. = FALSE)
  }
  if (!identical(as.integer(meta$years), as.integer(spec$years)))
    stop(what, " was made for different years than requested.", call. = FALSE)
  invisible(TRUE)
}

load_rda_object <- function(file, name) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  e <- new.env()
  load(file, envir = e)
  if (!exists(name, envir = e, inherits = FALSE))
    stop(file, " does not contain an object named '", name, "'.")
  get(name, envir = e)
}

# Number of workers: --workers, else SLURM_CPUS_PER_TASK, else 1.
n_workers <- function(args) {
  w <- if (!is.null(args$workers)) args$workers else
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")
  w <- suppressWarnings(as.integer(w))
  if (is.na(w) || w < 1) 1L else w
}

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "--", ..., "\n")

# Build a run spec from R (the scripts build it from command-line arguments).
# years = NULL means all years with a data file covering the window.
make_spec <- function(region, years = NULL, doy_start, doy_end, direction) {
  check_window(doy_start, doy_end)
  check_direction(direction)
  cfg <- get_region(region)
  if (is.null(years)) years <- data_years(cfg, doy_end)
  years <- sort(unique(as.integer(years)))
  list(region = region, years = years, doy_start = as.integer(doy_start),
       doy_end = as.integer(doy_end), direction = check_direction(direction),
       tag = run_tag(region, years, doy_start, doy_end, direction))
}
