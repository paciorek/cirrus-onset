# setup.R
#
# Load the region config and all shared functions.  Every script calls
#   source(file.path(base_dir, "R", "setup.R")); onset_setup(base_dir)
# where base_dir is the path of repository directory.  Relative paths
# in config/regions.R (data folders, gold standard files, shapefiles) are
# resolved against `base_dir`.

onset_setup <- function(base_dir, envir = globalenv()) {
  base_dir <- normalizePath(base_dir, mustWork = TRUE)
  options(onset.base_dir = base_dir)

  sys.source(file.path(base_dir, "config", "regions.R"), envir = envir)

  r_files <- list.files(file.path(base_dir, "R"), pattern = "\\.R$",
                        full.names = TRUE)
  r_files <- r_files[basename(r_files) != "setup.R"]
  for (f in r_files) sys.source(f, envir = envir)

  invisible(base_dir)
}

# Directory of the currently running Rscript (falls back to the working
# directory when run interactively).
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
