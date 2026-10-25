# setup.R
#
# Load the region config and all shared functions.  Every script calls
#   source(file.path(refactor_dir, "R", "setup.R")); onset_setup(refactor_dir)
# where refactor_dir is the path of the `refactor/` directory.  Relative paths
# in config/regions.R (data folders, gold standard files, shapefiles) are
# resolved against refactor_dir.

onset_setup <- function(refactor_dir, envir = globalenv()) {
  refactor_dir <- normalizePath(refactor_dir, mustWork = TRUE)
  options(onset.refactor_dir = refactor_dir)

  sys.source(file.path(refactor_dir, "config", "regions.R"), envir = envir)

  r_files <- list.files(file.path(refactor_dir, "R"), pattern = "\\.R$",
                        full.names = TRUE)
  r_files <- r_files[basename(r_files) != "setup.R"]
  for (f in r_files) sys.source(f, envir = envir)

  invisible(refactor_dir)
}

# Directory of the currently running Rscript (falls back to the working
# directory when run interactively).
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
