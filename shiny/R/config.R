# config.R
#
# Locate the directory and load the shared region config and
# functions (config/regions.R and R/*.R), so the app uses the same
# readers, map specs and gold standard handling as the pipeline.
#
# The root directory is the parent of this app's directory by default;
# set ONSET_DIR to override.

local({
  env <- Sys.getenv("ONSET_DIR", unset = "")
  candidates <- c(env, "..", "../..")
  candidates <- candidates[nzchar(candidates)]
  found <- NULL
  for (d in candidates)
    if (file.exists(file.path(d, "R", "setup.R"))) { found <- d; break }
  if (is.null(found))
    stop("Could not find the root directory (looked in: ",
         paste(candidates, collapse = ", "),
         "). Set ONSET_DIR to the path of the root dir.")
  source(file.path(found, "R", "setup.R"), local = TRUE)
  onset_setup(found, envir = globalenv())
})
