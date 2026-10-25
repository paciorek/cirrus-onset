# config.R
#
# Locate the refactor/ directory and load the shared region config and
# functions (config/regions.R and refactor/R/*.R), so the app uses the same
# readers, map specs and gold standard handling as the pipeline.
#
# The refactor directory is the parent of this app's directory by default;
# set ONSET_REFACTOR_DIR to override.

local({
  env <- Sys.getenv("ONSET_REFACTOR_DIR", unset = "")
  candidates <- c(env, "..", "~/research/boos/onset/refactor")
  candidates <- candidates[nzchar(candidates)]
  found <- NULL
  for (d in candidates)
    if (file.exists(file.path(d, "R", "setup.R"))) { found <- d; break }
  if (is.null(found))
    stop("Could not find the refactor directory (looked in: ",
         paste(candidates, collapse = ", "),
         "). Set ONSET_REFACTOR_DIR to the path of refactor/.")
  source(file.path(found, "R", "setup.R"), local = TRUE)
  onset_setup(found, envir = globalenv())
})
