# steps.R
#
# The pipeline steps as R functions, so they can be run from an R session as
# well as from the scripts (which are thin wrappers around these).
#
#   source("R/setup.R"); onset_setup(".")
#   spec <- make_spec("india", doy_start = 1, doy_end = 240, direction = "onset")
#   res   <- run_changepoints_step(spec)
#   stats <- run_stats_step(spec, metrics = c("pre_post", "dev_full"), workers = 8)
#   run_maps_step(spec)
#
# run_stats_step() is in R/stats_driver.R.

# Step 1: changepoints for every cell and year; saves changepoints.Rda and
# returns the `results` list invisibly.
run_changepoints_step <- function(spec,
                                  gamma_floor = analysis_defaults$gamma_floor,
                                  min_seg = analysis_defaults$min_seg) {
  # Check settings before reading any data.
  n_win <- spec$doy_end - spec$doy_start + 1
  if (length(min_seg) != 1 || !is.finite(min_seg) || min_seg < 1 ||
      min_seg != round(min_seg))
    stop("min_seg must be a whole number of at least 1 (got ", min_seg, ").",
         call. = FALSE)
  if (2 * min_seg > n_win)
    stop("min_seg = ", min_seg, " is too large for a ", n_win,
         "-day window (each segment needs min_seg days).", call. = FALSE)
  if (length(gamma_floor) != 1 || !is.finite(gamma_floor) || gamma_floor <= 0)
    stop("gamma_floor must be a positive number (got ", gamma_floor, ").",
         call. = FALSE)
  min_seg <- as.integer(min_seg)

  cfg <- get_region(spec$region)
  out_file <- file.path(run_dir(spec$region, spec$tag), "changepoints.Rda")

  msg("estimate_changepoints:", spec$tag)
  msg("gamma_floor =", gamma_floor, " min_seg =", min_seg)

  t0   <- Sys.time()
  data <- read_region(cfg, spec$years, spec$doy_start, spec$doy_end)
  msg("data read in", format(round(Sys.time() - t0, 1)))

  t1  <- Sys.time()
  chg <- estimate_changepoints_grid(data, gamma_floor = gamma_floor,
                                    min_seg = min_seg, direction = spec$direction)
  msg("changepoints computed in", format(round(Sys.time() - t1, 1)))

  results <- list(
    chpt     = chg$chpt,
    dev_expl = chg$dev_expl,
    lon      = data$lon,
    lat      = data$lat,
    years    = data$years,
    land     = apply(!is.na(chg$dev_expl), c(1, 2), any),
    meta     = make_meta(spec, gamma_floor = gamma_floor, min_seg = min_seg,
                         data_files = data$files)
  )
  save(results, file = out_file)
  msg("saved", out_file)
  invisible(results)
}

# Step 3: one PDF per map into output/<region>/<tag>/maps/.
#
# maps      : names of map specs to draw (NULL = all available for the run)
# summaries : names from year_summaries (default: all)
#
# Returns the paths of the PDFs written, invisibly.
run_maps_step <- function(spec, maps = NULL, summaries = names(year_summaries)) {
  run <- load_run(spec$region, spec$tag)
  check_meta(run$meta, spec, "changepoints.Rda")

  bad <- setdiff(summaries, names(year_summaries))
  if (length(bad)) stop("Unknown summaries: ", paste(bad, collapse = ", "))

  specs <- specs_for_run(run$r, run$cfg, spec$direction,
                         spec$doy_start, spec$doy_end)
  if (!is.null(maps)) {
    missing <- setdiff(maps, names(specs))
    if (length(missing))
      warning("Not available for this run (skipped): ",
              paste(missing, collapse = ", "))
    specs <- specs[intersect(maps, names(specs))]
  }

  overlay <- load_overlay(run$cfg)
  out_dir <- maps_dir(spec$region, spec$tag)
  yr_txt  <- sprintf("%d-%d", min(run$years), max(run$years))
  files   <- character(0)

  for (nm in names(specs)) {
    s <- specs[[nm]]
    f <- file.path(out_dir, paste0(nm, ".pdf"))
    pdf(f, width = 4 * length(summaries), height = 4)
    tryCatch(
      plot_map_summaries(
        s$fun(run$r), run$lon, run$lat,
        title = sprintf("%s  [%s, doy %d-%d, %s, %s]", s$label, spec$region,
                        spec$doy_start, spec$doy_end, spec$direction, yr_txt),
        zlim = s$zlim, summaries = summaries,
        map_regions = run$cfg$map_regions, overlay = overlay),
      finally = dev.off())
    msg("wrote", f)
    files <- c(files, f)
  }
  msg(length(specs), "maps written to", out_dir)
  invisible(files)
}
