# stats_driver.R
#
# compute_stats(): apply metric functions to every cell-year, in parallel over
# longitudes.  Each worker receives one longitude slice (lat x day x year) and
# handles all lat-year cases for it.
#
# Slices are sent to workers in batches (one longitude per worker per batch),
# so only a few slices are copied at a time even with a "multisession" plan.

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

# Compute all metrics for one longitude.
#
# slice     : lat x day x year precipitation (the window)
# chpt      : lat x year changepoints (absolute doy)
# dev_expl  : lat x year
# Returns list(values = list(<output> = lat x year matrix), errors = named int)
stats_one_lon <- function(slice, chpt, dev_expl, lon, lat, years, doy,
                          metrics, gamma_floor, direction) {
  n_lat  <- dim(slice)[1]
  n_year <- dim(slice)[3]
  outs   <- unlist(lapply(metrics, `[[`, "outputs"))
  values <- lapply(setNames(outs, outs), function(o)
    matrix(NA_real_, n_lat, n_year))
  errors <- setNames(integer(length(metrics)), names(metrics))
  warns  <- errors
  first_error <- setNames(rep(NA_character_, length(metrics)), names(metrics))

  for (j in seq_len(n_lat)) {
    for (k in seq_len(n_year)) {
      y <- slice[j, , k]
      if (anyNA(y)) next
      ctx <- list(chpt = chpt[j, k] - doy[1] + 1, chpt_doy = chpt[j, k],
                  dev_expl = dev_expl[j, k], doy = doy, direction = direction,
                  gamma_floor = gamma_floor, year = years[k],
                  lon = lon, lat = lat[j])
      for (m in names(metrics)) {
        mt <- metrics[[m]]
        if (mt$needs_chpt && is.na(ctx$chpt)) next
        r <- tryCatch(
          withCallingHandlers(do.call(mt$fun, c(list(y, ctx), mt$args)),
                              warning = function(w) {
                                warns[m] <<- warns[m] + 1L
                                invokeRestart("muffleWarning")
                              }),
          error = function(e) e)
        if (inherits(r, "error")) {
          errors[m] <- errors[m] + 1L
          if (is.na(first_error[m])) first_error[m] <- conditionMessage(r)
          next
        }
        if (length(r) != length(mt$outputs)) {
          errors[m] <- errors[m] + 1L
          if (is.na(first_error[m]))
            first_error[m] <- sprintf("returned %d values, expected %d",
                                      length(r), length(mt$outputs))
          next
        }
        for (q in seq_along(mt$outputs))
          values[[mt$outputs[q]]][j, k] <- as.numeric(r[q])
      }
    }
  }
  list(values = values, errors = errors, warnings = warns,
       first_error = first_error)
}

# Worker wrapper, defined at top level (not as a closure inside
# compute_stats) so that the full data array is not serialized with it.
stats_worker <- function(s, common) do.call(stats_one_lon, c(s, common))

# data     : list from read_region()
# chg      : changepoint results (list with chpt, dev_expl arrays)
# metrics  : anything normalize_metrics() accepts
# workers  : number of parallel workers (1 = sequential)
# plan     : "multicore" (forked; shares memory, Linux/macOS outside RStudio)
#            or "multisession"; NULL picks multicore when supported
#
# Returns a named list of lon x lat x year arrays, plus attribute "errors"
# (count of failed cell-years per metric; their values are NA).  Warnings
# from metric functions are counted and reported rather than printed.
compute_stats <- function(data, chg, metrics, gamma_floor = 1,
                          direction = "onset", workers = 1, plan = NULL,
                          verbose = TRUE) {
  metrics <- normalize_metrics(metrics)
  d <- dim(data$rain)
  n_lon <- d[1]; n_lat <- d[2]; n_year <- d[4]
  if (!identical(dim(chg$chpt), c(n_lon, n_lat, n_year)))
    stop("Changepoint arrays do not match the data dimensions.")

  packages <- unique(unlist(lapply(metrics, `[[`, "packages")))
  for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))

  args_common <- list(lat = data$lat, years = data$years, doy = data$doy,
                      metrics = metrics, gamma_floor = gamma_floor,
                      direction = direction)
  # One longitude's inputs, keeping dimensions even when n_lat or n_year is 1
  lon_inputs <- function(i) list(
    slice    = array(data$rain[i, , , , drop = FALSE], d[2:4]),
    chpt     = matrix(chg$chpt[i, , , drop = FALSE], n_lat, n_year),
    dev_expl = matrix(chg$dev_expl[i, , , drop = FALSE], n_lat, n_year),
    lon      = data$lon[i])
  run_one <- function(i) do.call(stats_one_lon, c(lon_inputs(i), args_common))

  if (workers > 1) {
    if (is.null(plan))
      plan <- if (supportsMulticore()) "multicore" else "multisession"
    old_plan <- future::plan(plan, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
    if (verbose) msg("compute_stats:", workers, "workers,", plan, "plan")
  }

  results <- vector("list", n_lon)
  batch_size <- if (workers > 1) 4L * workers else 1L
  batches <- split(seq_len(n_lon), ceiling(seq_len(n_lon) / batch_size))
  for (b in batches) {
    if (workers > 1) {
      slices <- lapply(b, lon_inputs)
      res <- future_lapply(slices, stats_worker, common = args_common,
                           future.packages = packages,
                           future.globals = list(stats_one_lon = stats_one_lon),
                           future.seed = TRUE, future.chunk.size = 1)
      results[b] <- res
    } else {
      results[b] <- lapply(b, run_one)
    }
    if (verbose) msg("compute_stats: lon", max(b), "of", n_lon)
  }

  dn <- dimnames(chg$chpt)
  outs <- unlist(lapply(metrics, `[[`, "outputs"))
  stats <- lapply(setNames(outs, outs), function(o) {
    a <- array(NA_real_, c(n_lon, n_lat, n_year), dimnames = dn)
    for (i in seq_len(n_lon)) a[i, , ] <- results[[i]]$values[[o]]
    a
  })
  errors <- Reduce(`+`, lapply(results, `[[`, "errors"))
  warns  <- Reduce(`+`, lapply(results, `[[`, "warnings"))
  first_error <- Reduce(function(a, b) ifelse(is.na(a), b, a),
                        lapply(results, `[[`, "first_error"))
  if (verbose && any(errors > 0)) {
    for (m in names(errors)[errors > 0])
      msg(sprintf("metric %s failed for %d cell-years (first error: %s)",
                  m, errors[m], first_error[m]))
  }
  if (verbose && any(warns > 0)) {
    for (m in names(warns)[warns > 0])
      msg(sprintf("metric %s gave warnings for %d cell-years (e.g. GAM convergence)",
                  m, warns[m]))
  }
  attr(stats, "errors") <- errors
  stats
}

# The whole compute_stats step for a run: load the changepoints, read the
# data, compute the metrics and save (merging with an existing stats.Rda).
# Can be called from R with your own metric functions, e.g.
#   spec <- make_spec("india", 1951:2025, 1, 240, "onset")
#   run_stats_step(spec, metrics = list(max_post = function(y, ctx)
#                    max(y[(ctx$chpt + 1):length(y)])), workers = 4)
run_stats_step <- function(spec, cfg = get_region(spec$region), metrics,
                           workers = 1, plan = NULL) {
  metrics <- normalize_metrics(metrics)
  res <- load_rda_object(changepoints_path(spec$region, spec$tag), "results")
  check_meta(res$meta, spec, "changepoints.Rda")

  msg("compute_stats:", spec$tag)
  msg("metrics:", paste(names(metrics), collapse = ", "), "| workers:", workers)

  t0   <- Sys.time()
  data <- read_region(cfg, spec$years, spec$doy_start, spec$doy_end)
  msg("data read in", format(round(Sys.time() - t0, 1)))

  t1 <- Sys.time()
  arrays <- compute_stats(data, res, metrics,
                          gamma_floor = res$meta$gamma_floor,
                          direction = spec$direction,
                          workers = workers, plan = plan)
  errors <- attr(arrays, "errors")
  attr(arrays, "errors") <- NULL
  msg("stats computed in", format(round(Sys.time() - t1, 1)))

  info <- lapply(metrics, function(m)
    list(args = m$args, outputs = m$outputs, needs_chpt = m$needs_chpt))

  out_file <- stats_path(spec$region, spec$tag)
  if (file.exists(out_file)) {
    old <- load_rda_object(out_file, "stats")
    check_meta(old$meta, spec, "existing stats.Rda")
    keep_metrics <- Filter(function(m) !any(m$outputs %in% names(arrays)),
                           old$meta$metrics[setdiff(names(old$meta$metrics),
                                                    names(info))])
    keep_arrays <- old$arrays[unlist(lapply(keep_metrics, `[[`, "outputs"))]
    arrays <- c(keep_arrays, arrays)
    info   <- c(keep_metrics, info)
    msg("merged with existing stats; kept:",
        if (length(keep_metrics)) paste(names(keep_metrics), collapse = ", ") else "none")
  }

  stats <- list(
    arrays = arrays,
    meta   = make_meta(spec, gamma_floor = res$meta$gamma_floor,
                       min_seg = res$meta$min_seg, metrics = info,
                       errors = errors,
                       changepoints_created = res$meta$created)
  )
  save(stats, file = out_file)
  msg("saved", out_file)
  invisible(stats)
}
