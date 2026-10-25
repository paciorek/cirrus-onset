# run_io.R
#
# Reading and writing a run's output files, and listing available runs.
#
# output/<region>/<tag>/changepoints.Rda  object `results`:
#     chpt, dev_expl   lon x lat x year arrays (chpt = absolute day of year)
#     lon, lat, years, land (lon x lat logical: cell has data), meta
# output/<region>/<tag>/stats.Rda         object `stats`:
#     arrays   named list of lon x lat x year arrays
#     meta     run metadata plus the metrics computed and their settings

changepoints_path <- function(region, tag) {
  file.path(run_dir(region, tag, create = FALSE), "changepoints.Rda")
}

stats_path <- function(region, tag) {
  file.path(run_dir(region, tag, create = FALSE), "stats.Rda")
}

maps_dir <- function(region, tag) {
  d <- file.path(run_dir(region, tag, create = FALSE), "maps")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Everything needed to map or browse a run.
#
# Returns list(r, lon, lat, years, land, meta, stats_meta, gold, cfg) where r
# holds all lon x lat x year arrays: chpt, dev_expl, any stats, and gold
# standard arrays named gold_<name>.
load_run <- function(region, tag, with_gold = TRUE) {
  res <- load_rda_object(changepoints_path(region, tag), "results")
  r <- res[c("chpt", "dev_expl")]

  stats_meta <- NULL
  sf <- stats_path(region, tag)
  if (file.exists(sf)) {
    st <- load_rda_object(sf, "stats")
    r <- c(r, st$arrays)
    stats_meta <- st$meta
  }

  cfg  <- get_region(region)
  gold <- character(0)
  if (with_gold) {
    g <- read_gold_standards(cfg, res$meta$direction, res$lon, res$lat, res$years)
    if (length(g)) {
      names(g) <- paste0("gold_", names(g))
      r <- c(r, g)
      gold <- names(g)
    }
  }

  list(r = r, lon = res$lon, lat = res$lat, years = res$years, land = res$land,
       meta = res$meta, stats_meta = stats_meta, gold = gold, cfg = cfg)
}

# Runs with a changepoints.Rda under output/, as a data frame with one row
# per run (region, tag, parsed tag fields, label).
list_runs <- function() {
  files <- Sys.glob(file.path(output_root(), "*", "*", "changepoints.Rda"))
  if (!length(files)) return(data.frame())
  tags <- basename(dirname(files))
  regs <- basename(dirname(dirname(files)))
  rows <- lapply(seq_along(tags), function(k) {
    p <- parse_run_tag(tags[k])
    if (is.null(p) || !regs[k] %in% names(regions)) return(NULL)
    data.frame(region = regs[k], tag = tags[k],
               year_start = p$year_start, year_end = p$year_end,
               doy_start = p$doy_start, doy_end = p$doy_end,
               direction = p$direction,
               # paste0 rather than sprintf: sprintf can mangle non-ASCII
               # characters in the region label under a non-UTF-8 locale
               label = paste0(regions[[regs[k]]]$label, " \u00b7 ",
                              p$year_start, "\u2013", p$year_end,
                              " \u00b7 doy ", p$doy_start, "\u2013", p$doy_end,
                              " \u00b7 ", p$direction),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$region, out$direction, out$doy_start, out$year_start), ]
}
