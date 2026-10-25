# estimate_changepoints.R
#
# Step 1: changepoint day and deviance explained for every cell and year.
#
# Usage (from any directory):
#   Rscript refactor/scripts/estimate_changepoints.R --region=india \
#       --years=1951:2025 --doy=1:240 --direction=onset \
#       [--gamma_floor=1] [--min_seg=5]
#
# --years is optional (here and in the other scripts): without it, every year
# with a data file covering the window is used.
#
# Output: refactor/output/<region>/<tag>/changepoints.Rda (object `results`;
# see R/run_io.R).  Runs serially; India for all years takes a few minutes of
# computation plus the time to read the files.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

usage <- paste(
  "Usage: Rscript scripts/estimate_changepoints.R --region=<region>",
  "[--years=<y1:y2>] --doy=<start:end> --direction=<onset|retreat>",
  "[--gamma_floor=1] [--min_seg=5]")

args <- parse_cli()
spec <- run_spec_from_args(args, usage)
cfg  <- get_region(spec$region)

gamma_floor <- if (!is.null(args$gamma_floor)) as.numeric(args$gamma_floor) else
  analysis_defaults$gamma_floor
min_seg <- if (!is.null(args$min_seg)) as.integer(args$min_seg) else
  analysis_defaults$min_seg

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
