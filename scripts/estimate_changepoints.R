# estimate_changepoints.R
#
# Step 1: changepoint day and deviance explained for every cell and year.
#
# Usage (from any directory):
#   Rscript scripts/estimate_changepoints.R --region=india \
#       --years=1951:2025 --doy=1:240 --direction=onset \
#       [--gamma_floor=1] [--min_seg=5]
#
# --years is optional (here and in the other scripts): without it, every year
# with a data file covering the window is used.
#
# Output: output/<region>/<tag>/changepoints.Rda (object `results`;
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
gamma_floor <- if (!is.null(args$gamma_floor)) as.numeric(args$gamma_floor) else
  analysis_defaults$gamma_floor
# (min_seg and gamma_floor are checked in run_changepoints_step())
min_seg <- if (!is.null(args$min_seg)) as.numeric(args$min_seg) else
  analysis_defaults$min_seg

run_changepoints_step(spec, gamma_floor = gamma_floor, min_seg = min_seg)
