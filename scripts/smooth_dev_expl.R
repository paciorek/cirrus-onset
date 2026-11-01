# smooth_dev_expl.R
#
# Fit a spatial smooth of mean dev_expl for a run and save it to
# output/<region>/<tag>/smooth/, where the Shiny app finds it (so slow fits,
# e.g. betar with large k, can be run offline or under Slurm).
#
# Usage:
#   Rscript scripts/smooth_dev_expl.R --region=india --doy=1:240 \
#       --direction=onset [--years=1901:2025] --k=500 --family=betar \
#       [--smooth_years=2000:2024] [--engine=gam|bam] [--refit=1] \
#       [--threshold=0.3]
#
# --years         identifies the run (as for the other scripts)
# --smooth_years  years to average dev_expl over (default: all run years).
#                 The Shiny app averages over the years it displays, so to
#                 precompute for the app use the same years (e.g. 2000:2024).
# --family        gaussian, betar or quasibinomial (default gaussian)
# --threshold     also save the region mask (smoothed mean dev_expl >=
#                 threshold) to output/<region>/<tag>/masks/

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

usage <- paste(
  "Usage: Rscript scripts/smooth_dev_expl.R --region=<region>",
  "[--years=<y1:y2>] --doy=<start:end> --direction=<onset|retreat>",
  "[--k=500] [--family=gaussian|betar|quasibinomial]",
  "[--smooth_years=<y1:y2>] [--engine=gam|bam] [--refit=1] [--threshold=0.3]")

args <- parse_cli()
spec <- run_spec_from_args(args, usage)

sm <- run_smooth_step(
  spec,
  k      = if (!is.null(args$k)) as.integer(args$k) else 500L,
  family = if (!is.null(args$family)) args$family else "gaussian",
  years  = if (!is.null(args$smooth_years)) parse_years(args$smooth_years) else NULL,
  engine = if (!is.null(args$engine)) args$engine else "gam",
  refit  = isTRUE(as.logical(as.integer(if (is.null(args$refit)) 0 else args$refit))))

if (!is.null(args$threshold)) {
  run <- load_run(spec$region, spec$tag, with_gold = FALSE)
  m <- make_region_mask(run, as.numeric(args$threshold), source = "smooth",
                        years = sm$years, smooth = sm)
  print(m)
  msg("mask saved to", save_region_mask(m))
}
