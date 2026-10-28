# make_maps.R
#
# Step 3: one PDF per map (mean and SD across years) for a run, written to
# output/<region>/<tag>/maps/.  Uses changepoints.Rda, stats.Rda (if
# present) and any gold standard dataset for the run's direction.
#
# Usage:
#   Rscript scripts/make_maps.R --region=india --years=1951:2025 \
#       --doy=1:240 --direction=onset [--summaries=mean,sd] [--maps=chpt,dev_expl]
#
# --summaries   names from year_summaries in R/maps.R (default: all)
# --maps        names of map specs to draw (default: all available)

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

usage <- paste(
  "Usage: Rscript scripts/make_maps.R --region=<region> [--years=<y1:y2>]",
  "--doy=<start:end> --direction=<onset|retreat> [--summaries=mean,sd]",
  "[--maps=a,b]")

args <- parse_cli()
spec <- run_spec_from_args(args, usage)

split_arg <- function(x) if (is.null(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]

run_maps_step(spec,
              maps      = split_arg(args$maps),
              summaries = if (is.null(args$summaries)) names(year_summaries) else
                            split_arg(args$summaries))
