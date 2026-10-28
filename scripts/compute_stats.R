# compute_stats.R
#
# Step 2: per cell-year summary statistics from metric functions, given the
# changepoints from step 1.  Runs in parallel over longitudes.
#
# Usage:
#   Rscript scripts/compute_stats.R --region=india --years=1951:2025 \
#       --doy=1:240 --direction=onset \
#       [--metrics=dev_full,pre_post,transition_sd,transition_edf] \
#       [--metric_args=dev_full.spline_k=20,transition_sd.window=7] \
#       [--workers=8] [--plan=multicore|multisession]
#
# --metrics     registry names from R/metrics.R (default: all)
# --workers     default: SLURM_CPUS_PER_TASK if set, else 1
#
# Output: output/<region>/<tag>/stats.Rda (object `stats`).  If the
# file exists, arrays for the metrics computed now are replaced and the others
# kept.
#
# To use your own metric functions, call compute_stats() from R instead (see
# README.md).

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

usage <- paste(
  "Usage: Rscript scripts/compute_stats.R --region=<region> [--years=<y1:y2>]",
  "--doy=<start:end> --direction=<onset|retreat> [--metrics=a,b]",
  "[--metric_args=m.arg=value,...] [--workers=N] [--plan=multicore|multisession]")

args <- parse_cli()
spec <- run_spec_from_args(args, usage)
cfg  <- get_region(spec$region)

metric_names <- if (!is.null(args$metrics))
  strsplit(args$metrics, ",", fixed = TRUE)[[1]] else names(metric_registry)
metrics <- normalize_metrics(metric_names)
metrics <- apply_metric_args(metrics, args$metric_args)
workers <- n_workers(args)

run_stats_step(spec, cfg, metrics, workers = workers, plan = args$plan)
