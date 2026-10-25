# make_maps.R
#
# Step 3: one PDF per map (mean and SD across years) for a run, written to
# refactor/output/<region>/<tag>/maps/.  Uses changepoints.Rda, stats.Rda (if
# present) and any gold standard dataset for the run's direction.
#
# Usage:
#   Rscript refactor/scripts/make_maps.R --region=india --years=1951:2025 \
#       --doy=1:240 --direction=onset [--summaries=mean,sd] [--maps=chpt,dev_rel]
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

run <- load_run(spec$region, spec$tag)
check_meta(run$meta, spec, "changepoints.Rda")

summaries <- if (!is.null(args$summaries))
  strsplit(args$summaries, ",", fixed = TRUE)[[1]] else names(year_summaries)
bad <- setdiff(summaries, names(year_summaries))
if (length(bad)) stop("Unknown summaries: ", paste(bad, collapse = ", "))

specs <- specs_for_run(run$r, run$cfg, spec$direction)
if (!is.null(args$maps)) {
  want <- strsplit(args$maps, ",", fixed = TRUE)[[1]]
  missing <- setdiff(want, names(specs))
  if (length(missing))
    warning("Not available for this run (skipped): ", paste(missing, collapse = ", "))
  specs <- specs[intersect(want, names(specs))]
}

overlay <- load_overlay(run$cfg)
out_dir <- maps_dir(spec$region, spec$tag)
yr_txt  <- sprintf("%d-%d", min(run$years), max(run$years))

for (nm in names(specs)) {
  s <- specs[[nm]]
  f <- file.path(out_dir, paste0(nm, ".pdf"))
  pdf(f, width = 4 * length(summaries), height = 4)
  plot_map_summaries(
    s$fun(run$r), run$lon, run$lat,
    title = sprintf("%s  [%s, doy %d-%d, %s, %s]", s$label, spec$region,
                    spec$doy_start, spec$doy_end, spec$direction, yr_txt),
    zlim = s$zlim, summaries = summaries,
    map_regions = run$cfg$map_regions, overlay = overlay)
  dev.off()
  msg("wrote", f)
}
msg(length(specs), "maps written to", out_dir)
