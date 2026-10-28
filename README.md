# Onset/retreat analysis pipeline

Region-agnostic version of the India and Ethiopia onset code. 

Every run is defined by four inputs: **region**, **years**, **window**
(`--doy=start:end`, within one calendar year, end at most 365) and
**direction** (`onset` = dry to wet, `retreat` = wet to dry). Outputs go to
`output/<region>/<tag>/`, where the tag is e.g. `india_1951-2025_doy1-240_onset`.

## Steps

Scripts can be run from any directory. `--years` is optional: without it,
every year with a data file in the region's data folder is used (years whose
file ends before the window does are skipped, with a message). Steps 2 and 3
look the years up again the same way, so if you add data files between steps,
pass `--years` explicitly to keep referring to the same run.

```sh
# 1. Changepoints (serial; India, all years: about a minute plus file reading)
Rscript scripts/estimate_changepoints.R --region=india \
    --years=1901:2025 --doy=1:240 --direction=onset

# 2. Summary statistics (parallel over longitudes)
Rscript scripts/compute_stats.R --region=india \
    --years=1901:2025 --doy=1:240 --direction=onset --workers=24

# 3. Maps -> output/<region>/<tag>/maps/*.pdf
Rscript scripts/make_maps.R --region=india \
    --years=1901:2025 --doy=1:240 --direction=onset

# 4. Browse runs
R -e 'shiny::runApp("shiny")'
```

### From an R session

Each step is also an R function (the scripts are thin wrappers around them).
A run is described by a `spec` from `make_spec()`; leave out `years` to use
all years with data files.

```r
source("R/setup.R"); onset_setup(".")   # run from the base directory

spec <- make_spec("india", years = 1951:2025, doy_start = 1, doy_end = 240,
                  direction = "onset")

res   <- run_changepoints_step(spec)                 # gamma_floor =, min_seg =
stats <- run_stats_step(spec, metrics = c("pre_post", "dev_full"), workers = 8)
run_maps_step(spec, maps = c("chpt", "dev_rel"))     # summaries = c("mean", "sd")

run <- load_run(spec$region, spec$tag)   # all arrays for a run, e.g. run$r$chpt
```

`run_stats_step()` takes registry names or your own functions (see below);
`metrics = names(metric_registry)` runs them all. After editing
`config/regions.R` or a file in `R/`, call `onset_setup(".")` again to
reload.

Optional arguments:

| Script | Argument | Default |
| --- | --- | --- |
| estimate_changepoints | `--gamma_floor`, `--min_seg` | 1 mm, 5 days |
| compute_stats | `--metrics=dev_full,pre_post,...` | all registry metrics |
| compute_stats | `--metric_args=dev_full.spline_k=30` | registry defaults |
| compute_stats | `--workers`, `--plan=multicore\|multisession` | `SLURM_CPUS_PER_TASK` or 1; multicore if supported |
| make_maps | `--maps=chpt,dev_rel`, `--summaries=mean,sd` | all available |

`compute_stats` replaces only the arrays for the metrics it computes, so
metrics can be added to an existing `stats.Rda` one at a time. `dev_full` and
`transition_edf` fit one GAM per cell-year and are the slow ones (roughly
20-40 ms each per cell-year per core).

## Changepoint model

Values below `gamma_floor` are raised to it; each candidate changepoint `c`
(from `min_seg` to `n - min_seg`) splits the window into two segments with
their own Gamma mean. The deviance is computed in closed form from the
segment means (see `R/changepoints.R`), so no model is fit. `chpt` is the last
day of the pre segment, as day of year. Any missing day in the window gives
`NA`; if no candidate has the right direction, `chpt` is `NA` and `dev_expl` 0.

`scripts/validate_changepoints.R` checks this against `glm()` on a random
sample of cell-years, and with `--old=<file>` compares against an older
results file:

```sh
Rscript scripts/validate_changepoints.R --region=india \
    --years=2020:2021 --doy=1:240 --direction=onset --n=100 \
    --old=scratch/india_changepoints_doy1-240.Rda
```

(`--old` reads the new changepoints from the run's `changepoints.Rda`, so run
step 1 for the same inputs first.)

## Your own metrics

A metric is `f(y, ctx, ...)`: `y` is the raw window for one cell-year, `ctx`
has `chpt` (index into `y`), `chpt_doy`, `dev_expl`, `doy`, `direction`,
`gamma_floor`, `year`, `lon`, `lat`. Return one number or a named vector. Add
it to `metric_registry` in `R/metrics.R` to use it from the command line, or
call the step from R:

```r
source("R/setup.R"); onset_setup(".")
spec <- make_spec("india", 1901:2025, 1, 240, "onset")
# or all years with data files:
# spec <- make_spec("india", doy_start = 1, doy_end = 240, direction = "onset")
run_stats_step(spec, workers = 8, metrics = list(
  max_pre = function(y, ctx) max(y[seq_len(ctx$chpt)]),
  wet_days = list(fun = function(y, ctx, thr) c(pre = mean(y[1:ctx$chpt] > thr),
                                                post = mean(y[-(1:ctx$chpt)] > thr)),
                  args = list(thr = 1), outputs = c("pre_wet", "post_wet"))))
```

Metric functions should be self-contained (package functions are fine;
list the packages in `packages`), because with a `multisession` plan they are
sent to worker processes without the rest of the code. To map a new
statistic, add an entry to `map_specs()` in `R/maps.R`; the Shiny app picks
it up too.

## Regions and gold standard dates

`config/regions.R` holds everything region-specific: data folder and file
names, NetCDF variable and dimension names, map borders, an optional
shapefile, color-limit overrides, Shiny defaults, and optional gold standard
datasets. A gold standard file is a NetCDF holding all years on the
precipitation grid; `origin` gives the date a value of 1 refers to, and day
of year is computed per year, so leap years are handled (O2: `"04-01"`, which
gives doy 91 in normal years and 92 in leap years). Datasets are used only for
runs with matching `direction`; the one with `overlay = TRUE` is drawn on the
Shiny time series.

### Map color limits

Both `make_maps` and the Shiny app resolve each map's color limits in this
order (later wins):

1. Defaults in `map_specs()` (`R/maps.R`). The changepoint map (and the gold
   standard maps) default to the run's own window, e.g. 150 to 365 for a
   `--doy=150:365` run, so they usually need no setting.
2. `zlim$<direction>` in the region config: applies to every window of that
   direction.
3. `zlim$windows$"<direction>_doy<start>-<end>"`: applies to that window only.

```r
zlim = list(
  onset   = list(dev_expl = list(mean = c(0, 1))),
  retreat = list(),
  windows = list(
    "onset_doy1-240"     = list(chpt = list(mean = c(0, 240))),
    "retreat_doy210-365" = list(chpt = list(mean = c(210, 365), sd = c(0, 40)))
  )
)
```

Give only the panels you want to change (`mean`, `sd`); `NULL` means
automatic. Map names are those of the PDFs in `maps/` (`chpt`, `dev_expl`,
`dev_full`, `dev_rel`, `pre_mean`, `post_mean`, `post_frac`, `transition_sd`,
`transition_edf`).

## Files

| Path | Contents |
| --- | --- |
| `config/regions.R` | Region settings and analysis defaults |
| `R/setup.R` | `onset_setup()`: loads config and all functions |
| `R/utils.R` | Run tags, argument parsing, metadata checks |
| `R/read_data.R` | `read_region()`: window of the yearly NetCDF files |
| `R/changepoints.R` | Closed-form changepoint model |
| `R/metrics.R` | Metric functions and registry |
| `R/stats_driver.R` | `compute_stats()`, `run_stats_step()` |
| `R/gold_standard.R` | Gold standard reader |
| `R/maps.R` | Map specs, year summaries, map plotting |
| `R/run_io.R` | `load_run()`, `list_runs()` |
| `scripts/` | The step scripts and the validation script |
| `shiny/` | Adapted copy of the original app (`scratch/shiny`, untouched) |
