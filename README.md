# Onset/retreat analysis pipeline

Region-agnostic version of the India and Ethiopia onset code. 

Every run is defined by four inputs: **region**, **years**, **window**
(`--doy=start:end`, within one calendar year, end at most 365) and
**direction** (`onset` = dry to wet, `retreat` = wet to dry). Outputs go to
`output/<region>/<tag>/`, where the tag is e.g. `india_1951-2025_doy1-240_onset`.

## Data

The raw precipitation data are not in the repository. Download them to:

| Region | Data | Folder (relative to this directory) | File names |
| --- | --- | --- | --- |
| India | IMD25x25 | `data/india_onset/IMD_25by25` | `data_YYYY.nc` |
| Ethiopia | CHIRPS | `data/ethiopia_onset/CHIRPS_IMERG` | `YYYY.nc` |

The scripts stop with a message saying where to put the files if they are
missing; the Shiny app still shows maps from existing runs but says so in
place of the time series.

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

## Smoothing dev_expl and choosing regions

To pick out regions whose onset (or retreat) may be predictable:

1. At each grid cell, average the yearly `dev_expl` values over the years
   used, giving one value per location (a spatial field).
2. Smooth that field with a thin-plate spline `s(lon, lat, k)`, which gives
   smoother region boundaries and fewer holes.
3. Threshold the smoothed (or raw) field to get a region mask.

Functions are in `R/smooth.R`:

```r
run <- load_run("india", spec$tag)
z   <- mean_dev_expl(run$r$dev_expl, run$years, years = 2000:2024)  # lon x lat
sm  <- smooth_dev_expl(z, run$lon, run$lat, k = 500, family = "betar")
m   <- make_region_mask(run, threshold = 0.3, source = "smooth",
                        years = 2000:2024, smooth = sm)   # or source = "raw"
m                       # summary: cells kept
save_region_mask(m)     # -> output/<region>/<tag>/masks/mask_...rds
m$keep                  # lon x lat logical (TRUE = in the region)
mask_cells(m)           # data frame: lon, lat, value, keep
m <- load_region_mask("output/india/<tag>/masks/<file>.rds")
```

`family` is `"gaussian"`, `"betar"` or `"quasibinomial"` (or any family
object); `engine = "bam"` uses `bam(discrete = TRUE)`, which helps for
gaussian with large k. Fitting time rises steeply with k, and betar and
quasibinomial are much slower than gaussian (India, about 5000 cells:
gaussian k = 500 about 20 s; betar k = 500 several minutes; betar k = 3000
hours). Slow fits can be run offline (e.g. under Slurm); they are saved to
`output/<region>/<tag>/smooth/` (fitted values only, not the model) and
picked up by the app. `--threshold` also saves a mask:

```sh
Rscript scripts/smooth_dev_expl.R --region=india --doy=1:240 --direction=onset \
    --k=3000 --family=betar --smooth_years=2000:2024 [--threshold=0.3]
```

From R, `run_smooth_step(spec, k =, family =, years =)` does the same, and
`fit_or_load_smooth(run, ...)` loads a saved fit or fits and saves one.

In the Shiny app, the panel under the map statistic holds the smoothing
settings (k, family, "Fit smooth"); the smooth uses the years displayed.

- Map statistic "Smoothed mean dev_expl" shows the smoothed field.
- "Mask cells with low mean dev_expl" leaves uncolored, on any map, the
  cells whose raw or smoothed mean dev_expl is below the threshold.
  "Save mask" writes the mask (as from `save_region_mask()`) to
  `output/<region>/<tag>/masks/`; "Download mask (CSV)" gives lon, lat,
  value and keep for each cell.
- Fits are kept in memory (model and fitted values, for all sessions of the
  app) and saved to disk, so changing the map, summary or threshold, or going
  back to earlier settings, does not refit; a saved fit for the current
  settings and years is used without pressing the button. To precompute for
  the app, use `--smooth_years` equal to the years the app displays.

## Regions and gold standard dates

`config/regions.R` holds everything region-specific: data folder and file
names, NetCDF variable and dimension names, map borders, an optional
shapefile, color-limit overrides, Shiny defaults, and optional gold standard
datasets. A gold standard file is a NetCDF holding all years on the
precipitation grid; `origin` gives the date a value of 1 refers to, and day
of year is computed per year, so leap years are handled (O2: `"04-01"`, which
gives doy 91 in normal years and 92 in leap years). Datasets are used only for
runs with matching `direction`. In the Shiny time series a "Gold standard
line" choice picks which one (or none) is drawn; the one with
`overlay = TRUE` is the default (for India: O2 unrestricted, or O2
restricted by MOK).

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
| `R/smooth.R` | Smoothing and thresholding dev_expl |
| `R/steps.R` | Step functions (`run_changepoints_step()`, `run_maps_step()`, `run_smooth_step()`) |
| `scripts/` | The step scripts, `smooth_dev_expl.R` and the validation script |
| `shiny/` | Adapted copy of the original app (`scratch/shiny`, untouched) |
