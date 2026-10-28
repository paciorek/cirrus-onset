# regions.R
#
# Region configuration.  Everything that differs between regions lives here,
# so adding a region (with yearly NetCDF precipitation files) means adding an
# entry to `regions`.
#
# Relative paths are resolved against the refactor/ directory.
#
# Fields:
#   label          display name
#   data_dir       folder with one NetCDF file per year
#   file_fmt       sprintf format giving a year's file name
#   var_name       precipitation variable in the NetCDF files
#   dim_names      names of the lon, lat and time dimensions in the files
#   map_regions    country outlines drawn on maps (maps::map("world", ...))
#   shapefile      optional polygon overlay drawn on maps (NULL = none)
#   gold_standard  NULL, or a named list of gold standard datasets.  Each:
#                    file       NetCDF file holding all years
#                    var_name   variable name
#                    dim_names  c(lon = , lat = , year = ) dimension names
#                    direction  "onset" or "retreat": used only for runs with
#                               this direction
#                    origin     "MM-DD": the date that a value of 1 refers to;
#                               day of year is computed per year, so Feb 29
#                               is handled correctly
#                    label      display name
#                    overlay    TRUE for the (one) dataset drawn as a vertical
#                               bar in the Shiny time series
#   zlim           optional colour-limit overrides for maps, by direction and
#                  map name: zlim$onset$chpt = list(mean = c(70, 280), sd = NULL)
#   shiny          defaults for the Shiny app

india_map_regions <- c("India", "Pakistan", "China", "Nepal",
                       "Bhutan", "Bangladesh", "Myanmar", "Sri Lanka")

ethiopia_map_regions <- c("Ethiopia", "Eritrea", "Djibouti", "Somalia",
                          "Somaliland", "Kenya", "South Sudan", "Sudan")

regions <- list(

  india = list(
    label       = "India (IMD 0.25\u00b0)",
    data_dir    = "../data/india_onset/IMD_25by25",
    file_fmt    = "data_%d.nc",
    var_name    = "RAINFALL",
    dim_names   = c(lon = "LONGITUDE", lat = "LATITUDE", time = "TIME"),
    map_regions = india_map_regions,
    shapefile   = NULL,
    gold_standard = list(
      o2_any = list(
        file      = "../O2_onset_full_timeseries_025x025_withoutMOK.nc4",
        var_name  = "O2_onset",
        dim_names = c(lon = "lon", lat = "lat", year = "year"),
        direction = "onset",
        origin    = "04-01",
        label     = "O2 onset (unrestricted)",
        overlay   = TRUE
      ),
      o2_postMOK = list(
        file      = "../O2_onset_full_timeseries_025x025.nc4",
        var_name  = "O2_onset",
        dim_names = c(lon = "lon", lat = "lat", year = "year"),
        direction = "onset",
        origin    = "04-01",
        label     = "O2 onset (post MOK)",
        overlay   = FALSE
      )
    ),
    zlim = list(
      onset   = list(chpt = list(mean = c(0,240), sd = NULL),
                     dev_expl = list(mean = c(0, 1), sd = NULL)),
      retreat = list(chpt = list(mean = c(210, 365), sd = NULL),
                     dev_expl = list(mean = c(0, 1), sd = NULL))
    ),
    shiny = list(default_lon = 77, default_lat = 25,
                 default_start = 2000, n_years_show = 25)
  ),

  ethiopia = list(
    label       = "Ethiopia (CHIRPS 0.25\u00b0)",
    data_dir    = "../data/ethiopia_onset/CHIRPS_IMERG",
    file_fmt    = "%d.nc",
    var_name    = "RAINFALL",
    dim_names   = c(lon = "LONGITUDE", lat = "LATITUDE", time = "TIME"),
    map_regions = ethiopia_map_regions,
    shapefile   = "../manual_zones_woredas.shp",
    gold_standard = NULL,
    zlim = list(onset = list(chpt = list(mean = c(0, 180), sd = NULL)),
                retreat = list()),
    shiny = list(default_lon = 38.5, default_lat = 9,
                 default_start = 2000, n_years_show = 25)
  )
)

# Model and analysis defaults (can be overridden on the command line).
analysis_defaults <- list(
  gamma_floor = 1,    # precip below this (mm) is raised to it before fitting
  min_seg     = 5,    # minimum days on each side of a changepoint
  spline_k    = 20    # basis size of the spline used for dev_full
)
