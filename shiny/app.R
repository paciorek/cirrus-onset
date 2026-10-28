# app.R -- Rainy-season onset/retreat explorer
#
# Adapted copy of scratch/shiny/app.R that reads the pipeline's
# outputs.  A dropdown lists every run found under output/
# (region, years, window, direction).
#
# Left panel: a map of a chosen statistic summarized over the selected years.
# Click the map to select a grid cell.
# Right panel: one panel per year of daily precipitation at that cell, with
# the changepoint and spline fits overlaid, and the gold standard date when
# the region has one for the run's direction.
#
# Run from R with
#   shiny::runApp("./shiny")
# Set ONSET_DIR if not in the root directory.

library(shiny)

# Shiny renders plots with ragg if it is installed.  A ragg binary built for
# an older R fails with "Graphics API version mismatch"; in that case fall
# back to base png() (and reinstall ragg to fix it properly).
local({
  ragg_ok <- requireNamespace("ragg", quietly = TRUE) && tryCatch({
    f <- tempfile(fileext = ".png")
    ragg::agg_png(f); plot.new(); dev.off(); unlink(f)
    TRUE
  }, error = function(e) {
    message("ragg graphics device failed (", conditionMessage(e),
            "); using base png() instead. ",
            "Fix with install.packages(\"ragg\").")
    FALSE
  })
  if (!ragg_ok) options(shiny.useragg = FALSE)
})

# Shiny auto-sources R/*.R for runApp(); source explicitly too so the app
# also works via shiny::shinyAppFile() or source("app.R").
if (!exists("read_cell_series")) {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
}

runs <- list_runs()
if (!nrow(runs))
  stop("No runs found under ", output_root(),
       ". Run scripts/estimate_changepoints.R first.")
run_choices <- setNames(runs$tag, runs$label)
message("Found ", nrow(runs), " run(s) in ", output_root())

# Shapefile overlays, read once per region on first use
overlay_cache <- new.env()
region_overlay <- function(cfg) {
  if (!exists(cfg$key, envir = overlay_cache, inherits = FALSE))
    assign(cfg$key, load_overlay(cfg), envir = overlay_cache)
  get(cfg$key, envir = overlay_cache)
}


# -- UI ----------------------------------------------------------------------

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .cell-info { font-family: monospace; font-size: 0.9em; margin-top: 6px; }
    .help-note { color: #666; font-size: 0.85em; }
  "))),
  titlePanel("Rainy-season onset/retreat explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 5,
      selectInput("run", "Run", choices = run_choices,
                  selected = run_choices[1], width = "100%"),
      fluidRow(
        column(6,
          numericInput("start_year", "Start year", value = NA, step = 1),
          div(class = "help-note", textOutput("years_text", inline = TRUE))),
        column(6, radioButtons("map_summary", "Map summary over years",
                               choices = c("Mean" = "mean", "Std dev" = "sd"),
                               selected = "mean", inline = TRUE))
      ),
      selectInput("map_stat", "Map statistic", choices = NULL, width = "100%"),
      plotOutput("map", height = "460px", click = "map_click"),
      div(class = "help-note",
          "Click the map, or enter a latitude/longitude below, to select a ",
          "grid cell (snaps to the nearest cell with data)."),
      div(class = "cell-info", textOutput("cell_info")),
      fluidRow(
        column(4, numericInput("lat_in", "Latitude", value = NA, step = 0.25)),
        column(4, numericInput("lon_in", "Longitude", value = NA, step = 0.25)),
        column(4, div(style = "margin-top: 25px;",
                      actionButton("go_latlon", "Go to lat/lon")))
      ),
      div(class = "help-note", textOutput("data_status")),
      hr(),
      fluidRow(
        column(4, numericInput("y_max", "Y-axis max (mm)", value = 50,
                               min = 1, step = 10)),
        column(4, checkboxInput("free_y", "Free y-axis", value = FALSE)),
        column(4, uiOutput("gold_ui"))
      ),
      downloadButton("download_ts", "Download time-series PDF")
    ),
    mainPanel(
      width = 7,
      plotOutput("ts", height = "820px"),
      div(class = "help-note", textOutput("ts_legend"))
    )
  )
)


# -- Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # -- The selected run --
  run <- reactive({
    req(input$run)
    row <- runs[runs$tag == input$run, ]
    withProgress(message = "Loading run...", {
      r <- load_run(row$region, row$tag)
    })
    cfg <- r$cfg
    dir <- r$meta$direction

    # Gold standard dataset to overlay on the time series (if any)
    gs <- gold_standards_for(cfg, dir)
    ov <- names(Filter(function(e) isTRUE(e$overlay), gs))
    r$gold_overlay <- if (length(ov) && paste0("gold_", ov[1]) %in% names(r$r))
      paste0("gold_", ov[1]) else NULL
    r$gold_label <- if (!is.null(r$gold_overlay)) gs[[ov[1]]]$label else NULL

    r$menu    <- specs_for_run(r$r, cfg, dir, r$meta$doy_start, r$meta$doy_end)
    r$n_show  <- min(cfg$shiny$n_years_show, length(r$years))
    r$start_min <- min(r$years)
    r$start_max <- r$years[length(r$years) - r$n_show + 1]
    r$spline_k <- {
      k <- r$stats_meta$metrics$dev_full$args$spline_k
      if (is.null(k)) analysis_defaults$spline_k else k
    }
    r$data_years <- available_data_years(cfg, r$years)
    r
  })

  cell <- reactiveVal(NULL)

  # When the run changes, reset the inputs that depend on it
  observeEvent(run(), {
    r <- run()
    cfg <- r$cfg
    start <- min(max(cfg$shiny$default_start, r$start_min), r$start_max)
    updateNumericInput(session, "start_year",
                       label = sprintf("Start year (%d-%d)", r$start_min, r$start_max),
                       value = start, min = r$start_min, max = r$start_max)
    sel <- if (isolate(input$map_stat) %in% names(r$menu)) isolate(input$map_stat) else "chpt"
    updateSelectInput(session, "map_stat",
                      choices = setNames(names(r$menu),
                                         vapply(r$menu, `[[`, "", "label")),
                      selected = sel)
    c0 <- nearest_cell(cfg$shiny$default_lon, cfg$shiny$default_lat,
                       r$lon, r$lat, valid = r$land)
    cell(c0)
    updateNumericInput(session, "lat_in", value = c0$lat,
                       min = min(r$lat), max = max(r$lat))
    updateNumericInput(session, "lon_in", value = c0$lon,
                       min = min(r$lon), max = max(r$lon))
  }, priority = 10)

  output$gold_ui <- renderUI({
    r <- run()
    if (is.null(r$gold_overlay)) return(NULL)
    checkboxInput("show_gold", paste("Show", r$gold_label), value = TRUE)
  })

  output$data_status <- renderText({
    r <- run()
    dd <- resolve_path(r$cfg$data_dir)
    if (length(r$data_years))
      sprintf("Precip data: %s (%d of %d years present)", dd,
              length(r$data_years), length(r$years)) else
      sprintf("Precip data: no yearly files found in %s", dd)
  })

  # Start year: debounced so typing "19.." doesn't trigger intermediate
  # redraws.  Out-of-range values are clamped.
  start_year <- debounce(reactive({
    r <- run()
    s <- suppressWarnings(as.integer(round(input$start_year)))
    if (length(s) != 1 || is.na(s)) s <- r$cfg$shiny$default_start
    min(max(s, r$start_min), r$start_max)
  }), 600)

  years_sel <- reactive({
    r <- run()
    yrs <- r$years[r$years >= start_year()]
    head(yrs, r$n_show)
  })

  output$years_text <- renderText(
    sprintf("Showing %d-%d (%d years)", min(years_sel()),
            max(years_sel()), length(years_sel())))
  year_idx <- reactive(match(years_sel(), run()$years))

  observeEvent(input$map_click, {
    r <- run()
    x <- input$map_click$x
    y <- input$map_click$y
    if (is.null(x) || is.null(y)) return()
    # Ignore clicks outside the grid (e.g. on the colour legend)
    if (x < min(r$lon) - 0.5 || x > max(r$lon) + 0.5 ||
        y < min(r$lat) - 0.5 || y > max(r$lat) + 0.5) return()
    select_location(x, y)
  })

  # Typed lat/lon override the map selection
  observeEvent(input$go_latlon, {
    r <- run()
    lat0 <- input$lat_in
    lon0 <- input$lon_in
    if (!is.numeric(lat0) || !is.numeric(lon0) || is.na(lat0) || is.na(lon0)) {
      showNotification("Enter a numeric latitude and longitude.", type = "error")
      return()
    }
    if (lon0 < min(r$lon) || lon0 > max(r$lon) ||
        lat0 < min(r$lat) || lat0 > max(r$lat)) {
      showNotification(sprintf("Location outside the grid (lon %g-%g, lat %g-%g).",
                               min(r$lon), max(r$lon), min(r$lat), max(r$lat)),
                       type = "error")
      return()
    }
    select_location(lon0, lat0, update_inputs = FALSE)
  })

  # Select the nearest grid cell with data; optionally sync the lat/lon boxes
  select_location <- function(lon0, lat0, update_inputs = TRUE) {
    r <- run()
    c0 <- nearest_cell(lon0, lat0, r$lon, r$lat, valid = r$land)
    cell(c0)
    if (update_inputs) {
      updateNumericInput(session, "lat_in", value = c0$lat)
      updateNumericInput(session, "lon_in", value = c0$lon)
    }
  }

  # -- Map --
  map_field <- reactive({
    r <- run()
    req(input$map_stat %in% names(r$menu))
    s <- r$menu[[input$map_stat]]
    arr <- s$fun(r$r)[, , year_idx(), drop = FALSE]
    list(z     = summarize_years(arr, input$map_summary),
         zlim  = s$zlim[[input$map_summary]],
         label = s$label)
  })

  output$map <- renderPlot({
    r  <- run()
    mf <- map_field()
    c0 <- cell()
    req(c0)
    yrs <- years_sel()
    par(mai = c(0.4, 0.4, 0.5, 0.2), mgp = c(1, 0.5, 0), cex = 0.85,
        cex.main = 1.0)
    plot_map_panel(
      mf$z, r$lon, r$lat,
      main = sprintf("%s\n%s, %d-%d", mf$label,
                     if (input$map_summary == "mean") "mean" else "std dev",
                     min(yrs), max(yrs)),
      zlim = mf$zlim, map_regions = r$cfg$map_regions,
      overlay = region_overlay(r$cfg), mark = c(c0$lon, c0$lat))
  })

  output$cell_info <- renderText({
    c0 <- cell()
    req(c0)
    mf <- map_field()
    sprintf("Grid cell: lon %.3f, lat %.3f  |  map value: %.3g",
            c0$lon, c0$lat, mf$z[c0$i_lon, c0$i_lat])
  })

  # -- Time series --
  cell_rain <- reactive({
    c0 <- cell()
    req(c0)
    read_cell_series(run()$cfg, years_sel(), c0$lon, c0$lat)
  })

  ts_args <- reactive({
    r   <- run()
    c0  <- cell()
    req(c0)
    yi  <- year_idx()
    gold <- if (!is.null(r$gold_overlay) && isTRUE(input$show_gold))
      r$r[[r$gold_overlay]][c0$i_lon, c0$i_lat, yi] else NULL
    list(rain        = cell_rain(),
         years       = years_sel(),
         chpt        = r$r$chpt[c0$i_lon, c0$i_lat, yi],
         dev_expl    = r$r$dev_expl[c0$i_lon, c0$i_lat, yi],
         dev_full    = if (!is.null(r$r$dev_full))
                         r$r$dev_full[c0$i_lon, c0$i_lat, yi] else NULL,
         doy_start   = r$meta$doy_start,
         doy_end     = r$meta$doy_end,
         gold        = gold,
         y_max       = if (isTRUE(input$free_y)) NULL else input$y_max,
         spline_k    = r$spline_k,
         gamma_floor = r$meta$gamma_floor,
         # Outer title only for the PDF download (the app shows the location
         # in the left panel)
         main        = sprintf("%s - lat %.3f, lon %.3f  |  doy %d-%d, %s",
                               r$cfg$label, c0$lat, c0$lon, r$meta$doy_start,
                               r$meta$doy_end, r$meta$direction))
  })

  output$ts <- renderPlot({
    args <- ts_args()
    validate(need(any(!is.na(args$rain)),
                  paste("No precipitation data files found for these years in",
                        resolve_path(run()$cfg$data_dir))))
    args$main <- NULL
    do.call(plot_location_panels, c(args, cex = 0.85))
  })

  output$ts_legend <- renderText({
    r <- run()
    sprintf(paste0("Gray: daily precip; shaded: doy %d-%d analysis window; ",
                   "red: spline fit; blue: changepoint fit (%s); dashed blue: ",
                   "changepoint DOY%s. Titles: year  changepoint dev.expl%s."),
            r$meta$doy_start, r$meta$doy_end, r$meta$direction,
            if (!is.null(r$gold_overlay) && isTRUE(input$show_gold))
              paste0("; green: ", r$gold_label) else "",
            if (!is.null(r$r$dev_full)) " : spline dev.expl" else "")
  })

  output$download_ts <- downloadHandler(
    filename = function() {
      c0 <- cell()
      sprintf("%s_lon%.2f_lat%.2f_%d-%d.pdf", run()$meta$tag, c0$lon, c0$lat,
              min(years_sel()), max(years_sel()))
    },
    content = function(file) {
      pdf(file, width = 11, height = 9)
      on.exit(dev.off())
      do.call(plot_location_panels, ts_args())
    }
  )
}

shinyApp(ui, server)
