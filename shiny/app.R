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


# Smoothed dev_expl fits, shared by all sessions of this app process and kept
# for the life of the process.  Keyed by run, family, k and years; holds the
# fitted model (when fitted here) and the fitted values.  At most
# smooth_cache_max fits are kept in memory; fits are also saved to disk by
# fit_or_load_smooth(), so evicted or earlier fits reload quickly.
smooth_cache <- new.env()
smooth_cache_order <- character(0)
smooth_cache_max <- 6
smooth_cache_put <- function(id, sm) {
  assign(id, sm, envir = smooth_cache)
  smooth_cache_order <<- c(setdiff(smooth_cache_order, id), id)
  while (length(smooth_cache_order) > smooth_cache_max) {
    rm(list = smooth_cache_order[1], envir = smooth_cache)
    smooth_cache_order <<- smooth_cache_order[-1]
  }
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
      uiOutput("data_status_ui"),
      wellPanel(
        style = "padding: 8px 12px; margin-top: 12px;",
        strong("Mean dev_expl: smoothing and region mask"),
        div(class = "help-note",
            "Yearly dev_expl is averaged over the displayed years at each cell, ",
            "and that field is smoothed with s(lon, lat, k). Show it with the map ",
            "statistic 'Smoothed mean dev_expl' (above), or use it to mask the map."),
        fluidRow(
          column(4, numericInput("smooth_k", "k (basis size)", value = 500,
                                 min = 10, step = 50)),
          column(4, selectInput("smooth_family", "GAM family",
                                choices = names(smooth_families),
                                selected = "gaussian")),
          column(4, div(style = "margin-top: 25px;",
                        actionButton("fit_smooth", "Fit smooth")))
        ),
        div(class = "help-note", textOutput("smooth_status")),
        checkboxInput("thresh_on", "Mask cells with low mean dev_expl",
                      value = FALSE),
        conditionalPanel(
          "input.thresh_on",
          fluidRow(
            column(6, radioButtons("thresh_src", "Threshold applied to",
                                   choices = c("Smoothed" = "smooth",
                                               "Raw" = "raw"),
                                   selected = "smooth", inline = TRUE)),
            column(6, numericInput("thresh", "Threshold",
                                   value = 0.3, min = 0, max = 1, step = 0.05))
          ),
          div(class = "help-note", textOutput("mask_status")),
          div(style = "margin-top: 6px;",
              actionButton("save_mask", "Save mask"),
              downloadButton("download_mask", "Download mask (CSV)"))
        )
      )
    ),
    mainPanel(
      width = 7,
      plotOutput("ts", height = "820px"),
      div(class = "help-note", textOutput("ts_legend")),
      hr(),
      fluidRow(
        column(3, numericInput("y_max", "Y-axis max (mm)", value = 50,
                               min = 1, step = 10)),
        column(2, div(style = "margin-top: 25px;",
                      checkboxInput("free_y", "Free y-axis", value = FALSE))),
        column(4, uiOutput("gold_ui")),
        column(3, div(style = "margin-top: 25px;",
                      downloadButton("download_ts", "Download time-series PDF")))
      )
    )
  )
)


# Map statistic menu: the run's map specs plus the smoothed mean dev_expl
# (a single field, not a per-year array), placed after dev_expl.
map_menu_choices <- function(menu) {
  ch <- setNames(names(menu), vapply(menu, `[[`, "", "label"))
  i  <- match("dev_expl", ch)
  sm <- c("Smoothed mean dev_expl (see settings below)" = "dev_expl_smooth")
  if (is.na(i)) c(ch, sm) else c(ch[seq_len(i)], sm, ch[-seq_len(i)])
}


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
    # Gold standard datasets available for the time-series line: a named
    # vector label -> array name, and the default (overlay = TRUE, else first)
    avail <- names(gs)[paste0("gold_", names(gs)) %in% names(r$r)]
    r$gold_choices <- if (length(avail))
      setNames(paste0("gold_", avail), vapply(gs[avail], `[[`, "", "label")) else
      character(0)
    r$gold_default <- if (length(intersect(ov, avail)))
      paste0("gold_", intersect(ov, avail)[1]) else unname(r$gold_choices[1])

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
    choices <- map_menu_choices(r$menu)
    sel <- if (isolate(input$map_stat) %in% choices) isolate(input$map_stat) else "chpt"
    updateSelectInput(session, "map_stat", choices = choices, selected = sel)
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
    if (!length(r$gold_choices)) return(NULL)
    radioButtons("gold_sel", "Gold standard line",
                 choices = c("None" = "none", r$gold_choices),
                 selected = r$gold_default)
  })

  # The selected gold standard array name, or NULL
  gold_sel <- reactive({
    r <- run()
    g <- input$gold_sel
    if (is.null(g) || !g %in% r$gold_choices) NULL else g
  })

  output$data_status_ui <- renderUI({
    r <- run()
    missing <- length(r$data_years) < length(r$years)
    div(class = "help-note",
        style = if (missing) "color: #b00020; white-space: pre-line;" else NULL,
        textOutput("data_status"))
  })

  output$data_status <- renderText({
    r <- run()
    dd <- resolve_path(r$cfg$data_dir)
    if (!length(r$data_years)) return(no_data_message(r$cfg))
    msg <- sprintf("Precip data: %s (%d of %d years present)", dd,
                   length(r$data_years), length(r$years))
    if (length(r$data_years) < length(r$years))
      msg <- paste(msg, data_download_hint(r$cfg))
    msg
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

  # -- dev_expl threshold mask (raw or smoothed) --
  # Mean dev_expl over the displayed years, per cell.
  raw_dev_mean <- reactive(mean_dev_expl(run()$r$dev_expl, run()$years,
                                         years_sel()))

  # Settings that identify a smooth; debounced so typing k doesn't trigger
  # lookups for intermediate values.
  smooth_key <- debounce(reactive({
    r <- run()
    k <- suppressWarnings(as.integer(input$smooth_k))
    req(length(k) == 1, !is.na(k), k >= 10, input$smooth_family)
    yrs <- years_sel()
    list(id = paste(r$meta$tag, input$smooth_family, k,
                    smooth_years_label(yrs), sep = "|"),
         k = k, family = input$smooth_family, years = yrs)
  }), 600)

  # The id of the smooth the user asked to fit (via the button)
  fit_request <- reactiveVal(NULL)
  observeEvent(input$fit_smooth, fit_request(smooth_key()$id))
  smooth_error <- reactiveVal(NULL)

  # The smooth for the current settings: from the in-memory cache, else from
  # disk, else fitted -- but only once "Fit smooth" has been pressed for these
  # settings.  NULL if not available yet.  Changing the map statistic, the
  # summary or the threshold does not touch this.
  smooth_fit <- reactive({
    key <- smooth_key()
    if (exists(key$id, envir = smooth_cache, inherits = FALSE))
      return(get(key$id, envir = smooth_cache))
    r <- run()
    do_fit <- identical(fit_request(), key$id)
    sm <- tryCatch({
      if (do_fit) {
        withProgress(
          message = sprintf("Fitting smooth (%s, k = %d)...", key$family, key$k),
          detail = "large k or non-gaussian families can take minutes",
          fit_or_load_smooth(r, k = key$k, family = key$family, years = key$years))
      } else {
        fit_or_load_smooth(r, k = key$k, family = key$family, years = key$years,
                           fit = FALSE)
      }
    }, error = function(e) {
      # Remember the error for these settings and don't retry automatically
      smooth_error(list(id = key$id, msg = conditionMessage(e)))
      fit_request(NULL)
      NULL
    })
    if (!is.null(sm)) smooth_cache_put(key$id, sm)
    sm
  })

  output$smooth_status <- renderText({
    key <- smooth_key()
    sm  <- smooth_fit()
    yrs <- years_sel()
    err <- smooth_error()
    if (!is.null(err) && identical(err$id, key$id) && is.null(sm))
      return(paste("Smoothing failed:", err$msg))
    if (is.null(sm))
      return(sprintf(paste0("No smooth yet for %s, k = %d, years %d-%d: press ",
                            "'Fit smooth' (or precompute with ",
                            "scripts/smooth_dev_expl.R)."),
                     key$family, key$k, min(yrs), max(yrs)))
    sprintf("Smooth: %s, k = %d, edf %.0f, years %d-%d (%s).",
            sm$family, sm$k, sm$edf, min(yrs), max(yrs),
            if (identical(sm$source, "loaded")) "loaded from saved fit" else
              sprintf("fitted in %.0f s", sm$seconds))
  })

  # The region mask (see make_region_mask()), or NULL for no masking.
  region_mask <- reactive({
    if (!isTRUE(input$thresh_on)) return(NULL)
    thr <- input$thresh
    req(is.numeric(thr), !is.na(thr))
    r <- run()
    if (input$thresh_src == "raw")
      return(make_region_mask(r, thr, source = "raw", years = years_sel()))
    sm <- smooth_fit()
    if (is.null(sm)) return(NULL)
    make_region_mask(r, thr, source = "smooth", years = years_sel(), smooth = sm)
  })
  dev_mask <- reactive({
    m <- region_mask()
    if (is.null(m)) NULL else m$keep
  })

  output$mask_status <- renderText({
    req(isTRUE(input$thresh_on))
    m <- region_mask()
    if (is.null(m))
      return("Fit a smooth (above) to mask on smoothed values. Map not masked.")
    sprintf("%s mean dev_expl >= %g: %d of %d cells kept.",
            if (m$source == "raw") "Raw" else "Smoothed", m$threshold,
            sum(m$keep), sum(!is.na(m$values)))
  })

  observeEvent(input$save_mask, {
    m <- region_mask()
    if (is.null(m)) {
      showNotification("No mask to save yet.", type = "warning")
      return()
    }
    f <- tryCatch(save_region_mask(m), error = function(e) e)
    if (inherits(f, "error"))
      showNotification(paste("Could not save mask:", conditionMessage(f)),
                       type = "error")
    else
      showNotification(paste("Mask saved to", f), duration = 10)
  })

  output$download_mask <- downloadHandler(
    filename = function() {
      m <- region_mask()
      if (is.null(m)) "mask.csv" else
        sub("\\.rds$", ".csv", basename(region_mask_file(m)))
    },
    content = function(file) {
      m <- region_mask()
      validate(need(!is.null(m), "No mask to download yet."))
      write.csv(mask_cells(m), file, row.names = FALSE)
    }
  )

  # -- Map --
  map_field <- reactive({
    r <- run()
    if (identical(input$map_stat, "dev_expl_smooth")) {
      sm <- smooth_fit()
      validate(need(!is.null(sm), paste(
        "No smoothed dev_expl for these settings and years yet:",
        "press 'Fit smooth' (below the location and y-axis controls).")))
      return(list(z = sm$fitted, zlim = r$menu$dev_expl$zlim$mean,
                  label = sprintf("Smoothed mean dev_expl (%s, k = %d)",
                                  sm$family, sm$k),
                  single = TRUE))
    }
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
    z <- mf$z
    m <- dev_mask()
    if (!is.null(m)) z[!m] <- NA
    plot_map_panel(
      z, r$lon, r$lat,
      main = sprintf("%s\n%s, %d-%d%s", mf$label,
                     if (isTRUE(mf$single)) "mean over years" else
                       if (input$map_summary == "mean") "mean" else "std dev",
                     min(yrs), max(yrs),
                     if (is.null(m)) "" else sprintf(
                       " | %s dev_expl >= %g",
                       if (input$thresh_src == "raw") "raw" else "smoothed",
                       input$thresh)),
      zlim = mf$zlim, map_regions = r$cfg$map_regions,
      overlay = region_overlay(r$cfg), mark = c(c0$lon, c0$lat))
  })

  output$cell_info <- renderText({
    c0 <- cell()
    req(c0)
    mf <- map_field()
    txt <- sprintf("Grid cell: lon %.3f, lat %.3f  |  map value: %.3g",
                   c0$lon, c0$lat, mf$z[c0$i_lon, c0$i_lat])
    txt <- paste0(txt, sprintf("  |  mean dev_expl: %.3g",
                               raw_dev_mean()[c0$i_lon, c0$i_lat]))
    sm <- smooth_fit()
    if (!is.null(sm))
      txt <- paste0(txt, sprintf(", smoothed: %.3g", sm$fitted[c0$i_lon, c0$i_lat]))
    txt
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
    g <- gold_sel()
    gold <- if (!is.null(g)) r$r[[g]][c0$i_lon, c0$i_lat, yi] else NULL
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
    r <- run()
    validate(need(length(r$data_years) > 0, no_data_message(r$cfg)))
    missing <- setdiff(years_sel(), r$data_years)
    validate(need(any(!is.na(args$rain)), paste0(
      "No precipitation data files for the displayed years (",
      paste(range(missing), collapse = "-"), ") in ",
      resolve_path(r$cfg$data_dir), ".\n", data_download_hint(r$cfg))))
    args$main <- NULL
    do.call(plot_location_panels, c(args, cex = 0.85))
  })

  output$ts_legend <- renderText({
    r <- run()
    sprintf(paste0("Gray: daily precip; shaded: doy %d-%d analysis window; ",
                   "red: spline fit; blue: changepoint fit (%s); dashed blue: ",
                   "changepoint DOY%s. Titles: year, deviance explained."),
            r$meta$doy_start, r$meta$doy_end, r$meta$direction,
            if (!is.null(gold_sel()))
              paste0("; green: ", names(r$gold_choices)[r$gold_choices == gold_sel()])
            else "")
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
