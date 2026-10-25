# metrics.R
#
# Metric functions for compute_stats, and the registry that names them.
#
# A metric function has the signature  f(y, ctx, ...)  where
#   y    the raw precipitation window for one cell-year (not floored)
#   ctx  list(chpt, chpt_doy, dev_expl, doy, direction, gamma_floor,
#             year, lon, lat)
#        chpt is the changepoint as an index into y (last day of the pre
#        segment); chpt_doy is the same as day of year.
#   ...  the metric's own settings
# and returns one number, or a named numeric vector of related quantities.
#
# The driver never calls a metric on a window with missing values (the result
# is NA), and, for metrics with needs_chpt = TRUE, never when the changepoint
# is NA.
#
# Metric functions must be self-contained (no calls to helper functions
# defined in this project), because with a "multisession" parallel plan they
# are sent to worker processes without the rest of the code.  Package
# functions are fine if the package is listed in `packages`.
#
# Registry entries:
#   fun         the function
#   args        default settings passed as ...
#   outputs     names of the values returned, in order (one array each)
#   needs_chpt  skip (NA) when the changepoint is NA
#   packages    packages the function uses

# Deviance explained by a smooth spline of time, fit to the floored series so
# it is on the same scale as the changepoint model's dev_expl.
metric_dev_full <- function(y, ctx, spline_k = 20) {
  y <- pmax(y, ctx$gamma_floor)
  t <- seq_along(y)
  fit <- mgcv::gam(y ~ s(t, k = spline_k), family = Gamma(link = "log"))
  (fit$null.deviance - fit$deviance) / fit$null.deviance
}

# Means, totals and SDs of daily precip before and after the changepoint
# (raw, unfloored values).
metric_pre_post <- function(y, ctx) {
  pre  <- y[seq_len(ctx$chpt)]
  post <- y[(ctx$chpt + 1):length(y)]
  c(pre_mean  = mean(pre),  post_mean  = mean(post),
    pre_total = sum(pre),   post_total = sum(post),
    pre_sd    = sd(pre),    post_sd    = sd(post))
}

# Transition sharpness: sqrt(var_pre + var_post) of the `window`-day running
# mean, using running-mean windows starting [pre_start, pre_end] days before
# and [post_start, post_end] days after the changepoint.
# (Ported from compute_onset_running_mean_sd().)
metric_transition_sd <- function(y, ctx, window = 7,
                                 pre_start = 36, pre_end = 6,
                                 post_start = 1, post_end = 30) {
  n <- length(y)
  n_windows <- n - window + 1L
  if (n_windows < 2) return(NA_real_)
  cs <- cumsum(c(0, y))
  rm <- (cs[(window + 1):(n + 1)] - cs[seq_len(n_windows)]) / window
  pre_idx  <- seq(ctx$chpt - pre_start, ctx$chpt - pre_end)
  pre_idx  <- pre_idx[pre_idx >= 1 & pre_idx <= n_windows]
  post_idx <- seq(ctx$chpt + post_start, ctx$chpt + post_end)
  post_idx <- post_idx[post_idx >= 1 & post_idx <= n_windows]
  if (length(pre_idx) < 2 || length(post_idx) < 2) return(NA_real_)
  sqrt(var(rm[pre_idx]) + var(rm[post_idx]))
}

# Transition complexity: EDF of a GAM smooth over a `window_len`-day window
# centred on the changepoint, with a pre/post step term.
# (Ported from compute_onset_smooth_edf().)
metric_transition_edf <- function(y, ctx, window_len = 60, smooth_k = 30) {
  n <- length(y)
  half_before <- floor((window_len - 1) / 2)
  half_after  <- window_len - half_before - 1L
  idx <- seq(ctx$chpt - half_before, ctx$chpt + half_after)
  if (min(idx) < 1 || max(idx) > n) return(NA_real_)
  yy   <- pmax(y[idx], ctx$gamma_floor)
  day  <- seq_along(idx)
  post <- as.numeric(idx > ctx$chpt)
  fit <- mgcv::gam(yy ~ post + s(day, k = smooth_k), family = Gamma(link = "log"))
  s_tab <- mgcv::summary.gam(fit)$s.table
  as.numeric(s_tab[1L, "edf"])
}

metric_registry <- list(
  dev_full = list(
    fun = metric_dev_full, args = list(spline_k = 20),
    outputs = "dev_full", needs_chpt = FALSE, packages = "mgcv"),
  pre_post = list(
    fun = metric_pre_post, args = list(),
    outputs = c("pre_mean", "post_mean", "pre_total", "post_total",
                "pre_sd", "post_sd"),
    needs_chpt = TRUE, packages = character(0)),
  transition_sd = list(
    fun = metric_transition_sd,
    args = list(window = 7, pre_start = 36, pre_end = 6,
                post_start = 1, post_end = 30),
    outputs = "transition_sd", needs_chpt = TRUE, packages = character(0)),
  transition_edf = list(
    fun = metric_transition_edf, args = list(window_len = 60, smooth_k = 30),
    outputs = "transition_edf", needs_chpt = TRUE, packages = "mgcv")
)

# Turn a user-supplied list (registry names, registry-style entries, or bare
# functions) into a list of complete registry entries.
#
#   normalize_metrics(c("dev_full", "pre_post"))
#   normalize_metrics(list(my_stat = function(y, ctx) max(y)))
#   normalize_metrics(list(my_pair = list(fun = f, outputs = c("a", "b"))))
normalize_metrics <- function(metrics) {
  if (is.character(metrics)) {
    unknown <- setdiff(metrics, names(metric_registry))
    if (length(unknown))
      stop("Unknown metric(s): ", paste(unknown, collapse = ", "),
           ". Available: ", paste(names(metric_registry), collapse = ", "))
    return(metric_registry[metrics])
  }
  if (is.null(names(metrics)) || any(names(metrics) == ""))
    stop("Metrics must be a named list.")
  out <- Map(function(m, nm) {
    if (is.function(m)) m <- list(fun = m)
    if (!is.function(m$fun)) stop("Metric '", nm, "' has no function.")
    if (is.null(m$args))       m$args <- list()
    if (is.null(m$outputs))    m$outputs <- nm
    if (is.null(m$needs_chpt)) m$needs_chpt <- TRUE
    if (is.null(m$packages))   m$packages <- character(0)
    m
  }, metrics, names(metrics))
  all_out <- unlist(lapply(out, `[[`, "outputs"))
  if (anyDuplicated(all_out))
    stop("Duplicate output names: ", paste(unique(all_out[duplicated(all_out)]),
                                           collapse = ", "))
  out
}

# Override registry defaults from the command line, e.g.
# --metric_args=dev_full.spline_k=30,transition_sd.window=5
apply_metric_args <- function(metrics, spec_string) {
  if (is.null(spec_string) || !nzchar(spec_string)) return(metrics)
  for (item in strsplit(spec_string, ",", fixed = TRUE)[[1]]) {
    kv  <- strsplit(item, "=", fixed = TRUE)[[1]]
    key <- strsplit(kv[1], ".", fixed = TRUE)[[1]]
    if (length(kv) != 2 || length(key) != 2 || !key[1] %in% names(metrics))
      stop("Bad --metric_args entry: ", item)
    val <- suppressWarnings(as.numeric(kv[2]))
    metrics[[key[1]]]$args[[key[2]]] <- if (is.na(val)) kv[2] else val
  }
  metrics
}
