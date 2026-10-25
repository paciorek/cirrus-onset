# changepoints.R
#
# Two-segment Gamma changepoint model, computed in closed form.
#
# Model: within the window (n days), values below gamma_floor are raised to
# gamma_floor.  A changepoint c splits the window into a pre segment (days
# 1..c) and a post segment (days c+1..n), each with its own Gamma mean.  The
# maximum-likelihood means are the segment sample means, and the Gamma
# deviance reduces to
#
#   D(c) = 2 * [ c log(mean_pre) + (n - c) log(mean_post) - sum(log y) ]
#
# (the (y - mu)/mu term of the Gamma deviance sums to zero within each
# segment).  The null deviance D0 uses the mean of the whole window, and
# deviance explained is 1 - D(c)/D0.  This is the same model as
# glm(y ~ seg - 1, family = Gamma(link = "log")) for each candidate c, so the
# chosen changepoint is the same; see scripts/validate_changepoints.R.
#
# D0 is the deviance of the intercept-only model, as in mgcv::gam(), so
# dev_expl matches the old gam()-based values and is on the same scale as
# dev_full.  (Plain glm() with a no-intercept formula would instead report a
# null deviance relative to mu = 1; the validation script uses glm(y ~ 1).)

# One window series.
#
# y           : numeric vector (the analysis window)
# gamma_floor : values below this are raised to it
# min_seg     : minimum number of days in each segment
# direction   : "onset" (post mean >= pre mean) or "retreat" (post <= pre)
#
# Returns c(chpt_rel, dev_expl):
#   chpt_rel  last day of the pre segment, as an index into y (NA if none)
#   dev_expl  deviance explained at that changepoint
# Any missing value in y, or a window with no variation after flooring,
# gives NA for both.  If no candidate has the right direction, chpt_rel is
# NA and dev_expl is 0.
changepoint_series <- function(y, gamma_floor = 1, min_seg = 5,
                               direction = "onset") {
  out <- c(chpt_rel = NA_real_, dev_expl = NA_real_)
  y <- as.numeric(y)              # drop names (they would leak into results)
  n <- length(y)
  if (n < 2 * min_seg || anyNA(y) || any(!is.finite(y))) return(out)

  y   <- pmax(y, gamma_floor)
  sly <- sum(log(y))
  d0  <- 2 * (n * log(mean(y)) - sly)
  if (!(d0 > 1e-10)) return(out)

  cs   <- cumsum(y)
  cand <- min_seg:(n - min_seg)
  m1   <- cs[cand] / cand
  m2   <- (cs[n] - cs[cand]) / (n - cand)
  dev  <- 2 * (cand * log(m1) + (n - cand) * log(m2) - sly)
  de   <- 1 - dev / d0

  ok <- if (direction == "onset") m2 >= m1 else m2 <= m1
  if (!any(ok)) {
    out["dev_expl"] <- 0
    return(out)
  }
  de[!ok] <- -Inf
  k <- which.max(de)              # ties go to the earliest day
  c(chpt_rel = cand[k], dev_expl = de[[k]])
}

# Fitted values of the changepoint model (segment means of the floored
# series), used for plotting.
changepoint_fitted <- function(y, chpt_rel, gamma_floor = 1) {
  y <- pmax(y, gamma_floor)
  n <- length(y)
  c(rep(mean(y[1:chpt_rel]), chpt_rel),
    rep(mean(y[(chpt_rel + 1):n]), n - chpt_rel))
}

# All cells and years.
#
# data : list from read_region()
#
# Returns list(chpt, dev_expl): lon x lat x year arrays; chpt is absolute day
# of year (last day of the pre segment).
estimate_changepoints_grid <- function(data, gamma_floor = 1, min_seg = 5,
                                       direction = "onset", verbose = TRUE) {
  d <- dim(data$rain)
  n_lon <- d[1]; n_lat <- d[2]; n_year <- d[4]
  doy_start <- data$doy[1]

  dn <- dimnames(data$rain)[c("lon", "lat", "year")]
  chpt     <- array(NA_real_, c(n_lon, n_lat, n_year), dimnames = dn)
  dev_expl <- array(NA_real_, c(n_lon, n_lat, n_year), dimnames = dn)

  for (i in seq_len(n_lon)) {
    for (j in seq_len(n_lat)) {
      for (k in seq_len(n_year)) {
        y <- data$rain[i, j, , k]
        if (anyNA(y)) next            # includes ocean / outside-mask cells
        r <- changepoint_series(y, gamma_floor, min_seg, direction)
        chpt[i, j, k]     <- r[["chpt_rel"]] + doy_start - 1
        dev_expl[i, j, k] <- r[["dev_expl"]]
      }
    }
    if (verbose && (i %% 10 == 0 || i == n_lon))
      msg("changepoints: lon", i, "of", n_lon)
  }
  list(chpt = chpt, dev_expl = dev_expl)
}
