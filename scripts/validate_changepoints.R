# validate_changepoints.R
#
# Check the closed-form changepoint calculation (R/changepoints.R) against
# explicit glm() fits, on a random sample of cell-years.
#
# For each sampled cell-year and every candidate changepoint, fits
#   glm(y ~ seg - 1, family = Gamma(link = "log"))
# on the floored window, applies the same rules (minimum segment length,
# direction), and compares the chosen day and deviance explained (relative to
# glm(y ~ 1)) with changepoint_series().
#
# Optionally (--old=<file>) also compares the new changepoints for the run's
# years with an older results file (e.g. ../india_changepoints_doy1-240.Rda),
# reading the new ones from this run's changepoints.Rda.
#
# Usage:
#   Rscript refactor/scripts/validate_changepoints.R --region=india \
#       --years=2000:2005 --doy=1:240 --direction=onset [--n=200] [--seed=1] \
#       [--gamma_floor=1] [--min_seg=5] [--old=../india_changepoints_doy1-240.Rda]

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f) == 1) dirname(normalizePath(f)) else getwd()
  source(file.path(d, "..", "R", "setup.R"))
  onset_setup(file.path(d, ".."))
})

usage <- paste(
  "Usage: Rscript scripts/validate_changepoints.R --region=<region>",
  "[--years=<y1:y2>] --doy=<start:end> --direction=<onset|retreat> [--n=200]",
  "[--seed=1] [--gamma_floor=1] [--min_seg=5] [--old=<file>]")

args <- parse_cli()
spec <- run_spec_from_args(args, usage)
cfg  <- get_region(spec$region)
n_sample    <- if (!is.null(args$n)) as.integer(args$n) else 200L
gamma_floor <- if (!is.null(args$gamma_floor)) as.numeric(args$gamma_floor) else
  analysis_defaults$gamma_floor
min_seg <- if (!is.null(args$min_seg)) as.integer(args$min_seg) else
  analysis_defaults$min_seg
set.seed(if (!is.null(args$seed)) as.integer(args$seed) else 1L)

# Reference: one glm per candidate changepoint.
changepoint_glm <- function(y, gamma_floor, min_seg, direction) {
  y <- pmax(y, gamma_floor)
  n <- length(y)
  null_dev <- glm(y ~ 1, family = Gamma(link = "log"))$deviance
  cand <- min_seg:(n - min_seg)
  de <- rep(-Inf, length(cand))
  for (q in seq_along(cand)) {
    seg <- factor(rep(c("pre", "post"), c(cand[q], n - cand[q])),
                  levels = c("pre", "post"))
    fit <- glm(y ~ seg - 1, family = Gamma(link = "log"),
               control = glm.control(epsilon = 1e-12, maxit = 100))
    b <- coef(fit)
    ok <- if (direction == "onset") b[2] >= b[1] else b[2] <= b[1]
    if (ok) de[q] <- 1 - fit$deviance / null_dev
  }
  if (all(!is.finite(de))) return(c(chpt_rel = NA, dev_expl = 0))
  k <- which.max(de)
  # report near-ties so they can be told apart from real disagreements
  second <- max(de[-k])
  c(chpt_rel = cand[k], dev_expl = de[k], gap = de[k] - second)
}

data <- read_region(cfg, spec$years, spec$doy_start, spec$doy_end)
d <- dim(data$rain)
has_data <- which(!apply(is.na(data$rain), c(1, 2, 4), any), arr.ind = TRUE)
if (!nrow(has_data)) stop("No complete cell-years in the data.")
pick <- has_data[sample(nrow(has_data), min(n_sample, nrow(has_data))), , drop = FALSE]

msg("validating", nrow(pick), "cell-years against glm()")
rows <- lapply(seq_len(nrow(pick)), function(q) {
  i <- pick[q, 1]; j <- pick[q, 2]; k <- pick[q, 3]
  y <- data$rain[i, j, , k]
  a <- changepoint_series(y, gamma_floor, min_seg, spec$direction)
  b <- changepoint_glm(y, gamma_floor, min_seg, spec$direction)
  data.frame(lon = data$lon[i], lat = data$lat[j], year = data$years[k],
             chpt_direct = a[["chpt_rel"]], chpt_glm = b[["chpt_rel"]],
             dev_direct = a[["dev_expl"]], dev_glm = b[["dev_expl"]],
             gap = if ("gap" %in% names(b)) b[["gap"]] else NA)
})
v <- do.call(rbind, rows)

same_chpt <- (is.na(v$chpt_direct) & is.na(v$chpt_glm)) |
  (!is.na(v$chpt_direct) & !is.na(v$chpt_glm) & v$chpt_direct == v$chpt_glm)
dev_diff <- abs(v$dev_direct - v$dev_glm)
near_tie <- !is.na(v$gap) & v$gap < 1e-6

cat("\n=== Direct calculation vs glm() ===\n")
cat("cell-years checked          :", nrow(v), "\n")
cat("same changepoint            :", sum(same_chpt), "\n")
cat("different, at near-ties     :", sum(!same_chpt & near_tie), "\n")
cat("different, not near-ties    :", sum(!same_chpt & !near_tie), "\n")
cat("max |dev_expl difference|   :", signif(max(dev_diff, na.rm = TRUE), 3), "\n")
if (any(!same_chpt & !near_tie)) {
  cat("\nDisagreements:\n")
  print(v[!same_chpt & !near_tie, ])
}

# -- Optional comparison with an older results file --------------------------
if (!is.null(args$old)) {
  old <- load_rda_object(resolve_path(args$old), "results")
  new <- load_rda_object(changepoints_path(spec$region, spec$tag), "results")
  yrs <- intersect(dimnames(old$chpt)$year, dimnames(new$chpt)$year)
  if (!length(yrs)) stop("No years in common with ", args$old)
  oc <- old$chpt[, , yrs, drop = FALSE]
  nc <- new$chpt[, , yrs, drop = FALSE]
  both <- !is.na(oc) & !is.na(nc)
  diff <- nc[both] - oc[both]
  near_edge <- pmin(oc[both] - spec$doy_start + 1, spec$doy_end - oc[both]) < min_seg
  cat("\n=== New vs old changepoints (", length(yrs), " years) ===\n", sep = "")
  cat("cells with both              :", sum(both), "\n")
  cat("identical                    :", sum(diff == 0), "\n")
  cat("different, old within", min_seg, "days of a window edge:",
      sum(diff != 0 & near_edge), "\n")
  cat("different, otherwise         :", sum(diff != 0 & !near_edge), "\n")
  cat("NA in old only / new only    :", sum(is.na(oc) & !is.na(nc)), "/",
      sum(!is.na(oc) & is.na(nc)), "\n")
  od <- old$dev_expl[, , yrs, drop = FALSE]
  nd <- new$dev_expl[, , yrs, drop = FALSE]
  same <- both & nc == oc
  cat("max |dev_expl diff| where the changepoint is identical:",
      signif(max(abs(nd[same] - od[same]), na.rm = TRUE), 3), "\n")
}
