#' Daily weather forcing for a simulation
#'
#' Checks a daily weather series and precomputes the temperature-driven model
#' terms for every day (model specification S1-S4, S7). The solvers hold these
#' values constant within each day.
#'
#' @param weather A data frame with one row per day and columns:
#'   \describe{
#'     \item{`date`}{Dates (class `Date`, or text in `"YYYY-MM-DD"` format),
#'       consecutive days without gaps or repeats.}
#'     \item{`tmean`}{Daily mean temperature (degrees Celsius). If absent,
#'       computed as `(tmax + tmin) / 2` from columns `tmax` and `tmin`.}
#'     \item{`rain`}{Optional daily rainfall (mm), used for the pupal
#'       moisture index.}
#'   }
#'   Gaps must be filled before use, e.g. with [sg_fill_weather()].
#' @param p An `sg_params` object. Only the temperature, development and
#'   mortality constants are used here; `sg_model()` (milestone M1, WP5) checks that the model is
#'   assembled with the same values.
#' @param start,end Optional first and last date of the simulation (`Date` or
#'   `"YYYY-MM-DD"`). Day 1 of a simulation is `start` (D-015). Defaults: the
#'   whole series.
#'
#' @details
#' For each day the forcing holds:
#' * `tmean` and `tmean5`, the running mean of the current and previous four
#'   days. Days before `start` in `weather` are used for the first values; at
#'   the start of the series the mean of the available days is used, as in the
#'   Gutierrez & Ponti code.
#' * `dx`: degree-days accumulated (S1), driving pupal and adult ageing.
#' * `dx_egg`: degree-days for eggs, equal to `dx` when
#'   `egg_temperature = "ambient"` and the constant `dx(T_host)` when
#'   `egg_temperature = "host"` (D-014).
#' * `mu_T`: pupal/adult temperature mortality rate per day (S2) and `p_A`,
#'   the corresponding daily proportion.
#' * `egg_viability` and `mu_E`, the egg hatch probability and mortality rate
#'   from the 5-day mean temperature (S3).
#' * `phi_ov`: the oviposition temperature scalar (S4).
#' * `moisture`: the moisture index I_W in \[0, 1\] (S7), equal to
#'   `min(1, annual rainfall / rain_threshold)` for the calendar year, with
#'   annual rainfall scaled up from incomplete years; 1 (no stress) when there
#'   is no `rain` column.
#'
#' @return An object of class `sg_forcing`: a list with `dates`, `day`
#'   (1, 2, ...) and the daily vectors above, plus `n_days`, `has_rain` and
#'   `param_values` (the parameter values used, so `sg_model()` can check
#'   consistency).
#' @seealso [sg_example_weather()], [sg_fill_weather()], [sg_cold_index()],
#'   [temperature_functions]
#' @export
#' @examples
#' p <- sg_params()
#' wx <- sg_forcing(sg_example_weather(), p, start = "1961-01-01", end = "1961-12-31")
#' wx
#' summary(wx$dx)
sg_forcing <- function(weather, p = sg_params(), start = NULL, end = NULL) {
  p <- .sg_as_params(p)
  w <- .sg_check_weather(weather)

  tmean5_all <- .sg_running_mean(w$tmean, 5L)
  keep <- .sg_date_window(w$date, start, end)
  w <- w[keep, , drop = FALSE]
  tmean5 <- tmean5_all[keep]

  dx <- sg_degree_days(w$tmean, p)
  dx_egg <- if (p$egg_temperature == "host") {
    rep(sg_degree_days(p$T_host, p), nrow(w))
  } else {
    dx
  }
  viab <- sg_egg_viability(tmean5, p)

  structure(
    list(
      dates = w$date,
      day = seq_len(nrow(w)),
      n_days = nrow(w),
      tmean = w$tmean,
      tmean5 = tmean5,
      dx = dx,
      dx_egg = dx_egg,
      mu_T = sg_temp_mortality(w$tmean, p, type = "hazard"),
      p_A = sg_temp_mortality(w$tmean, p, type = "proportion"),
      egg_viability = viab,
      mu_E = -log(pmax(p$eps_s, viab)),
      phi_ov = sg_ovip_scalar(w$tmean, p),
      has_rain = "rain" %in% names(w),
      rain = if ("rain" %in% names(w)) w$rain else NULL,
      moisture = .sg_moisture_index(w, p$rain_threshold),
      param_values = unclass(p)[.sg_forcing_param_names()]
    ),
    class = "sg_forcing"
  )
}

#' @export
print.sg_forcing <- function(x, ...) {
  cat(
    "<sg_forcing> ", x$n_days, " days, ", format(x$dates[1L]), " to ",
    format(x$dates[x$n_days]), "\n",
    sep = ""
  )
  rng <- function(v) {
    sprintf("mean %.2f, range [%.2f, %.2f]", mean(v), min(v), max(v))
  }
  cat("  daily mean temperature (degC): ", rng(x$tmean), "\n", sep = "")
  cat("  degree-days per day (dx):      ", rng(x$dx), "\n", sep = "")
  cat("  days without development:      ", sum(x$dx == 0), "\n", sep = "")
  cat("  adult/pupal mortality mu_T:    ", rng(x$mu_T), " per day\n", sep = "")
  cat("  egg viability:                 ", rng(x$egg_viability), "\n", sep = "")
  cat(
    "  eggs aged at:                  ", x$param_values$egg_temperature,
    " temperature\n",
    sep = ""
  )
  cat(
    "  rainfall:                      ",
    if (x$has_rain) {
      sprintf("present; moisture index %s", rng(x$moisture))
    } else {
      "absent (moisture index 1)"
    },
    "\n",
    sep = ""
  )
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# Parameters that sg_forcing() uses; sg_model() compares these with its own.
.sg_forcing_param_names <- function() {
  c(
    "T_L", "T_U", "a_dev", "b_dev", "dx_floor", "X_pre", "T_opt", "TU_ovip",
    "egg_temperature", "T_host", "c_A", "p_A0", "c_E", "p_E0", "eps_s",
    "rain_threshold"
  )
}

# Validate a weather data frame; return date, tmean (and rain) columns.
.sg_check_weather <- function(weather) {
  if (!is.data.frame(weather)) {
    .sg_stop("`weather` must be a data frame with columns date and tmean (or tmax and tmin).")
  }
  if (!"date" %in% names(weather)) {
    .sg_stop("`weather` needs a `date` column.")
  }
  date <- tryCatch(as.Date(weather$date), error = function(e) NULL)
  if (is.null(date) || anyNA(date)) {
    .sg_stop("`weather$date` must be dates (class Date or \"YYYY-MM-DD\" text) without NA.")
  }
  if (nrow(weather) == 0L) {
    .sg_stop("`weather` has no rows.")
  }
  if ("tmean" %in% names(weather)) {
    tmean <- weather$tmean
  } else if (all(c("tmax", "tmin") %in% names(weather))) {
    tmean <- (weather$tmax + weather$tmin) / 2
    if (any(weather$tmin > weather$tmax, na.rm = TRUE)) {
      .sg_warn("Some days have tmin > tmax; check the weather data.")
    }
  } else {
    .sg_stop("`weather` needs a `tmean` column, or both `tmax` and `tmin`.")
  }
  if (!is.numeric(tmean)) {
    .sg_stop("Temperatures in `weather` must be numeric (degrees Celsius).")
  }

  ord <- order(date)
  if (!identical(ord, seq_along(date))) {
    .sg_stop("`weather` must be sorted by date (earliest first).")
  }
  steps <- diff(as.numeric(date))
  if (any(steps == 0)) {
    .sg_stop("`weather` has repeated dates, e.g. ", format(date[which(steps == 0)[1L]]), ".")
  }
  if (any(steps != 1)) {
    i <- which(steps != 1)[1L]
    .sg_stop(
      "`weather` must have one row per consecutive day; there is a gap after ",
      format(date[i]), " (next date ", format(date[i + 1L]), "). ",
      "Fill gaps first, e.g. with sg_fill_weather()."
    )
  }
  if (anyNA(tmean)) {
    .sg_stop(
      "`weather` has missing temperatures, e.g. on ", format(date[which(is.na(tmean))[1L]]),
      ". Fill gaps first, e.g. with sg_fill_weather()."
    )
  }
  if (any(tmean < -40 | tmean > 55)) {
    i <- which(tmean < -40 | tmean > 55)[1L]
    .sg_stop(
      "Implausible daily mean temperature ", .sg_fmt(tmean[i]), " degC on ",
      format(date[i]), " (allowed -40 to 55). Are the units degrees Celsius?"
    )
  }
  out <- data.frame(date = date, tmean = as.numeric(tmean))
  if ("rain" %in% names(weather)) {
    rain <- weather$rain
    if (!is.numeric(rain) || anyNA(rain) || any(rain < 0)) {
      .sg_stop("`weather$rain` must be numeric, >= 0 and without NA (mm per day).")
    }
    out$rain <- as.numeric(rain)
  }
  out
}

# Logical index of the rows between start and end (inclusive).
.sg_date_window <- function(date, start, end) {
  to_date <- function(x, what) {
    d <- tryCatch(as.Date(x), error = function(e) NA)
    if (length(d) != 1L || is.na(d)) {
      .sg_stop("`", what, "` must be a single date (Date or \"YYYY-MM-DD\").")
    }
    d
  }
  s <- if (is.null(start)) min(date) else to_date(start, "start")
  e <- if (is.null(end)) max(date) else to_date(end, "end")
  if (s < min(date) || e > max(date)) {
    .sg_stop(
      "The period ", format(s), " to ", format(e), " is not covered by `weather` (",
      format(min(date)), " to ", format(max(date)), ")."
    )
  }
  if (e < s) {
    .sg_stop("`end` (", format(e), ") is before `start` (", format(s), ").")
  }
  date >= s & date <= e
}

# Moisture index per day: min(1, annual rain / threshold) for the calendar
# year, scaling incomplete years to a full year (S7).
.sg_moisture_index <- function(w, threshold) {
  if (!"rain" %in% names(w)) {
    return(rep(1, nrow(w)))
  }
  yr <- format(w$date, "%Y")
  total <- tapply(w$rain, yr, sum)
  n_days <- tapply(w$rain, yr, length)
  days_in_year <- ifelse(.sg_leap_year(as.integer(names(total))), 366, 365)
  annual <- total * days_in_year / n_days
  idx <- stats::setNames(pmin(1, as.numeric(annual) / threshold), names(total))
  unname(idx[yr])
}

.sg_leap_year <- function(y) (y %% 4L == 0L & y %% 100L != 0L) | y %% 400L == 0L
