#' Example weather: McAllen, Texas, 1960-1983
#'
#' Daily weather for McAllen, Texas (USA), from 1 January 1960 to 31 December
#' 1983, the reference location of Gutierrez & Ponti (2014) and Gutierrez et
#' al. (2019) in the south Texas transition zone of screwworm persistence.
#'
#' @details
#' Source: NOAA NCEI Global Historical Climatology Network - Daily (GHCN-Daily),
#' station `USC00415701` (McAllen). Missing days were filled with
#' [sg_fill_weather()] from station `USW00012959` (McAllen Miller International
#' Airport, less than 1 km away), then by the mean of the neighbouring days,
#' then by the mean of the preceding and following five days, following
#' Gutierrez et al. (2019, Appendix S1). Of 8766 days, 8612 come from the
#' primary station, 152 from the backup station and 2 were interpolated.
#'
#' @return A data frame with 8766 rows and columns `date` (Date), `tmax`,
#'   `tmin`, `tmean` (degrees Celsius; `tmean = (tmax + tmin) / 2`), `rain`
#'   (mm per day) and `source` (`"primary"`, `"backup"` or `"interpolated"`).
#' @references Menne, M.J. et al. (2012). An overview of the Global Historical
#'   Climatology Network-Daily database. *Journal of Atmospheric and Oceanic
#'   Technology* 29:897-910.
#' @seealso [sg_forcing()]
#' @export
#' @examples
#' wx <- sg_example_weather()
#' head(wx)
#' table(wx$source)
sg_example_weather <- function() {
  file <- system.file("extdata", "mcallen_tx_1960_1983.csv", package = "SGDrivE")
  w <- utils::read.csv(file, stringsAsFactors = FALSE)
  w$date <- as.Date(w$date)
  w[, c("date", "tmax", "tmin", "tmean", "rain", "source")]
}

#' Read a NOAA GHCN-Daily station file
#'
#' Reads a station file in the comma-separated "access" format of NOAA's
#' Global Historical Climatology Network - Daily
#' (<https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/>,
#' one file per station, e.g. `USC00415701.csv`) and converts it to daily
#' maximum and minimum temperature and rainfall. The file is read as is:
#' missing days stay missing (use [sg_fill_weather()] to fill them).
#'
#' @param file Path to the station file.
#' @return A data frame with columns `date` (Date), `tmax`, `tmin` (degrees
#'   Celsius) and `rain` (mm per day), with attributes `station` and `name`.
#'   GHCN-Daily stores these values in tenths of a unit; they are divided by 10.
#' @seealso [sg_fill_weather()]
#' @export
#' @examples
#' \dontrun{
#' # after downloading USC00415701.csv from NOAA
#' w <- sg_read_ghcn("USC00415701.csv")
#' }
sg_read_ghcn <- function(file) {
  if (!is.character(file) || length(file) != 1L || !file.exists(file)) {
    .sg_stop("`file` must be the path of an existing GHCN-Daily station file.")
  }
  raw <- utils::read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
  if (!all(c("DATE", "TMAX", "TMIN") %in% names(raw))) {
    .sg_stop(
      "`file` does not look like a GHCN-Daily access file: columns DATE, TMAX and ",
      "TMIN are required."
    )
  }
  num <- function(col) {
    if (col %in% names(raw)) suppressWarnings(as.numeric(raw[[col]])) / 10 else NA_real_
  }
  out <- data.frame(
    date = as.Date(raw$DATE),
    tmax = num("TMAX"),
    tmin = num("TMIN"),
    rain = num("PRCP")
  )
  attr(out, "station") <- if ("STATION" %in% names(raw)) raw$STATION[1L] else NA_character_
  attr(out, "name") <- if ("NAME" %in% names(raw)) raw$NAME[1L] else NA_character_
  out
}

#' Fill gaps in daily weather data
#'
#' Produces a complete daily series from a primary station, optionally using a
#' nearby backup station, following the procedure of Gutierrez et al. (2019,
#' Appendix S1) for McAllen, Texas.
#'
#' @param primary Data frame with columns `date`, `tmax`, `tmin` and `rain`
#'   (e.g. from [sg_read_ghcn()]); missing values and missing days allowed.
#' @param backup Optional data frame in the same format from a nearby station.
#' @param start,end First and last date of the output (`Date` or
#'   `"YYYY-MM-DD"`). Default: the range of `primary`.
#'
#' @details
#' Each variable (`tmax`, `tmin`, `rain`) is filled in three steps:
#' 1. values missing at the primary station are taken from the backup station;
#' 2. isolated missing days are set to the mean of the preceding and following
#'    day when both are available;
#' 3. any remaining gaps are set to the mean of the available values among the
#'    preceding and following five days, working forwards in time.
#'
#' Days whose `tmax` and `tmin` both came from the primary station are labelled
#' `"primary"`, those completed from the backup station `"backup"`, and the
#' rest `"interpolated"`. Values that cannot be filled (no data within five
#' days) stay `NA`.
#'
#' @return A data frame with one row per day from `start` to `end` and columns
#'   `date`, `tmax`, `tmin`, `rain`, `source` and `tmean = (tmax + tmin) / 2`.
#' @references Gutierrez, A.P., Ponti, L. & Arias, P.A. (2019). *Medical and
#'   Veterinary Entomology* 33:282-295. \doi{10.1111/mve.12362}
#' @seealso [sg_read_ghcn()], [sg_example_weather()]
#' @export
#' @examples
#' w <- data.frame(
#'   date = as.Date("2000-01-01") + c(0:2, 4:6),
#'   tmax = c(20, NA, 22, 23, 24, 25), tmin = c(10, 11, 12, 13, 14, 15),
#'   rain = 0
#' )
#' sg_fill_weather(w)
sg_fill_weather <- function(primary, backup = NULL, start = NULL, end = NULL) {
  pri <- .sg_check_station(primary, "primary")
  rng_start <- if (is.null(start)) min(pri$date) else as.Date(start)
  rng_end <- if (is.null(end)) max(pri$date) else as.Date(end)
  if (is.na(rng_start) || is.na(rng_end) || rng_end < rng_start) {
    .sg_stop("`start` and `end` must be dates with start <= end.")
  }
  w <- merge(data.frame(date = seq(rng_start, rng_end, by = "day")), pri, all.x = TRUE)
  vars <- c("tmax", "tmin", "rain")
  w$source <- ifelse(is.na(w$tmax) | is.na(w$tmin), NA_character_, "primary")

  if (!is.null(backup)) {
    bk <- .sg_check_station(backup, "backup")
    b <- merge(w["date"], bk, all.x = TRUE)
    for (v in vars) {
      fill <- is.na(w[[v]]) & !is.na(b[[v]])
      w[[v]][fill] <- b[[v]][fill]
    }
    w$source[is.na(w$source) & !is.na(w$tmax) & !is.na(w$tmin)] <- "backup"
  }
  for (v in vars) {
    w[[v]] <- .sg_fill_gaps(w[[v]])
  }
  w$source[is.na(w$source)] <- "interpolated"
  w$tmean <- (w$tmax + w$tmin) / 2
  # merge() can store dates as integers; keep the usual double representation
  w$date <- as.Date(as.numeric(w$date), origin = "1970-01-01")
  w
}

#' Cold- and heat-stress indices
#'
#' Computes the annual cold-stress index mu_cold of Gutierrez & Ponti (2014):
#' the sum of the daily pupal/adult mortality proportions p_A(T) over days with
#' mean temperature below `optimal_survival_temp` (27.2 degC). Values up to about 10 indicate
#' that the screwworm can persist year-round; higher values indicate winter
#' die-off. The heat-stress index mu_hot sums p_A(T) over days at or above
#' `optimal_survival_temp`.
#'
#' @param x An `sg_forcing` object, or a weather data frame accepted by
#'   [sg_forcing()].
#' @param params An `sg_params` object (used when `x` is a data frame; ignored for
#'   an `sg_forcing` object, which carries its own values).
#' @param period `"calendar"` (1 January to 31 December, as in the published
#'   index values) or `"jul_jun"` (1 July of the previous year to 30 June,
#'   reported for the year in which it ends).
#' @return A data frame with columns `year`, `mu_cold`, `mu_hot`, `n_days` and
#'   `complete` (`TRUE` when the period has all its days).
#' @references Gutierrez, A.P. & Ponti, L. (2014). *Agricultural and Forest
#'   Entomology* 16:158-173. \doi{10.1111/afe.12046}
#' @export
#' @examples
#' ci <- sg_cold_index(sg_example_weather())
#' ci[ci$year %in% c(1971, 1979), ]
sg_cold_index <- function(x, params = sg_params(), period = c("calendar", "jul_jun")) {
  period <- match.arg(period)
  if (inherits(x, "sg_forcing")) {
    dates <- x$dates
    p_a <- x$p_A
    tmean <- x$tmean
    t_opt <- x$param_values$optimal_survival_temp
  } else {
    f <- sg_forcing(x, .sg_as_params(params))
    dates <- f$dates
    p_a <- f$p_A
    tmean <- f$tmean
    t_opt <- params$optimal_survival_temp
  }
  yr <- as.integer(format(dates, "%Y"))
  if (period == "jul_jun") {
    yr <- yr + (as.integer(format(dates, "%m")) >= 7L)
  }
  cold <- tmean < t_opt
  years <- sort(unique(yr))
  mu_cold <- vapply(years, function(y) sum(p_a[yr == y & cold]), 0)
  mu_hot <- vapply(years, function(y) sum(p_a[yr == y & !cold]), 0)
  n_days <- vapply(years, function(y) sum(yr == y), 0L)
  # Both a calendar year y and the period Jul (y-1) - Jun (y) contain
  # 29 February of year y when y is a leap year.
  expected <- ifelse(.sg_leap_year(years), 366L, 365L)
  data.frame(
    year = years, mu_cold = mu_cold, mu_hot = mu_hot, n_days = n_days,
    complete = n_days == expected
  )
}

# ---- internal ------------------------------------------------------------------

.sg_check_station <- function(x, what) {
  if (!is.data.frame(x) || !all(c("date", "tmax", "tmin", "rain") %in% names(x))) {
    .sg_stop("`", what, "` must be a data frame with columns date, tmax, tmin and rain.")
  }
  x$date <- as.Date(x$date)
  if (anyNA(x$date) || anyDuplicated(x$date) > 0L) {
    .sg_stop("`", what, "$date` must be dates without NA or repeats.")
  }
  x[, c("date", "tmax", "tmin", "rain")]
}

# Steps 2 and 3 of sg_fill_weather() for one variable (same order and logic
# as the procedure used to build the bundled McAllen series).
.sg_fill_gaps <- function(x) {
  n <- length(x)
  i <- which(is.na(x))
  i <- i[i > 1L & i < n]
  ok <- !is.na(x[i - 1L]) & !is.na(x[i + 1L])
  x[i[ok]] <- (x[i[ok] - 1L] + x[i[ok] + 1L]) / 2
  for (j in which(is.na(x))) {
    # At the series ends these index ranges run past the data; the extra
    # entries are NA (or index 0, dropped) and are ignored by the mean.
    nb <- x[c(max(1L, j - 5L):(j - 1L), (j + 1L):min(n, j + 5L))]
    if (any(!is.na(nb))) x[j] <- mean(nb, na.rm = TRUE)
  }
  x
}
