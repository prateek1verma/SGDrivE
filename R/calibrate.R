# Calibration of the wound density H to a target population density (D-026).

#' Calibrate the wound density to a target population density
#'
#' Finds the wound density `H` (wounds per km^2) at which a summary of the
#' simulated population, by default the mean adult density over the model's
#' forcing period, equals a target.
#'
#' @param model An `sg_model` object ([sg_model()]); its forcing period is
#'   the reference period (for example one year).
#' @param target Target density, individuals per km^2 (> 0).
#' @param stat Summary over the reference period: `"mean"` (default),
#'   `"peak"` (maximum), `"min"` or `"quantile"`.
#' @param stages Stages counted (see [sg_totals()]); default adults.
#' @param period Optional two dates within the forcing: the reference period
#'   (default the whole forcing period).
#' @param prob For `stat = "quantile"`: the probability.
#' @param years,weather Burn-in before the reference period (see [sg_init()]
#'   with `method = "burn_in"`); without `weather` the model's first year of
#'   forcing is repeated.
#' @param tol Relative tolerance of the root search.
#'
#' @details Each evaluation runs the ODE from a burn-in. With no Allee effect
#'   (`h_M = 0`) and no releases every equation is homogeneous of degree one
#'   in the state and H, so after the burn-in the population is exactly
#'   proportional to H: one evaluation gives the calibrated H in closed form
#'   (`H * target / value`), checked by a second. Otherwise (Allee effect, or
#'   releases, which are usually left out of a calibration of the wild
#'   population) H is found by a one-dimensional root search on log H.
#'   Because H is per km^2, the result does not depend on the node area.
#' @return An object of class `sg_calibration`, a list with `H` (calibrated
#'   wounds per km^2), `params` (the model's parameters with that H),
#'   `model` (the model rebuilt with it), `target`, `achieved`, `stat`,
#'   `stages`, `method` (`"proportional"` or `"root search"`) and
#'   `evaluations`.
#' @seealso [check_model()], [sg_init()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   cal <- calibrate_H(mod, target = 100) # mean of 100 adults per km^2 in 1971
#'   cal
#'   # simulate with the calibrated model
#'   sim <- simulate_ode(cal$model, init = sg_init(cal$model, "burn_in"))
#' }
calibrate_H <- function(model, target, stat = c("mean", "peak", "min", "quantile"),
                        stages = "adults", period = NULL, prob = 0.5, years = 1L,
                        weather = NULL, tol = 1e-8) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target) || target <= 0) {
    .sg_stop("`target` must be a single positive density (individuals per km^2).")
  }
  stat <- .sg_match_choice(stat, c("mean", "peak", "min", "quantile"), "stat")
  if (stat == "quantile" && !(is.numeric(prob) && length(prob) == 1L && prob >= 0 && prob <= 1)) {
    .sg_stop("`prob` must be a single probability between 0 and 1.")
  }
  .sg_expand_stages(stages)
  days <- .sg_calib_days(model, period)
  if (model$counts$H <= 0) .sg_stop("The model has no wounds (H = 0); give a positive H to start from.")

  evaluations <- 0L
  value_at <- function(h) {
    evaluations <<- evaluations + 1L
    .sg_calib_value(.sg_model_with_h(model, h), stages, days, stat, prob, years, weather)
  }
  h0 <- model$params$H
  v0 <- value_at(h0)
  if (!is.finite(v0) || v0 < 1e-9) { # essentially extinct
    .sg_stop(
      "The population is zero over the reference period at H = ", format(h0),
      ", so H cannot be calibrated; check that the population can persist (check_model())."
    )
  }
  proportional <- model$params$h_M == 0 && nrow(model$releases) == 0L
  if (nrow(model$releases) > 0L) {
    .sg_warn(
      "The model has releases; calibrating with them included. Calibrations of the wild ",
      "population usually use a model without releases."
    )
  }
  if (proportional) {
    h <- h0 * target / v0
    method <- "proportional"
  } else {
    gap <- function(log_h) log(value_at(exp(log_h))) - log(target)
    guess <- log(h0 * target / v0)
    lo <- hi <- guess
    g_lo <- g_hi <- gap(guess)
    step <- log(2)
    while (g_lo > 0) {
      lo <- lo - step
      g_lo <- gap(lo)
      step <- 2 * step
      if (step > 80) .sg_stop("No wound density gives a value as low as the target.")
    }
    step <- log(2)
    while (g_hi < 0) {
      hi <- hi + step
      g_hi <- gap(hi)
      step <- 2 * step
      if (step > 80) .sg_stop("No wound density gives a value as high as the target.")
    }
    h <- if (lo == hi) exp(lo) else exp(stats::uniroot(gap, c(lo, hi), tol = tol)$root)
    method <- "root search"
  }
  calibrated <- .sg_model_with_h(model, h)
  achieved <- value_at(h)
  structure(
    list(
      H = h, params = calibrated$params, model = calibrated, target = target,
      achieved = achieved, stat = stat, stages = stages, method = method,
      evaluations = evaluations
    ),
    class = "sg_calibration"
  )
}

#' @export
print.sg_calibration <- function(x, ...) {
  cat("<sg_calibration>\n")
  cat(sprintf(
    "  H = %.6g wounds per km^2 (%s, %d model runs)\n", x$H, x$method, x$evaluations
  ))
  cat(sprintf(
    "  %s of %s per km^2: target %.6g, achieved %.6g (relative difference %.2g)\n",
    x$stat, paste(x$stages, collapse = ", "), x$target, x$achieved,
    abs(x$achieved / x$target - 1)
  ))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# The model with wound density h (per km^2), its hash updated.
.sg_model_with_h <- function(model, h) {
  m <- model
  m$params <- stats::update(model$params, H = h)
  m$counts$H <- h * model$area
  m$meta$hash <- .sg_model_hash(m)
  m
}

# Forcing days (1-based) of the reference period.
.sg_calib_days <- function(model, period) {
  dates <- model$forcing$dates
  if (is.null(period)) {
    return(seq_along(dates))
  }
  period <- tryCatch(as.Date(period), error = function(e) NULL)
  if (length(period) != 2L || anyNA(period) || period[1L] > period[2L]) {
    .sg_stop("`period` must be two dates (start, end).")
  }
  days <- which(dates >= period[1L] & dates <= period[2L])
  if (!length(days)) .sg_stop("`period` is outside the model's forcing.")
  days
}

# Summary of the ODE population (per km^2) over the reference days, after a
# burn-in.
.sg_calib_value <- function(m, stages, days, stat, prob, years, weather) {
  ok <- is.numeric(years) && length(years) == 1L && is.finite(years) && years >= 1 &&
    years == round(years)
  if (!ok) .sg_stop("`years` must be a whole number >= 1.")
  init <- suppressWarnings(.sg_init_burn_in(m, as.numeric(sg_init(m)), years, weather,
    check_start = FALSE
  ))$x
  sim <- simulate_ode(m, init = init)
  v <- sg_totals(sim, by = character(0), stages = stages, per = "km2")$value[-1L][days]
  switch(EXPR = stat,
    mean = mean(v),
    peak = max(v),
    min = min(v),
    quantile = unname(stats::quantile(v, prob))
  )
}
