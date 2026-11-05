# Calibration of a population-size parameter to a target density (D-026;
# general form, decision 9 of D-030).

#' Calibrate a parameter to a target population density
#'
#' Finds the value of a parameter at which a summary of the simulated
#' population, by default the mean adult density over the model's forcing
#' period, equals a target.
#'
#' @param model An `sg_model` object ([sg_model()]); its forcing period is
#'   the reference period (for example one year).
#' @param target Target density, individuals per km^2 (> 0).
#' @param parameter The parameter to calibrate: `"wound_density"` (wounds per
#'   km^2; default with `wound_model = "static"`), `"wounded_fraction"` (the
#'   fraction of hosts with an attractive wound without flies; default with
#'   `wound_model = "hosts"`, because host density usually comes from a
#'   census) or `"host_density"` (animals per km^2, host-wound model).
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
#' @details Each evaluation runs the ODE from a burn-in. Without an Allee
#'   effect (`male_density_half_mating = 0`) and releases every equation is
#'   homogeneous of degree one in the state and the wound density (static
#'   model) or the host density (host-wound model), so after the burn-in the
#'   population is exactly proportional to that parameter: one evaluation
#'   gives the calibrated value in closed form (`value * target / achieved`),
#'   checked by a second. Otherwise (Allee effect, releases, or
#'   `wounded_fraction`) the value is found by a one-dimensional root search,
#'   on the log scale (densities) or the logit scale (`wounded_fraction`). The
#'   densities are per km^2, so the result does not depend on the node area.
#' @return An object of class `sg_calibration`, a list with `parameter`,
#'   `value` (the calibrated value), an element named after the parameter
#'   with the same value (e.g. `wound_density`), `params` (the model's
#'   parameters with that value), `model` (the model rebuilt with it),
#'   `target`, `achieved`, `stat`, `stages`, `method` (`"proportional"` or
#'   `"root search"`) and `evaluations`.
#' @seealso [calibrate_wound_density()], [check_model()], [sg_init()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   calibrate_population(mod, target = 100) # wound density for 100 adults per km^2
#' }
calibrate_population <- function(model, target, parameter = NULL,
                                 stat = c("mean", "peak", "min", "quantile"),
                                 stages = "adults", period = NULL, prob = 0.5, years = 1L,
                                 weather = NULL, tol = 1e-8) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  .sg_one_node_only(model, "calibrate_population()")
  hosts <- !is.null(model$hosts)
  allowed <- if (hosts) c("wounded_fraction", "host_density") else "wound_density"
  if (is.null(parameter)) parameter <- allowed[1L]
  parameter <- .sg_match_choice(
    parameter, c("wound_density", "wounded_fraction", "host_density"), "parameter"
  )
  if (!parameter %in% allowed) {
    .sg_stop(
      "parameter = \"", parameter, "\" is not used with wound_model = \"",
      model$params$wound_model, "\"; use ", paste0("\"", allowed, "\"", collapse = " or "), "."
    )
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
  v_start <- model$params[[parameter]]
  if (!(v_start > 0)) {
    .sg_stop(
      "The model has ", if (parameter == "wound_density") "no wounds (" else "", parameter,
      " = ", format(v_start), if (parameter == "wound_density") ")" else "",
      "; give a positive ", parameter, " to start from."
    )
  }

  runs <- new.env(parent = emptyenv()) # counts model runs across value_at() calls
  runs$n <- 0L
  value_at <- function(v) {
    runs$n <- runs$n + 1L
    .sg_calib_value(
      .sg_model_with_value(model, parameter, v), stages, days, stat, prob, years,
      weather
    )
  }
  v0 <- value_at(v_start)
  if (!is.finite(v0) || v0 < 1e-9) { # essentially extinct
    .sg_stop(
      "The population is zero over the reference period at ", parameter, " = ", format(v_start),
      ", so ", parameter, " cannot be calibrated; check that the population can persist ",
      "(check_model())."
    )
  }
  if (nrow(model$releases) > 0L) {
    .sg_warn(
      "The model has releases; calibrating with them included. Calibrations of the wild ",
      "population usually use a model without releases."
    )
  }
  proportional <- parameter != "wounded_fraction" &&
    model$params$male_density_half_mating == 0 && nrow(model$releases) == 0L
  if (proportional) {
    v <- v_start * target / v0
    method <- "proportional"
  } else {
    # densities on the log scale, wounded_fraction on the logit scale
    to_x <- if (parameter == "wounded_fraction") stats::qlogis else log
    from_x <- if (parameter == "wounded_fraction") stats::plogis else exp
    gap <- function(x) log(value_at(from_x(x))) - log(target)
    guess <- if (parameter == "wounded_fraction") {
      to_x(min(max(v_start * target / v0, 1e-6), 0.999))
    } else {
      to_x(v_start * target / v0)
    }
    limit <- function(dir) {
      if (parameter == "wounded_fraction") {
        "No wounded fraction between 0 and 1 gives the target; check the host density."
      } else {
        paste0("No ", parameter, " gives a value as ", dir, " as the target.")
      }
    }
    # search range: a wounded fraction between 1e-6 and 0.999 (wounding rates
    # grow without bound as it approaches 1); densities over 80 e-folds
    bounds <- if (parameter == "wounded_fraction") stats::qlogis(c(1e-6, 0.999)) else c(-Inf, Inf)
    lo <- hi <- guess
    g_lo <- g_hi <- gap(guess)
    step <- log(2)
    while (g_lo > 0) {
      if (lo <= bounds[1L] || step > 80) .sg_stop(limit("low"))
      lo <- max(lo - step, bounds[1L])
      g_lo <- gap(lo)
      step <- 2 * step
    }
    step <- log(2)
    while (g_hi < 0) {
      if (hi >= bounds[2L] || step > 80) .sg_stop(limit("high"))
      hi <- min(hi + step, bounds[2L])
      g_hi <- gap(hi)
      step <- 2 * step
    }
    x <- if (lo == hi) lo else stats::uniroot(gap, c(lo, hi), tol = tol)$root
    v <- from_x(x)
    method <- "root search"
  }
  calibrated <- .sg_model_with_value(model, parameter, v)
  achieved <- value_at(v)
  out <- list(
    parameter = parameter, value = v, params = calibrated$params, model = calibrated,
    target = target, achieved = achieved, stat = stat, stages = stages, method = method,
    evaluations = runs$n
  )
  out[[parameter]] <- v
  structure(out, class = "sg_calibration")
}

#' Calibrate the wound density to a target population density
#'
#' Finds the wound density `wound_density` (wounds per km^2) at which a
#' summary of the simulated population, by default the mean adult density
#' over the model's forcing period, equals a target. A shortcut for
#' [calibrate_population()] with `parameter = "wound_density"`, for models
#' with `wound_model = "static"`.
#'
#' @inheritParams calibrate_population
#' @return An object of class `sg_calibration`; see [calibrate_population()].
#'   The calibrated value is also in `wound_density`.
#' @seealso [calibrate_population()], [check_model()], [sg_init()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   cal <- calibrate_wound_density(mod, target = 100) # mean of 100 adults per km^2 in 1971
#'   cal
#'   # simulate with the calibrated model
#'   sim <- simulate_ode(cal$model, init = sg_init(cal$model, "burn_in"))
#' }
calibrate_wound_density <- function(model, target, stat = c("mean", "peak", "min", "quantile"),
                                    stages = "adults", period = NULL, prob = 0.5, years = 1L,
                                    weather = NULL, tol = 1e-8) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  if (!is.null(model$hosts)) {
    .sg_stop(
      "wound_density is not used with wound_model = \"hosts\"; use calibrate_population(), ",
      "e.g. with parameter = \"wounded_fraction\"."
    )
  }
  calibrate_population(model, target,
    parameter = "wound_density", stat = stat, stages = stages,
    period = period, prob = prob, years = years, weather = weather, tol = tol
  )
}

#' @export
print.sg_calibration <- function(x, ...) {
  units <- switch(EXPR = x$parameter,
    wound_density = " wounds per km^2",
    host_density = " animals per km^2",
    ""
  )
  cat("<sg_calibration>\n")
  cat(sprintf(
    "  %s = %.6g%s (%s, %d model runs)\n", x$parameter, x$value, units, x$method,
    x$evaluations
  ))
  cat(sprintf(
    "  %s of %s per km^2: target %.6g, achieved %.6g (relative difference %.2g)\n",
    x$stat, paste(x$stages, collapse = ", "), x$target, x$achieved,
    abs(x$achieved / x$target - 1)
  ))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# The model with one population-size parameter changed, its derived counts
# and hash updated.
.sg_model_with_value <- function(model, parameter, value) {
  m <- model
  args <- stats::setNames(list(value), parameter)
  m$params <- suppressWarnings(do.call(stats::update, c(list(model$params), args)))
  if (is.null(model$hosts)) {
    m$counts$H <- m$params$wound_density * model$area
  } else {
    m$hosts <- .sg_host_setup(m$params, model$area)
  }
  m$meta$adult_females <- NULL # the area was chosen for the old value
  m$meta$hash <- .sg_model_hash(m)
  m
}

# The model with wound density h (per km^2); used by older code paths.
.sg_model_with_h <- function(model, h) .sg_model_with_value(model, "wound_density", h)

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
