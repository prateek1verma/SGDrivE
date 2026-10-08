# Temperature-driven functions of Gutierrez & Ponti (2014), model
# specification S1-S4. All functions are vectorised over temperature and take
# their constants from an sg_params object.

#' Temperature-driven life-history functions
#'
#' Functions of daily mean temperature from Gutierrez & Ponti (2014) that drive
#' development, mortality and oviposition (model specification S1-S4). They are
#' vectorised over temperature.
#'
#' @param temp Daily mean temperature (degrees Celsius), numeric vector.
#' @param temp5 Running 5-day mean temperature (degrees Celsius), numeric vector.
#' @param p An `sg_params` object supplying the constants (see
#'   [sg_param_table()], groups "temperature" and "mortality").
#' @param type What to return:
#'   * `sg_temp_mortality()`: `"hazard"` (default), the instantaneous rate
#'     mu_T = -log(max(eps_s, 1 - p_A(T))) per day used by the solvers; or
#'     `"proportion"`, the daily proportion dying p_A(T) clamped to \[0, 1\]
#'     (used for the cold-stress index).
#'   * `sg_egg_viability()`: `"viability"` (default), the hatch probability
#'     h = 1 - min(1, c_E (T5 - T_opt)^2 + p_E0); or `"hazard"`, the egg
#'     mortality rate mu_E = -log(max(eps_s, h)).
#'
#' @details
#' * **S1 physiological time.** `sg_degree_days()` returns the degree-days
#'   accumulated in one day,
#'   \deqn{dx(T) = X_pre max(dx_floor, 24 a (T - T_L) / (1 + b^(T - T_U))),}
#'   the development rate of Adams (1979) pre-oviposition data scaled by the
#'   61 dd pre-oviposition period. It is close to 0.95 (T - 14.5) up to about
#'   33 degC, peaks near 36 degC and declines at higher temperatures.
#' * **S2 pupal/adult mortality.** p_A(T) = c_A (T - T_opt)^2 + p_A0.
#' * **S3 egg viability**, from the 5-day mean temperature.
#' * **S4 oviposition scalar.**
#'   phi_ov(T) = max(0, 1 - ((T - T_L - psi) / psi)^2) with
#'   psi = (TU_ovip - T_L) / 2: 1 at about 29 degC, 0 outside
#'   \[T_L, TU_ovip\].
#'
#' @return Numeric vector of the same length as the temperature input: degree-days
#'   per day (`sg_degree_days`), a rate per day or a daily proportion
#'   (`sg_temp_mortality`, `sg_egg_viability`), or a scalar in \[0, 1\]
#'   (`sg_ovip_scalar`).
#' @references Gutierrez, A.P. & Ponti, L. (2014). *Agricultural and Forest
#'   Entomology* 16:158-173. \doi{10.1111/afe.12046}
#' @seealso [sg_forcing()], [sg_cold_index()]
#' @name temperature_functions
#' @examples
#' p <- sg_params()
#' temp <- c(10, 20, 27.2, 35, 40)
#' sg_degree_days(temp, p)
#' sg_temp_mortality(temp, p)
#' sg_temp_mortality(temp, p, type = "proportion")
#' sg_egg_viability(temp, p)
#' sg_ovip_scalar(temp, p)
NULL

#' @rdname temperature_functions
#' @export
sg_degree_days <- function(temp, p = sg_params()) {
  .sg_check_temperature_input(temp, "temp")
  p <- .sg_as_params(p)
  r_hour <- p$a_dev * (temp - p$T_L) / (1 + p$b_dev^(temp - p$T_U))
  p$X_pre * pmax(p$dx_floor, 24 * r_hour)
}

#' @rdname temperature_functions
#' @export
sg_temp_mortality <- function(temp, p = sg_params(), type = c("hazard", "proportion")) {
  .sg_check_temperature_input(temp, "temp")
  p <- .sg_as_params(p)
  type <- match.arg(type)
  p_a <- p$c_A * (temp - p$T_opt)^2 + p$p_A0
  if (type == "proportion") {
    return(pmin(1, pmax(0, p_a)))
  }
  -log(pmax(p$eps_s, 1 - p_a))
}

#' @rdname temperature_functions
#' @export
sg_egg_viability <- function(temp5, p = sg_params(), type = c("viability", "hazard")) {
  .sg_check_temperature_input(temp5, "temp5")
  p <- .sg_as_params(p)
  type <- match.arg(type)
  h <- 1 - pmin(1, p$c_E * (temp5 - p$T_opt)^2 + p$p_E0)
  if (type == "viability") {
    return(h)
  }
  -log(pmax(p$eps_s, h))
}

#' @rdname temperature_functions
#' @export
sg_ovip_scalar <- function(temp, p = sg_params()) {
  .sg_check_temperature_input(temp, "temp")
  p <- .sg_as_params(p)
  psi <- (p$TU_ovip - p$T_L) / 2
  pmax(0, 1 - ((temp - p$T_L - psi) / psi)^2)
}

# ---- internal ------------------------------------------------------------------

.sg_check_temperature_input <- function(x, what) {
  if (!is.numeric(x)) {
    .sg_stop("`", what, "` must be numeric (temperature in degrees Celsius).")
  }
  invisible(TRUE)
}

.sg_as_params <- function(p) {
  if (!inherits(p, "sg_params")) {
    .sg_stop("`p` must be an sg_params object; create one with sg_params().")
  }
  p
}

# Running mean over the current and previous (k - 1) values; at the start of
# the series the mean of the available values (Gutierrez & Ponti code, Aver5).
# A direct window mean (not a running sum), so long series accumulate no
# rounding error.
.sg_running_mean <- function(x, k = 5L) {
  vapply(seq_along(x), function(t) mean(x[max(1L, t - k + 1L):t]), numeric(1))
}
