#' Model parameters
#'
#' Creates a validated set of ScrewDrivE model parameters. Every parameter has a
#' default (see [sg_param_table()] for the full list with units, ranges and
#' sources); arguments given in `...` override the defaults.
#'
#' @param ... Named parameter values that override the defaults, e.g.
#'   `adult_substages = 8, eggs_per_day_max = 67`. Names must match
#'   `sg_param_table()$name`. To pass a
#'   list `x` of values, use `do.call(sg_params, x)`.
#' @param preset Which set of defaults to start from:
#'   \describe{
#'     \item{`"default"`}{ScrewDrivE defaults: Gutierrez & Ponti (2014) biology
#'       with speed-tuned sub-stage numbers (at most 10 per stage), mating at
#'       rate `nu`, and searching mated females as competitors for wounds.}
#'     \item{`"gp2014"`}{The configuration of the Gutierrez & Ponti (2014)
#'       Pascal code, for regression tests and benchmarking: 25 sub-stages per
#'       stage, egg laying from adult class 8 (`oviposition_age_dd` = 7/25 of
#'       `adult_lifespan_dd`), the
#'       0.61 dd per day floor on
#'       development, all adults competing for wounds, and instant mating (the
#'       2014 model has no mating step).}
#'   }
#'   Arguments in `...` are applied on top of the preset.
#'
#' @details
#' Each value is checked against its type and allowed range. Cross-checks
#' then catch combinations that cannot work, for example an age at first
#' oviposition `oviposition_age_dd` longer than the adult lifespan
#' `adult_lifespan_dd`. Settings that
#' are valid but have no effect (such as `larval_crowding_coef > 0` without larval density
#' dependence) give a warning.
#'
#' Units follow the conventions of the package: temperatures in degrees
#' Celsius, durations in degree-days (`dd`) or days (`d`), rates per day.
#' Parameters given per km^2 (`wound_density`, `male_density_half_mating`,
#' `host_density`) are converted to counts per node when the model is
#' assembled.
#'
#' **Wound model** (D-030). With `wound_model = "static"` (default) the number
#' of wounds is fixed (`wound_density`) and every fertile egg batch is a case.
#' With `wound_model = "hosts"` hosts move between having no wound, an
#' attractive wound and an infestation, and a case is a new infestation of a
#' host; `host_density` and `wounded_fraction` must then be given (they have
#' no defaults), and `wound_density` is not used.
#'
#' @return An object of class `sg_params`: a named list of parameter values in
#'   the order of [sg_param_table()], with attribute `preset`.
#' @seealso [sg_param_table()] for definitions; [update.sg_params()] to change
#'   values.
#' @export
#' @examples
#' p <- sg_params()
#' p
#'
#' # Change values
#' p2 <- sg_params(adult_substages = 8, mating_rate = 1)
#' p2$adult_substages
#'
#' # Configuration of the Gutierrez & Ponti (2014) code
#' gp <- sg_params(preset = "gp2014")
#' gp$adult_substages
sg_params <- function(..., preset = c("default", "gp2014")) {
  preset <- match.arg(preset)
  overrides <- list(...)
  values <- .sg_param_defaults()
  values[names(.sg_param_presets()[[preset]])] <- .sg_param_presets()[[preset]]
  values <- .sg_apply_overrides(values, overrides)
  values <- .sg_validate_params(values)
  structure(values, class = "sg_params", preset = preset)
}

# Values that each preset changes relative to the defaults.
.sg_param_presets <- function() {
  list(
    default = list(),
    gp2014 = list(
      egg_substages = 25L, larva_substages = 25L, pupa_substages = 25L, adult_substages = 25L,
      # laying from class 8 of 25: 7 classes of 202.8 / 25 dd before laying
      oviposition_age_dd = 7 / 25 * 202.8,
      dev_rate_floor = 0.01,
      wound_competitors = "all_adults",
      mating_mode = "instant"
    )
  )
}

# Former parameter names (Gutierrez & Ponti symbols) and their current names
# (D-027), for informative errors; the symbols remain in sg_param_table().
.sg_param_renames <- function() {
  c(
    T_L = "dev_temp_min",
    T_U = "dev_temp_decline",
    a_dev = "dev_rate_slope",
    b_dev = "dev_decline_steepness",
    dx_floor = "dev_rate_floor",
    X_pre = "dd_scale_preoviposition",
    T_opt = "optimal_survival_temp",
    TU_ovip = "oviposition_temp_max",
    egg_temperature = "egg_temp_source",
    T_host = "host_body_temp",
    X_E = "egg_duration_dd",
    tau_L = "larva_duration_days",
    X_P = "pupa_duration_dd",
    X_A = "adult_lifespan_dd",
    n_E = "egg_substages",
    n_L = "larva_substages",
    n_P = "pupa_substages",
    n_A = "adult_substages",
    X_ovi = "oviposition_age_dd",
    X_rec = "female_receptive_age_dd",
    X_mat = "male_mature_age_dd",
    p_A0 = "adult_mortality_min",
    c_A = "adult_mortality_curvature",
    p_E0 = "egg_failure_min",
    c_E = "egg_failure_curvature",
    s_L = "larva_survival_daily",
    mu_b_F = "female_background_mortality",
    mu_b_M = "male_background_mortality",
    larval_dd = "larval_crowding_mode",
    gamma_L = "larval_crowding_coef",
    mu_W = "pupa_drought_mortality",
    rain_threshold = "drought_free_rainfall",
    eps_s = "survival_floor",
    R = "eggs_per_day_max",
    H = "wound_density",
    a_s = "wound_search_rate",
    sigma0 = "wound_finding_max",
    search_denominator = "wound_competitors",
    oviposition = "oviposition_mode",
    nu = "mating_rate",
    h_M = "male_density_half_mating",
    mating = "mating_mode",
    ext_threshold = "extinction_threshold"
  )
}

# Check names of user overrides and apply them.
.sg_apply_overrides <- function(values, overrides) {
  if (length(overrides) == 0L) {
    return(values)
  }
  nms <- names(overrides)
  if (is.null(nms) || any(!nzchar(nms))) {
    .sg_stop(
      "All parameters must be named, e.g. sg_params(adult_substages = 8). ",
      "To pass a list x, use do.call(sg_params, x)."
    )
  }
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0L) {
    .sg_stop("Parameter(s) given more than once: ", paste(dup, collapse = ", "), ".")
  }
  unknown <- setdiff(nms, names(values))
  if (length(unknown) > 0L) {
    renamed <- .sg_param_renames()
    msgs <- vapply(unknown, function(u) {
      if (u %in% names(renamed)) {
        return(paste0("'", u, "' (renamed to '", renamed[[u]], "')"))
      }
      paste0("'", u, "'", .sg_suggest(u, names(values)))
    }, "")
    .sg_stop(
      "Unknown parameter(s): ", paste(msgs, collapse = "; "),
      " See sg_param_table() for valid names."
    )
  }
  for (nm in nms) {
    if (is.null(overrides[[nm]])) {
      .sg_stop("Parameter '", nm, "' is NULL; give a value or leave it out.")
    }
    values[[nm]] <- overrides[[nm]]
  }
  values
}

# Validate every value against the table, then cross-check combinations.
# Returns the values with integers stored as integer.
.sg_validate_params <- function(values) {
  rows <- .sg_param_rows()
  for (r in rows) {
    values[[r$name]] <- .sg_check_value(values[[r$name]], r)
  }
  .sg_cross_check_params(values)
  values
}

.sg_check_value <- function(x, r) {
  nm <- r$name
  if (length(x) != 1L) {
    .sg_stop("Parameter '", nm, "' must be a single value, not ", .sg_fmt(x), ".")
  }
  if (is.na(x)) {
    if (!r$allow_na) {
      .sg_stop("Parameter '", nm, "' cannot be NA.")
    }
    return(if (r$type == "integer") NA_integer_ else x)
  }
  if (r$type == "choice") {
    choices <- strsplit(r$choices, "|", fixed = TRUE)[[1L]]
    if (!is.character(x) || !x %in% choices) {
      .sg_stop(
        "Parameter '", nm, "' must be one of ",
        paste0("\"", choices, "\"", collapse = ", "), "; got ", .sg_fmt(x), ".",
        if (is.character(x)) .sg_suggest(x, choices) else ""
      )
    }
    return(x)
  }
  if (!is.numeric(x) || !is.finite(x)) {
    .sg_stop("Parameter '", nm, "' must be a finite number; got ", .sg_fmt(x), ".")
  }
  if (r$type == "integer") {
    if (abs(x - round(x)) > 1e-8) {
      .sg_stop("Parameter '", nm, "' must be a whole number; got ", .sg_fmt(x), ".")
    }
    x <- as.integer(round(x))
  }
  lo_ok <- if (substr(r$bounds, 1L, 1L) == "[") x >= r$lower else x > r$lower
  hi_ok <- if (substr(r$bounds, 2L, 2L) == "]") x <= r$upper else x < r$upper
  if (!lo_ok || !hi_ok) {
    .sg_stop(
      "Parameter '", nm, "' = ", .sg_fmt(x), " is outside its allowed range ",
      substr(r$bounds, 1L, 1L), .sg_fmt(r$lower), ", ", .sg_fmt(r$upper),
      substr(r$bounds, 2L, 2L), " (", r$units, "). ", r$description
    )
  }
  x
}

.sg_cross_check_params <- function(v) {
  if (v$dev_temp_decline <= v$dev_temp_min) {
    .sg_stop(
      "dev_temp_decline (", v$dev_temp_decline, ") must be greater than dev_temp_min (",
      v$dev_temp_min, ")."
    )
  }
  if (v$oviposition_temp_max <= v$dev_temp_min) {
    .sg_stop(
      "oviposition_temp_max (", v$oviposition_temp_max, ") must be greater than dev_temp_min (",
      v$dev_temp_min, ")."
    )
  }
  for (nm in c("oviposition_age_dd", "female_receptive_age_dd", "male_mature_age_dd")) {
    if (v[[nm]] >= v$adult_lifespan_dd) {
      .sg_stop(
        nm, " (", v[[nm]], " dd) must be shorter than the adult lifespan adult_lifespan_dd (",
        v$adult_lifespan_dd, " dd); otherwise no wild-type adult reaches that age."
      )
    }
  }
  if (v$larval_crowding_mode == "none" && v$larval_crowding_coef > 0) {
    .sg_warn(
      "larval_crowding_coef = ", v$larval_crowding_coef,
      " has no effect unless larval_crowding_mode is \"per_wound\" or ",
      "\"per_infested_host\"."
    )
  }
  .sg_cross_check_hosts(v)
  if (v$oviposition_mode == "eggs" && !is.na(v$clutch_size)) {
    .sg_warn("clutch_size has no effect unless oviposition_mode = \"batches\".")
  }
  if (v$mating_mode == "instant" && v$female_receptive_age_dd > 0) {
    .sg_warn(
      "With mating_mode = \"instant\" females mate at emergence; ",
      "female_receptive_age_dd is ignored."
    )
  }
  invisible(TRUE)
}

# Host-wound model (D-030): required values in host mode, crowding modes that
# need a given wound model, and host settings that have no effect.
.sg_cross_check_hosts <- function(v) {
  hosts <- v$wound_model == "hosts"
  if (hosts) {
    missing <- c("host_density", "wounded_fraction")[
      c(is.na(v$host_density), is.na(v$wounded_fraction))
    ]
    if (length(missing) > 0L) {
      .sg_stop(
        "wound_model = \"hosts\" needs ", paste(missing, collapse = " and "),
        " (no defaults). Set them from the region, e.g. a livestock census for ",
        "host_density; host_density * wounded_fraction = 5 per km^2 matches the ",
        "default wound_density."
      )
    }
    if (v$larval_crowding_mode == "per_wound") {
      .sg_stop(
        "larval_crowding_mode = \"per_wound\" uses the fixed wounds of wound_model = ",
        "\"static\"; with wound_model = \"hosts\" use \"per_infested_host\"."
      )
    }
    if (!identical(v$wound_density, .sg_param_defaults()$wound_density)) {
      .sg_warn(
        "wound_density has no effect with wound_model = \"hosts\"; the wounds come ",
        "from host_density and wounded_fraction."
      )
    }
    return(invisible(TRUE))
  }
  if (v$larval_crowding_mode == "per_infested_host") {
    .sg_stop("larval_crowding_mode = \"per_infested_host\" needs wound_model = \"hosts\".")
  }
  defaults <- .sg_param_defaults()
  host_par <- c(
    "host_density", "wounded_fraction", "wound_healing_days", "infestation_days",
    "infestation_substages", "infested_wound_attractiveness", "treatment_rate"
  )
  changed <- host_par[!vapply(host_par, function(nm) {
    isTRUE(all.equal(v[[nm]], defaults[[nm]])) || (is.na(v[[nm]]) && is.na(defaults[[nm]]))
  }, TRUE)]
  if (length(changed) > 0L) {
    .sg_warn(
      paste(changed, collapse = ", "), " has no effect with wound_model = \"static\"; ",
      "set wound_model = \"hosts\" to use the host-wound model."
    )
  }
  invisible(TRUE)
}

#' Change model parameters
#'
#' Returns a copy of an `sg_params` object with some values changed. The
#' result is validated again, so invalid combinations are caught.
#'
#' @param object An `sg_params` object.
#' @param ... Named parameter values to change.
#' @return An `sg_params` object with the same preset attribute.
#' @export
#' @examples
#' p <- sg_params()
#' p2 <- update(p, eggs_per_day_max = 67, adult_substages = 8)
#' c(p$eggs_per_day_max, p2$eggs_per_day_max)
update.sg_params <- function(object, ...) {
  values <- .sg_apply_overrides(unclass(object), list(...))
  attr(values, "preset") <- NULL
  values <- .sg_validate_params(values)
  structure(values, class = "sg_params", preset = attr(object, "preset"))
}

# `row.names` matches the signature of the base::as.data.frame() generic.
#' @export
# nolint start: object_name_linter.
as.data.frame.sg_params <- function(x, row.names = NULL, optional = FALSE, ...) {
  # nolint end
  tab <- sg_param_table()
  value <- vapply(tab$name, function(nm) .sg_fmt(x[[nm]]), "")
  out <- data.frame(
    name = tab$name,
    value = unname(value),
    default = tab$default,
    changed = unname(value) != tab$default,
    units = tab$units,
    group = tab$group,
    source = tab$source,
    stringsAsFactors = FALSE
  )
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}

#' @export
print.sg_params <- function(x, ...) {
  df <- as.data.frame(x)
  cat(
    "<sg_params> preset \"", attr(x, "preset"), "\": ", nrow(df),
    " parameters, ", sum(df$changed), " differ from the defaults (*)\n",
    sep = ""
  )
  width <- max(nchar(df$name)) # one name column across all groups
  for (g in .sg_param_groups) {
    d <- df[df$group == g, , drop = FALSE]
    if (nrow(d) == 0L) next
    cat("\n", g, "\n", sep = "")
    lines <- sprintf(
      "  %s %-*s %-20s %-21s %s",
      ifelse(d$changed, "*", " "), width, d$name, d$value, d$units, d$source
    )
    cat(lines, sep = "\n")
  }
  invisible(x)
}
