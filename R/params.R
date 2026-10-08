#' Model parameters
#'
#' Creates a validated set of SGDrivE model parameters. Every parameter has a
#' default (see [sg_param_table()] for the full list with units, ranges and
#' sources); arguments given in `...` override the defaults.
#'
#' @param ... Named parameter values that override the defaults, e.g.
#'   `n_A = 8, R = 67`. Names must match `sg_param_table()$name`. To pass a
#'   list `x` of values, use `do.call(sg_params, x)`.
#' @param preset Which set of defaults to start from:
#'   \describe{
#'     \item{`"default"`}{SGDrivE defaults: Gutierrez & Ponti (2014) biology
#'       with speed-tuned sub-stage numbers (at most 10 per stage), mating at
#'       rate `nu`, and searching mated females as competitors for wounds.}
#'     \item{`"gp2014"`}{The configuration of the Gutierrez & Ponti (2014)
#'       Pascal code, for regression tests and benchmarking: 25 sub-stages per
#'       stage, egg laying from adult class 8 (`X_ovi` = 7/25 of `X_A`), the
#'       0.61 dd per day floor on
#'       development, all adults competing for wounds, and instant mating (the
#'       2014 model has no mating step).}
#'   }
#'   Arguments in `...` are applied on top of the preset.
#'
#' @details
#' Each value is checked against its type and allowed range. Cross-checks
#' then catch combinations that cannot work, for example an age at first
#' oviposition `X_ovi` longer than the adult lifespan `X_A`. Settings that
#' are valid but have no effect (such as `gamma_L > 0` without larval density
#' dependence) give a warning.
#'
#' Units follow the conventions of the package: temperatures in degrees
#' Celsius, durations in degree-days (`dd`) or days (`d`), rates per day.
#' Parameters given per km^2 (`H`, `h_M`) are converted to counts per node
#' when the model is assembled.
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
#' p2 <- sg_params(n_A = 8, nu = 1)
#' p2$n_A
#'
#' # Configuration of the Gutierrez & Ponti (2014) code
#' gp <- sg_params(preset = "gp2014")
#' gp$n_A
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
      n_E = 25L, n_L = 25L, n_P = 25L, n_A = 25L,
      # laying from class 8 of 25: 7 classes of 202.8 / 25 dd before laying
      X_ovi = 7 / 25 * 202.8,
      dx_floor = 0.01,
      search_denominator = "all_adults",
      mating = "instant"
    )
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
      "All parameters must be named, e.g. sg_params(n_A = 8). ",
      "To pass a list x, use do.call(sg_params, x)."
    )
  }
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0L) {
    .sg_stop("Parameter(s) given more than once: ", paste(dup, collapse = ", "), ".")
  }
  unknown <- setdiff(nms, names(values))
  if (length(unknown) > 0L) {
    msgs <- vapply(unknown, function(u) {
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
  if (v$T_U <= v$T_L) {
    .sg_stop("T_U (", v$T_U, ") must be greater than T_L (", v$T_L, ").")
  }
  if (v$TU_ovip <= v$T_L) {
    .sg_stop("TU_ovip (", v$TU_ovip, ") must be greater than T_L (", v$T_L, ").")
  }
  for (nm in c("X_ovi", "X_rec", "X_mat")) {
    if (v[[nm]] >= v$X_A) {
      .sg_stop(
        nm, " (", v[[nm]], " dd) must be shorter than the adult lifespan X_A (",
        v$X_A, " dd); otherwise no wild-type adult reaches that age."
      )
    }
  }
  if (v$larval_dd == "none" && v$gamma_L > 0) {
    .sg_warn("gamma_L = ", v$gamma_L, " has no effect unless larval_dd = \"per_wound\".")
  }
  if (v$mating == "instant" && v$X_rec > 0) {
    .sg_warn("With mating = \"instant\" females mate at emergence; X_rec is ignored.")
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
#' p2 <- update(p, R = 67, n_A = 8)
#' c(p$R, p2$R)
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
  for (g in .sg_param_groups) {
    d <- df[df$group == g, , drop = FALSE]
    if (nrow(d) == 0L) next
    cat("\n", g, "\n", sep = "")
    lines <- sprintf(
      "  %s %-19s %-20s %-21s %s",
      ifelse(d$changed, "*", " "), d$name, d$value, d$units, d$source
    )
    cat(lines, sep = "\n")
  }
  invisible(x)
}
