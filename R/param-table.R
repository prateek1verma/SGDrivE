# Parameter table: the single source of truth for every model parameter
# (D-017). Defaults, allowed ranges, units, sources and the equation of
# dev/MODEL_SPEC.md that uses each parameter are defined here once; sg_params(),
# its validation, print() and the documentation all read this table.
#
# Strings are ASCII only (R CMD check portability): "degC" for degrees Celsius,
# "dd" for degree-days, "d^-1" for per day, "km^-2" for per square kilometre.
#
# Source tags: GP2014 = Gutierrez & Ponti (2014) paper or Pascal code;
# G2019 = Gutierrez, Ponti & Arias (2019); MGDrivE2 = MGDrivE2 lifecycle model;
# adapted = new for SGDrivE (decision number given in the description).

# One row of the table. `default` keeps its exact R value (e.g. log(2)).
.prow <- function(name, default, group, units, source, spec, description,
                  type = "real", lower = -Inf, upper = Inf, bounds = "[]",
                  choices = NA_character_, allow_na = FALSE, symbol = name) {
  list(
    name = name, symbol = symbol, group = group, type = type,
    default = default, lower = lower, upper = upper, bounds = bounds,
    choices = choices, allow_na = allow_na, units = units, source = source,
    spec = spec, description = description
  )
}

# Display order of groups in print() and in the documentation.
.sg_param_groups <- c(
  "temperature", "development", "stages", "thresholds", "mortality",
  "reproduction", "mating", "numerical"
)

.sg_param_rows <- function() {
  list(
    # ---- temperature: Gutierrez & Ponti (2014) temperature functions ------
    .prow("T_L", 14.5, "temperature", "degC", "GP2014", "S1, S4",
      "Lower developmental threshold temperature.",
      lower = -10, upper = 40
    ),
    .prow("T_U", 40.5, "temperature", "degC", "GP2014", "S1",
      paste(
        "Temperature at which the high-temperature decline of the",
        "development rate is centred (40 in the paper, 40.5 in the code)."
      ),
      lower = 20, upper = 60
    ),
    .prow("a_dev", 0.00065, "temperature", "h^-1 degC^-1", "GP2014", "S1",
      "Initial slope of the development-rate curve r(T).",
      lower = 0, bounds = "(]", symbol = "a"
    ),
    .prow("b_dev", 1.75, "temperature", "-", "GP2014", "S1",
      "Base of the high-temperature decline term of r(T).",
      lower = 1, bounds = "()", symbol = "b"
    ),
    .prow("dx_floor", 0, "temperature", "-", "GP2014", "S1",
      paste(
        "Floor on the daily development fraction 24 r(T). 0 stops",
        "development below T_L; 0.01 reproduces the Pascal code's minimum of",
        "0.61 dd per day."
      ),
      lower = 0, upper = 1
    ),
    .prow("T_opt", 27.2, "temperature", "degC", "GP2014", "S2, S3",
      "Temperature of minimum pupal/adult mortality and egg hatch failure.",
      lower = 0, upper = 50
    ),
    .prow("TU_ovip", 43.5, "temperature", "degC", "GP2014", "S4",
      "Upper temperature limit for oviposition.",
      lower = 20, upper = 60
    ),
    .prow("egg_temperature", "ambient", "temperature", "-", "adapted", "S5",
      paste(
        "Temperature used to age eggs: \"ambient\" (Gutierrez & Ponti code)",
        "or \"host\" (constant T_host). D-014."
      ),
      type = "choice", choices = "ambient|host"
    ),
    .prow("T_host", 38.5, "temperature", "degC", "GP2014", "S5",
      "Host body temperature, used for eggs when egg_temperature = \"host\".",
      lower = 20, upper = 50
    ),

    # ---- development: mean stage durations ---------------------------------
    .prow("X_pre", 61, "development", "dd", "GP2014", "S1, S6",
      paste(
        "Pre-oviposition period. Scales the degree-day increment and sets",
        "the age at which females start laying eggs."
      ),
      lower = 0, bounds = "(]"
    ),
    .prow("X_E", 10, "development", "dd", "GP2014", "S5",
      "Mean egg stage duration (Pascal code; 8.64 dd in the paper).",
      lower = 0, bounds = "(]"
    ),
    .prow("tau_L", 6, "development", "d", "GP2014", "S5",
      "Mean larval stage duration in calendar days (larvae develop in hosts).",
      lower = 0, bounds = "(]", symbol = "tau_L"
    ),
    .prow("X_P", 95, "development", "dd", "GP2014", "S5",
      "Mean pupal stage duration (6-10 d, mean 7.5 d at 27.5 degC).",
      lower = 0, bounds = "(]"
    ),
    .prow("X_A", 202.8, "development", "dd", "GP2014", "S5, S6",
      paste(
        "Mean adult physiological lifespan of the wild type; adults leaving",
        "the last age class die of old age. Scaled by l_g per genotype (D-002)."
      ),
      lower = 0, bounds = "(]"
    ),

    # ---- stages: Erlang sub-stage numbers (D-003: defaults <= 10) ----------
    .prow("n_E", 3L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for eggs. D-003.",
      type = "integer", lower = 1, upper = 1000
    ),
    .prow("n_L", 6L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for larvae. D-003.",
      type = "integer", lower = 1, upper = 1000
    ),
    .prow("n_P", 10L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for pupae. D-003.",
      type = "integer", lower = 1, upper = 1000
    ),
    .prow("n_A", 10L, "stages", "-", "adapted", "S5, S6",
      "Number of adult age classes (Erlang sub-stages). D-003.",
      type = "integer", lower = 1, upper = 1000
    ),

    # ---- thresholds: maturity ages in absolute degree-days (D-002, D-004) --
    .prow("X_rec", 0, "thresholds", "dd", "adapted", "S6, S12",
      "Female receptivity age: virgins younger than this cannot mate. D-004.",
      lower = 0
    ),
    .prow("X_mat", 0, "thresholds", "dd", "adapted", "S6, S12",
      "Male maturity age: younger males are not available as mates. D-004.",
      lower = 0
    ),
    .prow("i0", NA_integer_, "thresholds", "-", "adapted", "S6",
      paste(
        "Override: first egg-laying adult class for every genotype. NA",
        "derives it from X_pre per genotype."
      ),
      type = "integer", lower = 1, upper = 1000, allow_na = TRUE
    ),
    .prow("i_m", NA_integer_, "thresholds", "-", "adapted", "S6",
      "Override: first receptive female class. NA derives it from X_rec.",
      type = "integer", lower = 1, upper = 1000, allow_na = TRUE
    ),
    .prow("i_mM", NA_integer_, "thresholds", "-", "adapted", "S6",
      "Override: first mature male class. NA derives it from X_mat.",
      type = "integer", lower = 1, upper = 1000, allow_na = TRUE
    ),

    # ---- mortality ----------------------------------------------------------
    .prow("c_A", 0.00036, "mortality", "degC^-2", "GP2014", "S2",
      "Curvature of the daily pupal/adult mortality curve p_A(T).",
      lower = 0
    ),
    .prow("p_A0", 0.0035, "mortality", "-", "GP2014", "S2",
      "Minimum daily pupal/adult mortality (at T_opt).",
      lower = 0, upper = 1
    ),
    .prow("c_E", 0.00375, "mortality", "degC^-2", "GP2014", "S3",
      "Curvature of the egg hatch-failure curve.",
      lower = 0
    ),
    .prow("p_E0", 0.3, "mortality", "-", "GP2014", "S3",
      "Minimum egg hatch failure (at T_opt).",
      lower = 0, upper = 1
    ),
    .prow("s_L", 0.95, "mortality", "-", "GP2014", "S7",
      "Daily larval survival (density independent).",
      lower = 0, upper = 1, bounds = "(]"
    ),
    .prow("mu_b_F", 0, "mortality", "d^-1", "adapted", "S7",
      "Background (field) mortality of adult females, added to mu_T.",
      lower = 0, symbol = "mu_b,F"
    ),
    .prow("mu_b_M", 0, "mortality", "d^-1", "adapted", "S7",
      "Background (field) mortality of adult males, added to mu_T.",
      lower = 0, symbol = "mu_b,M"
    ),
    .prow("larval_dd", "none", "mortality", "-", "adapted", "S7",
      paste(
        "Larval density dependence: \"none\" or \"per_wound\"",
        "(mu_L + gamma_L * L / H)."
      ),
      type = "choice", choices = "none|per_wound"
    ),
    .prow("gamma_L", 0, "mortality", "d^-1 wound larva^-1", "adapted", "S7",
      "Larval crowding coefficient, used when larval_dd = \"per_wound\".",
      lower = 0
    ),
    .prow("mu_W", 0, "mortality", "d^-1", "adapted", "S7",
      "Extra pupal mortality under full moisture stress (rain index 0).",
      lower = 0
    ),
    .prow("rain_threshold", 450, "mortality", "mm yr^-1", "GP2014", "S7",
      "Annual rainfall above which there is no moisture stress.",
      lower = 0, bounds = "(]"
    ),
    .prow("eps_s", 1e-5, "mortality", "-", "GP2014", "S2, S3",
      "Lower clamp on daily survival before converting to a hazard.",
      lower = 0, upper = 1, bounds = "()", symbol = "epsilon_s"
    ),

    # ---- reproduction -------------------------------------------------------
    .prow("R", 90, "reproduction", "eggs female^-1 d^-1", "GP2014", "S9",
      "Maximum eggs per egg-laying female per day (at about 29 degC).",
      lower = 0
    ),
    .prow("a_s", 0.001, "reproduction", "dd^-1", "GP2014", "S8",
      "Wound-search constant.",
      lower = 0
    ),
    .prow("H", 5, "reproduction", "wounds km^-2", "GP2014", "S8",
      paste(
        "Density of wounds available for oviposition; converted to a count",
        "per node with the node area. The default is the relative value of",
        "Gutierrez & Ponti (2014); calibrate for absolute densities."
      ),
      lower = 0
    ),
    .prow("sigma0", 1, "reproduction", "-", "adapted", "S8",
      "Baseline wound-finding efficiency at low density.",
      lower = 0, upper = 1, symbol = "sigma_0"
    ),
    .prow("search_denominator", "searching_females", "reproduction", "-",
      "adapted", "S8",
      paste(
        "Who competes for wounds: \"searching_females\" (egg-laying mated",
        "females, default) or \"all_adults\" (Gutierrez & Ponti code). D-006."
      ),
      type = "choice", choices = "searching_females|all_adults"
    ),

    # ---- mating -------------------------------------------------------------
    .prow("nu", log(2), "mating", "d^-1", "G2019", "S12",
      "Mating rate of receptive virgins (half mate per day). D-004.",
      lower = 0
    ),
    .prow("h_M", 0, "mating", "males km^-2", "adapted", "S12",
      "Half-saturation male density for mate finding; 0 = no Allee effect.",
      lower = 0
    ),
    .prow("mating", "rate", "mating", "-", "MGDrivE2", "S12, S14",
      paste(
        "\"rate\": virgins mate at rate nu; \"instant\": females mate at",
        "emergence when males are present (MGDrivE2)."
      ),
      type = "choice", choices = "rate|instant"
    ),

    # ---- numerical ----------------------------------------------------------
    .prow("ext_threshold", 1, "numerical", "individuals", "adapted", "S10",
      "Extinction threshold for all wild stages combined.",
      lower = 0
    )
  )
}

# Named list of exact default values, in table order.
.sg_param_defaults <- function() {
  rows <- .sg_param_rows()
  stats::setNames(lapply(rows, `[[`, "default"), vapply(rows, `[[`, "", "name"))
}

#' Parameter reference table
#'
#' Returns the table that defines every SGDrivE model parameter: its name,
#' symbol, default value, allowed range, units, source and the equation of the
#' model specification that uses it. This table is the single source for the
#' defaults and checks in [sg_params()].
#'
#' @return A data frame with one row per parameter and columns:
#'   \describe{
#'     \item{`name`}{Argument name used in [sg_params()].}
#'     \item{`symbol`}{Mathematical symbol in the model specification.}
#'     \item{`group`}{Parameter group (temperature, development, stages,
#'       thresholds, mortality, reproduction, mating, numerical).}
#'     \item{`type`}{`"real"`, `"integer"` or `"choice"`.}
#'     \item{`default`}{Default value, as text.}
#'     \item{`lower`, `upper`}{Allowed range (numeric parameters).}
#'     \item{`bounds`}{Whether the range ends are included: `"[]"`, `"(]"`,
#'       `"[)"` or `"()"`.}
#'     \item{`choices`}{Allowed values of choice parameters, separated by
#'       `"|"`.}
#'     \item{`allow_na`}{Whether `NA` is allowed (meaning "derive
#'       automatically").}
#'     \item{`units`}{Units: `degC` degrees Celsius, `dd` degree-days, `d`
#'       days, `d^-1` per day, `km^-2` per square kilometre.}
#'     \item{`source`}{`GP2014` (Gutierrez & Ponti 2014), `G2019`
#'       (Gutierrez et al. 2019), `MGDrivE2`, or `adapted` (new in SGDrivE).}
#'     \item{`spec`}{Equation(s) of the model specification using the
#'       parameter.}
#'     \item{`description`}{What the parameter does.}
#'   }
#' @seealso [sg_params()]
#' @export
#' @examples
#' tab <- sg_param_table()
#' tab[tab$group == "mating", c("name", "default", "units", "description")]
sg_param_table <- function() {
  rows <- .sg_param_rows()
  col <- function(field) vapply(rows, function(r) r[[field]], rows[[1L]][[field]])
  data.frame(
    name = col("name"),
    symbol = col("symbol"),
    group = col("group"),
    type = col("type"),
    default = vapply(rows, function(r) .sg_fmt(r$default), ""),
    lower = col("lower"),
    upper = col("upper"),
    bounds = col("bounds"),
    choices = col("choices"),
    allow_na = col("allow_na"),
    units = col("units"),
    source = col("source"),
    spec = col("spec"),
    description = col("description"),
    stringsAsFactors = FALSE
  )
}
