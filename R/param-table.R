# Parameter table: the single source of truth for every model parameter
# (D-017). Defaults, allowed ranges, units, sources and the equation of
# dev/MODEL_SPEC.md that uses each parameter are defined here once; sg_params(),
# its validation, print() and the documentation all read this table.
#
# Strings are ASCII only (R CMD check portability): "degC" for degrees Celsius,
# "dd" for degree-days, "d^-1" for per day, "km^-2" for per square kilometre.
#
# Spec tags H1-H3 refer to the host-wound model (dev/HOST_WOUND_DESIGN.md,
# section 3, D-030): H1 host states and equations, H2 coupling to the flies,
# H3 chain-binomial host step. Tags N1-N5 refer to several nodes and adult
# movement (MODEL_SPEC section 8a, D-032).
#
# Source tags: GP2014 = Gutierrez & Ponti (2014) paper or Pascal code;
# G2019 = Gutierrez, Ponti & Arias (2019); MGDrivE2 = MGDrivE2 lifecycle model;
# adapted = new for ScrewDrivE (decision number given in the description).

# One row of the table. `default` keeps its exact R value (e.g. log(2)).
# `allow_inf` lets a real parameter take the value Inf (its range must then
# include Inf, i.e. upper = Inf with a closed upper bound).
.prow <- function(name, default, group, units, source, spec, description,
                  type = "real", lower = -Inf, upper = Inf, bounds = "[]",
                  choices = NA_character_, allow_na = FALSE, symbol = name,
                  allow_inf = FALSE) {
  list(
    name = name, symbol = symbol, group = group, type = type,
    default = default, lower = lower, upper = upper, bounds = bounds,
    choices = choices, allow_na = allow_na, allow_inf = allow_inf, units = units,
    source = source, spec = spec, description = description
  )
}

# Display order of groups in print() and in the documentation.
.sg_param_groups <- c(
  "temperature", "development", "stages", "thresholds", "mortality",
  "reproduction", "mating", "hosts", "movement", "numerical"
)

.sg_param_rows <- function() {
  list(
    # ---- temperature: Gutierrez & Ponti (2014) temperature functions ------
    .prow("dev_temp_min", 14.5, "temperature", "degC", "GP2014", "S1, S4",
      "Lower developmental threshold temperature.",
      lower = -10, upper = 40, symbol = "T_L"
    ),
    .prow("dev_temp_decline", 40.5, "temperature", "degC", "GP2014", "S1",
      paste(
        "Temperature at which the high-temperature decline of the",
        "development rate is centred (40 in the paper, 40.5 in the code)."
      ),
      lower = 20, upper = 60, symbol = "T_U"
    ),
    .prow("dev_rate_slope", 0.00065, "temperature", "h^-1 degC^-1", "GP2014", "S1",
      "Initial slope of the development-rate curve r(T).",
      lower = 0, bounds = "(]", symbol = "a"
    ),
    .prow("dev_decline_steepness", 1.75, "temperature", "-", "GP2014", "S1",
      "Base of the high-temperature decline term of r(T).",
      lower = 1, bounds = "()", symbol = "b"
    ),
    .prow("dev_rate_floor", 0, "temperature", "-", "GP2014", "S1",
      paste(
        "Floor on the daily development fraction 24 r(T). 0 stops",
        "development below dev_temp_min; 0.01 reproduces the Pascal code's minimum of",
        "0.61 dd per day."
      ),
      lower = 0, upper = 1, symbol = "dx_floor"
    ),
    .prow("optimal_survival_temp", 27.2, "temperature", "degC", "GP2014", "S2, S3",
      "Temperature of minimum pupal/adult mortality and egg hatch failure.",
      lower = 0, upper = 50, symbol = "T_opt"
    ),
    .prow("oviposition_temp_max", 43.5, "temperature", "degC", "GP2014", "S4",
      "Upper temperature limit for oviposition.",
      lower = 20, upper = 60, symbol = "TU_ovip"
    ),
    .prow("egg_temp_source", "ambient", "temperature", "-", "adapted", "S5",
      paste(
        "Temperature used to age eggs: \"ambient\" (Gutierrez & Ponti code)",
        "or \"host\" (constant host_body_temp). D-014."
      ),
      type = "choice", choices = "ambient|host", symbol = "egg_temperature"
    ),
    .prow("host_body_temp", 38.5, "temperature", "degC", "GP2014", "S5",
      "Host body temperature, used for eggs when egg_temp_source = \"host\".",
      lower = 20, upper = 50, symbol = "T_host"
    ),

    # ---- development: mean stage durations ---------------------------------
    .prow("dd_scale_preoviposition", 61, "development", "dd", "GP2014", "S1",
      paste(
        "Pre-oviposition period of Adams (1979), used to scale the daily",
        "degree-day increment dx. The egg-laying age itself is oviposition_age_dd."
      ),
      lower = 0, bounds = "(]", symbol = "X_pre"
    ),
    .prow("egg_duration_dd", 10, "development", "dd", "GP2014", "S5",
      "Mean egg stage duration (Pascal code; 8.64 dd in the paper).",
      lower = 0, bounds = "(]", symbol = "X_E"
    ),
    .prow("larva_duration_days", 6, "development", "d", "GP2014", "S5",
      "Mean larval stage duration in calendar days (larvae develop in hosts).",
      lower = 0, bounds = "(]", symbol = "tau_L"
    ),
    .prow("pupa_duration_dd", 95, "development", "dd", "GP2014", "S5",
      "Mean pupal stage duration (6-10 d, mean 7.5 d at 27.5 degC).",
      lower = 0, bounds = "(]", symbol = "X_P"
    ),
    .prow("adult_lifespan_dd", 202.8, "development", "dd", "GP2014", "S5, S6",
      paste(
        "Mean adult physiological lifespan of the wild type; adults leaving",
        "the last age class die of old age. Scaled by l_g per genotype (D-002)."
      ),
      lower = 0, bounds = "(]", symbol = "X_A"
    ),

    # ---- stages: Erlang sub-stage numbers (D-003: defaults <= 10) ----------
    .prow("egg_substages", 3L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for eggs. D-003.",
      type = "integer", lower = 1, upper = 1000, symbol = "n_E"
    ),
    .prow("larva_substages", 6L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for larvae. D-003.",
      type = "integer", lower = 1, upper = 1000, symbol = "n_L"
    ),
    .prow("pupa_substages", 10L, "stages", "-", "adapted", "S5",
      "Number of Erlang sub-stages for pupae. D-003.",
      type = "integer", lower = 1, upper = 1000, symbol = "n_P"
    ),
    .prow("adult_substages", 10L, "stages", "-", "adapted", "S5, S6",
      paste(
        "Number of adult sub-stages per sex, allocated to the segments of the",
        "adult chain (S6). D-003, D-018."
      ),
      type = "integer", lower = 1, upper = 1000, symbol = "n_A"
    ),

    # ---- thresholds: maturity ages in absolute degree-days (D-002, D-004) --
    .prow("oviposition_age_dd", 61, "thresholds", "dd", "GP2014", "S6, S8",
      paste(
        "Age at first oviposition: mated females older than this lay eggs",
        "and search for wounds (segment boundary of the female chain). D-018."
      ),
      lower = 0, bounds = "(]", symbol = "X_ovi"
    ),
    .prow("female_receptive_age_dd", 0, "thresholds", "dd", "adapted", "S6, S12",
      "Female receptivity age: virgins younger than this cannot mate. D-004, D-018.",
      lower = 0, symbol = "X_rec"
    ),
    .prow("male_mature_age_dd", 0, "thresholds", "dd", "adapted", "S6, S12",
      "Male maturity age: younger males are not available as mates. D-004, D-018.",
      lower = 0, symbol = "X_mat"
    ),

    # ---- mortality ----------------------------------------------------------
    .prow("adult_mortality_curvature", 0.00036, "mortality", "degC^-2", "GP2014", "S2",
      "Curvature of the daily pupal/adult mortality curve p_A(T).",
      lower = 0, symbol = "c_A"
    ),
    .prow("adult_mortality_min", 0.0035, "mortality", "-", "GP2014", "S2",
      "Minimum daily pupal/adult mortality (at optimal_survival_temp).",
      lower = 0, upper = 1, symbol = "p_A0"
    ),
    .prow("egg_failure_curvature", 0.00375, "mortality", "degC^-2", "GP2014", "S3",
      "Curvature of the egg hatch-failure curve.",
      lower = 0, symbol = "c_E"
    ),
    .prow("egg_failure_min", 0.3, "mortality", "-", "GP2014", "S3",
      "Minimum egg hatch failure (at optimal_survival_temp).",
      lower = 0, upper = 1, symbol = "p_E0"
    ),
    .prow("larva_survival_daily", 0.95, "mortality", "-", "GP2014", "S7",
      "Daily larval survival (density independent).",
      lower = 0, upper = 1, bounds = "(]", symbol = "s_L"
    ),
    .prow("female_background_mortality", 0, "mortality", "d^-1", "adapted", "S7",
      "Background (field) mortality of adult females, added to mu_T.",
      lower = 0, symbol = "mu_b,F"
    ),
    .prow("male_background_mortality", 0, "mortality", "d^-1", "adapted", "S7",
      "Background (field) mortality of adult males, added to mu_T.",
      lower = 0, symbol = "mu_b,M"
    ),
    .prow("larval_crowding_mode", "none", "mortality", "-", "adapted", "S7",
      paste(
        "Larval density dependence: \"none\", \"per_wound\" (larval hazard +",
        "larval_crowding_coef * larvae / wounds) or \"per_infested_host\" (+",
        "larval_crowding_coef * larvae / infested hosts; wound_model = \"hosts\" only,",
        "D-030)."
      ),
      type = "choice", choices = "none|per_wound|per_infested_host", symbol = "larval_dd"
    ),
    .prow("larval_crowding_coef", 0, "mortality", "d^-1 wound larva^-1", "adapted", "S7, H2",
      "Larval crowding coefficient, used when larval_crowding_mode is not \"none\".",
      lower = 0, symbol = "gamma_L"
    ),
    .prow("pupa_drought_mortality", 0, "mortality", "d^-1", "adapted", "S7",
      "Extra pupal mortality under full moisture stress (rain index 0).",
      lower = 0, symbol = "mu_W"
    ),
    .prow("drought_free_rainfall", 450, "mortality", "mm yr^-1", "GP2014", "S7",
      "Annual rainfall above which there is no moisture stress.",
      lower = 0, bounds = "(]", symbol = "rain_threshold"
    ),
    .prow("survival_floor", 1e-5, "mortality", "-", "GP2014", "S2, S3",
      "Lower clamp on daily survival before converting to a hazard.",
      lower = 0, upper = 1, bounds = "()", symbol = "epsilon_s"
    ),

    # ---- reproduction -------------------------------------------------------
    .prow("eggs_per_day_max", 90, "reproduction", "eggs female^-1 d^-1", "GP2014", "S9",
      "Maximum eggs per egg-laying female per day (at about 29 degC).",
      lower = 0, symbol = "R"
    ),
    .prow("wound_search_rate", 0.001, "reproduction", "dd^-1", "GP2014", "S8",
      "Wound-search constant.",
      lower = 0, symbol = "a_s"
    ),
    .prow("wound_density", 5, "reproduction", "wounds km^-2", "GP2014", "S8",
      paste(
        "Density of wounds available for oviposition; converted to a count",
        "per node with the node area. The default is the relative value of",
        "Gutierrez & Ponti (2014); calibrate for absolute densities."
      ),
      lower = 0, symbol = "H"
    ),
    .prow("wound_finding_max", 1, "reproduction", "-", "adapted", "S8",
      "Baseline wound-finding efficiency at low density.",
      lower = 0, upper = 1, symbol = "sigma_0"
    ),
    .prow("wound_competitors", "searching_females", "reproduction", "-",
      "adapted", "S8",
      paste(
        "Who competes for wounds: \"searching_females\" (egg-laying mated",
        "females, default) or \"all_adults\" (Gutierrez & Ponti code). D-006."
      ),
      type = "choice", choices = "searching_females|all_adults", symbol = "search_denominator"
    ),
    .prow("oviposition_mode", "batches", "reproduction", "-", "adapted", "S17",
      paste(
        "Egg laying in the chain-binomial solver: \"batches\" (default;",
        "Poisson numbers of oviposition events, each of clutch_size eggs on",
        "average) or \"eggs\" (Poisson numbers of single eggs). Same mean;",
        "batches give the larger variance of real egg masses. The ODE is not",
        "affected. D-024."
      ),
      type = "choice", choices = "batches|eggs", symbol = "oviposition"
    ),
    .prow("clutch_size", NA_real_, "reproduction", "eggs", "adapted", "S17",
      paste(
        "Mean eggs per oviposition event. NA (default): eggs_per_day_max * phi_ov,",
        "the eggs a female lays in one day at the day's temperature (one wound found",
        "per day, as in Gutierrez & Ponti). When set, events occur at",
        "eggs_per_day_max * phi_ov * Psi / clutch_size per egg-laying female per day, so",
        "mean egg production is unchanged. D-024."
      ),
      lower = 0, bounds = "(]", allow_na = TRUE
    ),

    # ---- mating -------------------------------------------------------------
    .prow("mating_rate", log(2), "mating", "d^-1", "G2019", "S12",
      "Mating rate of receptive virgins (half mate per day). D-004.",
      lower = 0, symbol = "nu"
    ),
    .prow("male_density_half_mating", 0, "mating", "males km^-2", "adapted", "S12",
      "Half-saturation male density for mate finding; 0 = no Allee effect.",
      lower = 0, symbol = "h_M"
    ),
    .prow("mating_mode", "rate", "mating", "-", "MGDrivE2", "S12, S14",
      paste(
        "\"rate\": virgins mate at rate nu; \"instant\": females mate at",
        "emergence when males are present (MGDrivE2)."
      ),
      type = "choice", choices = "rate|instant", symbol = "mating"
    ),

    # ---- hosts: host-wound model (D-030) -------------------------------------
    .prow("wound_model", "static", "hosts", "-", "adapted", "S8, H1",
      paste(
        "\"static\": a fixed number of wounds (wound_density), every fertile egg",
        "batch counts as a case (Gutierrez & Ponti). \"hosts\": hosts without a",
        "wound, wounded and infested, with cases as new infestations of hosts, so",
        "at most host_density * area hosts are infested at once. D-030."
      ),
      type = "choice", choices = "static|hosts"
    ),
    .prow("host_density", NA_real_, "hosts", "animals km^-2", "adapted", "H1",
      paste(
        "Density of hosts (animals). No default: required when wound_model =",
        "\"hosts\", e.g. from a livestock census."
      ),
      lower = 0, bounds = "(]", allow_na = TRUE, symbol = "N_H"
    ),
    .prow("wounded_fraction", NA_real_, "hosts", "-", "adapted", "H1",
      paste(
        "Fraction of hosts with an attractive wound when there are no flies.",
        "No default: required when wound_model = \"hosts\";",
        "host_density * wounded_fraction plays the role of wound_density."
      ),
      lower = 0, upper = 1, bounds = "()", allow_na = TRUE, symbol = "p_0"
    ),
    .prow("wound_healing_days", 7, "hosts", "d", "adapted", "H1",
      paste(
        "Mean time an uninfested wound stays attractive before it heals.",
        "Assumption: no published value found."
      ),
      lower = 0, bounds = "(]", symbol = "1/gamma_W"
    ),
    .prow("infestation_days", 6, "hosts", "d", "adapted", "H1",
      paste(
        "Mean duration of an untreated infestation: eggs hatch within about a",
        "day and larvae feed for 5-7 days (WOAH; CFSPH; FAO)."
      ),
      lower = 0, bounds = "(]", symbol = "n_I/q_I"
    ),
    .prow("infestation_substages", 6L, "hosts", "-", "adapted", "H1",
      "Erlang sub-stages of an infestation (coefficient of variation 1/sqrt(n)).",
      type = "integer", lower = 1, upper = 10, symbol = "n_I"
    ),
    .prow("infested_wound_attractiveness", 1, "hosts", "-", "adapted", "H2",
      paste(
        "Attractiveness of an infested wound to egg-laying females relative to",
        "an uninfested one; 1 = neutral. Infested wounds attract further females",
        "(WOAH; CFSPH), but no value is known."
      ),
      lower = 0, symbol = "alpha"
    ),
    .prow("treatment_rate", 0, "hosts", "d^-1", "adapted", "H1, H2",
      paste(
        "Rate at which infested hosts are treated; treatment ends the",
        "infestation and kills its larvae."
      ),
      lower = 0, symbol = "kappa"
    ),

    # ---- movement: several nodes (MODEL_SPEC section 8a, D-032) -------------
    .prow("movement_t_ref", 27, "movement", "degC", "adapted", "N3",
      paste(
        "Reference temperature of the movement rates: leave rates are given per",
        "day at this temperature, and the physiological clock scales them by",
        "dx(T) / dx(movement_t_ref) (D-032). Near the survival optimum and the",
        "temperatures at which field leave rates were measured; not a biological",
        "limit."
      ),
      lower = 0, upper = 45, symbol = "T_ref"
    ),
    .prow("movement_kappa_max", 1, "movement", "-", "adapted", "N3",
      paste(
        "Cap on the physiological movement clock kappa = min(kappa_max,",
        "dx(T) / dx(T_ref)): 1 = the given rate at and above T_ref (capped),",
        "Inf = uncapped degree-day scaling (D-032). No study supports either;",
        "compare both."
      ),
      lower = 0, bounds = "(]", symbol = "kappa_max", allow_inf = TRUE
    ),

    # ---- numerical ----------------------------------------------------------
    .prow("extinction_threshold", 1, "numerical", "individuals", "adapted", "S10",
      "Extinction threshold for all wild stages combined.",
      lower = 0, symbol = "ext_threshold"
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
#' Returns the table that defines every ScrewDrivE model parameter: its name,
#' symbol, default value, allowed range, units, source and the equation of the
#' model specification that uses it. This table is the single source for the
#' defaults and checks in [sg_params()].
#'
#' @return A data frame with one row per parameter and columns:
#'   \describe{
#'     \item{`name`}{Argument name used in [sg_params()].}
#'     \item{`symbol`}{Mathematical symbol in the model specification.}
#'     \item{`group`}{Parameter group (temperature, development, stages,
#'       thresholds, mortality, reproduction, mating, hosts, movement,
#'       numerical).}
#'     \item{`type`}{`"real"`, `"integer"` or `"choice"`.}
#'     \item{`default`}{Default value, as text.}
#'     \item{`lower`, `upper`}{Allowed range (numeric parameters).}
#'     \item{`bounds`}{Whether the range ends are included: `"[]"`, `"(]"`,
#'       `"[)"` or `"()"`.}
#'     \item{`choices`}{Allowed values of choice parameters, separated by
#'       `"|"`.}
#'     \item{`allow_na`}{Whether `NA` is allowed (meaning "derive
#'       automatically", or for the host parameters "not set").}
#'     \item{`allow_inf`}{Whether `Inf` is allowed (only
#'       `movement_kappa_max`, where it means uncapped).}
#'     \item{`units`}{Units: `degC` degrees Celsius, `dd` degree-days, `d`
#'       days, `d^-1` per day, `km^-2` per square kilometre.}
#'     \item{`source`}{`GP2014` (Gutierrez & Ponti 2014), `G2019`
#'       (Gutierrez et al. 2019), `MGDrivE2`, or `adapted` (new in ScrewDrivE).}
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
    allow_inf = col("allow_inf"),
    units = col("units"),
    source = col("source"),
    spec = col("spec"),
    description = col("description"),
    stringsAsFactors = FALSE
  )
}
