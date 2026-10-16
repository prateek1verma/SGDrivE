# Model checks (D-026): R0 by temperature, flow balance of simulations, and
# check_model().

#' Eggs per egg at low density (R0)
#'
#' The basic reproduction number of a wild-type population at constant
#' temperature: the expected number of eggs laid by the female descendants of
#' one egg when the population is small (wound search at its low-density
#' value, no larval crowding, mating as soon as males are present). The
#' population can grow at a temperature only where R0 > 1.
#'
#' @param model An `sg_model` object ([sg_model()]) with a single wild-type
#'   genotype that produces both sexes and breeds true.
#' @param temperature Constant temperature(s) in degC.
#' @return A numeric vector of R0 values, one per temperature.
#' @details The same quantity defines the constant-temperature equilibrium
#'   ([sg_init()] with `method = "equilibrium"`, model specification section
#'   11). The Allee effect (`male_density_half_mating > 0`) is not included: with it, mating fails
#'   at very low densities.
#' @seealso [check_model()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-01-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   temps <- seq(14, 40, by = 2)
#'   data.frame(temperature = temps, R0 = signif(sg_r0(mod, temps), 3))
#' }
sg_r0 <- function(model, temperature) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  if (!is.numeric(temperature) || length(temperature) < 1L || any(!is.finite(temperature))) {
    .sg_stop("`temperature` must be finite numbers (degC).")
  }
  w <- .sg_wild_setup(model, "sg_r0()")
  vapply(temperature, function(temp) {
    .sg_r0_of(w, .sg_constant_rates(w$s, .sg_constant_forcing(model$params, temp)))
  }, 0)
}

#' Flow balance of a simulation
#'
#' Checks that every change in the total number of individuals is explained
#' by the recorded flows: releases plus eggs laid minus deaths (by hazard, old
#' age and at emergence), between consecutive output days.
#'
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @return A data frame with columns `day` (end of each interval between
#'   output days), `replicate` (stochastic runs only), `change` (change in
#'   the total), `releases`, `eggs_laid`, `deaths` and `residual` (change
#'   minus releases and eggs plus deaths). The residual is about 1e-12 of
#'   the population for the ODE and exactly zero for stochastic runs.
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   max(abs(sg_flow_balance(sim)$residual))
#' }
sg_flow_balance <- function(sim) {
  .sg_sim_info(sim)
  a <- .sg_sim_array(sim)$array
  stochastic <- inherits(sim, "sg_sim_cb")
  k <- sim$counters
  rel <- sim$releases
  if (!stochastic) {
    dim(k) <- c(dim(k), 1L)
    dimnames(k) <- list(NULL, colnames(sim$counters), NULL)
    dim(rel) <- c(dim(rel), 1L)
  }
  cn <- dimnames(k)[[2L]]
  out_days <- sim$day[-1L]
  interval <- findInterval(seq_len(dim(k)[1L]), c(0, out_days) + 1L)
  rows <- lapply(seq_len(dim(a)[3L]), function(r) {
    total <- rowSums(matrix(a[, , r], dim(a)[1L]))
    per_day <- function(cols) rowSums(matrix(k[, cols, r], dim(k)[1L]))
    sum_by <- function(v) as.vector(tapply(v, interval, sum))
    eggs <- sum_by(per_day(grepl("^eggs_laid", cn)))
    deaths <- sum_by(per_day(grepl("^deaths_|^senescence_|^emergence_deaths$", cn)))
    released <- sum_by(rowSums(matrix(rel[, , r], dim(rel)[1L])))
    change <- diff(total)
    d <- data.frame(
      day = out_days, change = change, releases = released, eggs_laid = eggs,
      deaths = deaths, residual = change - released - eggs + deaths
    )
    if (stochastic) d <- cbind(d[1L], replicate = r, d[-1L])
    d
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Check a model before simulating
#'
#' Summarises a model's size, options and fastest processes, the season in
#' which a wild population can grow, possible configuration problems and,
#' for stochastic simulation, whether the population is large enough to
#' persist.
#'
#' @param model An `sg_model` object ([sg_model()]).
#' @param dt Time step of the chain-binomial solver (days), for the per-step
#'   mating probability and memory estimates.
#' @param scale If `TRUE` (default), estimate the smallest wild-type
#'   population over the forcing period with the ODE (from a one-year
#'   burn-in; needs a year of forcing), to judge stochastic extinction risk.
#' @return An object of class `sg_check` (printed as a report), a list with
#'   `size`, `options`, `rates`, `r0`, `scale`, `memory` and `warnings` (a
#'   character vector of the problems found).
#' @details
#' **Rates.** The largest total rate of any sub-stage (transfer plus
#' mortality, without larval crowding) over the forcing, and separately the
#' mating rate of receptive virgins, nu (O-007, D-026). Rates above about
#' 1000 per day can make the ODE stiff on some days (handled by the `lsoda`
#' fallback); the chain-binomial solver is not affected.
#'
#' **R0.** For a single true-breeding wild type, [sg_r0()] at each day's
#' forcing (as if it lasted): its range, the number of days with R0 > 1 and
#' the longest run of such days.
#'
#' **Warnings.** Genotypes that never lay eggs, never become receptive or
#' never mature as mates (their lifespan ends first, D-018); released
#' genotypes whose crosses with wild-type females give no viable eggs are
#' listed as sterile (expected for sterile insect releases); populations
#' below about ten wild-type individuals per node at their seasonal low.
#' @seealso [sg_r0()], [sg_flow_balance()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   check_model(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#' }
check_model <- function(model, dt = 0.1, scale = TRUE) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  .sg_require_static(model, "check_model()")
  .sg_cb_inputs(model, NULL, dt, stochastic = FALSE, n_days = 1L) # checks dt
  s <- .sg_reference_setup(model)
  p <- s$p
  f <- model$forcing
  cb <- model$cube
  warn <- character(0)

  rel <- model$releases
  size <- list(
    genotypes = s$n_g, states = s$n_states, days = f$n_days,
    period = range(f$dates), area = model$area, wounds = s$counts$H,
    release_events = nrow(rel), released = sum(rel$count),
    released_genotypes = unique(rel$genotype)
  )
  options <- p[c(
    "mating_mode", "oviposition_mode", "clutch_size", "wound_competitors",
    "larval_crowding_mode", "egg_temp_source", "male_density_half_mating"
  )]

  # ---- fastest processes (without crowding; mating separately) ------------
  day_max <- vapply(seq_len(f$n_days), function(d) {
    juv <- c(
      s$n_E * f$dx_egg[d] / p$egg_duration_dd + f$mu_E[d],
      s$n_L / p$larva_duration_days - log(p$larva_survival_daily),
      s$n_P * f$dx[d] / p$pupa_duration_dd + f$mu_T[d] +
        p$pupa_drought_mortality * (1 - f$moisture[d])
    )
    q_f <- (f$dx[d] / s$w_f)[is.finite(s$w_f)]
    q_m <- (f$dx[d] / s$w_m)[is.finite(s$w_m)]
    max(
      juv, max(q_f) + f$mu_T[d] + p$female_background_mortality,
      max(q_m) + f$mu_T[d] + p$male_background_mortality
    )
  }, 0)
  rates <- list(
    max_rate = max(day_max), max_rate_date = f$dates[which.max(day_max)],
    mating_rate = p$mating_rate, mating_probability_per_step = 1 - exp(-p$mating_rate * dt)
  )
  if (p$mating_rate >= 1000 || rates$max_rate >= 1000) {
    warn <- c(warn, sprintf(
      "Rates up to %.3g per day: the ODE may be stiff on some days (%s).",
      max(p$mating_rate, rates$max_rate), "handled by the lsoda fallback"
    ))
  }

  # ---- R0 over the forcing --------------------------------------------------
  r0 <- tryCatch(
    {
      w <- .sg_wild_setup(model, "R0")
      daily <- vapply(seq_len(f$n_days), function(d) {
        .sg_r0_of(w, .sg_constant_rates(s, list(
          dx = f$dx[d], dx_egg = f$dx_egg[d], mu_E = f$mu_E[d], mu_T = f$mu_T[d],
          phi_ov = f$phi_ov[d], moisture = f$moisture[d]
        )))
      }, 0)
      runs <- rle(daily > 1)
      list(
        available = TRUE, daily = daily, range = range(daily), days_above_1 = sum(daily > 1),
        longest_run = if (any(runs$values)) max(runs$lengths[runs$values]) else 0L,
        at_mean_temperature = sg_r0(model, mean(f$tmean))
      )
    },
    error = function(e) list(available = FALSE, note = conditionMessage(e))
  )
  if (isTRUE(r0$available) && r0$days_above_1 == 0) {
    warn <- c(warn, "R0 never exceeds 1 over the forcing: a wild population cannot grow.")
  }

  # ---- genotypes that cannot take part ---------------------------------------
  first_class <- function(flags) if (any(flags)) which(flags)[1L] else Inf
  females <- cb$genotypes[cb$phi > 0]
  males <- cb$genotypes[cb$phi < 1]
  never_lay <- females[s$last_f[females] < first_class(s$laying)]
  never_rec <- if (p$mating_mode == "rate") {
    females[s$last_f[females] < first_class(s$receptive)]
  } else {
    character(0)
  }
  never_mature <- males[s$last_m[males] < first_class(s$mature)]
  if (length(never_lay)) {
    warn <- c(warn, paste0(
      "Genotype(s) ", paste(never_lay, collapse = ", "),
      ": females die before they can lay eggs (lifespan shorter than oviposition_age_dd)."
    ))
  }
  if (length(never_rec)) {
    warn <- c(warn, paste0(
      "Genotype(s) ", paste(never_rec, collapse = ", "),
      ": females die before they become receptive (female_receptive_age_dd)."
    ))
  }
  if (length(never_mature)) {
    warn <- c(warn, paste0(
      "Genotype(s) ", paste(never_mature, collapse = ", "),
      ": males die before they mature as mates (male_mature_age_dd)."
    ))
  }
  wild_f <- intersect(cb$wild_type, females)
  sterile <- vapply(size$released_genotypes, function(g) {
    all(cb$egg_total[paste0(wild_f, "|", g)] == 0)
  }, TRUE)
  size$sterile_releases <- size$released_genotypes[sterile]

  # ---- population scale for stochastic simulation ----------------------------
  scale_info <- list(available = FALSE)
  if (scale && f$n_days >= 365L) {
    init <- suppressWarnings(sg_init(model, "burn_in"))
    sim <- simulate_ode(model, init = init)
    geno <- sub("^.\\[([^],|]+).*$", "\\1", colnames(sim$state))
    wild_total <- rowSums(sim$state[, geno %in% cb$wild_type, drop = FALSE])
    scale_info <- list(
      available = TRUE, min_wild = min(wild_total), max_wild = max(wild_total),
      min_wild_date = sim$date[which.min(wild_total)]
    )
    if (scale_info$min_wild < 10) {
      warn <- c(warn, sprintf(
        "The wild population falls to %.3g individuals per node (%s); stochastic runs %s",
        scale_info$min_wild, format(scale_info$min_wild_date),
        "may go extinct by chance. Use a larger node area."
      ))
    }
  } else if (scale) {
    scale_info$note <- "needs at least 365 days of forcing"
  }

  memory <- list(
    state_bytes_per_replicate = 8 * (f$n_days + 1) * s$n_states,
    totals_bytes_per_replicate = 8 * (f$n_days + 1) * 6 * s$n_g,
    counters_bytes_per_replicate = 8 * f$n_days * length(s$counter_names)
  )
  structure(
    list(
      size = size, options = options, rates = rates, r0 = r0, scale = scale_info,
      memory = memory, warnings = warn
    ),
    class = "sg_check"
  )
}

#' @export
print.sg_check <- function(x, ...) {
  s <- x$size
  cat("<sg_check>\n")
  cat(sprintf(
    "  model: %d genotype(s), %d states, %d days (%s to %s), node %g km^2, %g wounds\n",
    s$genotypes, s$states, s$days, format(s$period[1L]), format(s$period[2L]), s$area, s$wounds
  ))
  cat(sprintf(
    "  releases: %d event(s), %g individuals%s\n", s$release_events, s$released,
    if (length(s$sterile_releases)) {
      paste0("; sterile genotype(s): ", paste(s$sterile_releases, collapse = ", "))
    } else {
      ""
    }
  ))
  o <- x$options
  cat(sprintf(
    paste(
      "  options: mating \"%s\", eggs in %s, search \"%s\", larval crowding \"%s\",",
      "eggs at %s temperature\n"
    ),
    o$mating_mode, o$oviposition_mode, o$wound_competitors,
    o$larval_crowding_mode, o$egg_temp_source
  ))
  r <- x$rates
  cat(sprintf(
    "  fastest sub-stage: %.3g per day (%s); mating rate %.3g per day (%.3g per step)\n",
    r$max_rate, format(r$max_rate_date), r$mating_rate, r$mating_probability_per_step
  ))
  if (isTRUE(x$r0$available)) {
    cat(sprintf(
      paste(
        "  R0 (eggs per egg at low density): %.3g to %.3g; above 1 on %d days",
        "(longest run %d); %.3g at the mean temperature\n"
      ),
      x$r0$range[1L], x$r0$range[2L], x$r0$days_above_1, x$r0$longest_run, x$r0$at_mean_temperature
    ))
  } else {
    cat("  R0: not available (", x$r0$note, ")\n", sep = "")
  }
  if (isTRUE(x$scale$available)) {
    cat(sprintf(
      "  wild population (ODE after a one-year burn-in): %.3g to %.3g individuals per node\n",
      x$scale$min_wild, x$scale$max_wild
    ))
  }
  cat(sprintf(
    "  memory per stochastic replicate: %.3g MB (states) or %.3g MB (totals)\n",
    (x$memory$state_bytes_per_replicate + x$memory$counters_bytes_per_replicate) / 1e6,
    (x$memory$totals_bytes_per_replicate + x$memory$counters_bytes_per_replicate) / 1e6
  ))
  if (length(x$warnings)) {
    cat("  warnings:\n")
    for (w in x$warnings) cat("   - ", w, "\n", sep = "")
  } else {
    cat("  no problems found\n")
  }
  invisible(x)
}
