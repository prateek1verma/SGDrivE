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
#' @return A numeric vector of R0 values, one per temperature. R0 at low
#'   density does not depend on a node's area, wounds or hosts (wound search
#'   is at its low-density value), so for a model with a landscape it applies
#'   to every node at that temperature (D-032).
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
#' age and at emergence), between consecutive output days. With several
#' nodes (D-032) the balance is node by node and includes movement
#' (immigrants minus emigrants, emigrants including those that leave the
#' landscape); with the node extinction rule ([simulate_ode()]) it includes
#' the removed flies.
#'
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @return A data frame with columns `day` (end of each interval between
#'   output days), `replicate` (stochastic runs only), `change` (change in
#'   the total), `releases`, `eggs_laid`, `deaths` and `residual` (change
#'   minus releases and eggs plus deaths). The residual is about 1e-12 of
#'   the population for the ODE and exactly zero for stochastic runs. Host
#'   states are not counted. With the host-wound model (D-030) a column
#'   `host_residual` gives the largest mismatch between the change of a host
#'   block (HS, HW, HI) and its host flows (round-off for the ODE, zero for
#'   stochastic runs). With several nodes a column `node` after `day`, one row
#'   per node and interval; with movement a column `movement` (net
#'   immigration); with the node extinction rule a column `removed`; the
#'   residual includes them.
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   max(abs(sg_flow_balance(sim)$residual))
#' }
sg_flow_balance <- function(sim) {
  sel <- .sg_sim_nodes(sim, default = "all")
  .sg_need_counters(sim, "sg_flow_balance()")
  if (!sel$multi) {
    return(.sg_flow_balance_unit(sim, .sg_sim_array(sim)$array, sim$counters, sim$releases))
  }
  # several nodes: each node's states, counters and releases (D-032)
  a_all <- .sg_sim_array(sim, default = "all")
  k_node <- .sg_name_node(colnames(sim$counters))
  r_node <- .sg_name_node(colnames(sim$releases))
  rows <- lapply(sel$ids, function(id) {
    a <- a_all$array[, a_all$node == id, , drop = FALSE]
    k <- sim$counters[, k_node == id, drop = FALSE]
    colnames(k) <- sub("@.*$", "", colnames(k))
    rel <- sim$releases[, r_node == id, drop = FALSE]
    d <- .sg_flow_balance_unit(sim, a, k, rel)
    cbind(d[1L], node = id, d[-1L])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Flow balance of one node (or of a one-node simulation): states a [time,
# column, replicate], counters k and releases rel without node suffixes.
.sg_flow_balance_unit <- function(sim, a, k, rel) {
  stochastic <- inherits(sim, "sg_sim_cb")
  if (!stochastic) {
    cn0 <- colnames(k)
    dim(k) <- c(dim(k), 1L)
    dimnames(k) <- list(NULL, cn0, NULL)
    dim(rel) <- c(dim(rel), 1L)
  }
  cn <- dimnames(k)[[2L]]
  out_days <- sim$day[-1L]
  interval <- findInterval(seq_len(dim(k)[1L]), c(0, out_days) + 1L)
  host <- .sg_is_host_name(dimnames(a)[[2L]])
  moving <- any(grepl("^emigrants_", cn))
  removing <- any(grepl("^removed_", cn))
  rows <- lapply(seq_len(dim(a)[3L]), function(r) {
    total <- rowSums(matrix(a[, !host, r], dim(a)[1L]))
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
    if (moving) {
      d$movement <- sum_by(per_day(grepl("^immigrants_", cn))) -
        sum_by(per_day(grepl("^emigrants_", cn)))
      d$residual <- d$residual - d$movement
    }
    if (removing) {
      d$removed <- sum_by(per_day(grepl("^removed_", cn)))
      d$residual <- d$residual + d$removed
    }
    if (any(host)) {
      # host conservation (D-030): the change of each host block equals its
      # flows; the largest absolute mismatch of HS, HW and HI
      block <- substr(dimnames(a)[[2L]], 1L, 2L)
      hb <- function(b) diff(rowSums(matrix(a[, block == b, r], dim(a)[1L])))
      flow <- function(name) sum_by(per_day(cn == name))
      to_s <- -flow("hosts_wounded") + flow("wounds_healed") + flow("infestations_ended") +
        flow("infestations_treated")
      to_w <- flow("hosts_wounded") - flow("wounds_healed") - flow("infestations")
      to_i <- flow("infestations") - flow("infestations_ended") - flow("infestations_treated")
      d$host_residual <- pmax(abs(hb("HS") - to_s), abs(hb("HW") - to_w), abs(hb("HI") - to_i))
    }
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
#'
#' **Several nodes** (D-032). Every node is checked as a one-node model (its
#' forcing, area, wounds or hosts), giving a table by node (`nodes`) with the
#' fastest rate, the R0 range and days with R0 > 1 and, with `scale`, the
#' smallest wild population in the node from a burn-in and run of the whole
#' landscape with movement. `movement` summarises the movement groups
#' (largest leave rate by class at the reference temperature, nodes with open
#' edges). Warnings are those of the nodes and of the model; with more than
#' one node a note recalls that ODE runs without `node_extinction` can show
#' recolonisation by fractions of a fly.
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
  if (!is.null(model$landscape)) {
    return(.sg_check_landscape(model, dt, scale))
  }
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
    "larval_crowding_mode", "egg_temp_source", "male_density_half_mating", "wound_model"
  )]
  h <- model$hosts
  if (!is.null(h)) { # host-wound model (D-030)
    # a wound receives at most wound_finding_max * wound_search_rate * dx
    # oviposition events per day, however many flies there are (S8)
    strike <- p$wound_finding_max * p$wound_search_rate * f$dx
    size$hosts <- list(
      n = h$n, wounded_fraction = p$wounded_fraction, wounding = h$wounding,
      healing = h$healing, infestation_days = p$infestation_days, n_I = h$n_I,
      treatment = h$treatment, attractiveness = h$attractiveness,
      max_strike_rate = max(strike), max_strike_date = f$dates[which.max(strike)]
    )
  }

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
  if (!is.null(h)) day_max <- pmax(day_max, h$wounding, h$healing, h$advance + h$treatment)
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
    totals_bytes_per_replicate = 8 * (f$n_days + 1) * (6 * s$n_g + if (is.null(h)) 0 else 3),
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

# check_model() of a model with a landscape (D-032): each node as a one-node
# model (without the population-scale run), the movement groups, and the
# population scale from a burn-in and run of the whole landscape.
.sg_check_landscape <- function(model, dt, scale) {
  ids <- .sg_node_ids(model)
  per_node <- lapply(seq_along(ids), function(i) {
    check_model(.sg_node_model(model, i), dt = dt, scale = FALSE)
  })
  nd <- model$landscape$nodes
  r0_end <- function(ck, k) if (isTRUE(ck$r0$available)) ck$r0$range[k] else NA_real_
  nodes <- data.frame(
    node = ids, role = nd$role, area = nd$area, forcing_series = nd$forcing_series,
    wounds = if (is.null(model$hosts)) model$counts$H else NA_real_,
    hosts = if (is.null(model$hosts)) NA_real_ else model$hosts$n,
    max_rate = vapply(per_node, function(ck) ck$rates$max_rate, 0),
    r0_min = vapply(per_node, function(ck) r0_end(ck, 1L), 0),
    r0_max = vapply(per_node, function(ck) r0_end(ck, 2L), 0),
    days_r0_above_1 = vapply(per_node, function(ck) {
      if (isTRUE(ck$r0$available)) as.numeric(ck$r0$days_above_1) else NA_real_
    }, 0),
    stringsAsFactors = FALSE
  )
  # warnings: those of the model (the same in every node) once, node ones by node
  node_warn <- unlist(lapply(seq_along(ids), function(i) per_node[[i]]$warnings))
  warn <- character(0)
  for (w in unique(node_warn)) {
    hit <- ids[vapply(per_node, function(ck) w %in% ck$warnings, TRUE)]
    warn <- c(warn, if (length(hit) == length(ids)) {
      w
    } else {
      paste0("node(s) ", paste(hit, collapse = ", "), ": ", w)
    })
  }
  if (length(ids) > 1L) {
    warn <- c(warn, paste(
      "With several nodes, ODE runs without node_extinction = TRUE can show",
      "recolonisation by fractions of a fly (MODEL_SPEC N5)."
    ))
  }
  mv <- model$movement
  movement <- if (is.null(mv)) {
    NULL
  } else {
    list(
      t_ref = mv$t_ref, kappa_max = mv$kappa_max,
      groups = lapply(mv$groups, function(g) {
        list(
          clock = g$clock, edges = length(g$col),
          max_rate = apply(g$rate, 1L, max),
          open_nodes = ids[apply(g$loss, 2L, function(v) any(v > 0))]
        )
      })
    )
  }
  scale_info <- list(available = FALSE)
  f <- model$forcing
  if (scale && f$n_days >= 365L) {
    init <- suppressWarnings(sg_init(model, "burn_in"))
    sim <- simulate_ode(model, init = init)
    nm <- colnames(sim$state)
    geno <- sub("^.\\[([^],|]+).*$", "\\1", nm)
    wild <- geno %in% model$cube$wild_type
    node_of <- .sg_name_node(nm)
    if (length(ids) == 1L) node_of <- rep(ids, length(nm))
    min_wild <- vapply(ids, function(id) {
      min(rowSums(sim$state[, wild & node_of == id, drop = FALSE]))
    }, 0)
    nodes$min_wild <- unname(min_wild)
    scale_info <- list(available = TRUE, min_wild = min(min_wild))
    small <- ids[min_wild < 10]
    if (length(small) > 0L) {
      warn <- c(warn, paste0(
        "The wild population falls below 10 individuals in node(s) ", paste(small, collapse = ", "),
        " (ODE with movement after a one-year burn-in); stochastic runs may lose them by chance."
      ))
    }
  } else if (scale) {
    scale_info$note <- "needs at least 365 days of forcing"
  }
  first <- per_node[[1L]]
  rel <- model$releases
  size <- first$size
  size$states <- model$n_states
  size$area <- sum(nd$area)
  size$n_nodes <- length(ids)
  size$release_events <- nrow(rel)
  size$released <- sum(rel$count)
  size$released_genotypes <- unique(rel$genotype)
  wild_f <- intersect(model$cube$wild_type, model$cube$genotypes[model$cube$phi > 0])
  size$sterile_releases <- size$released_genotypes[vapply(size$released_genotypes, function(g) {
    all(model$cube$egg_total[paste0(wild_f, "|", g)] == 0)
  }, TRUE)]
  n_counters <- length(.sg_model_counter_names(model))
  memory <- list(
    state_bytes_per_replicate = 8 * (f$n_days + 1) * model$n_states,
    totals_bytes_per_replicate = 8 * (f$n_days + 1) * length(ids) *
      (6 * model$dims$n_genotypes + if (is.null(model$hosts)) 0 else 3),
    counters_bytes_per_replicate = 8 * f$n_days * n_counters
  )
  structure(
    list(
      size = size, options = first$options, rates = list(
        max_rate = max(nodes$max_rate), mating_rate = first$rates$mating_rate,
        mating_probability_per_step = first$rates$mating_probability_per_step
      ),
      r0 = list(available = any(!is.na(nodes$r0_max))), scale = scale_info,
      memory = memory, warnings = warn, nodes = nodes, movement = movement
    ),
    class = "sg_check"
  )
}

#' @export
print.sg_check <- function(x, ...) {
  if (!is.null(x$nodes)) {
    return(.sg_print_check_landscape(x))
  }
  s <- x$size
  cat("<sg_check>\n")
  wounds <- if (is.null(s$hosts)) sprintf("%g wounds", s$wounds) else sprintf("%g hosts", s$hosts$n)
  cat(sprintf(
    "  model: %d genotype(s), %d states, %d days (%s to %s), node %g km^2, %s\n",
    s$genotypes, s$states, s$days, format(s$period[1L]), format(s$period[2L]), s$area, wounds
  ))
  if (!is.null(s$hosts)) {
    hh <- s$hosts
    cat(sprintf(
      paste0(
        "  hosts: %.3g wounded without flies (wounds heal in %.3g d); infestations %.3g d ",
        "in %d sub-stages; treatment %.3g per day; infested wounds x%.3g attractive; ",
        "at most %.3g oviposition events per wound per day (%s)\n"
      ),
      hh$wounded_fraction, 1 / hh$healing, hh$infestation_days, hh$n_I, hh$treatment,
      hh$attractiveness, hh$max_strike_rate, format(hh$max_strike_date)
    ))
  }
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

.sg_print_check_landscape <- function(x) {
  s <- x$size
  cat("<sg_check>\n")
  cat(sprintf(
    "  model: %d genotype(s), %d states in %d nodes (%s km^2), %d days (%s to %s)\n",
    s$genotypes, s$states, s$n_nodes, format(s$area, big.mark = ","), s$days,
    format(s$period[1L]), format(s$period[2L])
  ))
  cat(sprintf(
    "  releases: %d event(s), %g individuals%s\n", s$release_events, s$released,
    if (length(s$sterile_releases)) {
      paste0("; sterile genotype(s): ", paste(s$sterile_releases, collapse = ", "))
    } else {
      ""
    }
  ))
  cat("  nodes:\n")
  print(x$nodes, row.names = FALSE, digits = 4L)
  if (!is.null(x$movement)) {
    for (g in x$movement$groups) {
      cat(sprintf(
        paste0(
          "  movement (%s clock): %d edges; largest leave rate at %g degC: male %.3g, ",
          "virgin %.3g, mated %.3g per day; open edges in %s\n"
        ),
        g$clock, g$edges, x$movement$t_ref, g$max_rate[1L], g$max_rate[2L], g$max_rate[3L],
        if (length(g$open_nodes)) paste(g$open_nodes, collapse = ", ") else "no node"
      ))
    }
  } else {
    cat("  movement: none\n")
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
