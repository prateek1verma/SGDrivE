# Initial states from the model's own dynamics (D-025): a burn-in run of the
# ODE, and the closed-form equilibrium at a constant temperature.

# ---- burn-in -----------------------------------------------------------------

# Run the ODE (no releases) for `years` before the model's start and return
# the state at the start, with a record of how it was obtained and how much
# it still depends on the starting state. The burn-in
# weather is either taken from `weather` (the calendar years before the start
# date) or, without it, the model's first 365 days of forcing repeated (the
# same season, one year later).
.sg_init_burn_in <- function(model, x_start, years, weather, check_start = TRUE) {
  ok <- is.numeric(years) && length(years) == 1L && is.finite(years) && years >= 1 &&
    years == round(years)
  if (!ok) .sg_stop("`years` must be a whole number >= 1.")
  start <- model$forcing$dates[1L]
  if (is.null(weather)) {
    if (model$forcing$n_days < 365L) {
      .sg_stop(
        "A burn-in without `weather` repeats the model's first year of forcing, but the ",
        "model has only ", model$forcing$n_days, " days. Give the weather of the years ",
        "before the start in `weather`."
      )
    }
    idx <- rep_len(seq_len(365L), 365L * years)
    forcing <- .sg_forcing_subset(model$forcing, idx, start - rev(seq_along(idx)))
    source <- "first year of the model's forcing, repeated"
  } else {
    burn_start <- seq(start, by = "-1 year", length.out = years + 1L)[years + 1L]
    forcing <- sg_forcing(weather, model$params, start = burn_start, end = start - 1)
    source <- "weather before the start date"
  }
  burn <- model
  burn$forcing <- forcing
  burn$releases <- model$releases[0L, , drop = FALSE]
  run <- function(x0) {
    x <- unname(simulate_ode(burn, init = as.numeric(x0), output_every = forcing$n_days)$state)
    x <- x[nrow(x), ]
    x[x < 0] <- 0 # round-off near extinction
    x
  }
  x <- run(x_start)
  # how much the result still depends on the starting state: the same
  # burn-in from a start ten times larger (skipped by internal callers)
  effect <- NA_real_
  if (check_start) {
    fly <- setdiff(seq_along(x_start), .sg_host_index(model)) # hosts keep their number
    x10 <- x_start
    x10[fly] <- 10 * x_start[fly]
    x10 <- run(x10)
    effect <- max(abs(x10 - x)) / max(max(abs(x)), .Machine$double.xmin)
  }
  if (sum(x) == 0) {
    .sg_warn("The population died out during the burn-in; the initial state is empty.")
  } else if (!is.na(effect) && effect > 1e-3) {
    .sg_warn(sprintf(
      "After %d year(s) of burn-in the state still depends on the starting state (%.2g); %s",
      as.integer(years), effect, "use more years."
    ))
  }
  list(x = x, info = list(
    years = as.integer(years), source = source, days = forcing$n_days,
    start_effect = effect
  ))
}

# A forcing object with the days `idx` of `forcing`, on new dates.
.sg_forcing_subset <- function(forcing, idx, dates) {
  out <- forcing
  daily <- c(
    "tmean", "tmean5", "dx", "dx_egg", "mu_T", "p_A", "egg_viability", "mu_E",
    "phi_ov", "rain", "moisture"
  )
  for (nm in intersect(daily, names(forcing))) out[[nm]] <- forcing[[nm]][idx]
  out$dates <- dates
  out$day <- seq_along(idx)
  out$n_days <- length(idx)
  out
}

# ---- equilibria at a constant temperature -------------------------------------

# All equilibria of a wild-type population at constant temperature (model
# specification section 11, D-025, D-029, D-030). All chains are linear, so
# for an egg input B (eggs per day) every compartment is proportional to B;
# only wound search Psi (and, if enabled, larval crowding, the Allee effect
# and the host-wound model) is nonlinear. B is an equilibrium where the eggs
# laid by the population it sustains, R phi_ov T_ww Psi W(B), equal B.
#
# * Static wounds without crowding or Allee effect: one positive equilibrium
#   in closed form (R phi_ov T_ww Psi* W(1) = 1 gives Psi*, S8 solved for N_s
#   gives N_s*, B* = N_s* / N_s(1)).
# * Otherwise (decision 10 of D-030): the balance gap(log B) = log(eggs laid)
#   - log B is evaluated on a log grid of B (20 points per decade over 30
#   decades below an upper bound where egg laying falls short), every sign
#   change is refined with uniroot(), and each root is classified: a crossing
#   from above to below zero is stable, from below to above unstable (a
#   threshold). Roots closer together than the grid spacing can be missed.
# * With the host-wound model, egg balance fixes the fertile-batch rate for
#   each B, F = B * batch factor * fertile share / (R phi_ov T_ww), so the host
#   steady state (and crowding per infested host) is solved before the fly
#   chains (.sg_host_steady()).
# The extinction equilibrium (B = 0) is listed too: stable when egg laying
# per egg at low density is below one (R0 <= 1, or an Allee effect).
.sg_equilibria_core <- function(model, temperature, what) {
  ok <- is.numeric(temperature) && length(temperature) == 1L && is.finite(temperature)
  if (!ok) .sg_stop(what, " needs a single `temperature` (degC).")
  w <- .sg_wild_setup(model, what)
  s <- w$s
  p <- s$p
  gi <- w$gi
  h <- s$hosts
  rates <- .sg_constant_rates(s, .sg_constant_forcing(model$params, temperature))
  per_b <- .sg_eq_given_b(1, s, gi, rates, crowd = FALSE)
  egg_balance <- p$eggs_per_day_max * rates$phi_ov * w$t_ww # eggs per (Psi * W)
  r0 <- egg_balance * per_b$W * p$wound_finding_max # eggs per egg at low density
  wounds_max <- if (is.null(h)) s$counts$H else max(1, h$attractiveness) * h$n
  pair <- (gi - 1L) * s$n_g + gi
  # fertile batches per egg laid: s_w batches per searching female, each
  # fertile with the share of its batch viability (D-026, D-031)
  fertile_per_egg <- .sg_batch_factor(p, rates$phi_ov) * s$search_weight[gi] *
    .sg_fertile_share(s$batch_viability[pair], .sg_clutch(p, rates$phi_ov)) / egg_balance

  psi_of <- function(n_s, wounds) {
    if (n_s > 0) {
      p$wound_finding_max * (1 - exp(-p$wound_search_rate * rates$dx * wounds / n_s))
    } else {
      p$wound_finding_max
    }
  }
  state_at <- function(b, crowd) {
    hs <- if (!is.null(h)) .sg_host_steady(h, b * fertile_per_egg)
    e <- .sg_eq_given_b(b, s, gi, rates, crowd, infested = hs$infested)
    if (!is.null(h)) {
      e$x[s$idx$HS] <- hs$hs
      e$x[s$idx$HW] <- hs$hw
      e$x[s$idx$HI] <- hs$hi
    }
    wounds <- if (is.null(h)) s$counts$H else hs$sites
    c(e, list(psi = psi_of(e$N_s, wounds), hosts = hs))
  }
  describe <- function(b, st, stable) {
    females <- sum(st$x[c(s$idx$U, s$idx$F)])
    data.frame(
      temperature = temperature, eggs_per_day = b, adult_females = females,
      adult_males = sum(st$x[s$idx$M]),
      infested_hosts = if (is.null(h)) NA_real_ else st$hosts$infested,
      prevalence = if (is.null(h)) NA_real_ else st$hosts$infested / h$n,
      wound_search = st$psi, stable = stable
    )
  }
  extinction_state <- function() {
    st <- state_at(0, crowd = FALSE)
    st$psi <- p$wound_finding_max
    st
  }

  possible <- is.finite(r0) && r0 > 1 && rates$dx > 0 && wounds_max > 0
  allee <- p$male_density_half_mating > 0 && p$mating_mode == "rate"
  crowd_wound <- p$larval_crowding_mode == "per_wound" && isTRUE(s$counts$H > 0)
  crowd_host <- p$larval_crowding_mode == "per_infested_host"
  crowded <- (crowd_wound || crowd_host) && p$larval_crowding_coef > 0
  nonlinear <- crowded || allee || !is.null(h)
  roots <- numeric(0)
  stable <- logical(0)
  extinction_stable <- !possible || allee
  if (possible && !nonlinear) {
    psi <- 1 / (egg_balance * per_b$W)
    n_s <- p$wound_search_rate * rates$dx * s$counts$H / (-log(1 - psi / p$wound_finding_max))
    roots <- n_s / per_b$N_s
    stable <- TRUE
  } else if (possible) {
    gap <- function(log_b) {
      st <- state_at(exp(log_b), crowd = TRUE)
      v <- log(egg_balance * st$psi * st$W) - log_b
      if (is.finite(v)) v else -1e6 # no egg laying at all
    }
    # upper bound: above it egg laying falls short of the input (gap < 0)
    psi0 <- 1 / (egg_balance * per_b$W)
    n_s0 <- p$wound_search_rate * rates$dx * wounds_max / (-log(1 - psi0 / p$wound_finding_max))
    hi <- log(n_s0 / per_b$N_s)
    tries <- 0L
    while (gap(hi) > 0 && tries < 60L) {
      hi <- hi + log(10)
      tries <- tries + 1L
    }
    grid <- seq(hi - 30 * log(10), hi, length.out = 601L)
    g <- vapply(grid, gap, 0)
    extinction_stable <- g[1L] < 0
    up <- g > 0
    cross <- which(up[-1L] != up[-length(up)])
    for (i in cross) {
      r <- stats::uniroot(gap, c(grid[i], grid[i + 1L]), tol = 1e-12)$root
      roots <- c(roots, exp(r))
      stable <- c(stable, up[i]) # from above to below zero: stable
    }
  }
  states <- c(list(extinction_state()), lapply(roots, state_at, crowd = TRUE))
  tab <- do.call(rbind, c(
    list(describe(0, states[[1L]], extinction_stable)),
    Map(describe, roots, states[-1L], stable)
  ))
  rownames(tab) <- NULL
  list(table = tab, states = states, r0 = r0, temperature = temperature, nonlinear = nonlinear)
}

# Host steady state for a constant fertile-batch rate F (design note 3.7,
# D-030): with b the infestation hazard of a wounded host, I = b W D,
# S = (gamma_W + b) W / omega, W = N_H / (1 + (gamma_W + b) / omega + b D),
# and b W(b) (1 + alpha b D) = F, where D is the mean time a host stays
# infested. Without attraction to infested wounds (alpha = 0) F can exceed
# what the wounded hosts can absorb; then no wounded host stays uninfested.
.sg_host_steady <- function(h, fertile) {
  q <- h$advance
  k <- h$treatment
  n <- h$n_I
  d_inf <- if (k > 0) -expm1(n * log(q / (q + k))) / k else n / q
  w_of <- function(z) h$n / (1 + (h$healing + z) / h$wounding + z * d_inf)
  lhs <- function(z) z * w_of(z) * (1 + h$attractiveness * z * d_inf)
  cap <- if (h$attractiveness > 0) Inf else h$n / (1 / h$wounding + d_inf)
  if (fertile <= 0) {
    hw <- w_of(0)
    flux <- 0
    hs <- h$n - hw
  } else if (fertile >= cap) {
    hw <- 0
    flux <- cap
    hs <- h$n - flux * d_inf
  } else {
    hi <- 1
    while (lhs(hi) < fertile) hi <- 2 * hi
    z <- stats::uniroot(function(z) lhs(z) - fertile, c(0, hi), tol = 1e-15 * hi)$root
    hw <- w_of(z)
    flux <- z * hw
    hs <- (h$healing + z) * hw / h$wounding
  }
  hi_vec <- flux / (q + k) * (q / (q + k))^(0:(n - 1L))
  list(
    hs = hs, hw = hw, hi = hi_vec, infested = sum(hi_vec),
    sites = hw + h$attractiveness * sum(hi_vec)
  )
}

# Equilibrium state for sg_init() and adult_females (D-025, D-029, D-030):
# root = "upper" the highest stable equilibrium, root = "lower" the unstable
# equilibrium just below it (the edge of its basin, e.g. the Allee
# threshold). Returns the state and the information stored with it.
.sg_init_equilibrium <- function(model, temperature, what = "method = \"equilibrium\"",
                                 root = c("upper", "lower")) {
  root <- match.arg(root)
  eq <- .sg_equilibria_core(model, temperature, what)
  tab <- eq$table
  pos <- which(tab$eggs_per_day > 0)
  stable_pos <- pos[tab$stable[pos]]
  info <- function(i) {
    list(
      temperature = temperature, R0 = eq$r0,
      psi = if (i == 1L) NA_real_ else tab$wound_search[i],
      eggs_per_day = tab$eggs_per_day[i], adult_females = tab$adult_females[i],
      equilibria = tab
    )
  }
  zero <- list(x = eq$states[[1L]]$x, info = info(1L))
  if (length(stable_pos) == 0L) {
    if (!(eq$r0 > 1)) {
      .sg_message(sprintf(
        "No positive equilibrium at %g degC: eggs per egg at low density R0 = %.3g.",
        temperature, eq$r0
      ))
    } else {
      .sg_message(sprintf(
        "No positive equilibrium at %g degC with larval crowding, the Allee effect or the hosts.",
        temperature
      ))
    }
    return(zero)
  }
  top <- max(stable_pos)
  if (length(stable_pos) > 1L) {
    .sg_message(sprintf(
      paste(
        "%d stable equilibria at %g degC (%s adult females); using the upper one;",
        "see sg_equilibria()."
      ),
      length(stable_pos), temperature,
      paste(format(tab$adult_females[rev(stable_pos)], digits = 4L), collapse = " and ")
    ))
  }
  if (root == "upper") {
    return(list(x = eq$states[[top]]$x, info = info(top)))
  }
  below <- pos[!tab$stable[pos] & pos < top]
  if (length(below) == 0L) {
    .sg_message(
      "No Allee threshold: without an Allee effect (male_density_half_mating > 0 with ",
      "mating_mode = \"rate\") or another unstable equilibrium a population of any size ",
      "grows when R0 > 1."
    )
    return(zero)
  }
  i <- max(below)
  list(x = eq$states[[i]]$x, info = info(i))
}

#' Equilibria at a constant temperature
#'
#' Finds every equilibrium of a wild-type population at a constant
#' temperature and whether it is stable.
#'
#' @param model An `sg_model` object ([sg_model()]) with a single wild-type
#'   genotype that breeds true.
#' @param temperature Constant temperature (degC). Default: the temperature
#'   of the model's forcing, when it is constant (e.g. from
#'   [sg_constant_weather()]).
#' @details An equilibrium is an egg input B (eggs per day) at which the eggs
#'   laid by the population it sustains equal B. Without larval crowding, an
#'   Allee effect or the host-wound model there is one positive equilibrium,
#'   found in closed form. Otherwise the balance is scanned on a fine
#'   logarithmic grid of B (20 points per decade) and every crossing is
#'   refined: a crossing where egg laying goes from above to below the input
#'   is stable, the opposite unstable (a threshold, D-029). With the
#'   host-wound model (`wound_model = "hosts"`) and infested wounds more
#'   attractive than clean ones there can be several stable equilibria
#'   (D-030). Roots closer together than the grid spacing can be missed. The
#'   extinction equilibrium (B = 0) is always listed; it is stable when the
#'   population cannot grow from low numbers (R0 <= 1, or an Allee effect).
#' @return A data frame with one row per equilibrium, ordered by egg input:
#'   `temperature`, `eggs_per_day`, `adult_females` and `adult_males` (in the
#'   node), `infested_hosts` and `prevalence` (infested / hosts; `NA` without
#'   the host-wound model), `wound_search` (Psi) and `stable`. Attribute `R0`:
#'   eggs per egg at low density.
#' @seealso [sg_init()] with `method = "equilibrium"` or `"threshold"`
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params(male_density_half_mating = 5)
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()),
#'     sg_forcing(sg_constant_weather(25, days = 30), p),
#'     area = 100
#'   )
#'   sg_equilibria(mod) # extinction (stable), Allee threshold, upper equilibrium
#' }
sg_equilibria <- function(model, temperature = NULL) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  if (is.null(temperature)) {
    temperature <- .sg_forcing_constant_temp(model$forcing, model$params)
    if (is.null(temperature)) {
      .sg_stop(
        "sg_equilibria() needs a `temperature` (degC): the model's forcing does not ",
        "have one constant temperature."
      )
    }
  }
  eq <- .sg_equilibria_core(model, temperature, "sg_equilibria()")
  structure(eq$table, R0 = eq$r0)
}

# The wild-type genotype of a model for equilibrium and R0 calculations: one
# wild type that produces both sexes and breeds true.
.sg_wild_setup <- function(model, what) {
  cb <- model$cube
  w <- cb$wild_type
  if (length(w) != 1L || cb$ih[w, w, w] < 1 - 1e-12 || cb$phi[[w]] <= 0 || cb$phi[[w]] >= 1) {
    .sg_stop(
      what, " needs a single wild-type genotype that produces both sexes ",
      "and breeds true (wild x wild gives only wild type); this cube has wild type(s) ",
      .sg_wrap_names(w), "."
    )
  }
  s <- .sg_reference_setup(model)
  gi <- match(w, s$geno)
  list(s = s, gi = gi, t_ww = s$egg_tensor[(gi - 1L) * s$n_g + gi, gi])
}

# The temperature of a forcing whose every day equals one constant-temperature
# day (tmean and its 5-day mean constant, and no moisture stress that the
# equilibrium would miss), or NULL.
.sg_forcing_constant_temp <- function(forcing, params) {
  t0 <- forcing$tmean[1L]
  flat <- function(v) all(abs(v - t0) <= 1e-9 * max(1, abs(t0)))
  dry <- params$pupa_drought_mortality > 0 && any(forcing$moisture < 1)
  if (!flat(forcing$tmean) || !flat(forcing$tmean5) || dry) {
    return(NULL)
  }
  t0
}

# Node area (km^2) at which the wild-type equilibrium at the forcing's constant
# temperature holds `n` adult females (D-028). The model is homogeneous of
# degree one in the state and the area, so the area is n divided by the
# equilibrium density.
.sg_area_for_females <- function(model, n) {
  what <- "`adult_females`"
  temperature <- .sg_forcing_constant_temp(model$forcing, model$params)
  if (is.null(temperature)) {
    .sg_stop(
      what, " needs a forcing with one constant temperature, e.g. ",
      "sg_forcing(sg_constant_weather(20), params); for varying weather set `area`, ",
      "or calibrate the density with calibrate_population()."
    )
  }
  eq <- suppressMessages(.sg_init_equilibrium(model, temperature, what))
  females <- eq$info$adult_females
  if (!(females > 0)) {
    .sg_stop(sprintf(
      "%s cannot be reached: no wild-type population persists at %g degC (R0 = %.3g).",
      what, temperature, eq$info$R0
    ))
  }
  density <- females / model$area
  list(adult_females = n, temperature = temperature, density = density, area = n / density)
}

# Rates of the equilibrium calculation from one day's forcing values.
.sg_constant_rates <- function(s, f) {
  p <- s$p
  list(
    dx = f$dx, q_e = s$n_E * f$dx_egg / p$egg_duration_dd, mu_e = f$mu_E,
    q_l = s$n_L / p$larva_duration_days, q_p = s$n_P * f$dx / p$pupa_duration_dd,
    mu_p = f$mu_T + p$pupa_drought_mortality * (1 - f$moisture),
    mu_f = f$mu_T + p$female_background_mortality,
    mu_m = f$mu_T + p$male_background_mortality,
    phi_ov = f$phi_ov
  )
}

# Eggs per egg at low density (R0) of the wild type for given rates: the
# lifetime egg production of an egg's female descendants with wound search at
# its low-density value sigma0 (no crowding; mating as soon as males exist).
.sg_r0_of <- function(w, rates) {
  p <- w$s$p
  per_b <- .sg_eq_given_b(1, w$s, w$gi, rates, crowd = FALSE)
  p$eggs_per_day_max * rates$phi_ov * w$t_ww * per_b$W * p$wound_finding_max
}

# Forcing values at a constant temperature (after five days, so that the
# five-day mean used for egg viability equals the temperature), computed by
# sg_forcing() from constant weather without rainfall.
.sg_constant_forcing <- function(params, temperature) {
  dates <- as.Date("2001-01-01") + 0:9
  wx <- data.frame(date = dates, tmax = temperature, tmin = temperature, tmean = temperature)
  f <- sg_forcing(wx, params)
  list(
    dx = f$dx[10L], dx_egg = f$dx_egg[10L], mu_T = f$mu_T[10L], mu_E = f$mu_E[10L],
    phi_ov = f$phi_ov[10L], moisture = f$moisture[10L]
  )
}

# Steady state of all chains for an egg input b of the wild genotype gi
# (eggs per day), with the rates of a constant temperature. With crowd =
# TRUE the larval hazard includes crowding (solved for) and mating includes
# the Allee effect.
# With the host-wound model the larval hazard includes the treatment rate,
# and with crowding per infested host `infested` (the infested hosts of the
# host steady state) is the crowding denominator (D-030).
.sg_eq_given_b <- function(b, s, gi, rates, crowd, infested = NULL) {
  p <- s$p
  n_a <- s$n_A
  # Erlang chain with rate q and hazard mu, inflow `inflow`: steady state
  chain <- function(inflow, n, q, mu) inflow / (q + mu) * (q / (q + mu))^(0:(n - 1L))
  eggs <- chain(b, s$n_E, rates$q_e, rates$mu_e)
  into_l <- rates$q_e * eggs[s$n_E]
  # larval mortality without crowding
  mu_l0 <- -log(p$larva_survival_daily) + (if (is.null(s$hosts)) 0 else s$hosts$treatment)
  mu_l <- mu_l0
  crowd_coef <- p$larval_crowding_coef
  denom <- if (p$larval_crowding_mode == "per_wound") {
    s$counts$H
  } else if (p$larval_crowding_mode == "per_infested_host" && !is.null(infested)) {
    infested + .sg_crowd_floor * s$hosts$n
  } else {
    0
  }
  if (crowd && p$larval_crowding_mode != "none" && isTRUE(denom > 0) && crowd_coef > 0) {
    l_tot <- function(mu) sum(chain(into_l, s$n_L, rates$q_l, mu))
    fixed <- function(mu) mu - mu_l0 - crowd_coef * l_tot(mu) / denom
    hi <- mu_l0 + crowd_coef * l_tot(mu_l0) / denom
    if (hi > mu_l0) {
      f_hi <- fixed(hi)
      # at very low density the crowding term is below round-off: hi is the root
      mu_l <- if (f_hi <= 0) {
        hi
      } else {
        stats::uniroot(fixed, c(mu_l0, hi),
          f.lower = fixed(mu_l0), f.upper = f_hi, tol = 1e-14
        )$root
      }
    }
  }
  lar <- chain(into_l, s$n_L, rates$q_l, mu_l)
  pup <- chain(rates$q_l * lar[s$n_L], s$n_P, rates$q_p, rates$mu_p)
  emerge <- rates$q_p * pup[s$n_P]
  # segmented adult chains (S6): rate dx / width per class, classes after the
  # last reached one stay empty
  adult <- function(inflow, w, last, mu, extra = rep(0, n_a)) {
    q <- rates$dx / w
    x <- numeric(n_a)
    for (i in seq_len(last)) {
      x[i] <- (if (i == 1L) inflow else x[i - 1L] * q[i - 1L]) / (q[i] + mu + extra[i])
    }
    x
  }
  mal <- adult((1 - s$phi[gi]) * s$xiM[gi] * emerge, s$w_m[, gi], s$last_m[gi], rates$mu_m)
  m_bar <- sum(mal[s$mature])
  mate_ok <- m_bar > 0 && s$eta[gi, gi] > 0
  g_find <- if (crowd && p$male_density_half_mating > 0) {
    m_bar / (s$counts$h_M + m_bar)
  } else {
    as.numeric(m_bar > 0)
  }
  into_f <- s$phi[gi] * s$xiF[gi] * emerge
  all_f <- adult(into_f, s$w_f[, gi], s$last_f[gi], rates$mu_f)
  unm <- if (p$mating_mode == "instant") {
    if (mate_ok) numeric(n_a) else all_f
  } else {
    nu <- if (mate_ok) p$mating_rate * s$receptive * g_find else rep(0, n_a)
    adult(into_f, s$w_f[, gi], s$last_f[gi], rates$mu_f, extra = nu)
  }
  fem <- all_f - unm
  w_lay <- sum(fem[s$laying])
  n_s <- if (p$wound_competitors == "searching_females") {
    s$search_weight[gi] * w_lay # weighted by fecundity (D-031)
  } else {
    sum(mal, unm, fem)
  }

  x <- numeric(s$n_states)
  idx <- function(block, offset_len) s$idx[[block]][(gi - 1L) * offset_len + seq_len(offset_len)]
  x[idx("E", s$n_E)] <- eggs
  x[idx("L", s$n_L)] <- lar
  x[idx("P", s$n_P)] <- pup
  x[idx("M", n_a)] <- mal
  x[idx("U", n_a)] <- unm
  # F block: [genotype, mate, substage]; the mate is the wild type itself
  x[s$idx$F[((gi - 1L) * s$n_g + (gi - 1L)) * n_a + seq_len(n_a)]] <- fem
  list(x = x, W = w_lay, N_s = n_s)
}
