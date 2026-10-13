# Initial states from the model's own dynamics (D-025): a burn-in run of the
# ODE, and the closed-form equilibrium at a constant temperature.

# ---- burn-in -----------------------------------------------------------------

# Run the ODE (no releases) for `years` before the model's start and return
# the state at the start, with a record of how it was obtained and how much
# it still depends on the starting state. The burn-in
# weather is either taken from `weather` (the calendar years before the start
# date) or, without it, the model's first 365 days of forcing repeated (the
# same season, one year later).
.sg_init_burn_in <- function(model, x_start, years, weather) {
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
  # burn-in from a start ten times larger
  x10 <- run(10 * x_start)
  effect <- max(abs(x10 - x)) / max(max(abs(x)), .Machine$double.xmin)
  if (sum(x) == 0) {
    .sg_warn("The population died out during the burn-in; the initial state is empty.")
  } else if (effect > 1e-3) {
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

# ---- equilibrium at a constant temperature -----------------------------------

# Equilibrium of a wild-type population at constant temperature (model
# specification section 11). All chains are linear, so for an egg input B
# (eggs per day) every compartment is proportional to B; only wound search
# Psi (and, if enabled, larval crowding and the Allee effect) is nonlinear.
# Without crowding and Allee effect the equilibrium is in closed form:
# R phi_ov T_ww Psi* W(1) = 1 gives Psi*, Psi* = sigma0 [1 - exp(-a_s dx H /
# N_s*)] gives N_s*, and B* = N_s* / N_s(1). Otherwise B* is found by a
# one-dimensional root search (the upper, stable root).
.sg_init_equilibrium <- function(model, temperature) {
  ok <- is.numeric(temperature) && length(temperature) == 1L && is.finite(temperature)
  if (!ok) .sg_stop("method = \"equilibrium\" needs a single `temperature` (degC).")
  w <- .sg_wild_setup(model, "method = \"equilibrium\"")
  s <- w$s
  p <- s$p
  gi <- w$gi
  t_ww <- w$t_ww
  rates <- .sg_constant_rates(s, .sg_constant_forcing(model$params, temperature))
  per_b <- .sg_eq_given_b(1, s, gi, rates, crowd = FALSE)
  r0 <- p$R * rates$phi_ov * t_ww * per_b$W * p$sigma0 # eggs per egg at low density
  nonlinear <- (p$larval_dd == "per_wound" && s$counts$H > 0 && p$gamma_L > 0) || p$h_M > 0
  zero <- list(
    x = numeric(s$n_states),
    info = list(temperature = temperature, R0 = r0, psi = NA_real_, eggs_per_day = 0)
  )
  if (!is.finite(r0) || r0 <= 1 || rates$dx <= 0 || s$counts$H <= 0) {
    .sg_message(sprintf(
      "No positive equilibrium at %g degC: eggs per egg at low density R0 = %.3g.",
      temperature, r0
    ))
    return(zero)
  }
  psi_of <- function(n_s) {
    if (n_s > 0) p$sigma0 * (1 - exp(-p$a_s * rates$dx * s$counts$H / n_s)) else p$sigma0
  }
  if (!nonlinear) {
    psi <- 1 / (p$R * rates$phi_ov * t_ww * per_b$W)
    n_s <- p$a_s * rates$dx * s$counts$H / (-log(1 - psi / p$sigma0))
    b <- n_s / per_b$N_s
    eq <- .sg_eq_given_b(b, s, gi, rates, crowd = FALSE)
  } else {
    gap <- function(log_b) {
      e <- .sg_eq_given_b(exp(log_b), s, gi, rates, crowd = TRUE)
      log(p$R * rates$phi_ov * t_ww * psi_of(e$N_s) * e$W) - log_b
    }
    # bracket the upper root: above it egg production falls short (gap < 0)
    psi0 <- 1 / (p$R * rates$phi_ov * t_ww * per_b$W)
    b_hi <- log(p$a_s * rates$dx * s$counts$H / (-log(1 - psi0 / p$sigma0)) / per_b$N_s)
    while (gap(b_hi) > 0) b_hi <- b_hi + log(10)
    b_lo <- b_hi
    tries <- 0L
    while (gap(b_lo) <= 0 && tries < 60L) {
      b_lo <- b_lo - log(10)
      tries <- tries + 1L
    }
    if (gap(b_lo) <= 0) {
      .sg_message(sprintf(
        "No positive equilibrium at %g degC with larval crowding or the Allee effect.",
        temperature
      ))
      return(zero)
    }
    b <- exp(stats::uniroot(gap, c(b_lo, b_hi), tol = 1e-12)$root)
    eq <- .sg_eq_given_b(b, s, gi, rates, crowd = TRUE)
    psi <- psi_of(eq$N_s)
  }
  list(x = eq$x, info = list(temperature = temperature, R0 = r0, psi = psi, eggs_per_day = b))
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

# Rates of the equilibrium calculation from one day's forcing values.
.sg_constant_rates <- function(s, f) {
  p <- s$p
  list(
    dx = f$dx, q_e = s$n_E * f$dx_egg / p$X_E, mu_e = f$mu_E,
    q_l = s$n_L / p$tau_L, q_p = s$n_P * f$dx / p$X_P,
    mu_p = f$mu_T + p$mu_W * (1 - f$moisture),
    mu_f = f$mu_T + p$mu_b_F, mu_m = f$mu_T + p$mu_b_M, phi_ov = f$phi_ov
  )
}

# Eggs per egg at low density (R0) of the wild type for given rates: the
# lifetime egg production of an egg's female descendants with wound search at
# its low-density value sigma0 (no crowding; mating as soon as males exist).
.sg_r0_of <- function(w, rates) {
  p <- w$s$p
  per_b <- .sg_eq_given_b(1, w$s, w$gi, rates, crowd = FALSE)
  p$R * rates$phi_ov * w$t_ww * per_b$W * p$sigma0
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
.sg_eq_given_b <- function(b, s, gi, rates, crowd) {
  p <- s$p
  n_a <- s$n_A
  # Erlang chain with rate q and hazard mu, inflow `inflow`: steady state
  chain <- function(inflow, n, q, mu) inflow / (q + mu) * (q / (q + mu))^(0:(n - 1L))
  eggs <- chain(b, s$n_E, rates$q_e, rates$mu_e)
  into_l <- rates$q_e * eggs[s$n_E]
  mu_l <- -log(p$s_L)
  if (crowd && p$larval_dd == "per_wound" && s$counts$H > 0 && p$gamma_L > 0) {
    l_tot <- function(mu) sum(chain(into_l, s$n_L, rates$q_l, mu))
    fixed <- function(mu) mu - (-log(p$s_L)) - p$gamma_L * l_tot(mu) / s$counts$H
    hi <- -log(p$s_L) + p$gamma_L * l_tot(-log(p$s_L)) / s$counts$H
    if (hi > -log(p$s_L)) mu_l <- stats::uniroot(fixed, c(-log(p$s_L), hi), tol = 1e-14)$root
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
  g_find <- if (crowd && p$h_M > 0) m_bar / (s$counts$h_M + m_bar) else as.numeric(m_bar > 0)
  into_f <- s$phi[gi] * s$xiF[gi] * emerge
  all_f <- adult(into_f, s$w_f[, gi], s$last_f[gi], rates$mu_f)
  unm <- if (p$mating == "instant") {
    if (mate_ok) numeric(n_a) else all_f
  } else {
    nu <- if (mate_ok) p$nu * s$receptive * g_find else rep(0, n_a)
    adult(into_f, s$w_f[, gi], s$last_f[gi], rates$mu_f, extra = nu)
  }
  fem <- all_f - unm
  w_lay <- sum(fem[s$laying])
  n_s <- if (p$search_denominator == "searching_females") w_lay else sum(mal, unm, fem)

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
