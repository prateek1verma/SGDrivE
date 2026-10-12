# Reference chain-binomial solver (pure R) of the SGDrivE v0.1 model.
#
# A direct translation of dev/MODEL_SPEC.md S15-S19 (D-022), the "truth"
# against which the C++ stepper (milestone M4) is tested. In stochastic mode
# it makes its random draws (stats::rbinom, stats::rpois) in the canonical
# order of the specification, so the C++ stepper can be compared draw for
# draw; in deterministic mode every draw is replaced by its expectation,
# which gives the deterministic map (the exact conditional mean of a step).
# Speed is not a goal. Internal; not a user-facing backend.

# Unbiased rounding to whole numbers: floor(x) + Bernoulli(x - floor(x)).
# Random numbers are drawn only for entries with a fractional part.
.sg_cb_round <- function(x) {
  fl <- floor(x)
  frac <- x - fl
  up <- frac > 0
  fl[up] <- fl[up] + (stats::runif(sum(up)) < frac[up])
  fl
}

# Draw functions: random (stochastic = TRUE) or expectations.
.sg_cb_draws <- function(stochastic) {
  if (!stochastic) {
    return(list(
      binom = function(n, p) n * p,
      pois = function(lambda) lambda,
      multinom = function(n, prob) n * prob
    ))
  }
  binom <- function(n, p) as.numeric(stats::rbinom(1L, n, p))
  list(
    binom = binom,
    pois = function(lambda) as.numeric(stats::rpois(1L, lambda)),
    # sequential conditional binomials in category order (D-022); the sums
    # of the remaining probabilities are accumulated in double precision from
    # the last category, as in C++
    multinom = function(n, prob) {
      k <- length(prob)
      out <- numeric(k)
      rest <- prob
      if (k > 1L) for (j in (k - 1L):1L) rest[j] <- prob[j] + rest[j + 1L]
      left <- n
      for (j in seq_len(k - 1L)) {
        if (left == 0) break
        pj <- if (rest[j] > 0) min(1, prob[j] / rest[j]) else 0
        out[j] <- binom(left, pj)
        left <- left - out[j]
      }
      out[k] <- left
      out
    }
  )
}

# A chain prepared for steps of length dt (S16), from sub-stage rates q,
# last reached sub-stage and hazard mu: for each starting sub-stage the
# probability of dying in the chain (die) and the distribution of survivors
# over sub-stages and "left the chain" (dest, rows summing to 1); for
# individuals arriving during a step (at a uniform time), the probability of
# dying before its end (arr_die) and the distribution of survivors over
# sub-stages (arr; the small probability of leaving the chain again within
# the step is kept in the last reached sub-stage). Sub-stages a genotype
# never reaches keep anyone placed there by hand, who only die (as in the
# ODE, rate 0).
.sg_cb_chain <- function(q, last, mu, dt) {
  n <- length(q)
  a <- .sg_chain_advance_r(q, last, dt, mu)
  keep <- rowSums(a)
  die <- if (mu > 0) pmin(pmax(1 - keep, 0), 1) else numeric(n)
  dest <- a / ifelse(keep > 0, keep, 1)
  if (last < n) {
    for (i in (last + 1L):n) {
      dest[i, ] <- 0
      dest[i, i] <- 1
      die[i] <- 1 - exp(-mu * dt)
    }
  }
  first <- .sg_chain_average_r(q, last, dt, mu)[1L, ]
  arr <- first[seq_len(n)]
  arr[last] <- arr[last] + first[n + 1L]
  list(
    dest = dest, die = die, arr = arr / sum(arr),
    arr_die = if (mu > 0) min(max(1 - sum(first), 0), 1) else 0
  )
}

# Forcing-dependent quantities of one day for steps of length dt: hazards and
# the prepared chains (S16) of eggs, pupae and adults; larvae are prepared in
# each step because their hazard depends on crowding.
.sg_cb_day_r <- function(day, s, dt) {
  f <- s$forcing
  p <- s$p
  dx <- f$dx[day]
  mu <- c(
    E = f$mu_E[day], P = f$mu_T[day] + p$mu_W * (1 - f$moisture[day]),
    F = f$mu_T[day] + p$mu_b_F, M = f$mu_T[day] + p$mu_b_M
  )
  q_e <- s$n_E * f$dx_egg[day] / p$X_E
  q_p <- s$n_P * dx / p$X_P
  list(
    day = day, phi_ov = f$phi_ov[day], mu = mu, q_l = s$n_L / p$tau_L,
    ch_e = .sg_cb_chain(rep(q_e, s$n_E), s$n_E, mu[["E"]], dt),
    ch_p = .sg_cb_chain(rep(q_p, s$n_P), s$n_P, mu[["P"]], dt),
    ch_m = lapply(seq_len(s$n_g), function(g) {
      .sg_cb_chain(dx / s$w_m[, g], s$last_m[g], mu[["M"]], dt)
    }),
    ch_f = lapply(seq_len(s$n_g), function(g) {
      .sg_cb_chain(dx / s$w_f[, g], s$last_f[g], mu[["F"]], dt)
    })
  )
}

# One step of length dt from state x (S15-S19). `info` is .sg_cb_day_r() of
# the step's day, `draw` is .sg_cb_draws(). Returns the new state and the
# flows of the step (counter order of .sg_counter_names()).
.sg_cb_step_r <- function(x, info, s, dt, draw) {
  p <- s$p
  n_g <- s$n_g
  n_a <- s$n_A

  # ---- aggregates (S15): from the state at the start of the step ----------
  a <- .sg_rates_r(info$day, x, s, aux = TRUE)$aux
  mate_ok <- rowSums(a$pi) > 0

  egg <- matrix(x[s$idx$E], s$n_E, n_g)
  lar <- matrix(x[s$idx$L], s$n_L, n_g)
  pup <- matrix(x[s$idx$P], s$n_P, n_g)
  mal <- matrix(x[s$idx$M], n_a, n_g)
  unm <- matrix(x[s$idx$U], n_a, n_g)
  fem <- array(x[s$idx$F], c(n_a, n_g, n_g)) # [substage, mate, genotype]

  # ---- eggs (S17): Poisson per female x mate pair, split over offspring ---
  eggs <- numeric(n_g)
  for (g in seq_len(n_g)) {
    for (m in seq_len(n_g)) {
      tgm <- s$egg_tensor[(g - 1L) * n_g + m, ]
      tot <- sum(tgm)
      if (tot <= 0) next
      n_eggs <- draw$pois(p$R * info$phi_ov * a$psi * a$w_lay[m, g] * tot * dt)
      eggs <- eggs + draw$multinom(n_eggs, tgm / tot)
    }
  }

  # ---- mating (S18), before ageing ------------------------------------------
  nu_i <- if (p$mating == "instant") {
    rep(p$nu * a$g_find, n_a)
  } else {
    p$nu * s$receptive * a$g_find
  }
  matings <- matrix(0, n_g, n_g) # [g, m]
  for (g in seq_len(n_g)) {
    if (!mate_ok[g]) next
    for (i in seq_len(n_a)) {
      y <- draw$binom(unm[i, g], 1 - exp(-nu_i[i] * dt))
      if (y == 0) next
      ym <- draw$multinom(y, a$pi[g, ])
      unm[i, g] <- unm[i, g] - y
      fem[i, , g] <- fem[i, , g] + ym
      matings[g, ] <- matings[g, ] + ym
    }
  }

  # ---- chain moves (S16): death in the chain, then destination --------------
  ch_l <- .sg_cb_chain(rep(info$q_l, s$n_L), s$n_L, a$mu[["L"]], dt)
  move <- function(xv, ch) {
    n <- length(xv)
    stay <- numeric(n)
    out <- 0
    deaths <- 0
    for (i in seq_len(n)) {
      if (xv[i] == 0) next
      d <- draw$binom(xv[i], ch$die[i])
      dest <- draw$multinom(xv[i] - d, ch$dest[i, i:(n + 1L)])
      stay[i:n] <- stay[i:n] + dest[-length(dest)]
      out <- out + dest[length(dest)]
      deaths <- deaths + d
    }
    list(x = stay, out = out, deaths = deaths)
  }
  # each block in state order: genotype (then mate) slowest, sub-stage fastest
  block_moves <- function(x, chain_of) {
    res <- lapply(seq_len(ncol(x)), function(j) move(x[, j], chain_of(j)))
    list(
      x = vapply(res, `[[`, numeric(nrow(x)), "x"),
      out = vapply(res, `[[`, 0, "out"),
      deaths = sum(vapply(res, `[[`, 0, "deaths"))
    )
  }
  as_mat <- function(v, nr) matrix(v, nr) # keep a matrix when n_g = 1
  me <- block_moves(egg, function(j) info$ch_e)
  ml <- block_moves(lar, function(j) ch_l)
  mp <- block_moves(pup, function(j) info$ch_p)
  mm <- block_moves(mal, function(j) info$ch_m[[j]])
  mu_ <- block_moves(unm, function(j) info$ch_f[[j]])
  fem_mat <- matrix(fem, n_a) # columns: (mate, genotype), mate fastest
  mf <- block_moves(fem_mat, function(j) info$ch_f[[(j - 1L) %/% n_g + 1L]])

  new_e <- as_mat(me$x, s$n_E)
  new_l <- as_mat(ml$x, s$n_L)
  new_p <- as_mat(mp$x, s$n_P)
  new_m <- as_mat(mm$x, n_a)
  new_u <- as_mat(mu_$x, n_a)
  new_f <- array(mf$x, c(n_a, n_g, n_g))

  # ---- emergence (S19): sex, then survival ------------------------------------
  em_m <- numeric(n_g)
  em_f <- numeric(n_g)
  for (g in seq_len(n_g)) {
    z <- mp$out[g]
    zf <- draw$binom(z, s$phi[g])
    em_f[g] <- draw$binom(zf, s$xiF[g])
    em_m[g] <- draw$binom(z - zf, s$xiM[g])
  }

  # ---- arrivals (S16): they entered at times spread over the step; death in
  # the new chain for the rest of the step, then placement at the expected
  # position (D-022); in state order
  arr_deaths <- c(E = 0, L = 0, P = 0, M = 0, U = 0, F = 0)
  # arrivals into a block, one column per genotype: deaths, then placement
  place_block <- function(x, arrivals, chain_of) {
    deaths <- 0
    for (g in seq_along(arrivals)) {
      if (arrivals[g] == 0) next
      ch <- chain_of(g)
      d <- draw$binom(arrivals[g], ch$arr_die)
      deaths <- deaths + d
      x[, g] <- x[, g] + draw$multinom(arrivals[g] - d, ch$arr)
    }
    list(x = x, deaths = deaths)
  }
  pe <- place_block(new_e, eggs, function(g) info$ch_e)
  pl <- place_block(new_l, me$out, function(g) ch_l)
  pp <- place_block(new_p, ml$out, function(g) info$ch_p)
  pm <- place_block(new_m, em_m, function(g) info$ch_m[[g]])
  new_e <- pe$x
  new_l <- pl$x
  new_p <- pp$x
  new_m <- pm$x
  arr_deaths[c("E", "L", "P", "M")] <- c(pe$deaths, pl$deaths, pp$deaths, pm$deaths)
  # Emerged females (S18, S19): death for the rest of the step, then mating:
  # at emergence with mating = "instant", otherwise at rate nu_1 for the rest
  # of the step (probability 1 - (1 - exp(-nu_1 dt)) / (nu_1 dt) for a
  # uniform arrival time); unmated survivors go to U, mated ones to F, both
  # placed like other arrivals.
  instant <- p$mating == "instant"
  x_arr <- nu_i[1L] * dt
  p_arr_mate <- if (x_arr > 0) 1 + expm1(-x_arr) / x_arr else 0
  for (g in seq_len(n_g)) {
    ch <- info$ch_f[[g]]
    ym <- numeric(n_g)
    unmated <- 0
    if (em_f[g] > 0) {
      mated_now <- instant && mate_ok[g]
      d <- draw$binom(em_f[g], ch$arr_die)
      arr_deaths[[if (mated_now) "F" else "U"]] <- arr_deaths[[if (mated_now) "F" else "U"]] + d
      surv <- em_f[g] - d
      y <- if (mated_now) {
        surv
      } else if (!instant && mate_ok[g]) {
        draw$binom(surv, p_arr_mate)
      } else {
        0
      }
      if (y > 0) ym <- draw$multinom(y, a$pi[g, ])
      unmated <- surv - y
      matings[g, ] <- matings[g, ] + ym
    }
    if (unmated > 0) new_u[, g] <- new_u[, g] + draw$multinom(unmated, ch$arr)
    for (m in seq_len(n_g)) {
      if (ym[m] > 0) new_f[, m, g] <- new_f[, m, g] + draw$multinom(ym[m], ch$arr)
    }
  }

  out <- numeric(s$n_states)
  out[s$idx$E] <- as.vector(new_e)
  out[s$idx$L] <- as.vector(new_l)
  out[s$idx$P] <- as.vector(new_p)
  out[s$idx$M] <- as.vector(new_m)
  out[s$idx$U] <- as.vector(new_u)
  out[s$idx$F] <- as.vector(new_f)

  flows <- c(
    eggs,
    a$psi * sum(a$w_lay) * dt, # egg_batches: expected value (D-022)
    sum(me$out), sum(ml$out), sum(em_m), sum(em_f),
    sum(mp$out) - sum(em_m) - sum(em_f),
    as.vector(t(matings)),
    c(me$deaths, ml$deaths, mp$deaths, mm$deaths, mu_$deaths, mf$deaths) + arr_deaths,
    sum(mm$out), sum(mu_$out), sum(mf$out)
  )
  names(flows) <- s$counter_names
  list(x = out, flows = flows)
}

# Day loop of the reference chain-binomial solver (tests only): releases at
# the start of each day, then 1/dt steps with the day's forcing. In
# stochastic mode the initial state and the releases are first rounded
# without bias (initial state, then releases in event order).
.sg_cb_run_r <- function(model, init = NULL, dt = 0.1, stochastic = TRUE, n_days = NULL) {
  s <- .sg_reference_setup(model)
  steps <- round(1 / dt)
  if (steps < 1L || abs(steps * dt - 1) > 1e-9) {
    .sg_stop("`dt` must divide one day into a whole number of steps.")
  }
  if (is.null(n_days)) n_days <- model$forcing$n_days
  x <- .sg_check_init(model, init)
  rel <- model$releases
  if (stochastic) {
    x <- .sg_cb_round(x)
    rel$count <- .sg_cb_round(rel$count)
  }
  draw <- .sg_cb_draws(stochastic)
  state <- matrix(NA_real_, n_days + 1L, s$n_states,
    dimnames = list(NULL, sg_state_names(model))
  )
  state[1L, ] <- x
  counters <- matrix(0, n_days, length(s$counter_names),
    dimnames = list(NULL, s$counter_names)
  )
  blocks <- c("E", "P", "M", "U")
  released <- matrix(0, n_days, 4L, dimnames = list(NULL, blocks))
  rel_block <- substr(sg_state_names(model)[rel$index], 1L, 1L)
  for (d in seq_len(n_days)) {
    for (j in which(rel$day == d)) {
      x[rel$index[j]] <- x[rel$index[j]] + rel$count[j]
      released[d, rel_block[j]] <- released[d, rel_block[j]] + rel$count[j]
    }
    info <- .sg_cb_day_r(d, s, dt)
    for (k in seq_len(steps)) {
      r <- .sg_cb_step_r(x, info, s, dt, draw)
      x <- r$x
      counters[d, ] <- counters[d, ] + r$flows
    }
    state[d + 1L, ] <- x
  }
  list(
    day = 0:n_days, state = state, counters = counters, releases = released,
    dt = dt, stochastic = stochastic
  )
}
