# Reference chain-binomial solver (pure R) of the ScrewDrivE v0.1 model.
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
      multinom = function(n, prob) n * prob,
      ztp_sum = function(n, lambda) n * lambda / -expm1(-lambda),
      occupancy = function(n, w) if (w > 0) w * -expm1(-n / w) else 0
    ))
  }
  binom <- function(n, p) as.numeric(stats::rbinom(1L, n, p))
  pois <- function(lambda) as.numeric(stats::rpois(1L, lambda))
  list(
    binom = binom,
    pois = pois,
    # Total of n zero-truncated Poisson(lambda) numbers: the eggs of n
    # batches known to hold at least one egg (D-026). Above lambda = 30 the
    # truncation changes probabilities by less than 1e-13, so one Poisson
    # draw is used; otherwise each batch by inversion (as in C++).
    ztp_sum = function(n, lambda) {
      if (lambda > 30) {
        return(pois(n * lambda))
      }
      k_max <- ceiling(lambda + 12 * sqrt(lambda) + 30)
      total <- 0
      for (i in seq_len(n)) {
        u <- stats::runif(1L) * -expm1(-lambda)
        p <- exp(-lambda) * lambda
        cum <- p
        k <- 1
        while (u > cum && k < k_max) {
          k <- k + 1
          p <- p * lambda / k
          cum <- cum + p
        }
        total <- total + k
      }
      total
    },
    # Number of distinct hosts hit when n batches fall at random on w wounded
    # hosts (host-wound model, D-030, design note 3.6): batch by batch, a
    # batch hits a host not hit before with probability (w - k) / w; stops
    # once every host is hit. Never more than w.
    occupancy = function(n, w) {
      k <- 0
      for (i in seq_len(n)) {
        if (k >= w) break
        if (stats::runif(1L) < (w - k) / w) k <- k + 1
      }
      k
    },
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
  if (!is.null(s$node_setups)) { # several nodes (D-032): the chains of each node
    return(list(day = day, nodes = lapply(s$node_setups, function(ns) .sg_cb_day_r(day, ns, dt))))
  }
  f <- s$forcing
  p <- s$p
  dx <- f$dx[day]
  mu <- c(
    E = f$mu_E[day], P = f$mu_T[day] + p$pupa_drought_mortality * (1 - f$moisture[day]),
    F = f$mu_T[day] + p$female_background_mortality, M = f$mu_T[day] + p$male_background_mortality
  )
  q_e <- s$n_E * f$dx_egg[day] / p$egg_duration_dd
  q_p <- s$n_P * dx / p$pupa_duration_dd
  list(
    day = day, phi_ov = f$phi_ov[day], mu = mu, q_l = s$n_L / p$larva_duration_days,
    ch_e = .sg_cb_chain(rep(q_e, s$n_E), s$n_E, mu[["E"]], dt),
    ch_p = .sg_cb_chain(rep(q_p, s$n_P), s$n_P, mu[["P"]], dt),
    ch_m = lapply(seq_len(s$n_g), function(g) {
      .sg_cb_chain(dx / s$w_m[, g], s$last_m[g], mu[["M"]], dt)
    }),
    ch_f = lapply(seq_len(s$n_g), function(g) {
      .sg_cb_chain(dx / s$w_f[, g], s$last_f[g], mu[["F"]], dt)
    }),
    # infestation chain of the host-wound model (D-030): treatment in the role
    # of the hazard; both exits return the host to HS
    ch_i = if (!is.null(s$hosts)) {
      h <- s$hosts
      .sg_cb_chain(rep(h$advance, h$n_I), h$n_I, h$treatment, dt)
    }
  )
}

# One step of length dt from state x (S15-S19). `info` is .sg_cb_day_r() of
# the step's day, `draw` is .sg_cb_draws(). Returns the new state and the
# flows of the step (counter order of .sg_model_counter_names()).
#
# With several nodes (D-032, MODEL_SPEC section 10) a step has three parts:
# part 1 (aggregates, eggs, hosts, mating, chain moves) of each node in index
# order, movement of all nodes (N4), then part 3 (emergence, arrivals) of
# each node in index order. With one node: part 1, then part 3.
.sg_cb_step_r <- function(x, info, s, dt, draw) {
  if (is.null(s$node_setups)) {
    return(.sg_cb_part3(.sg_cb_part1(x, info, s, dt, draw)))
  }
  n_s <- s$n_states_node
  at <- function(i) (i - 1L) * n_s + seq_len(n_s)
  work <- lapply(seq_len(s$n_node), function(i) {
    .sg_cb_part1(x[at(i)], info$nodes[[i]], s$node_setups[[i]], dt, draw)
  })
  new_x <- numeric(length(x))
  node_flows <- vector("list", s$n_node)
  for (i in seq_len(s$n_node)) {
    r <- .sg_cb_part3(work[[i]])
    new_x[at(i)] <- r$x
    node_flows[[i]] <- unname(r$flows)
  }
  all_flows <- unlist(node_flows)
  n_move <- length(s$counter_names) - length(all_flows)
  all_flows <- c(all_flows, numeric(n_move)) # movement counters: phase 2 WP3
  names(all_flows) <- s$counter_names
  list(x = new_x, flows = all_flows)
}

# Part 1 of a step of one node (S15, S17, H3, S18, S16): returns everything
# part 3 needs (its local variables).
.sg_cb_part1 <- function(x, info, s, dt, draw) {
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

  # chain moves of one chain (S16): death in the chain, then destination
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

  # ---- eggs (S17), per female x mate pair, split over offspring genotypes:
  # oviposition events (batches, D-024), those with at least one viable egg
  # (fertile, the case proxy; D-026), their viable eggs; or single eggs ----
  batches <- p$oviposition_mode == "batches"
  factor <- .sg_batch_factor(p, info$phi_ov) # events per wound found
  clutch <- .sg_clutch(p, info$phi_ov)
  eggs <- numeric(n_g)
  n_batches <- 0
  n_fertile <- 0
  expected_fertile <- 0
  for (g in seq_len(n_g)) {
    for (m in seq_len(n_g)) {
      pair <- (g - 1L) * n_g + m
      tgm <- s$egg_tensor[pair, ]
      tot <- sum(tgm)
      # females lay batches (and search) in proportion to s_g; viable eggs
      # per egg of a batch without s_g (D-031)
      w_search <- s$search_weight[g] * a$w_lay[m, g]
      viab <- s$batch_viability[pair]
      share <- .sg_fertile_share(viab, clutch)
      if (batches) {
        nb <- draw$pois(a$psi * factor * w_search * dt)
        n_batches <- n_batches + nb
        if (nb == 0 || share <= 0) next
        nf <- draw$binom(nb, share)
        n_fertile <- n_fertile + nf
        if (nf == 0) next
        n_eggs <- draw$ztp_sum(nf, clutch * viab)
      } else {
        expected_fertile <- expected_fertile + w_search * share
        if (tot <= 0) next
        n_eggs <- draw$pois(p$eggs_per_day_max * info$phi_ov * a$psi * a$w_lay[m, g] * tot * dt)
      }
      eggs <- eggs + draw$multinom(n_eggs, tgm / tot)
    }
  }

  # ---- hosts (H3, D-030): fertile batches on uninfested wounds infest the
  # distinct hosts they hit (cases); the rest are superinfestations; then
  # healing, the infestation chain (end or treatment) and wounding -------------
  h <- s$hosts
  if (!is.null(h)) {
    if (!batches) { # single eggs: fertile oviposition events drawn for the hosts
      n_fertile <- draw$pois(a$psi * factor * expected_fertile * dt)
    }
    h_s <- x[s$idx$HS]
    h_w <- x[s$idx$HW]
    h_i <- x[s$idx$HI]
    sites <- h_w + h$attractiveness * sum(h_i)
    n_on_w <- if (sites > 0 && n_fertile > 0) draw$binom(n_fertile, h_w / sites) else 0
    infested <- draw$occupancy(n_on_w, h_w)
    healed <- draw$binom(h_w - infested, 1 - exp(-h$healing * dt))
    mi <- move(h_i, info$ch_i)
    wounded <- draw$binom(h_s, 1 - exp(-h$wounding * dt))
    new_hi <- mi$x
    new_hi[1L] <- new_hi[1L] + infested
    # used in part 3
    host_flows <- c( # nolint: object_usage_linter.
      wounded, healed, infested, mi$out, mi$deaths, n_fertile - infested
    )
  }

  # ---- mating (S18), before ageing ------------------------------------------
  nu_i <- if (p$mating_mode == "instant") {
    rep(p$mating_rate * a$g_find, n_a)
  } else {
    p$mating_rate * s$receptive * a$g_find
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

  new_e <- as_mat(me$x, s$n_E) # nolint: object_usage_linter. (part 3)
  new_l <- as_mat(ml$x, s$n_L) # nolint: object_usage_linter. (part 3)
  new_p <- as_mat(mp$x, s$n_P) # nolint: object_usage_linter. (part 3)
  new_m <- as_mat(mm$x, n_a) # nolint: object_usage_linter. (part 3)
  new_u <- as_mat(mu_$x, n_a) # nolint: object_usage_linter. (part 3)
  new_f <- array(mf$x, c(n_a, n_g, n_g)) # nolint: object_usage_linter. (part 3)
  as.list(environment())
}

# Part 3 of a step of one node (S19, arrivals) from the result of part 1:
# the node's new state and flows. The variables of part 1 are restored into
# this function's environment, which the linter cannot see.
# nolint start: object_usage_linter.
.sg_cb_part3 <- function(work) {
  list2env(work, environment())

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
  instant <- p$mating_mode == "instant"
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
  if (!is.null(h)) {
    out[s$idx$HS] <- h_s - wounded + healed + mi$out + mi$deaths
    out[s$idx$HW] <- h_w - infested - healed + wounded
    out[s$idx$HI] <- new_hi
  }

  flows <- c(
    eggs,
    if (batches) { # egg_batches
      n_batches
    } else {
      a$psi * sum(sweep(a$w_lay, 2L, s$search_weight, `*`)) * factor * dt
    },
    if (batches || !is.null(h)) { # fertile_batches
      n_fertile
    } else {
      a$psi * factor * expected_fertile * dt
    },
    sum(me$out), sum(ml$out), sum(em_m), sum(em_f),
    sum(mp$out) - sum(em_m) - sum(em_f),
    as.vector(t(matings)),
    c(me$deaths, ml$deaths, mp$deaths, mm$deaths, mu_$deaths, mf$deaths) + arr_deaths,
    sum(mm$out), sum(mu_$out), sum(mf$out),
    if (!is.null(h)) host_flows
  )
  names(flows) <- s$counter_names
  list(x = out, flows = flows)
}
# nolint end

# Inputs of a chain-binomial run, shared by the R and C++ drivers so that
# both consume random numbers identically: steps per day, initial state and
# release counts (rounded without bias when stochastic: initial state, then
# releases in event order) and the matrix of individuals released per day
# and block.
.sg_cb_inputs <- function(model, init, dt, stochastic, n_days) {
  dt_ok <- is.numeric(dt) && length(dt) == 1L && dt > 0
  steps <- if (dt_ok) round(1 / dt) else 0
  if (!dt_ok || steps < 1L || abs(steps * dt - 1) > 1e-9) {
    .sg_stop("`dt` must divide one day into a whole number of steps.")
  }
  if (is.null(n_days)) n_days <- model$forcing$n_days
  x <- .sg_check_init(model, init)
  rel <- model$releases
  if (stochastic) {
    x <- .sg_cb_round(x)
    rel$count <- .sg_cb_round(rel$count)
  }
  released <- matrix(0, n_days, 4L, dimnames = list(NULL, c("E", "P", "M", "U")))
  rel_block <- substr(sg_state_names(model)[rel$index], 1L, 1L)
  for (j in which(rel$day <= n_days)) {
    released[rel$day[j], rel_block[j]] <- released[rel$day[j], rel_block[j]] + rel$count[j]
  }
  list(steps = steps, n_days = n_days, x = x, rel = rel, released = released)
}

# Day loop of the reference chain-binomial solver (tests only): releases at
# the start of each day, then 1/dt steps with the day's forcing.
.sg_cb_run_r <- function(model, init = NULL, dt = 0.1, stochastic = TRUE, n_days = NULL) {
  s <- .sg_reference_setup(model)
  inp <- .sg_cb_inputs(model, init, dt, stochastic, n_days)
  n_days <- inp$n_days
  x <- inp$x
  rel <- inp$rel
  draw <- .sg_cb_draws(stochastic)
  state <- matrix(NA_real_, n_days + 1L, s$n_states,
    dimnames = list(NULL, sg_state_names(model))
  )
  state[1L, ] <- x
  counters <- matrix(0, n_days, length(s$counter_names),
    dimnames = list(NULL, s$counter_names)
  )
  for (d in seq_len(n_days)) {
    for (j in which(rel$day == d)) x[rel$index[j]] <- x[rel$index[j]] + rel$count[j]
    info <- .sg_cb_day_r(d, s, dt)
    for (k in seq_len(inp$steps)) {
      r <- .sg_cb_step_r(x, info, s, dt, draw)
      x <- r$x
      counters[d, ] <- counters[d, ] + r$flows
    }
    state[d + 1L, ] <- x
  }
  list(
    day = 0:n_days, state = state, counters = counters, releases = inp$released,
    dt = dt, stochastic = stochastic
  )
}

# The same run with the C++ step and day loop (src/cb_step.cpp), returning
# the same structure (states of every day) plus the extinction day.
.sg_cb_run_cpp <- function(model, init = NULL, dt = 0.1, stochastic = TRUE, n_days = NULL) {
  inp <- .sg_cb_inputs(model, init, dt, stochastic, n_days)
  view <- .sg_model_view(model)
  res <- sg_cb_run_cpp(
    view, inp$x, as.numeric(inp$rel$count), dt, stochastic, inp$n_days,
    rep(TRUE, inp$n_days)
  )
  state <- rbind(inp$x, res$states)
  dimnames(state) <- list(NULL, sg_state_names(model))
  colnames(res$counters) <- attr(view, "counter_names")
  list(
    day = 0:inp$n_days, state = state, counters = res$counters, releases = inp$released,
    dt = dt, stochastic = stochastic, extinct_day = res$extinct_day,
    wild_extinct_day = res$wild_extinct_day, days_run = res$days_run
  )
}
