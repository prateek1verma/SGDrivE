# Reference rate function (pure R) of the ScrewDrivE model.
#
# A direct, readable translation of dev/MODEL_SPEC.md S7-S14 (one node) and
# section 8a N1-N3 (several nodes and adult movement, D-032), used as the
# "truth" against which the C++ kernel is tested. Speed is not a goal here;
# each block cites the specification equation it implements.
#
# With several nodes the state vector holds one per-node vector after another
# (node slowest); each node's rates are those of the one-node model with the
# node's forcing, wounds or hosts and mate-finding half-saturation, and the
# movement terms are added afterwards.
#
# State layout (sg_model(), CONVENTIONS.md): blocks E, L, P, M, U, F; within a
# block the sub-stage varies fastest, then mate genotype (F only), then
# genotype. Unpacked as matrices X[substage, genotype] and, for F, an array
# F[substage, mate, genotype]. With the host-wound model (D-030; design note
# dev/HOST_WOUND_DESIGN.md section 3, spec tags H1-H2) host blocks HS, HW and
# HI (infestation sub-stages) follow.

# Static, model-dependent quantities, computed once per simulation.
.sg_reference_setup <- function(model) {
  d <- model$dims
  p <- model$params
  cb <- model$cube
  geno <- cb$genotypes
  n_g <- d$n_genotypes
  lay <- model$layout
  idx <- lapply(stats::setNames(lay$block, lay$block), function(b) {
    i <- lay$block == b
    seq.int(lay$start[i], length.out = lay$length[i])
  })
  # adult sub-stage widths as [substage, genotype]; classes never reached get
  # width Inf (rate 0)
  widths <- function(ch) {
    w <- t(ch$width)
    w[is.na(w)] <- Inf
    w
  }
  # advance[s, g]: outflow of class s moves to s + 1 (TRUE) or is death (FALSE)
  advance <- function(ch) {
    outer(seq_len(d$n_A), ch$last_class, `<`)
  }
  fch <- model$adult$female
  mch <- model$adult$male
  flag_or_all <- function(ch, name) {
    if (name %in% names(ch$flags)) ch$flags[[name]] else rep(TRUE, d$n_A)
  }
  s <- list(
    geno = geno, n_g = n_g, n_E = d$n_E, n_L = d$n_L, n_P = d$n_P, n_A = d$n_A,
    idx = idx, n_states = model$n_states,
    p = p, counts = model$counts,
    egg_tensor = cb$egg_tensor, eta = cb$eta, phi = cb$phi, xiF = cb$xiF, xiM = cb$xiM,
    # D-031: females search and lay batches in proportion to s_g; viable eggs
    # per egg of a batch, per cross (cubes made before D-031 lack it)
    search_weight = as.numeric(cb$s),
    batch_viability = as.numeric(if (is.null(cb$batch_viability)) {
      .sg_batch_viability(cb$ih, cb$tau, geno)
    } else {
      cb$batch_viability
    }),
    w_f = widths(fch), w_m = widths(mch),
    adv_f = advance(fch), adv_m = advance(mch),
    last_f = fch$last_class, last_m = mch$last_class,
    laying = flag_or_all(fch, "X_ovi"),
    receptive = flag_or_all(fch, "X_rec"),
    mature = flag_or_all(mch, "X_mat"),
    forcing = model$forcing,
    hosts = model$hosts,
    counter_names = .sg_counter_names(geno, hosts = !is.null(model$hosts))
  )
  if (d$n_node == 1L && is.null(model$movement)) {
    return(s) # one node: the setup of the one-node model, unchanged
  }
  # several nodes or movement: one one-node setup per node, with the node's
  # forcing series, wounds, mate finding and hosts (MODEL_SPEC 8a)
  n_node_states <- .sg_n_states_node(model)
  series <- model$landscape$forcing
  series_of <- model$landscape$nodes$forcing_series
  s$node_setups <- lapply(seq_len(d$n_node), function(i) {
    si <- s
    si$forcing <- series[[series_of[i]]]
    si$counts$H <- model$counts$H[i]
    si$counts$h_M <- model$counts$h_M[i]
    if (!is.null(si$hosts)) si$hosts$n <- model$hosts$n[i]
    si$n_states <- n_node_states
    si
  })
  s$n_node <- d$n_node
  s$n_states_node <- n_node_states
  s$movement <- model$movement
  s$counter_names <- .sg_model_counter_names(model)
  s
}

# Counter names of a model: the one-node counters of every node (with an
# "@id" suffix when there are several nodes), then, with movement, the
# emigrants, immigrants and lost flies of each node by class (N2).
.sg_model_counter_names <- function(model) {
  base <- .sg_counter_names(model$cube$genotypes, hosts = !is.null(model$hosts))
  ids <- .sg_node_ids(model)
  per_node <- function(x) {
    if (length(ids) == 1L) x else paste0(rep(x, length(ids)), "@", rep(ids, each = length(x)))
  }
  out <- per_node(base)
  if (!is.null(model$movement)) out <- c(out, per_node(.sg_movement_counter_names))
  out
}

.sg_movement_counter_names <- paste0(
  rep(c("emigrants_", "immigrants_", "lost_"), each = 3L), c("M", "U", "F")
)

# Crowding floor of the per-infested-host larval hazard (D-030, design note
# 3.5): epsilon * N_H keeps the term finite when larvae outlive the
# infestations in the mean field and keeps the model homogeneous in the area.
.sg_crowd_floor <- 1e-6

# Oviposition events per wound found (D-024): 1 when a batch is the day's
# eggs (clutch_size = NA), otherwise R * phi_ov / clutch_size.
.sg_batch_factor <- function(p, phi_ov) {
  if (is.na(p$clutch_size)) 1 else p$eggs_per_day_max * phi_ov / p$clutch_size
}

# Mean eggs per oviposition event (D-024): R * phi_ov unless clutch_size is set.
.sg_clutch <- function(p, phi_ov) {
  if (is.na(p$clutch_size)) p$eggs_per_day_max * phi_ov else p$clutch_size
}

# Share of oviposition events with at least one viable egg (D-026): a batch
# of `clutch` eggs, each viable with probability `viability` (sum over
# offspring of ih tau for the female x mate pair, D-031).
.sg_fertile_share <- function(viability, clutch) -expm1(-clutch * viability)

# Names of the cumulative counters (rates per day in the rate function); the
# host-wound model adds one counter per host flow (D-030).
.sg_counter_names <- function(geno, hosts = FALSE) {
  pairs <- paste0(rep(geno, each = length(geno)), "|", rep(geno, times = length(geno)))
  c(
    paste0("eggs_laid[", geno, "]"),
    "egg_batches", "fertile_batches", "larvae_in", "pupae_in", "emerged_male", "emerged_female",
    "emergence_deaths",
    paste0("matings[", pairs, "]"),
    "deaths_E", "deaths_L", "deaths_P", "deaths_M", "deaths_U", "deaths_F",
    "senescence_M", "senescence_U", "senescence_F",
    if (hosts) .sg_host_counter_names
  )
}

.sg_host_counter_names <- c(
  "hosts_wounded", "wounds_healed", "infestations", "infestations_ended",
  "infestations_treated", "superinfestations"
)

# Erlang chain with one rate q and hazard mu for all sub-stages
# (X[substage, genotype]; inflow per genotype into sub-stage 1).
.sg_chain_const <- function(x, inflow, q, mu) {
  n <- nrow(x)
  dx <- -(q + mu) * x
  dx[1L, ] <- dx[1L, ] + inflow
  if (n > 1L) dx[-1L, ] <- dx[-1L, ] + q * x[-n, , drop = FALSE]
  dx
}

#' Reference rates of the model (pure R)
#'
#' Computes the time derivative of every state and the rates of the cumulative
#' counters on simulation day `day` (model specification S7-S14). Internal:
#' used by the reference solver and as the test oracle for the C++ kernel.
#'
#' @param day Simulation day (1 = first forcing day); the day's forcing is used.
#' @param state Numeric state vector (counts per node).
#' @param setup Result of `.sg_reference_setup(model)`.
#' @param aux If `TRUE`, also return intermediate quantities (W, N_s, Psi,
#'   g_find, pi, rates) for testing.
#' @return A list with `deriv` (length `n_states`), `counters` (named rates
#'   per day) and, if requested, `aux`.
#' @noRd
.sg_rates_r <- function(day, state, setup, aux = FALSE) {
  if (is.null(setup$node_setups)) {
    return(.sg_rates_node(day, state, setup, aux))
  }
  # several nodes (MODEL_SPEC 8a): each node's own rates, then movement (N1, N2)
  n_s <- setup$n_states_node
  deriv <- numeric(length(state))
  counters <- vector("list", setup$n_node)
  node_aux <- vector("list", setup$n_node)
  for (i in seq_len(setup$n_node)) {
    at <- (i - 1L) * n_s + seq_len(n_s)
    r <- .sg_rates_node(day, state[at], setup$node_setups[[i]], aux)
    deriv[at] <- r$deriv
    counters[[i]] <- unname(r$counters)
    node_aux[i] <- list(r$aux)
  }
  counters <- unlist(counters)
  mv <- .sg_movement_r(day, state, setup)
  if (!is.null(mv)) {
    deriv <- deriv + mv$deriv
    counters <- c(counters, mv$counters)
  }
  names(counters) <- setup$counter_names
  out <- list(deriv = deriv, counters = counters)
  if (aux) out$aux <- list(nodes = node_aux, lambda = if (is.null(mv)) NULL else mv$lambda)
  out
}

# Movement terms of all nodes (MODEL_SPEC 8a, N1-N3): derivatives (length of
# the state), counters (emigrants, immigrants, lost by class M, U, F for each
# node, node slowest) and the total leave rate lambda[class, node] of the day.
# NULL without movement.
.sg_movement_r <- function(day, state, setup) {
  mv <- setup$movement
  if (is.null(mv)) {
    return(NULL)
  }
  n_node <- setup$n_node
  n_s <- setup$n_states_node
  blocks <- setup$idx[c("M", "U", "F")] # positions within a node
  deriv <- numeric(length(state))
  emig <- immig <- lost <- lambda <- matrix(0, 3L, n_node)
  dx <- vapply(setup$node_setups, function(ns) ns$forcing$dx[day], 0)
  for (g in mv$groups) {
    # (N3) movement clock of each node on this day
    kappa <- if (g$clock == "calendar") rep(1, n_node) else pmin(mv$kappa_max, dx / mv$dx_ref)
    for (s in 1:3) {
      for (i in seq_len(n_node)) {
        lam <- kappa[i] * g$rate[s, i] # lambda^G_{s,i}
        if (lam == 0) next
        lambda[s, i] <- lambda[s, i] + lam
        from <- (i - 1L) * n_s + blocks[[s]]
        x <- state[from]
        # (N1) leaving node i ...
        deriv[from] <- deriv[from] - lam * x
        emig[s, i] <- emig[s, i] + lam * sum(x)
        # ... arriving in its destinations j
        edges <- seq.int(g$row_ptr[i], length.out = g$row_ptr[i + 1L] - g$row_ptr[i])
        for (e in edges) {
          j <- g$col[e]
          w <- lam * g$prob[s, e]
          to <- (j - 1L) * n_s + blocks[[s]]
          deriv[to] <- deriv[to] + w * x
          immig[s, j] <- immig[s, j] + w * sum(x)
        }
        # (N2) leaving the landscape
        lost[s, i] <- lost[s, i] + lam * g$loss[s, i] * sum(x)
      }
    }
  }
  list(deriv = deriv, counters = as.vector(rbind(emig, immig, lost)), lambda = lambda)
}

# Rates of one node (the one-node model, S7-S14).
.sg_rates_node <- function(day, state, setup, aux = FALSE) {
  s <- setup
  p <- s$p
  f <- s$forcing
  n_g <- s$n_g
  n_a <- s$n_A
  if (day < 1L || day > f$n_days) {
    .sg_stop("Day ", day, " is outside the forcing (1 to ", f$n_days, ").")
  }

  # ---- unpack the state -------------------------------------------------------
  egg <- matrix(state[s$idx$E], s$n_E, n_g)
  lar <- matrix(state[s$idx$L], s$n_L, n_g)
  pup <- matrix(state[s$idx$P], s$n_P, n_g)
  mal <- matrix(state[s$idx$M], n_a, n_g)
  unm <- matrix(state[s$idx$U], n_a, n_g)
  fem <- array(state[s$idx$F], c(n_a, n_g, n_g)) # [substage, mate, genotype]
  h <- s$hosts
  if (!is.null(h)) { # host-wound model (D-030)
    h_s <- state[s$idx$HS]
    h_w <- state[s$idx$HW]
    h_i <- state[s$idx$HI] # [infestation sub-stage]
    i_tot <- sum(h_i)
    sites <- h_w + h$attractiveness * i_tot # attractive sites A_s (H2)
  }

  # ---- daily forcing (S1-S4, S7) ---------------------------------------------
  dx <- f$dx[day]
  dx_egg <- f$dx_egg[day]
  mu_t <- f$mu_T[day]
  mu_e <- f$mu_E[day]
  phi_ov <- f$phi_ov[day]
  moisture <- f$moisture[day]

  # ---- transfer rates (S5, S6) -------------------------------------------------
  q_e <- s$n_E * dx_egg / p$egg_duration_dd
  q_l <- s$n_L / p$larva_duration_days
  q_p <- s$n_P * dx / p$pupa_duration_dd
  q_f <- dx / s$w_f # [substage, genotype], 0 for classes never reached
  q_m <- dx / s$w_m

  # ---- hazards (S7) ------------------------------------------------------------
  l_tot <- sum(lar)
  mu_l <- -log(p$larva_survival_daily)
  if (p$larval_crowding_mode == "per_wound" && s$counts$H > 0) {
    mu_l <- mu_l + p$larval_crowding_coef * l_tot / s$counts$H
  }
  if (!is.null(h)) {
    # treated hosts lose their larvae (H2); crowding per infested host (H2)
    mu_l <- mu_l + h$treatment
    if (p$larval_crowding_mode == "per_infested_host") {
      mu_l <- mu_l + p$larval_crowding_coef * l_tot / (i_tot + .sg_crowd_floor * h$n)
    }
  }
  mu_p <- mu_t + p$pupa_drought_mortality * (1 - moisture)
  mu_f <- mu_t + p$female_background_mortality
  mu_m <- mu_t + p$male_background_mortality

  # ---- wound search and egg production (S8, S9) -------------------------------
  w_lay <- apply(fem[s$laying, , , drop = FALSE], c(2L, 3L), sum) # [mate, genotype]
  # searching females: egg-laying females weighted by their fecundity s_g, so
  # females that produce no eggs do not compete for wounds (D-031)
  w_search <- sweep(w_lay, 2L, s$search_weight, `*`)
  n_s <- if (p$wound_competitors == "searching_females") {
    sum(w_search)
  } else {
    sum(mal) + sum(unm) + sum(fem)
  }
  wounds <- if (is.null(h)) s$counts$H else sites # fixed H, or A_s (H2)
  psi <- if (n_s > 0) {
    p$wound_finding_max * (1 - exp(-p$wound_search_rate * dx * wounds / n_s))
  } else {
    p$wound_finding_max
  }
  # Rows of the egg tensor run over crosses with the mate genotype varying
  # fastest, which is also the order of the flattened w_lay matrix.
  births <- p$eggs_per_day_max * phi_ov * psi * as.vector(crossprod(s$egg_tensor, as.vector(w_lay)))
  names(births) <- s$geno

  # ---- male availability and mating (S12) ------------------------------------
  m_bar <- stats::setNames(colSums(mal[s$mature, , drop = FALSE]), s$geno)
  m_tot <- sum(m_bar)
  g_find <- if (p$male_density_half_mating > 0) {
    m_tot / (s$counts$h_M + m_tot)
  } else {
    as.numeric(m_tot > 0)
  }
  # preference of females of genotype g for the mature males of genotype m
  pref <- sweep(s$eta, 2L, m_bar, `*`)
  pref_sum <- rowSums(pref)
  mate_ok <- pref_sum > 0
  pi_gm <- pref
  pi_gm[mate_ok, ] <- pref[mate_ok, , drop = FALSE] / pref_sum[mate_ok]
  pi_gm[!mate_ok, ] <- 0

  # ---- emergence (S10, S11) ------------------------------------------------------
  emerge <- q_p * pup[s$n_P, ]
  in_m <- (1 - s$phi) * s$xiM * emerge
  in_f <- s$phi * s$xiF * emerge
  emergence_deaths <- sum(emerge) - sum(in_m) - sum(in_f)

  # Mating flux per virgin [substage, genotype] (S12, S14).
  if (p$mating_mode == "instant") {
    # emerging females mate immediately when a compatible male is present;
    # virgins (from releases or times without males) mate at rate nu
    in_u <- ifelse(mate_ok, 0, in_f)
    in_f_mated <- ifelse(mate_ok, in_f, 0)
    nu_rate <- matrix(p$mating_rate * g_find, n_a, n_g)
  } else {
    in_u <- in_f
    in_f_mated <- rep(0, n_g)
    nu_rate <- outer(p$mating_rate * s$receptive * g_find, rep(1, n_g))
  }
  nu_rate[, !mate_ok] <- 0
  mating_flux <- nu_rate * unm # [substage, genotype]

  # ---- derivatives (S10, S13) ----------------------------------------------------
  d_e <- .sg_chain_const(egg, births, q_e, mu_e)
  d_l <- .sg_chain_const(lar, q_e * egg[s$n_E, ], q_l, mu_l)
  d_p <- .sg_chain_const(pup, q_l * lar[s$n_L, ], q_p, mu_p)

  adult_chain <- function(x, q, adv, mu) {
    out <- q * x
    dxa <- -out - mu * x
    if (n_a > 1L) {
      dxa[-1L, ] <- dxa[-1L, ] + (out * adv)[-n_a, , drop = FALSE]
    }
    list(d = dxa, senescence = sum(out * !adv))
  }
  am <- adult_chain(mal, q_m, s$adv_m, mu_m)
  d_m <- am$d
  d_m[1L, ] <- d_m[1L, ] + in_m

  au <- adult_chain(unm, q_f, s$adv_f, mu_f)
  d_u <- au$d - mating_flux
  d_u[1L, ] <- d_u[1L, ] + in_u

  d_f <- array(0, c(n_a, n_g, n_g))
  sen_f <- 0
  for (g in seq_len(n_g)) {
    x <- fem[, , g, drop = FALSE]
    dim(x) <- c(n_a, n_g) # [substage, mate]
    out <- q_f[, g] * x
    adv <- s$adv_f[, g]
    dfg <- -out - mu_f * x
    if (n_a > 1L) dfg[-1L, ] <- dfg[-1L, ] + (out * adv)[-n_a, , drop = FALSE]
    sen_f <- sen_f + sum(out[!adv, , drop = FALSE])
    # matings of virgins of genotype g with males of genotype m keep their age
    dfg <- dfg + outer(mating_flux[, g], pi_gm[g, ])
    dfg[1L, ] <- dfg[1L, ] + in_f_mated[g] * pi_gm[g, ]
    d_f[, , g] <- dfg
  }

  deriv <- numeric(s$n_states)
  deriv[s$idx$E] <- as.vector(d_e)
  deriv[s$idx$L] <- as.vector(d_l)
  deriv[s$idx$P] <- as.vector(d_p)
  deriv[s$idx$M] <- as.vector(d_m)
  deriv[s$idx$U] <- as.vector(d_u)
  deriv[s$idx$F] <- as.vector(d_f)

  # searching females weighted by the share of their batches with a viable
  # egg (D-026, D-031)
  share <- .sg_fertile_share(s$batch_viability, .sg_clutch(p, phi_ov))
  fertile_w <- sum(as.vector(w_search) * share)
  fertile_rate <- psi * .sg_batch_factor(p, phi_ov) * fertile_w # fertile batches per day

  # ---- hosts (H1, H2) ----------------------------------------------------------
  if (!is.null(h)) {
    # fertile batches are shared among wounds in proportion to attractiveness:
    # on an uninfested wound a batch starts an infestation (a case), on an
    # infested one it only adds larvae (superinfestation)
    new_inf <- if (sites > 0) fertile_rate * h_w / sites else 0
    superinf <- if (sites > 0) fertile_rate * h$attractiveness * i_tot / sites else 0
    n_i <- length(h_i)
    wounded <- h$wounding * h_s
    healed <- h$healing * h_w
    ended <- h$advance * h_i[n_i]
    treated <- h$treatment * i_tot
    d_i <- -(h$advance + h$treatment) * h_i
    d_i[1L] <- d_i[1L] + new_inf
    if (n_i > 1L) d_i[-1L] <- d_i[-1L] + h$advance * h_i[-n_i]
    deriv[s$idx$HS] <- -wounded + healed + ended + treated
    deriv[s$idx$HW] <- wounded - healed - new_inf
    deriv[s$idx$HI] <- d_i
  }

  # ---- counters (rates per day) -------------------------------------------------
  # matings[g, m]: virgins of genotype g mating males of genotype m
  matings <- (colSums(mating_flux) + in_f_mated) * pi_gm
  counters <- c(
    unname(births),
    psi * sum(w_search) * .sg_batch_factor(p, phi_ov), # egg_batches
    fertile_rate, # fertile_batches
    q_e * sum(egg[s$n_E, ]),
    q_l * sum(lar[s$n_L, ]),
    sum(in_m), sum(in_f), emergence_deaths,
    as.vector(t(matings)), # [g, m] with m fastest
    mu_e * sum(egg), mu_l * l_tot, mu_p * sum(pup),
    mu_m * sum(mal), mu_f * sum(unm), mu_f * sum(fem),
    am$senescence, au$senescence, sen_f,
    if (!is.null(h)) c(wounded, healed, new_inf, ended, treated, superinf)
  )
  names(counters) <- s$counter_names

  out <- list(deriv = deriv, counters = counters)
  if (aux) {
    out$aux <- list(
      w_lay = w_lay, n_s = n_s, psi = psi, births = births, m_bar = m_bar,
      g_find = g_find, pi = pi_gm, mating_flux = mating_flux,
      q = c(E = q_e, L = q_l, P = q_p), q_f = q_f, q_m = q_m,
      mu = c(E = mu_e, L = mu_l, P = mu_p, F = mu_f, M = mu_m),
      sites = if (is.null(h)) NA_real_ else sites,
      infestation_hazard = if (is.null(h) || h_w <= 0) 0 else new_inf / h_w
    )
  }
  out
}
