# Reference rate function (pure R) of the SGDrivE v0.1 model.
#
# A direct, readable translation of dev/MODEL_SPEC.md S7-S14, used as the
# "truth" against which the C++ kernel (milestone M3) is tested. Speed is not
# a goal here; each block cites the specification equation it implements.
#
# State layout (sg_model(), CONVENTIONS.md): blocks E, L, P, M, U, F; within a
# block the sub-stage varies fastest, then mate genotype (F only), then
# genotype. Unpacked as matrices X[substage, genotype] and, for F, an array
# F[substage, mate, genotype].

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
  list(
    geno = geno, n_g = n_g, n_E = d$n_E, n_L = d$n_L, n_P = d$n_P, n_A = d$n_A,
    idx = idx, n_states = model$n_states,
    p = p, counts = model$counts,
    egg_tensor = cb$egg_tensor, eta = cb$eta, phi = cb$phi, xiF = cb$xiF, xiM = cb$xiM,
    w_f = widths(fch), w_m = widths(mch),
    adv_f = advance(fch), adv_m = advance(mch),
    last_f = fch$last_class, last_m = mch$last_class,
    laying = flag_or_all(fch, "X_ovi"),
    receptive = flag_or_all(fch, "X_rec"),
    mature = flag_or_all(mch, "X_mat"),
    forcing = model$forcing,
    counter_names = .sg_counter_names(geno)
  )
}

# Oviposition events per wound found (D-024): 1 when a batch is the day's
# eggs (clutch_size = NA), otherwise R * phi_ov / clutch_size.
.sg_batch_factor <- function(p, phi_ov) {
  if (is.na(p$clutch_size)) 1 else p$R * phi_ov / p$clutch_size
}

# Mean eggs per oviposition event (D-024): R * phi_ov unless clutch_size is set.
.sg_clutch <- function(p, phi_ov) {
  if (is.na(p$clutch_size)) p$R * phi_ov else p$clutch_size
}

# Share of oviposition events with at least one viable egg (D-026): a batch
# of `clutch` eggs, each viable with probability `viability` (sum over
# offspring of s ih tau for the female x mate pair).
.sg_fertile_share <- function(viability, clutch) -expm1(-clutch * viability)

# Names of the cumulative counters (rates per day in the rate function).
.sg_counter_names <- function(geno) {
  pairs <- paste0(rep(geno, each = length(geno)), "|", rep(geno, times = length(geno)))
  c(
    paste0("eggs_laid[", geno, "]"),
    "egg_batches", "fertile_batches", "larvae_in", "pupae_in", "emerged_male", "emerged_female",
    "emergence_deaths",
    paste0("matings[", pairs, "]"),
    "deaths_E", "deaths_L", "deaths_P", "deaths_M", "deaths_U", "deaths_F",
    "senescence_M", "senescence_U", "senescence_F"
  )
}

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

  # ---- daily forcing (S1-S4, S7) ---------------------------------------------
  dx <- f$dx[day]
  dx_egg <- f$dx_egg[day]
  mu_t <- f$mu_T[day]
  mu_e <- f$mu_E[day]
  phi_ov <- f$phi_ov[day]
  moisture <- f$moisture[day]

  # ---- transfer rates (S5, S6) -------------------------------------------------
  q_e <- s$n_E * dx_egg / p$X_E
  q_l <- s$n_L / p$tau_L
  q_p <- s$n_P * dx / p$X_P
  q_f <- dx / s$w_f # [substage, genotype], 0 for classes never reached
  q_m <- dx / s$w_m

  # ---- hazards (S7) ------------------------------------------------------------
  l_tot <- sum(lar)
  mu_l <- -log(p$s_L)
  if (p$larval_dd == "per_wound" && s$counts$H > 0) {
    mu_l <- mu_l + p$gamma_L * l_tot / s$counts$H
  }
  mu_p <- mu_t + p$mu_W * (1 - moisture)
  mu_f <- mu_t + p$mu_b_F
  mu_m <- mu_t + p$mu_b_M

  # ---- wound search and egg production (S8, S9) -------------------------------
  w_lay <- apply(fem[s$laying, , , drop = FALSE], c(2L, 3L), sum) # [mate, genotype]
  n_s <- if (p$search_denominator == "searching_females") {
    sum(w_lay)
  } else {
    sum(mal) + sum(unm) + sum(fem)
  }
  psi <- if (n_s > 0) {
    p$sigma0 * (1 - exp(-p$a_s * dx * s$counts$H / n_s))
  } else {
    p$sigma0
  }
  # Rows of the egg tensor run over crosses with the mate genotype varying
  # fastest, which is also the order of the flattened w_lay matrix.
  births <- p$R * phi_ov * psi * as.vector(crossprod(s$egg_tensor, as.vector(w_lay)))
  names(births) <- s$geno

  # ---- male availability and mating (S12) ------------------------------------
  m_bar <- stats::setNames(colSums(mal[s$mature, , drop = FALSE]), s$geno)
  m_tot <- sum(m_bar)
  g_find <- if (p$h_M > 0) m_tot / (s$counts$h_M + m_tot) else as.numeric(m_tot > 0)
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
  if (p$mating == "instant") {
    # emerging females mate immediately when a compatible male is present;
    # virgins (from releases or times without males) mate at rate nu
    in_u <- ifelse(mate_ok, 0, in_f)
    in_f_mated <- ifelse(mate_ok, in_f, 0)
    nu_rate <- matrix(p$nu * g_find, n_a, n_g)
  } else {
    in_u <- in_f
    in_f_mated <- rep(0, n_g)
    nu_rate <- outer(p$nu * s$receptive * g_find, rep(1, n_g))
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

  # ---- counters (rates per day) -------------------------------------------------
  # egg-laying females weighted by the share of their batches with a viable
  # egg (D-026)
  share <- .sg_fertile_share(rowSums(s$egg_tensor), .sg_clutch(p, phi_ov))
  fertile_w <- sum(as.vector(w_lay) * share)
  # matings[g, m]: virgins of genotype g mating males of genotype m
  matings <- (colSums(mating_flux) + in_f_mated) * pi_gm
  counters <- c(
    unname(births),
    psi * sum(w_lay) * .sg_batch_factor(p, phi_ov), # egg_batches
    psi * .sg_batch_factor(p, phi_ov) * fertile_w, # fertile_batches
    q_e * sum(egg[s$n_E, ]),
    q_l * sum(lar[s$n_L, ]),
    sum(in_m), sum(in_f), emergence_deaths,
    as.vector(t(matings)), # [g, m] with m fastest
    mu_e * sum(egg), mu_l * l_tot, mu_p * sum(pup),
    mu_m * sum(mal), mu_f * sum(unm), mu_f * sum(fem),
    am$senescence, au$senescence, sen_f
  )
  names(counters) <- s$counter_names

  out <- list(deriv = deriv, counters = counters)
  if (aux) {
    out$aux <- list(
      w_lay = w_lay, n_s = n_s, psi = psi, births = births, m_bar = m_bar,
      g_find = g_find, pi = pi_gm, mating_flux = mating_flux,
      q = c(E = q_e, L = q_l, P = q_p), q_f = q_f, q_m = q_m,
      mu = c(E = mu_e, L = mu_l, P = mu_p, F = mu_f, M = mu_m)
    )
  }
  out
}
