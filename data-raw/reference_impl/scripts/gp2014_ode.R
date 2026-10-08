# =============================================================================
# MGDrivE2-style ODE version of the Gutierrez & Ponti (2014) screwworm PBDM
#
# Continuous-time Erlang-chain (linear chain trick) counterpart of the discrete
# daily-step distributed-delay model in `gutierrez_ponti_2014_pbdm.R`. Same
# functions and constants; only the time formulation changes:
#   * distributed delay (A = k*dx/X per day)  -> Erlang rate q = k*dx(T)/X
#   * daily survival s applied after ageing   -> hazard mu = -ln(s)
#   * pooled adults x sex ratio               -> female/male chains, phi = 0.5
#   * egg laying from adult classes 8..k      -> sum over F_i, i >= 8
#   * host search denominator = all adults (current, not previous day)
#   * optional floors on stage inputs mimic the Pascal max(1e-5, .) clamps
#   * optional larval-mortality correction reproduces the discrete scheme's
#     extra day of larval mortality (see match_discrete_larvae)
#
# States (relative density): E_1..k, L_1..k, P_1..k, F_1..k, M_1..k, plus
# cumulative counters for eggs laid, egg batches, pupal and adult inputs.
# Time is in calendar days; temperature is held constant within each day.
# =============================================================================

library(deSolve)
source("scripts/gutierrez_ponti_2014_pbdm.R")   # gp_params() and submodels

gp_ode_simulate <- function(weather, p = gp_params(), floors = TRUE,
                            match_discrete_larvae = FALSE,
                            rtol = 1e-6, atol = 1e-9) {
  k <- p$k; n <- nrow(weather)
  tm <- weather$tmean
  aver5 <- vapply(seq_len(n), function(t) mean(tm[max(1, t - 4):t]), numeric(1))

  # Daily forcing (constant within day d, i.e. on [d-1, d))
  dx   <- gp_degree_days(tm, p)
  muE  <- -log(pmax(0.00001, gp_egg_viability(aver5, p)))
  muA  <- -log(gp_temp_survival(tm, p))           # pupae and adults
  phiO <- gp_ff_temperature(tm, p$TL, p$TU_ovip)
  muL  <- -log(p$larval_surv)
  # The discrete scheme applies the daily larval survival on ~7 days for a
  # 6-day stage (stage survival 0.95^7 = 0.698, vs 0.95^6 here). Option to match it.
  if (match_discrete_larvae) muL <- muL * (p$larva_del + 1) / p$larva_del
  qL   <- k / p$larva_del
  i8   <- p$first_ovip_class:k

  iE <- 1:k; iL <- iE + k; iP <- iL + k; iF <- iP + k; iM <- iF + k
  iC <- 5 * k + 1:4      # cumulative: eggs laid, egg batches, pupae in, adults in

  chain <- function(X, inflow, q, mu) c(inflow, q * X[-k]) - (q + mu) * X

  derivs <- function(t, y, parms) {
    d <- min(n, floor(t) + 1)
    qE <- k * dx[d] / p$egg_del
    qP <- k * dx[d] / p$pupa_del
    qA <- k * dx[d] / p$adult_del
    E <- y[iE]; L <- y[iL]; P <- y[iP]; F <- y[iF]; M <- y[iM]

    Ntot <- sum(F) + sum(M)
    succ <- if (Ntot > 0) 1 - exp(-p$search_const * dx[d] * p$hosts / Ntot) else 1
    ovipF <- sum(F[i8])
    B <- p$fecundity * phiO[d] * ovipF * succ
    if (floors) B <- max(0.00001, B)

    pupIn <- qL * L[k]; adIn <- qP * P[k]
    if (floors) { pupIn <- max(0.0000001, pupIn); adIn <- max(0.000001, adIn) }

    list(c(chain(E, B, qE, muE[d]),
           chain(L, qE * E[k], qL, muL),
           chain(P, pupIn, qP, muA[d]),
           chain(F, p$sex_ratio * adIn, qA, muA[d]),
           chain(M, (1 - p$sex_ratio) * adIn, qA, muA[d]),
           B, ovipF * succ, pupIn, adIn))
  }

  y0 <- c(rep(p$init[["E"]], k), rep(p$init[["L"]], k), rep(p$init[["P"]], k),
          rep(p$init[["A"]] * p$sex_ratio, k), rep(p$init[["A"]] * (1 - p$sex_ratio), k),
          rep(0, 4))
  sol <- ode(y0, times = 0:n, func = derivs, parms = NULL, method = "lsoda",
             rtol = rtol, atol = atol, hmax = 0.25)

  y <- sol[-1, -1]                                 # states at end of each day
  cumd <- function(j) diff(c(0, y[, iC[j]]))       # daily increments
  data.frame(date = weather$date, tmean = tm, dd = dx,
             eggs = rowSums(y[, iE]), larvae = rowSums(y[, iL]),
             pupae = rowSums(y[, iP]),
             females = rowSums(y[, iF]), adults = rowSums(y[, iF]) + rowSums(y[, iM]),
             ovip_females = rowSums(y[, iF[i8]]),
             eggs_laid = cumd(1), egg_batches = cumd(2),
             pupae_in = cumd(3), adults_in = cumd(4))
}
