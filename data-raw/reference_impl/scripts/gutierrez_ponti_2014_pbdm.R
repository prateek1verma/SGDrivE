# =============================================================================
# Gutierrez & Ponti (2014) screwworm PBDM -- R translation
#
# Source: Gutierrez, A.P. & Ponti, L. (2014) The new world screwworm: prospective
#   distribution and role of weather in eradication. Agricultural and Forest
#   Entomology 16:158-173 (Appendix Eq. A1; Eqs. 1-6) and supplementary Doc S9
#   (Delphi Pascal unit `ScrewWormPop`, procedures SWSetup, SWmod, SWdynamics,
#   DelayNoPlr).
#
# This is a line-by-line translation of the published Pascal code. Where the
# code and paper disagree, the CODE is followed (see `gp_params()`); the paper
# values are given in comments. Two pieces are not in the published code and
# were filled in from the paper:
#   * ffTemperature(T, 14.5, 43.5): paper Eq. 6, phi = 1 - ((T - TL - psi)/psi)^2,
#     psi = (TU - TL)/2, truncated to [0, 1].
#   * Aver5 (running 5-day mean temperature): mean of the available days for the
#     first 4 days of the run.
#
# Model summary (all densities relative, k = 25 age classes per stage):
#   eggs (10 dd, ambient dd)  -> larvae (6 calendar days) -> pupae (95 dd)
#   -> adults (202.8 dd; exit from last class = senescence death).
#   Mortality is applied once per day after ageing (survival multipliers).
#   Eggs are laid by adult classes 8..25 (pre-oviposition), 90 eggs/female/day
#   x temperature scalar x sex ratio, reduced by a ratio-dependent host-search
#   term with H = 5 hosts and search constant 0.001 per dd.
# =============================================================================


# ---- Parameters (as in the Pascal code) --------------------------------------
gp_params <- function() {
  list(
    k          = 25,      # age classes per stage (KSWegg..KSWAdlt)
    egg_del    = 10,      # dd   (SwEggDel; paper 8.64 dd)
    larva_del  = 6,       # days (literal in DelayNoPLR call; paper 136.8 dd)
    larva_delk_base = 5.7,# SwLarvaeDel used for rate<->number conversion
    pupa_del   = 95,      # dd   (SwPupaeDel)
    adult_del  = 202.8,   # dd   (SwAdltDel)
    first_ovip_class = 8, # adult classes 8..k lay eggs
    fecundity  = 90,      # eggs per female per day at ~29 C
    sex_ratio  = 0.5,
    larval_surv = 0.95,   # daily larval survival
    TL = 14.5, TU_dev = 40.5, TU_ovip = 43.5,   # thresholds (paper Eq.1 uses TU = 40)
    a_dev = 0.00065, b_dev = 1.75, pre_ovip_dd = 61,
    dx_floor = 0.01,      # max(0.01, 24*R(T)) -> minimum 0.61 dd/day
    T_opt = 27.2,
    cA = 0.00036, pA0 = 0.0035,   # pupal/adult mortality (Eq. 2/3 in code)
    cE = 0.00375, pE0 = 0.3,      # egg hatch failure (Eq. 3/4 in code)
    search_const = 0.001, hosts = 5,
    init = c(E = 0.25, L = 0.25, P = 0.25, A = 0.25),   # SWE, SWL, SWP, SWA
    rain_threshold = 450          # mm/yr, post-hoc I_rain (not in dynamics)
  )
}


# ---- Submodels ---------------------------------------------------------------

# Daily physiological time (dd/day) -- SWmod: DDASw
gp_degree_days <- function(Tmean, p = gp_params()) {
  r_day <- 24 * (p$a_dev * (Tmean - p$TL) / (1 + p$b_dev^(Tmean - p$TU_dev)))
  p$pre_ovip_dd * pmax(p$dx_floor, r_day)
}

# Daily proportion of pupae/adults dying (SWTempMort), clamped to [0, 1]
gp_temp_mortality <- function(Tmean, p = gp_params()) {
  pmin(1, pmax(0, p$cA * (Tmean - p$T_opt)^2 + p$pA0))
}

# Daily pupal/adult survival (Tempsurvival)
gp_temp_survival <- function(Tmean, p = gp_params()) {
  pmax(0.00001, 1 - (p$cA * (Tmean - p$T_opt)^2 + p$pA0))
}

# Egg viability on 5-day mean temperature (Eggviability); applied as max(1e-5, .)
gp_egg_viability <- function(Taver5, p = gp_params()) {
  pmin(1, 1 - (p$cE * (Taver5 - p$T_opt)^2 + p$pE0))
}

# Temperature scalar for oviposition: ffTemperature(T, 14.5, 43.5), paper Eq. 6
gp_ff_temperature <- function(Tmean, TL = 14.5, TU = 43.5) {
  psi <- (TU - TL) / 2
  pmin(1, pmax(0, 1 - ((Tmean - TL - psi) / psi)^2))
}


# ---- Distributed delay without attrition: DelayNoPlr --------------------------
# r   : vector of k flow-rate densities (r[1] = Pascal r[0])
# vin : input increment during the time step (numbers)
# del : mean transit time; dt : time to process (dd or days)
# Returns list(r = updated vector, vout = outflow in rate units; multiply by
# del/k to convert to numbers, as the Pascal code does).
gp_delay_no_plr <- function(vin, r, del, dt, k = length(r)) {
  vout <- 0
  if (dt > 0) {
    idt <- trunc(1 + 2 * dt * (k / del))
    a   <- (dt / (del / k)) / idt
    for (j in seq_len(idt)) {
      vout <- vout + a * r[k]
      if (k > 1) r[k:2] <- r[k:2] + a * (r[(k - 1):1] - r[k:2])
      r[1] <- r[1] + a * ((vin / dt) - r[1])
    }
  }
  list(r = r, vout = vout)
}
# Note: the Pascal loop updates r[i] from i = k down to 1 using the OLD r[i-1]
# (descending order), which the vectorised update above reproduces exactly.


# ---- Main simulation ---------------------------------------------------------
# weather: data.frame with columns `date` (Date) and `tmean` (deg C); optional
#          `rain` (mm/day) for the post-hoc I_rain summary.
gp_simulate <- function(weather, p = gp_params()) {
  stopifnot(all(c("date", "tmean") %in% names(weather)))
  n <- nrow(weather); k <- p$k

  # SWSetup: widths of age classes
  egg_delk   <- p$egg_del / k
  larv_delk  <- p$larva_delk_base / k
  pup_delk   <- p$pupa_del / k
  adlt_delk  <- p$adult_del / k

  # Initial values apportioned to every age class (as rates)
  r_egg   <- rep(p$init[["E"]] / egg_delk,  k)
  r_larva <- rep(p$init[["L"]] / larv_delk, k)
  r_pupa  <- rep(p$init[["P"]] / pup_delk,  k)
  r_adult <- rep(p$init[["A"]] / adlt_delk, k)

  # Running 5-day mean (Aver5); first days use available values
  tm <- weather$tmean
  aver5 <- vapply(seq_len(n), function(t) mean(tm[max(1, t - 4):t]), numeric(1))

  tot_adult_prev <- 0   # TotSwAdlt is a global, 0 before the first SWmod call
  out <- matrix(NA_real_, n, 13,
                dimnames = list(NULL, c("tmean", "aver5", "dd", "temp_mort",
                                        "temp_surv", "egg_viab", "phi_temp",
                                        "eggs_laid", "larvae_in", "pupae_in",
                                        "adults_in", "adults_out", "dummy")))
  # Derived (not in the Pascal output): females, egg-laying females, and the
  # number of females that find a host each day (= egg batches laid = proxy
  # for new myiasis infestations). The 2014 model is unsexed and has no
  # mating, so every female is implicitly mated.
  derived <- matrix(NA_real_, n, 4, dimnames = list(NULL, c(
    "females", "ovip_females", "search_success", "egg_batches")))
  tot <- matrix(NA_real_, n, 4, dimnames = list(NULL, c("eggs", "larvae", "pupae", "adults")))

  for (t in seq_len(n)) {
    Tm <- tm[t]

    # --- SWmod: temperature indices and degree-days
    temp_mort <- gp_temp_mortality(Tm, p)
    dd <- gp_degree_days(Tm, p)

    # --- SWdynamics
    egg_viab  <- gp_egg_viability(aver5[t], p)
    temp_surv <- gp_temp_survival(Tm, p)
    phi_temp  <- gp_ff_temperature(Tm, p$TL, p$TU_ovip)

    # Reproduction (adult classes first_ovip_class..k)
    ovip_fem <- sum(r_adult[p$first_ovip_class:k]) * adlt_delk * p$sex_ratio
    succ <- 1
    ova <- sum(p$fecundity * r_adult[p$first_ovip_class:k])
    ova <- max(0.000001, ova * adlt_delk * phi_temp * p$sex_ratio)
    # Random search for hosts (Eq. 4/5); uses previous day's total adults
    if (tot_adult_prev > 0) {
      succ <- 1 - exp(-p$search_const * dd * p$hosts / tot_adult_prev)
      ova <- max(0.00001, ova * succ)
    }

    # Eggs (aged with ambient degree-days, as in the code)
    d <- gp_delay_no_plr(ova, r_egg, p$egg_del, dd, k)
    r_egg <- d$r * max(0.00001, egg_viab)
    egg_out <- d$vout

    # Larvae (calendar time: del = 6 days, dt = 1 day)
    vin_l <- egg_out * egg_delk
    d <- gp_delay_no_plr(vin_l, r_larva, p$larva_del, 1, k)
    r_larva <- d$r * p$larval_surv
    larv_out <- d$vout

    # Pupae
    vin_p <- max(0.0000001, larv_out * larv_delk)
    d <- gp_delay_no_plr(vin_p, r_pupa, p$pupa_del, dd, k)
    r_pupa <- d$r * temp_surv
    pup_out <- d$vout * p$pupa_del / k

    # Adults (outflow from last class = senescence death)
    vin_a <- max(0.000001, pup_out)
    d <- gp_delay_no_plr(vin_a, r_adult, p$adult_del, dd, k)
    r_adult <- d$r * temp_surv
    adlt_out <- d$vout * adlt_delk

    # --- SWmod: summarise population (numbers)
    tot[t, ] <- c(sum(r_egg) * egg_delk, sum(r_larva) * larv_delk,
                  sum(r_pupa) * pup_delk, sum(r_adult) * adlt_delk)
    tot_adult_prev <- tot[t, "adults"]

    derived[t, ] <- c(0.5 * tot[t, "adults"], ovip_fem, succ, ovip_fem * succ)
    out[t, ] <- c(Tm, aver5[t], dd, temp_mort, temp_surv, egg_viab, phi_temp,
                  ova, vin_l, vin_p, vin_a, adlt_out, NA)
  }

  daily <- data.frame(date = weather$date, out[, -13], tot, derived)
  if ("rain" %in% names(weather)) daily$rain <- weather$rain
  daily
}


# ---- Annual summaries (GIS-style outputs) ------------------------------------
# mu_cold / mu_hot: sums of daily temp_mort for tmean < 27.2 / >= 27.2.
#   Calendar-year sums as accumulated in the code, plus the paper's
#   1 July (y-1) - 30 June (y) version of mu_cold, reported for year y.
gp_annual <- function(daily, p = gp_params()) {
  daily$year <- as.integer(format(daily$date, "%Y"))
  m <- as.integer(format(daily$date, "%m"))
  daily$season_year <- daily$year + (m >= 7)   # Jul-Dec counted toward next year
  cold <- daily$tmean < p$T_opt

  ann <- aggregate(cbind(eggs_laid, larvae_in, pupae_in, adults_in) ~ year, daily, sum)
  names(ann)[-1] <- c("eggs_per_yr", "larvae_per_yr", "pupae_per_yr", "adults_per_yr")
  ann$mu_cold_cal <- tapply(daily$temp_mort * cold,  daily$year, sum)[as.character(ann$year)]
  ann$mu_hot      <- tapply(daily$temp_mort * !cold, daily$year, sum)[as.character(ann$year)]
  mc <- tapply(daily$temp_mort * cold, daily$season_year, sum)
  ann$mu_cold_julJun <- mc[as.character(ann$year)]
  if ("rain" %in% names(daily)) {
    ann$rain <- tapply(daily$rain, daily$year, sum, na.rm = TRUE)[as.character(ann$year)]
    ann$I_rain <- pmin(1, ann$rain / p$rain_threshold)
    ann$pupae_rain_corrected <- ann$pupae_per_yr * ann$I_rain
  }
  # First year is burn-in in the paper; first mu_cold_julJun year is incomplete
  ann$burn_in <- ann$year == min(ann$year)
  ann
}


# ---- Weather helpers ---------------------------------------------------------

# Synthetic daily temperatures (for testing only; NOT real station data).
gp_synthetic_weather <- function(start = "1961-01-01", years = 5,
                                 mean_T = 23.5, amp = 6.5, peak_doy = 200,
                                 sd_noise = 2, seed = 1) {
  set.seed(seed)
  date <- seq(as.Date(start), by = "day", length.out = round(365.25 * years))
  doy <- as.integer(format(date, "%j"))
  tmean <- mean_T + amp * cos(2 * pi * (doy - peak_doy) / 365.25) +
    as.numeric(stats::filter(rnorm(length(date), 0, sd_noise), 0.6, "recursive")) * 0.8
  data.frame(date = date, tmean = tmean, rain = rgamma(length(date), 0.3, 0.3))
}

# Read a NOAA GHCN-Daily "access" CSV (e.g. McAllen USC00415701) downloaded from
# https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/
# TMAX/TMIN/PRCP are stored in tenths (deg C, mm); check before use.
gp_read_ghcn <- function(file) {
  w <- utils::read.csv(file, stringsAsFactors = FALSE)
  w$date <- as.Date(w$DATE)
  tmax <- w$TMAX / 10; tmin <- w$TMIN / 10
  data.frame(date = w$date, tmean = (tmax + tmin) / 2,
             rain = ifelse(is.na(w$PRCP), 0, w$PRCP / 10))
}


# Combine a primary and backup GHCN station following Gutierrez et al. (2019,
# Appendix S1): (1) fill primary gaps from the backup station; (2) fill single
# gaps with the mean of the preceding and following day when both are valid;
# (3) fill anything left with the mean of the preceding and following 5 days.
gp_prepare_ghcn <- function(primary_file, backup_file = NULL,
                            start = NULL, end = NULL) {
  rd <- function(f) {
    w <- utils::read.csv(f, stringsAsFactors = FALSE)
    data.frame(date = as.Date(w$DATE),
               tmax = as.numeric(w$TMAX) / 10, tmin = as.numeric(w$TMIN) / 10,
               rain = as.numeric(w$PRCP) / 10)
  }
  pri <- rd(primary_file)
  if (is.null(start)) start <- min(pri$date)
  if (is.null(end))   end   <- max(pri$date)
  w <- merge(data.frame(date = seq(as.Date(start), as.Date(end), by = "day")),
             pri, all.x = TRUE)
  vars <- c("tmax", "tmin", "rain")
  w$source <- ifelse(is.na(w$tmax) | is.na(w$tmin), NA, "primary")

  if (!is.null(backup_file)) {
    bk <- rd(backup_file)
    b <- merge(w["date"], bk, all.x = TRUE)
    for (v in vars) {
      fill <- is.na(w[[v]]) & !is.na(b[[v]])
      w[[v]][fill] <- b[[v]][fill]
    }
    w$source[is.na(w$source) & !is.na(w$tmax) & !is.na(w$tmin)] <- "backup"
  }
  for (v in vars) {
    x <- w[[v]]; n <- length(x)
    i <- which(is.na(x)); i <- i[i > 1 & i < n]
    ok <- !is.na(x[i - 1]) & !is.na(x[i + 1])
    x[i[ok]] <- (x[i[ok] - 1] + x[i[ok] + 1]) / 2
    for (j in which(is.na(x))) {
      nb <- x[c(max(1, j - 5):(j - 1), (j + 1):min(n, j + 5))]
      if (any(!is.na(nb))) x[j] <- mean(nb, na.rm = TRUE)
    }
    w[[v]] <- x
  }
  w$source[is.na(w$source)] <- "interpolated"
  w$tmean <- (w$tmax + w$tmin) / 2
  w
}


# ---- Demo --------------------------------------------------------------------
if (sys.nframe() == 0) {
  p <- gp_params()

  # 1. Check the delay engine: a pulse with constant dt should exit with
  #    Erlang mean = del and SD ~ del/sqrt(k).
  r <- rep(0, p$k); outflow <- numeric(400)
  r <- gp_delay_no_plr(1, r, del = 95, dt = 1)$r            # pulse of 1 at step 1
  for (s in 2:400) { d <- gp_delay_no_plr(0, r, 95, 1); r <- d$r; outflow[s] <- d$vout * 95 / p$k }
  tt <- seq_along(outflow)
  m <- sum(tt * outflow) / sum(outflow)
  cat(sprintf("Delay check: total out = %.4f, mean transit = %.1f (target ~95), SD = %.1f (target ~%.1f)\n",
              sum(outflow), m, sqrt(sum((tt - m)^2 * outflow) / sum(outflow)), 95 / sqrt(p$k)))

  # 2. Run on synthetic weather (warm subtropical vs cooler temperate)
  warm <- gp_synthetic_weather(mean_T = 23.5, amp = 6.5)   # McAllen-like
  cool <- gp_synthetic_weather(mean_T = 18.5, amp = 9.5)   # Uvalde-like
  for (nm in c("warm", "cool")) {
    daily <- gp_simulate(get(nm), p)
    cat("\n", nm, "(synthetic weather)\n")
    print(round(gp_annual(daily, p)[, c("year", "pupae_per_yr", "mu_cold_cal",
                                        "mu_cold_julJun", "mu_hot")], 2))
  }
}
