# Scenario 20: calibration to Uvalde County myiasis, January-June 1962 (DRAFT)
#
# Fits the wild-type population of one county, Uvalde, Texas, to the weekly
# screwworm infestations recorded there in weeks 1-26 of 1962 (1 January to
# 1 July). The sterile-fly programme started on 6 July 1962, so the fit period
# has no releases. Weeks 27-52 are held out for the release test.
#
# Two parameters are fitted: the wound density (wounds per km^2) and a
# development-rate scale, which multiplies the initial slope of the
# development-rate curve (dev_rate_slope) for all stages.
#
# DRAFT: decisions still open are marked DECISION. Each has a default, so the
# script runs, but the defaults need checking before any result is reported.
#
# Data (not shipped with the package; see data-raw/texas_nws/):
# * weather: weather_1961_1962/station_daily_1961_1962_filled.csv (gap-filled)
# * cases:   uvalde_weekly_1962.csv (county weekly counts, texasScrewwormdat1962.xlsx)
#
# Run time: about a minute per fit.

library(ScrewDrivE)

# 0. Paths ---------------------------------------------------------------------
# Run from the repository root, or change `data_dir`.
data_dir <- "data-raw/texas_nws"
weather_file <- file.path(data_dir, "weather_1961_1962", "station_daily_1961_1962_filled.csv")
cases_file <- file.path(data_dir, "uvalde_weekly_1962.csv")

# 1. Node and weather ----------------------------------------------------------
# Uvalde County: total area 1,559 sq mi = 4,040 km^2 (Wikipedia, county entry).
# Station USC00419265 (UVALDE), 2.1 km from the county centroid.
node_area_km2 <- 4040
station <- "USC00419265"

wx_all <- read.csv(weather_file, stringsAsFactors = FALSE)
wx <- subset(wx_all, station_id == station)
wx <- data.frame(date = as.Date(wx$date), tmax = wx$tmax, tmin = wx$tmin, rain = wx$prcp)
wx <- wx[order(wx$date), ]
stopifnot(!anyNA(wx$tmax), !anyNA(wx$tmin))

# The run covers 1 January to 1 July 1962 (week 26 ends on 1 July). The
# burn-in uses the year before, so late 1961 sets the overwintering population.
burn_in_weather <- wx[wx$date < as.Date("1962-01-01"), ]
run_weather <- wx[wx$date >= as.Date("1962-01-01") & wx$date <= as.Date("1962-07-01"), ]

# 2. Observed weekly counts ----------------------------------------------------
# DECISION 1 (made): fit to the weekly county series, not the monthly
# statewide totals. Weekly data give 26 points, not 6, and the weeks line up
# with sg_cases(period = "week") because the run starts on 1 January.
obs <- read.csv(cases_file, stringsAsFactors = FALSE)
obs <- subset(obs, week <= 26)

# 3. Parameters, cube and model ------------------------------------------------
# DECISION 2: parameter preset. "gp2014" is the Gutierrez & Ponti (2014)
# configuration (scenario 09); "default" is the package's current defaults.
preset <- "gp2014" # switch to "default" to compare
p <- sg_params(preset = preset)

cube <- sg_cube(list( # one wild genotype, as in scenario 09
  genotypesID = "W", wildType = "W",
  ih = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  tau = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  eta = matrix(1, 1, 1, dimnames = list("W", "W")),
  phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1), omega = c(W = 1)
))

# No sterile releases before 6 July 1962.
releases <- NULL

# DECISION 3: reporting scale. The default is 1 (observed cases taken as given).
# A free factor would absorb any mismatch in size, so it is held at 1 here.
reporting_scale <- 1

# 4. Simulated weekly county counts --------------------------------------------
# Development rates enter the degree-day forcing, so the forcing is rebuilt
# for each value of the development scale.
simulate_weekly <- function(wound_density, dev_scale = 1) {
  p_w <- p
  p_w$wound_density <- wound_density # wounds per km^2 are a parameter
  p_w$dev_rate_slope <- p$dev_rate_slope * dev_scale # development-rate scale
  forcing <- sg_forcing(run_weather, p_w)
  m <- sg_model(p_w, cube, forcing, releases = releases, area = node_area_km2)
  init <- sg_init(m, "burn_in", years = 1L, weather = burn_in_weather) # DECISION 4: 1-year burn-in
  sim <- simulate_ode(m, init)
  cases <- sg_cases(sim, period = "week", per = "node") # county counts
  cases$cases[seq_len(nrow(obs))] * reporting_scale
}

# 5. Fit -----------------------------------------------------------------------
# DECISION 5: objective. Log(x + 1) squared error, since the weekly counts span
# several orders of magnitude. A Poisson or negative binomial likelihood is an
# alternative.
# DECISION 6: bounds. Wound density 0.001 to 10,000 per km^2; development scale
# 0.2 to 20 times the GP2014 slope (the fit with an upper bound of 5 stopped at
# the bound). Both are on the log scale. A scale near 6 is far from the
# literature value, so check it against the other parameters before use.
# fit_dev_scale = FALSE holds the development scale at 1 (the GP2014 or default
# slope) and fits the wound density alone, as a test of whether the scale is needed.
fit_dev_scale <- TRUE
objective <- function(par) {
  dev_scale <- if (fit_dev_scale) exp(par[2]) else 1
  sim_cases <- simulate_weekly(exp(par[1]), dev_scale)
  if (length(sim_cases) != nrow(obs) || anyNA(sim_cases)) return(1e6)
  sum((log(sim_cases + 1) - log(obs$cases + 1))^2)
}

if (fit_dev_scale) {
  fit <- optim(
    par = c(log(4), 0), fn = objective, method = "L-BFGS-B",
    lower = c(log(1e-3), log(0.2)), upper = c(log(1e4), log(20))
  )
  best_wound_density <- exp(fit$par[1])
  best_dev_scale <- exp(fit$par[2])
} else {
  fit1 <- optimize(function(x) objective(c(x, 0)), interval = log(c(1e-3, 1e4)))
  fit <- list(value = fit1$objective, convergence = 0L)
  best_wound_density <- exp(fit1$minimum)
  best_dev_scale <- 1
}
cat(sprintf(
  "Fitted wound density: %.4g per km^2; development scale: %.3g (objective %.4g, convergence %d)\n",
  best_wound_density, best_dev_scale, fit$value, fit$convergence
))

# 6. Comparison ----------------------------------------------------------------
sim_best <- simulate_weekly(best_wound_density, best_dev_scale)
comparison <- data.frame(
  week = obs$week, observed = obs$cases, simulated = round(sim_best, 2)
)
print(comparison)
cat(sprintf("Jan-Jun observed %d, simulated %.0f\n", sum(obs$cases), sum(sim_best)))

# DECISION 7: outputs. The draft prints the weekly table only. Plots, and the
# held-out weeks 27-52 (not used in the fit), come after the decisions.
