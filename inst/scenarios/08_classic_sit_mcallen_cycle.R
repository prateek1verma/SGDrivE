# Scenario 08: classic SIT under the mean annual temperature cycle of McAllen
#
# The classic SIT releases of scenario 07 (mgdrive-nws analysis 01, Classic
# SIT, parameter set "main"), but in a seasonal climate instead of a constant
# temperature. The temperature on each calendar day is the mean of that day
# over the 24 years of McAllen, Texas weather in the package (1960-1983), and
# the same yearly cycle repeats for the two simulated years.
#
# * The wild population follows the seasons: few adults in winter, many in
#   summer. The region is sized so that adult females average
#   wild_adults / 2 over a year of the cycle.
# * There is no equilibrium in a changing climate, so the population starts
#   from five years of the cycle without releases (a burn-in), which settles
#   on the yearly pattern.
# * Releases are sized from the annual mean of wild_adults, as in mgdrive-nws,
#   and start on 1 January, when the wild population is near its low.
#
# The plots show the temperature cycle, then the wild-type adult females of
# the stochastic replicates with the deterministic (ODE) runs with and without
# releases.
#
# Run time: a few seconds.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Settings -----------------------------------------------------------------
days <- 730 # simulated days (two years from 1 January)
start <- as.Date("2001-01-01") # first simulated day (any non-leap year)
wild_adults <- 20000 # mean adults over a year (mgdrive-nws AdPopEQ)
burn_in_years <- 5 # years of the cycle before the start, without releases
n_rep <- 20 # stochastic replicates
seed <- 20260811 # random seed (as in mgdrive-nws)

# Release strategy (mgdrive-nws R/parameters.R, STRATEGIES_MAIN$classic_sit)
eta_sterile <- 0.4 # mating competitiveness of sterile males
male_ratio <- 0.1 # nominal released males per wild adult, per release
female_percent <- 100 # released females, % of the nominal male release
mating_capacity <- 2 # females a sterile male can mate with (sequestration)
n_releases <- 26 # weekly releases
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)

# 2. Inheritance cube: wild type W and sterile S (as in scenario 07) ----------
genotypes <- c("W", "S")
dims <- list(genotypes, genotypes, genotypes)
ih <- array(0, dim = c(2, 2, 2), dimnames = dims)
ih[, , "W"] <- 1 # only W x W crosses hatch, and they give W
tau <- array(1, dim = c(2, 2, 2), dimnames = dims)
tau[, "S", ] <- 0 # eggs sired by sterile males do not hatch
eta <- matrix(1, 2, 2, dimnames = list(genotypes, genotypes))
eta[, "S"] <- eta_sterile # females prefer wild males (rows: female, cols: male)

cube <- sg_cube(list(
  genotypesID = genotypes, ih = ih, tau = tau, eta = eta,
  phi = c(W = 0.5, S = 0.5), # half of emerging adults are female
  xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), # emergence success
  s = c(W = 1, S = 0), # sterile females lay no viable eggs
  omega = c(W = 1, S = 1), # no lifespan cost of irradiation
  wildType = "W", releaseType = "S"
))

# 3. Mean annual temperature cycle at McAllen ---------------------------------
# For every calendar day ("01-01" ... "12-31"), the mean daily temperature
# over 1960-1983, and the 10% and 90% quantiles to show the spread between
# years. 29 February only occurs in leap years and is not used here.
mcallen <- sg_example_weather()
calendar_day <- format(mcallen$date, "%m-%d")
cycle <- tapply(mcallen$tmean, calendar_day, mean)
low <- tapply(mcallen$tmean, calendar_day, quantile, probs = 0.1)
high <- tapply(mcallen$tmean, calendar_day, quantile, probs = 0.9)

dates <- start + seq_len(days) - 1
weather <- data.frame(date = dates, tmean = as.numeric(cycle[format(dates, "%m-%d")]))
p <- sg_params()
forcing <- sg_forcing(weather, p)

# 4. Region size --------------------------------------------------------------
# Without releases the population is proportional to the area, so: run the
# cycle in a region of 1 km^2 (after the burn-in), then scale the area so that
# adult females average wild_adults / 2 over the first year.
unit <- sg_model(p, cube, forcing, area = 1)
unit_females <- sg_females(simulate_ode(unit, sg_init(unit, "burn_in", years = burn_in_years)))
area <- (wild_adults / 2) / mean(unit_females$females[unit_females$day %in% 1:365])
area # region area in km^2

# 5. Releases (sizes as in scenario 07) ---------------------------------------
female_ratio <- male_ratio * female_percent / 100
effective_male_ratio <- max(0, male_ratio - female_ratio / mating_capacity)
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
releases <- c(
  sg_releases("S", stage = "adult", sex = "male", times = release_days,
              number = effective_male_ratio * wild_adults),
  sg_releases("S", stage = "adult", sex = "female", times = release_days,
              number = female_ratio * wild_adults)
)

# 6. Models with and without releases, started after the burn-in --------------
model <- sg_model(p, cube, forcing, releases = releases, area = area)
no_release <- sg_model(p, cube, forcing, area = area)
init <- sg_init(model, "burn_in", years = burn_in_years)

# 7. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic, with releases
ode_wild <- simulate_ode(no_release, sg_init(no_release, "burn_in", years = burn_in_years))
cb <- simulate_cb(model, init, n_rep = n_rep, seed = seed) # stochastic

wt_cb <- sg_females(cb)
wt_ode <- sg_females(ode)
wt_wild <- sg_females(ode_wild)

final <- wt_cb[wt_cb$day == max(wt_cb$day), ]
n_eliminated <- sum(final$wild_type < 1) # replicates without wild females
n_eliminated

# 8. Plot: the temperature cycle ----------------------------------------------
temps <- data.frame(
  day = seq_len(days),
  mean = as.numeric(cycle[format(dates, "%m-%d")]),
  low = as.numeric(low[format(dates, "%m-%d")]),
  high = as.numeric(high[format(dates, "%m-%d")])
)
ggplot(temps, aes(day)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "#D55E00", alpha = 0.15) +
  geom_line(aes(y = mean), colour = "#D55E00", linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, days, by = 90)) +
  coord_cartesian(xlim = c(0, days)) +
  labs(
    x = "days since 1 January", y = "daily mean temperature (degC)",
    title = "McAllen, Texas: mean annual temperature cycle (1960-1983)",
    subtitle = "line: mean of each calendar day over 24 years; band: 10-90% of years"
  ) +
  theme_minimal()

# 9. Plot: wild-type adult females --------------------------------------------
# A release on day d is added at the start of that day, i.e. at time d - 1.
release_window <- c(min(release_days) - 1, max(release_days) - 1 + release_interval)
ggplot(wt_cb, aes(day, wild_type, group = replicate)) +
  annotate("rect",
    xmin = release_window[1], xmax = release_window[2],
    ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.12
  ) +
  geom_line(data = wt_wild, aes(day, wild_type), colour = "grey45",
            linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  geom_line(colour = "#0072B2", alpha = 0.3, linewidth = 0.4) +
  geom_line(data = wt_ode, aes(day, wild_type), colour = "black",
            linewidth = 0.8, inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(0, days, by = 90)) +
  # pseudo-log scale: shows the seasonal peaks and the decline towards zero
  scale_y_continuous(
    transform = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = c(0, 10^(0:5)), labels = c("0", "1", "10", "100", "1,000", "10,000", "100,000")
  ) +
  coord_cartesian(xlim = c(0, days)) +
  labs(
    x = "days since 1 January",
    y = "wild-type adult females in the region (pseudo-log scale)",
    title = "Classic SIT under the McAllen mean annual temperature cycle",
    subtitle = paste0(
      n_rep, " stochastic replicates (blue), ODE with releases (black), ",
      "without releases (dashed); shaded: releases; eliminated: ",
      n_eliminated, " of ", n_rep
    )
  ) +
  theme_minimal()
