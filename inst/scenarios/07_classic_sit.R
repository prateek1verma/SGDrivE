# Scenario 07: classic sterile insect technique (SIT) at a constant 25 degC
#
# Irradiated sterile males and females are released every week into a wild
# screwworm population at equilibrium. The release follows the "Classic SIT"
# strategy of the mgdrive-nws project (analysis 01, parameter set "main"):
#
# * 26 weekly releases from the first day, followed for two years;
# * released males: 5 per wild adult, less the males taken up by the released
#   females (mating capacity 2 per sterile male), so 2.5 per wild adult;
# * released females: as many as the nominal male release, 5 per wild adult;
# * sterile males mate with competitiveness eta = 0.4 relative to wild males,
#   and eggs they sire never hatch; sterile females lay no viable eggs;
# * no lifespan cost of irradiation.
#
# mgdrive-nws uses a six-genotype sex-explicit cube. Because sterile females
# have no viable eggs and eggs sired by sterile males never hatch, no
# heterozygote can ever be born, so two genotypes are enough here: W (wild)
# and S (sterile, released). Everything else (life cycle, temperature, wound
# search, mating) is ScrewDrivE's model, so the population responds differently
# from the MGDrivE version (see the notes at the end).
#
# Run time: a few seconds.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Settings -----------------------------------------------------------------
temperature <- 25 # constant temperature (degC)
days <- 730 # simulated days (two years)
wild_adults <- 20000 # adults at equilibrium (mgdrive-nws AdPopEQ)
n_rep <- 20 # stochastic replicates
seed <- 20260811 # random seed (as in mgdrive-nws)

# Release strategy (mgdrive-nws R/parameters.R, STRATEGIES_MAIN$classic_sit)
eta_sterile <- 0.4 # mating competitiveness of sterile males
male_ratio <- 1 # nominal released males per wild adult, per release
female_percent <- 100 # released females, % of the nominal male release
mating_capacity <- 2 # females a sterile male can mate with (sequestration)
n_releases <- 26 # weekly releases
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)

# 2. Inheritance cube: wild type W and sterile S ------------------------------
# ih[g, m, k]: offspring genotype k of a female g mated to a male m. Only
# W x W crosses have viable offspring, and they are W; the other crosses are
# set to W too but never hatch (tau = 0 for every cross with an S male, and
# s = 0 for S females).
genotypes <- c("W", "S")
dims <- list(genotypes, genotypes, genotypes)
ih <- array(0, dim = c(2, 2, 2), dimnames = dims)
ih[, , "W"] <- 1
tau <- array(1, dim = c(2, 2, 2), dimnames = dims)
tau[, "S", ] <- 0 # eggs sired by sterile males do not hatch
# tau["S", , ] <- 0 # eggs sired by sterile females do not hatch
eta <- matrix(1, 2, 2, dimnames = list(genotypes, genotypes))
eta[, "S"] <- eta_sterile # females prefer wild males (rows: female, cols: male)

cube <- sg_cube(list(
  genotypesID = genotypes, ih = ih, tau = tau, eta = eta,
  phi = c(W = 0.5, S = 0.5), # half of emerging adults are female
  xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), # emergence success
  s = c(W = 1, S = 0), # sterile females lay no viable eggs
  omega = c(W = 1, S = 1), # no lifespan cost of irradiation
  wildType = "W", releaseType = "S"
), lifespan_fraction = c(W = 1, S = 1))
cube

# 3. Weather and the wild population ------------------------------------------
# A region holding wild_adults adults (half of them female) at the
# equilibrium of a constant 25 degC, without releases.
p <- sg_params(female_background_mortality=0.1,male_background_mortality=0.12)
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
wild

# 4. Releases -----------------------------------------------------------------
# Release sizes as in mgdrive-nws: released females take up sterile-male
# mating capacity, which reduces the effective male release
# (effective ratio = male ratio - female ratio / mating capacity).
female_ratio <- male_ratio * female_percent / 100
effective_male_ratio <- max(0, male_ratio - female_ratio / mating_capacity)
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
releases <- c(
  sg_releases("S", stage = "adult", sex = "male", times = release_days,
              number = effective_male_ratio * wild_adults),
  sg_releases("S", stage = "adult", sex = "female", times = release_days,
              number = female_ratio * wild_adults)
)
releases

# 5. Model with releases, started at the wild equilibrium ---------------------
model <- sg_model(p, cube, forcing, releases = releases, area = wild$area)
init <- sg_init(model, "equilibrium")

# 6. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic reference
cb <- simulate_cb(model, init, n_rep = n_rep, seed = seed) # stochastic

# 7. Wild-type adult females --------------------------------------------------
# wild_type counts all adult females of genotype W, unmated or mated (also
# those mated to sterile males), as in the mgdrive-nws figures.
wt_cb <- sg_females(cb)
wt_ode <- sg_females(ode)

# Replicates with no wild-type adult female left at the end
final <- wt_cb[wt_cb$day == max(wt_cb$day), ]
n_eliminated <- sum(final$wild_type < 1)
n_eliminated

# Time in days since the start (day 0 = initial state). A release on day d is
# added at the start of that day, i.e. at time d - 1.
release_window <- c(min(release_days) - 1, max(release_days) - 1 + release_interval)
ggplot(wt_cb, aes(day, wild_type, group = replicate)) +
  annotate("rect",
    xmin = release_window[1], xmax = release_window[2],
    ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.12
  ) +
  geom_line(colour = "#0072B2", alpha = 0.45, linewidth = 0.4) +
  geom_line(data = wt_ode, aes(day, wild_type), colour = "black",
            linewidth = 0.8, inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(0, days, by = 90)) +
  coord_cartesian(xlim = c(0, days), ylim = c(0, 0.5 * wild_adults + 1000)) +
  labs(
    x = "days", y = "wild-type adult females in the region",
    title = "Classic SIT: 26 weekly releases of sterile males and females",
    subtitle = paste0(
      temperature, " degC; ", n_rep, " stochastic replicates (blue), ODE ",
      "(black); shaded: release period; eliminated: ", n_eliminated, " of ",
      n_rep
    )
  ) +
  theme_minimal()

# Notes: differences from the MGDrivE (mgdrive-nws) version
# * Density dependence is wound search, not larval carrying capacity:
#   females compete for wounds, so the growth rate at low density (R0, see
#   sg_r0()) is much higher than at equilibrium.
# * Development and mortality follow the temperature-driven model, so stage
#   durations and adult survival come from the 25 degC temperature, not from
#   MGDrivE's fixed durations and daily mortality.
# * Females mated to sterile males, and sterile females, still search for
#   wounds and lay (non-viable) eggs, competing with fertile females (D-006).
# * As in MGDrivE, a male can mate any number of times, which is why the
#   sequestration correction of the male release is kept.
