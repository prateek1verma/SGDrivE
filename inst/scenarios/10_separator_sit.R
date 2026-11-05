# Scenario 10: separator SIT at a constant 25 degC
#
# Sterile males are released every week into a wild screwworm population at
# equilibrium, together with a small contamination of sterile females that
# the sex separation (a separator machine or sorting step) fails to remove.
# The strategy follows "Separator SIT" of the mgdrive-nws project (analysis
# 01), with the values set in step 1:
#
# * 26 weekly releases from the first day, followed for two years;
# * released males: male_ratio (here 1) per wild adult per release;
# * released females: female_percent (here 0.1%) of the male release, i.e.
#   0.001 per wild adult or 20 per release;
# * sterile males mate with competitiveness eta = 0.6 relative to wild males,
#   and eggs they sire never hatch; sterile females produce no eggs;
# * no lifespan cost of irradiation.
#
# For comparison, mgdrive-nws uses 5 males per wild adult, 0.01% females
# (0.1% for its CGSS strategy) and eta = 0.4 for separator SIT.
#
# As in scenario 07, two genotypes are enough: W (wild) and S (sterile,
# released). Everything else (life cycle, temperature, wound search, mating)
# is ScrewDrivE's model, so the population responds differently from the
# MGDrivE version (see the notes at the end).
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

# Release strategy (cf. mgdrive-nws R/parameters.R, STRATEGIES_MAIN$separator_sit)
eta_sterile <- 0.5 # mating competitiveness of sterile males
male_ratio <- 1 # nominal released males per wild adult, per release
female_percent <- 0.01 # released females (contamination), % of the male release
mating_capacity <- 1 # not used here: no sequestration correction (step 4)
n_releases <- 26 # weekly releases
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)

# 2. Inheritance cube: wild type W and sterile S ------------------------------
# ih[g, m, k]: offspring genotype k of a female g mated to a male m. Only
# W x W crosses have viable offspring, and they are W; the other crosses are
# set to W too but never hatch (tau = 0 for every cross with an S male, and
# s = 0 for S females, which produce no eggs).
genotypes <- c("W", "S")
dims <- list(genotypes, genotypes, genotypes)
ih <- array(0, dim = c(2, 2, 2), dimnames = dims)
ih[, , "W"] <- 1
tau <- array(1, dim = c(2, 2, 2), dimnames = dims)
tau[, "S", ] <- 0 # eggs sired by sterile males do not hatch
# Alternative: S females lay eggs that do not hatch (set the "S" female
# slice of tau to 0).
eta <- matrix(1, 2, 2, dimnames = list(genotypes, genotypes))
eta[, "S"] <- eta_sterile # females prefer wild males (rows: female, cols: male)

cube <- sg_cube(list(
  genotypesID = genotypes, ih = ih, tau = tau, eta = eta,
  phi = c(W = 0.5, S = 0.5), # half of emerging adults are female
  xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), # emergence success
  s = c(W = 1, S = 0), # sterile females produce no eggs and do not search (D-031)
  omega = c(W = 1, S = 1), # no lifespan cost of irradiation
  wildType = "W", releaseType = "S"
), lifespan_fraction = c(W = 1, S = 0.6))
cube

# 3. Weather and the wild population ------------------------------------------
# A region holding wild_adults adults (half of them female) at the
# equilibrium of the constant temperature, without releases.
# Background mortalities: raise them for field mortality.
p <- sg_params(female_background_mortality = 0.1, male_background_mortality = 0.12)
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
wild

# 4. Releases -----------------------------------------------------------------
# Males: male_ratio per wild adult. Females: the contamination,
# female_percent of the male release. The released females are so few that
# the males they mate with are not subtracted from the male release (in
# mgdrive-nws this sequestration correction changes it by less than 0.01%).
female_ratio <- male_ratio * female_percent / 100
effective_male_ratio <- max(0, male_ratio)
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
    title = paste0(
      "Separator SIT: ", n_releases, " weekly releases of sterile males with ",
      female_percent, "% females"
    ),
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
# * Wild females mated to sterile males still search for wounds and lay
#   non-viable eggs, competing with fertile females (D-006); released sterile
#   females produce no eggs and do not search (D-031).
# * As in MGDrivE, a male can mate any number of times, so released males are
#   not used up by matings.
