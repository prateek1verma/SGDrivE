# Scenario 11: precision-guided SIT (pgSIT) at a constant 25 degC
#
# pgSIT males are released every week into a wild screwworm population at
# equilibrium. In pgSIT a genetic cross makes the released insects sterile
# males (females die during development), so the release is male-only and no
# irradiation is needed. The strategy follows "pgSIT" of the mgdrive-nws
# project (analysis 01; the two-genotype cube of analysis 03), with the
# values set in step 1:
#
# * 26 weekly releases from the first day, followed for two years;
# * released males: male_ratio (here 1) per wild adult per release, no
#   females;
# * pgSIT males mate with competitiveness eta = 0.95 relative to wild males
#   (not irradiated, so close to wild males);
# * offspring they sire survive with probability 1 - penetrance (here 0:
#   complete sterility);
# * no lifespan cost (no irradiation).
#
# For comparison, mgdrive-nws releases 5 males per wild adult. Release
# settings other than eta, lifespan and the male-only release are as in
# scenario 10 (separator SIT), so the two can be compared.
#
# Two genotypes are enough, as in the mgdrive-nws analysis 03 cube: W (wild)
# and S (released pgSIT male). Offspring of S males are recorded as W and
# survive with probability 1 - penetrance.
#
# Run time: a few seconds.

library(SGDrivE)
library(ggplot2) # plots

# 1. Settings -----------------------------------------------------------------
temperature <- 25 # constant temperature (degC)
days <- 730 # simulated days (two years)
wild_adults <- 20000 # adults at equilibrium (mgdrive-nws AdPopEQ)
n_rep <- 20 # stochastic replicates
seed <- 20260811 # random seed (as in mgdrive-nws)

# Release strategy (cf. mgdrive-nws R/parameters.R, STRATEGIES_MAIN$pgsit)
eta_pgsit <- 0.95 # mating competitiveness of pgSIT males
penetrance <- 1 # share of offspring of pgSIT males that die (1 = all)
lifespan_pgsit <- 1 # fraction of the wild adult lifespan kept (1 = no cost)
male_ratio <- 1 # released males per wild adult, per release
n_releases <- 26 # weekly releases
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)

# 2. Inheritance cube: wild type W and pgSIT S --------------------------------
# ih[g, m, k]: offspring genotype k of a female g mated to a male m. All
# offspring are recorded as W; those sired by S males survive with
# probability 1 - penetrance (tau). No S females are released or born.
genotypes <- c("W", "S")
dims <- list(genotypes, genotypes, genotypes)
ih <- array(0, dim = c(2, 2, 2), dimnames = dims)
ih[, , "W"] <- 1
tau <- array(1, dim = c(2, 2, 2), dimnames = dims)
tau[, "S", ] <- 1 - penetrance # offspring sired by pgSIT males
eta <- matrix(1, 2, 2, dimnames = list(genotypes, genotypes))
eta[, "S"] <- eta_pgsit # rows: female, cols: male

cube <- sg_cube(list(
  genotypesID = genotypes, ih = ih, tau = tau, eta = eta,
  phi = c(W = 0.5, S = 0.5), # half of emerging adults are female
  xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), # emergence success
  s = c(W = 1, S = 1), # fecundity (S females never occur)
  omega = c(W = 1, S = 1),
  wildType = "W", releaseType = "S"
), lifespan_fraction = c(W = 1, S = lifespan_pgsit))
cube

# 3. Weather and the wild population ------------------------------------------
# A region holding wild_adults adults (half of them female) at the
# equilibrium of the constant temperature, without releases.
p <- sg_params(female_background_mortality=0.1,male_background_mortality=0.12)
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
wild

# 4. Releases: pgSIT males only -----------------------------------------------
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
releases <- sg_releases("S", stage = "adult", sex = "male", times = release_days,
                        number = male_ratio * wild_adults)
releases

# 5. Model with releases, started at the wild equilibrium ---------------------
model <- sg_model(p, cube, forcing, releases = releases, area = wild$area)
init <- sg_init(model, "equilibrium")

# 6. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic reference
cb <- simulate_cb(model, init, n_rep = n_rep, seed = seed) # stochastic

# 7. Wild-type adult females --------------------------------------------------
# wild_type counts all adult females of genotype W, unmated or mated (also
# those mated to pgSIT males), as in the mgdrive-nws figures.
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
    title = paste0("pgSIT: ", n_releases, " weekly releases of pgSIT males (no females)"),
    subtitle = paste0(
      temperature, " degC; eta = ", eta_pgsit, ", penetrance = ", penetrance, "; ",
      n_rep, " stochastic replicates (blue), ODE (black); shaded: release period; ",
      "eliminated: ", n_eliminated, " of ", n_rep
    )
  ) +
  theme_minimal()

# Notes
# * pgSIT males are not irradiated, so their competitiveness (0.95) is much
#   higher than that of irradiated males (scenarios 07 and 10), and they have
#   no lifespan cost.
# * Wild females mated to pgSIT males still search for wounds and lay eggs
#   that do not hatch, competing with fertile females (D-006).
# * With penetrance < 1 some offspring of pgSIT males survive; they are
#   recorded as wild type, as in the mgdrive-nws analysis 03 cube.
