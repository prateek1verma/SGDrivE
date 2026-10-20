# Scenario 12: female-specific RIDL (fsRIDL) at a constant 25 degC
#
# fsRIDL males are released every week into a wild screwworm population at
# equilibrium. They carry a dominant transgene (allele H) that kills
# daughters as embryos but not sons, so the release is male-only, the males
# are fertile and not irradiated, and their sons carry the transgene on: a
# self-limiting effect that lasts a few generations after releases stop. The
# strategy follows "fsRIDL, not irradiated" of the mgdrive-nws project
# (analysis 01), with the values set in step 1:
#
# * 26 weekly releases from the first day, followed for two years;
# * released males: homozygous HH, male_ratio (here 1) per wild adult per
#   release, no females;
# * males carrying H mate with competitiveness eta = 0.9 relative to wild
#   males (0.8 in mgdrive-nws for irradiated fsRIDL males);
# * daughters carrying H die with probability female_lethality (here 1);
# * no lifespan cost and full fertility.
#
# For comparison, mgdrive-nws releases 5 males per wild adult. Other release
# settings are as in scenario 10 (separator SIT) and 11 (pgSIT).
#
# mgdrive-nws uses a six-genotype sex-explicit cube (sex in the genotype
# name). Here sex comes from phi, so three genotypes are enough: WW (wild),
# WH and HH. A daughter carrying H dies as an embryo, so for an offspring
# genotype carrying H a share phi * female_lethality of the eggs is not
# viable (tau) and the survivors are mostly sons (phi of that genotype).
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

# Release strategy (cf. mgdrive-nws R/parameters.R, STRATEGIES_MAIN$fsridl_norad)
eta_ridl <- 0.9 # mating competitiveness of males carrying H
female_lethality <- 1 # share of daughters carrying H that die (1 = all)
lifespan_ridl <- 1 # fraction of the wild adult lifespan kept by carriers
male_ratio <- 1 # released HH males per wild adult, per release
n_releases <- 26 # weekly releases
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)

# 2. Inheritance cube: WW, WH and HH ------------------------------------------
# Mendelian inheritance of the transgene H. phi0 is the share of daughters
# among offspring; for offspring genotypes carrying H, daughters die with
# probability female_lethality: tau removes those eggs, and phi is the share
# of daughters among the survivors (0 when all carrier daughters die).
genotypes <- c("WW", "WH", "HH")
dims <- list(genotypes, genotypes, genotypes)
gametes <- rbind(WW = c(W = 1, H = 0), WH = c(W = 0.5, H = 0.5), HH = c(W = 0, H = 1))
ih <- array(0, dim = c(3, 3, 3), dimnames = dims)
for (f in genotypes) {
  for (m in genotypes) {
    gf <- gametes[f, ]
    gm <- gametes[m, ]
    ih[f, m, ] <- c(gf[["W"]] * gm[["W"]], gf[["W"]] * gm[["H"]] + gf[["H"]] * gm[["W"]],
                    gf[["H"]] * gm[["H"]])
  }
}
phi0 <- 0.5
carriers <- c("WH", "HH")
tau <- array(1, dim = c(3, 3, 3), dimnames = dims)
tau[, , carriers] <- 1 - phi0 * female_lethality # carrier daughters die as embryos
phi <- c(WW = phi0, WH = 0, HH = 0)
phi[carriers] <- phi0 * (1 - female_lethality) / (1 - phi0 * female_lethality)
eta <- matrix(1, 3, 3, dimnames = list(genotypes, genotypes))
eta[, carriers] <- eta_ridl # rows: female, cols: male

cube <- sg_cube(list(
  genotypesID = genotypes, ih = ih, tau = tau, eta = eta, phi = phi,
  xiF = c(WW = 1, WH = 1, HH = 1), xiM = c(WW = 1, WH = 1, HH = 1), # emergence
  s = c(WW = 1, WH = 1, HH = 1), # full fertility
  omega = c(WW = 1, WH = 1, HH = 1),
  wildType = "WW", releaseType = "HH"
), lifespan_fraction = c(WW = 1, WH = lifespan_ridl, HH = lifespan_ridl))
cube

# 3. Weather and the wild population ------------------------------------------
# A region holding wild_adults adults (half of them female) at the
# equilibrium of the constant temperature, without releases.
p <- sg_params(female_background_mortality=0.1,male_background_mortality=0.12)
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
wild

# 4. Releases: HH males only --------------------------------------------------
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
releases <- sg_releases("HH", stage = "adult", sex = "male", times = release_days,
                        number = male_ratio * wild_adults)
releases

# 5. Model with releases, started at the wild equilibrium ---------------------
model <- sg_model(p, cube, forcing, releases = releases, area = wild$area)
init <- sg_init(model, "equilibrium")

# 6. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic reference
cb <- simulate_cb(model, init, n_rep = n_rep, seed = seed) # stochastic

# 7. Wild-type adult females --------------------------------------------------
# wild_type counts adult females of genotype WW (with complete female
# lethality there are no others), as in the mgdrive-nws figures.
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
    title = paste0("fsRIDL: ", n_releases, " weekly releases of HH males (no females)"),
    subtitle = paste0(
      temperature, " degC; eta = ", eta_ridl, ", female lethality = ", female_lethality,
      "; ", n_rep, " stochastic replicates (blue), ODE (black); shaded: release ",
      "period; eliminated: ", n_eliminated, " of ", n_rep
    )
  ) +
  theme_minimal()

# 8. The transgene among adult males ------------------------------------------
# Frequency of allele H among adult males: high during releases, then
# carried on by sons for a few generations and lost (self-limiting).

# h_cb <- sg_allele_freq(cb, stages = "males")
# h_cb <- h_cb[h_cb$allele == "H", ]
# h_ode <- sg_allele_freq(ode, stages = "males")
# h_ode <- h_ode[h_ode$allele == "H", ]
# ggplot(h_cb, aes(day, freq, group = replicate)) +
#   annotate("rect",
#     xmin = release_window[1], xmax = release_window[2],
#     ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.12
#   ) +
#   geom_line(colour = "#009E73", alpha = 0.45, linewidth = 0.4) +
#   geom_line(data = h_ode, aes(day, freq), colour = "black", linewidth = 0.8,
#             inherit.aes = FALSE) +
#   scale_x_continuous(breaks = seq(0, days, by = 90)) +
#   coord_cartesian(xlim = c(0, days), ylim = c(0, 1)) +
#   labs(
#     x = "days", y = "frequency of H among adult males",
#     title = "fsRIDL transgene among adult males",
#     subtitle = "stochastic replicates (green), ODE (black); shaded: release period"
#   ) +
#   theme_minimal()

# Notes
# * Unlike SIT and pgSIT males, fsRIDL males are fertile: their sons survive
#   and carry H, so the effect continues for a few generations after the
#   releases stop.
# * Daughters carrying H die as embryos, so their eggs neither hatch nor
#   become larvae; their mothers still search for wounds and lay (D-006).
# * With female_lethality < 1 some carrier daughters survive and pass H on.
