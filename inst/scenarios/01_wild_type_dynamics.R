# Scenario 01: wild-type screwworm population dynamics at McAllen, Texas
#
# A wild screwworm population (no genetic control) driven by the daily
# weather of McAllen, Texas, in 1971-1972. The script builds the model,
# checks it, starts the population from a one-year burn-in, runs the
# deterministic (ODE) and stochastic (chain-binomial) solvers, and plots the
# seasonal dynamics of life stages, females and cases, the agreement of the
# two solvers, and the thermal window of population growth (R0).
#
# Run time: well under a minute.

library(SGDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Parameters ---------------------------------------------------------------
# The defaults follow Gutierrez & Ponti (2014); sg_param_table() lists every
# parameter with its units, source and allowed range.
p <- sg_params()

# 2. Inheritance ----------------------------------------------------------------
# A one-locus Mendelian cube from MGDrivE. Only the wild type (AA) is present
# here, so inheritance plays no role yet; later scenarios release other
# genotypes into the same model.
cube <- sg_cube(MGDrivE::cubeMendelian())

# 3. Weather ---------------------------------------------------------------------
# Daily temperatures of McAllen, Texas, 1960-1983, are bundled with the
# package. sg_forcing() turns them into the daily development, mortality and
# oviposition rates of the model.
weather <- sg_example_weather()
forcing <- sg_forcing(weather, p, start = "1971-01-01", end = "1972-12-31")

# 4. Model -----------------------------------------------------------------------
# Counts are whole individuals per node, so the node area sets the population
# size. With 1000 km^2 the population stays large through the winter low,
# which the stochastic solver needs (check_model() warns otherwise).
model <- sg_model(p, cube, forcing, area = 1000)
check_model(model)

# 5. Initial population ------------------------------------------------------
# The ODE forgets its starting state within months, so the population is
# started from a one-year burn-in over the weather of 1970.
init <- sg_init(model, "burn_in", weather = weather)

# 6. Deterministic simulation ----------------------------------------------
ode <- simulate_ode(model, init)
ode

plot(ode) # all life stages, per km^2
plot(ode, "females") # unmated, mated and egg-laying females
plot(ode, "cases", period = "month") # wounds newly infested with viable eggs

# 7. Stochastic simulation -------------------------------------------------
# Ten replicates of the chain-binomial solver; each replicate has its own
# random-number stream, so the result is reproducible with the same seed.
cb <- simulate_cb(model, init, n_rep = 10, seed = 1)
cb

plot(cb, stages = "adults") # each replicate is one line

# 8. Deterministic and stochastic results together -------------------------
# The mean of the stochastic replicates follows the ODE; at 1000 km^2 the
# replicates stay close to it.
adults_ode <- sg_totals(ode, by = character(0), stages = "adults", per = "km2")
adults_cb <- sg_totals(cb, by = character(0), stages = "adults", per = "km2")
ggplot(adults_cb, aes(date, value, group = replicate)) +
  geom_line(colour = "#56B4E9", alpha = 0.6) +
  geom_line(
    data = adults_ode, aes(date, value), colour = "black", linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  labs(
    x = NULL, y = "adults per km^2",
    title = "Adults: ODE (black) and 10 stochastic replicates (blue)"
  ) +
  theme_minimal()

# 9. Weather and population ---------------------------------------------------
# Mean daily temperature against adult density: the population collapses in
# the cold months and recovers in spring.
temperature <- data.frame(date = forcing$dates, tmean = forcing$tmean)
both <- merge(adults_ode, temperature, by = "date")
scale <- max(both$value) / max(both$tmean)
ggplot(both, aes(date)) +
  geom_line(aes(y = tmean * scale), colour = "#E69F00") +
  geom_line(aes(y = value), colour = "black") +
  scale_y_continuous(
    name = "adults per km^2",
    sec.axis = sec_axis(~ . / scale, name = "mean temperature (degC)")
  ) +
  labs(x = NULL, title = "Adults (black) and daily mean temperature (orange)") +
  theme_minimal()



# 11. Tables ---------------------------------------------------------------------
# The same results as tables, for further analysis.
head(sg_totals(ode, by = "stage", per = "km2"))
sg_cases(ode, period = "month")
