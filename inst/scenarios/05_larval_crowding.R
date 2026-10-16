# Scenario 05: larval crowding
#
# A wild screwworm population at a constant 25 degC, in a region of 1000 km^2,
# with and without larval crowding. With crowding, larvae die faster the more
# larvae share a wound: the larval hazard grows by
# larval_crowding_coef * (larvae / wounds). Both populations start from the
# same small founding population. While it is small, crowding hardly matters
# and the two grow alike; as larvae per wound rise, crowding slows growth and
# the population levels off lower.
#
# Run time: a few seconds.

library(SGDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Parameters: without and with larval crowding ----------------------------
p_none <- sg_params() # default: no larval crowding
p_crowd <- sg_params(
  larval_crowding_mode = "per_wound", # crowding by larvae per wound
  larval_crowding_coef = 0.02 # extra larval hazard per day for each larva per wound
)
cube <- sg_cube(MGDrivE::cubeMendelian()) # only the wild type (AA) is present

# 2. Weather and models: 25 degC for one year, region of 1000 km^2 -----------
weather <- sg_constant_weather(temperature = 25, days = 365)
model_none <- sg_model(p_none, cube, sg_forcing(weather, p_none), area = 1000)
model_crowd <- sg_model(p_crowd, cube, sg_forcing(weather, p_crowd), area = 1000)

# 3. The same founding population for both: 1% of the equilibrium without
# crowding (all stages in their equilibrium proportions) -----------------------
start <- as.numeric(sg_init(model_none, "equilibrium")) * 0.01

# 4. Simulate -----------------------------------------------------------------
sim_none <- simulate_ode(model_none, start)
sim_crowd <- simulate_ode(model_crowd, start)

# 5. Plot adults per km^2 -------------------------------------------------------
adults <- rbind(
  cbind(sg_totals(sim_none, by = character(0), stages = "adults", per = "km2"),
    scenario = "no larval crowding"
  ),
  cbind(sg_totals(sim_crowd, by = character(0), stages = "adults", per = "km2"),
    scenario = "larval crowding"
  )
)
ggplot(adults, aes(date, value, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("no larval crowding" = "#0072B2", "larval crowding" = "#D55E00")) +
  expand_limits(y = 0) +
  labs(
    x = NULL, y = "adults per km^2", colour = NULL,
    title = "Growth from a small population, with and without larval crowding (25 degC)"
  ) +
  theme_minimal()

# 6. Where the two populations settle -------------------------------------------
# sg_equilibria() gives the equilibria at 25 degC (stable = TRUE is where the
# population settles); crowding lowers it.
sg_equilibria(model_none)
sg_equilibria(model_crowd)
