# Scenario 06: rainfall and pupal survival
#
# A wild screwworm population at a constant 25 degC, in a region of 1000 km^2,
# in a wet climate (600 mm of rain per year) and a dry one (300 mm per year).
# Pupae in dry soil die faster: their hazard grows by
# pupa_drought_mortality * (1 - moisture), where the moisture index is
# min(1, annual rainfall / 450 mm). With 600 mm the index is 1 (no stress),
# with 300 mm it is 0.67. Both populations start at the equilibrium of the
# wet climate; in the dry climate the population falls to a lower level.
#
# Run time: a few seconds.

library(ScrewDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Parameters ---------------------------------------------------------------
# pupa_drought_mortality is 0 by default (no effect of rainfall); 0.1 extra
# deaths per pupa per day at full drought is an assumed value for this
# demonstration.
p <- sg_params(
  pupa_drought_mortality = 0.1, # extra pupal hazard per day at full drought
  drought_free_rainfall = 450 # annual rainfall (mm) above which there is no stress
)
cube <- sg_cube(MGDrivE::cubeMendelian()) # only the wild type (AA) is present

# 2. Weather: 25 degC every day for one year, wet or dry -----------------------
wet <- sg_constant_weather(temperature = 25, days = 365)
wet$rain <- 600 / 365 # mm per day: 600 mm per year
dry <- wet
dry$rain <- 300 / 365 # 300 mm per year

forcing_wet <- sg_forcing(wet, p)
forcing_dry <- sg_forcing(dry, p)
c(wet = forcing_wet$moisture[1], dry = forcing_dry$moisture[1]) # moisture index

# 3. Models: region of 1000 km^2 -------------------------------------------------
model_wet <- sg_model(p, cube, forcing_wet, area = 1000)
model_dry <- sg_model(p, cube, forcing_dry, area = 1000)

# 4. Both start at the equilibrium of the wet climate ----------------------------
start <- as.numeric(sg_init(model_wet, "equilibrium"))

# 5. Simulate -----------------------------------------------------------------
sim_wet <- simulate_ode(model_wet, start)
sim_dry <- simulate_ode(model_dry, start)

# 6. Plot adults per km^2 -------------------------------------------------------
adults <- rbind(
  cbind(sg_totals(sim_wet, by = character(0), stages = "adults", per = "km2"),
    scenario = "wet (600 mm per year)"
  ),
  cbind(sg_totals(sim_dry, by = character(0), stages = "adults", per = "km2"),
    scenario = "dry (300 mm per year)"
  )
)
ggplot(adults, aes(date, value, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c("wet (600 mm per year)" = "#0072B2", "dry (300 mm per year)" = "#E69F00")
  ) +
  expand_limits(y = 0) +
  labs(
    x = NULL, y = "adults per km^2", colour = NULL,
    title = "Adults in a wet and a dry climate (25 degC)"
  ) +
  theme_minimal()
