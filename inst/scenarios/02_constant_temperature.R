# Scenario 02: wild-type screwworm population at a constant 20 degC
#
# A wild screwworm population (no genetic control) at a constant temperature
# of 20 degC, in a region just large enough to hold 1000 adult females. The
# population starts at its equilibrium, so the deterministic run stays flat
# and the stochastic runs fluctuate around it.
#
# Run time: a few seconds.

library(ScrewDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Parameters and inheritance -----------------------------------------------
# The defaults follow Gutierrez & Ponti (2014). Any parameter can be changed by
# name in sg_params(); the values below are the defaults, so edit them to
# explore. sg_param_table() lists every parameter with its units and range.
p <- sg_params(
  eggs_per_day_max = 90, # maximum eggs per egg-laying female per day
  wound_finding_max = 1, # wound-finding efficiency at low density (0 to 1)
  female_background_mortality = 0, # extra mortality of adult females (per day)
  male_background_mortality = 0 # extra mortality of adult males (per day)
)
cube <- sg_cube(MGDrivE::cubeMendelian()) # only the wild type (AA) is present

# 2. Weather: 20 degC every day for one year ----------------------------------
weather <- sg_constant_weather(temperature = 20, days = 365)
forcing <- sg_forcing(weather, p)

# 3. Model: a region holding 1000 adult females -------------------------------
# The parameters set the screwworm density per km^2. With adult_females = 1000
# the model chooses the area of the region (and so the number of hosts in it)
# that holds 1000 adult females at equilibrium; print(model) shows the area.
# Changing the parameters above changes the density, so the area changes and
# the number of adult females stays 1000.
model <- sg_model(p, cube, forcing, adult_females = 1000)
model

# 4. Start at the equilibrium -------------------------------------------------
init <- sg_init(model, "equilibrium")
init

# 5. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic
cb <- simulate_cb(model, init, n_rep = 10, seed = 1) # 10 stochastic replicates

# 6. Plots --------------------------------------------------------------------
plot(ode, per = "node") # all life stages, numbers in the region

females_ode <- sg_females(ode)
females_cb <- sg_females(cb)
ggplot(females_cb, aes(date, females, group = replicate)) +
  geom_line(colour = "#56B4E9", alpha = 0.6) +
  geom_line(data = females_ode, aes(date, females), colour = "black", inherit.aes = FALSE) +
  labs(
    x = NULL, y = "adult females in the region",
    title = "Adult females: ODE (black) and 10 stochastic replicates (blue)"
  ) +
  theme_minimal()
