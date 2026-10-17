# Scenario 04: myiasis cases limited by the number of hosts
#
# A wild screwworm population at McAllen, Texas (1971 weather), in a region of
# 1000 km^2 with 30 animals per km^2 (30,000 hosts), simulated with the
# host-wound model: animals get wounded, wounds heal or become infested, and
# a case is a new infestation of an animal. At most 30,000 animals can be
# infested at once, but the same animal can be infested again and again over
# the year. The script compares these cases with the egg batches holding
# viable eggs that the static wound model counts as cases.
#
# Run time: a few seconds.

library(SGDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Parameters: the host-wound model -----------------------------------------
# host_density and wounded_fraction have no defaults in the host-wound model.
# wound_search_rate sets how often a wound is found: at most
# wound_search_rate * (degree-days per day) egg batches per wound per day. The
# Gutierrez & Ponti value (0.001) was set for a relative wound density; with
# real numbers of animals it lets a wound be found only about once in 100
# days, so a larger value is used here. It still needs calibrating against
# case data.
p <- sg_params(
  wound_model = "hosts",
  host_density = 30, # animals per km^2
  wounded_fraction = 0.05, # share of animals with an attractive wound, without flies
  wound_healing_days = 7, # time an uninfested wound stays attractive
  infestation_days = 6, # eggs hatch within a day, larvae feed for 5-7 days
  wound_search_rate = 0.02
)
cube <- sg_cube(MGDrivE::cubeMendelian()) # only the wild type (AA) is present

# 2. Weather: McAllen, Texas, 1971 --------------------------------------------
weather <- sg_example_weather()
forcing <- sg_forcing(weather, p, start = "1971-01-01", end = "1971-12-31")

# 3. Model: a region of 1000 km^2 (30,000 animals) ---------------------------
model <- sg_model(p, cube, forcing, area = 1000)
check_model(model)

# 4. Start after a one-year burn-in over the weather of 1970 ----------------
init <- sg_init(model, "burn_in", weather = weather)

# 5. Simulate -----------------------------------------------------------------
ode <- simulate_ode(model, init) # deterministic
cb <- simulate_cb(model, init, n_rep = 10, seed = 1) # 10 stochastic replicates
ode
cb

# 6. Cases: animals newly infested, and egg batches with viable eggs ----------
# Every egg batch with a viable egg laid on an uninfested wound infests an
# animal (a case); batches laid on animals already infested add larvae but no
# case (superinfestations). The static wound model counts every such batch as
# a case.
cases <- sg_cases(ode, period = "month", per = "node")
cases
comparison <- rbind(
  data.frame(month = cases$period, value = cases$cases, count = "animals newly infested"),
  data.frame(
    month = cases$period, value = cases$cases + cases$superinfestations,
    count = "egg batches with viable eggs"
  )
)
comparison$count <- factor(comparison$count, levels = unique(comparison$count))
ggplot(comparison, aes(month, value, fill = count)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#0072B2", "#E69F00")) +
  labs(
    x = NULL, y = "per month in the region", fill = NULL,
    title = "Cases: animals newly infested and egg batches with viable eggs"
  ) +
  theme_minimal()

# 7. Infested animals through the year -----------------------------------------
hosts_ode <- sg_hosts(ode)
hosts_cb <- sg_hosts(cb)
ggplot(hosts_cb, aes(date, infested, group = replicate)) +
  geom_line(colour = "#56B4E9", alpha = 0.6) +
  geom_line(data = hosts_ode, aes(date, infested), colour = "black", inherit.aes = FALSE) +
  labs(
    x = NULL, y = "infested animals (of 30,000)",
    title = "Infested animals: ODE (black) and 10 stochastic replicates (blue)"
  ) +
  theme_minimal()
