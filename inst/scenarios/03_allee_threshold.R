# Scenario 03: the Allee threshold of a wild-type screwworm population
#
# With an Allee effect, unmated females find mates more slowly when males are
# scarce. The half-saturation male density male_density_half_mating (h_M,
# males per km^2) sets how strong this is: at h_M males per km^2 females mate
# at half the full rate. A very sparse population then cannot replace itself:
# below a threshold size it declines to extinction, even at 20 degC where R0
# is large (about 71), and above it grows to the equilibrium.
#
# The threshold is the lower, unstable equilibrium, which sg_init(model,
# "threshold") computes: a sparse population in which most females have not
# yet found a mate. For several values of h_M, the script starts populations
# just below and just above it (0.5, 0.9, 1.1 and 2 times the threshold) in a
# region that holds 1000 adult females at equilibrium. Without an Allee effect
# (h_M = 0) there is no threshold.
#
# Run time: a few seconds.

library(SGDrivE)
library(MGDrivE)
library(ggplot2) # plots

# 1. Settings -----------------------------------------------------------------
h_values <- c(1, 5, 20, 50) # half-saturation male density (males per km^2)
multiples <- c(0.5, 0.9, 1.1, 2) # starting population, times the threshold
cube <- sg_cube(MGDrivE::cubeMendelian()) # only the wild type (AA) is present
weather <- sg_constant_weather(temperature = 20, days = 730) # two years

# 2. Simulate every h_M and starting population -------------------------------
# For each h_M: a region holding 1000 adult females at equilibrium, its Allee
# threshold, and one deterministic run from each multiple of the threshold.
results <- NULL
thresholds <- NULL
for (h in h_values) {
  p <- sg_params(male_density_half_mating = h)
  model <- sg_model(p, cube, sg_forcing(weather, p), adult_females = 1000)

  threshold <- sg_init(model, "threshold")
  thresholds <- rbind(thresholds, data.frame(
    h_M = h, threshold = attr(threshold, "threshold")$adult_females
  ))

  for (k in multiples) {
    sim <- simulate_ode(model, init = threshold * k)
    females <- sg_females(sim)
    results <- rbind(results, data.frame(
      h_M = h, start = k, date = females$date, females = females$females
    ))
  }
}
thresholds # adult females at the Allee threshold, for each h_M

# 3. Plot ---------------------------------------------------------------------
# One panel per h_M; the dashed line is the threshold. Populations starting
# below it decline, those above it grow to the equilibrium of 1000. The
# threshold rises in proportion to h_M.
panel <- function(h) paste("h_M =", h, "males per km^2")
results$panel <- factor(panel(results$h_M), levels = panel(h_values))
thresholds$panel <- factor(panel(thresholds$h_M), levels = panel(h_values))

ggplot(results, aes(date, females, colour = factor(start))) +
  geom_line(linewidth = 0.8) +
  geom_hline(data = thresholds, aes(yintercept = threshold), linetype = "dashed") +
  facet_wrap(~panel) +
  scale_y_log10(breaks = 10^(-3:3), labels = as.character(10^(-3:3))) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  coord_cartesian(ylim = c(0.001, 2000)) +
  scale_colour_viridis_d(end = 0.9) +
  labs(
    x = NULL, y = "adult females in the region (log scale)",
    colour = "start\n(times the\nthreshold)",
    title = "Allee threshold: small populations die out when mates are hard to find",
    subtitle = "20 degC; dashed line: threshold (unstable equilibrium)"
  ) +
  theme_minimal()
