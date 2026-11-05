# Scenario 09: the Gutierrez & Ponti (2014) model at McAllen, Texas, 1960-1983
#
# Reproduces the configuration of Gutierrez & Ponti (2014) and their McAllen,
# Texas run with daily weather from 1960 to 1983:
#
# * sg_params(preset = "gp2014"): 25 sub-stages in every stage, first egg
#   laying at 56.784 dd, the Pascal code's development-rate floor, all adults
#   competing for wounds and mating at emergence;
# * one wild genotype (no genetic control);
# * the package's McAllen weather (NOAA, gap-filled) and the default start of
#   their code (sg_init(model)).
#
# With these settings ScrewDrivE matches the R translation of their Pascal code
# to a relative error of 1e-6 in every daily total over the 24 years
# (tests/testthat/test-regression-gp2014.R). The script checks the published
# result that can be compared numerically, the annual cold-stress index
# mu_cold (Gutierrez & Ponti 2014; Gutierrez et al. 2019), and plots the
# simulated seasonal dynamics. As in the paper, the wound density is a
# relative value, so the population is in relative units, not flies per km^2.
#
# Run time: a few seconds.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Parameters, cube and weather ---------------------------------------------
p <- sg_params(preset = "gp2014")
# For field mortality, raise female_background_mortality to 0.10 and
# male_background_mortality to 0.12.
cube <- sg_cube(list( # one wild genotype that breeds true
  genotypesID = "W", wildType = "W",
  ih = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  tau = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  eta = matrix(1, 1, 1, dimnames = list("W", "W")),
  phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1), omega = c(W = 1)
))
weather <- sg_example_weather() # McAllen, Texas, 1960-1983
forcing <- sg_forcing(weather, p)

# 2. Model, start and simulation ----------------------------------------------
model <- sg_model(p, cube, forcing)
model
sim <- simulate_ode(model, sg_init(model)) # starting state of the GP2014 code

# 3. Published result: the cold-stress index mu_cold --------------------------
# mu_cold sums the daily adult mortality p_A(T) over the days below 27.2 degC
# in each calendar year; up to about 10 the screwworm can persist through the
# winter. Published values for McAllen (Gutierrez & Ponti 2014; Gutierrez et
# al. 2019, text and Fig. 5B; for 1974 and 1975 the figure values are used,
# as the text gives them in swapped order).
published <- data.frame(
  year = c(1961, 1971, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981),
  mu_cold_published = c(8.31, 6.67, 8.12, 7.56, 10.3, 9.1, 10.3, 9.80, 9.2, 7.49)
)
cold <- sg_cold_index(forcing)
comparison <- merge(published, cold[, c("year", "mu_cold")])
comparison$difference <- comparison$mu_cold - comparison$mu_cold_published
comparison # differences of at most about 0.1

ggplot(cold, aes(year, mu_cold)) +
  geom_col(fill = "#56B4E9") +
  geom_point(data = published, aes(year, mu_cold_published), colour = "#D55E00", size = 2.5) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_x_continuous(breaks = seq(1960, 1983, by = 2)) +
  labs(
    x = NULL, y = expression(mu[cold]),
    title = "Cold-stress index at McAllen, Texas",
    subtitle = paste(
      "bars: ScrewDrivE; points: published values (Gutierrez & Ponti 2014;",
      "Gutierrez et al. 2019); dashed: about the limit for overwintering"
    )
  ) +
  theme_minimal()

# 4. Simulated dynamics -------------------------------------------------------
# Adult females and daily egg batches (oviposition events, the case proxy of
# the paper): low in winter, high from spring to autumn, with deeper winter
# declines after cold winters (high mu_cold).
females <- sg_females(sim)
batches <- data.frame(date = sim$date[-1], egg_batches = sim$counters[, "egg_batches"])
daily <- merge(females[, c("date", "females")], batches)
daily <- rbind(
  data.frame(date = daily$date, value = daily$females, what = "adult females"),
  data.frame(date = daily$date, value = daily$egg_batches, what = "egg batches per day")
)
ggplot(daily, aes(date, value)) +
  geom_line(colour = "#0072B2", linewidth = 0.4) +
  facet_wrap(~what, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    x = NULL, y = "relative units (relative wound density, as in the paper)",
    title = "Gutierrez & Ponti (2014) model at McAllen, Texas, 1960-1983",
    subtitle = 'sg_params(preset = "gp2014"), daily NOAA weather'
  ) +
  theme_minimal()

# 5. Cold winters and the population low --------------------------------------
# For each winter (1 July to 30 June, labelled by the year in which it ends),
# the lowest number of adult females against mu_cold of the same period: the
# colder the winter, the deeper the population falls.
winter <- sg_cold_index(forcing, period = "jul_jun")
fem <- females[females$day > 0, ]
month <- as.integer(format(fem$date, "%m"))
fem$year <- as.integer(format(fem$date, "%Y")) + (month >= 7)
lows <- aggregate(females ~ year, data = fem, FUN = min)
lows <- merge(lows, winter[winter$complete, c("year", "mu_cold")])
cor(lows$mu_cold, log(lows$females)) # about -0.95

ggplot(lows, aes(mu_cold, females)) +
  geom_point(colour = "#0072B2", size = 2) +
  geom_text(aes(label = year), nudge_y = 0.08, size = 3, check_overlap = TRUE) +
  scale_y_log10() +
  labs(
    x = expression(mu[cold] ~ "of the winter (July to June)"),
    y = "lowest adult females in the winter (log scale)",
    title = "Colder winters push the population lower"
  ) +
  theme_minimal()
