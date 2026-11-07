# Initial states from the model's dynamics (D-025): burn-in and the
# closed-form equilibrium at constant temperature.

skip_if_not_installed("deSolve")

ns <- asNamespace("ScrewDrivE")
wx_all <- sg_example_weather()
cube3 <- sg_cube(sg_test_cube_mendelian())
const_weather <- function(temp, days = 40) {
  data.frame(
    date = as.Date("2001-01-01") + seq_len(days) - 1, tmax = temp, tmin = temp, tmean = temp
  )
}

test_that("the equilibrium is a fixed point of the ODE at that temperature", {
  settings <- list(
    sg_params(),
    sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.05),
    sg_params(male_density_half_mating = 3),
    sg_params(mating_mode = "instant", wound_competitors = "all_adults"),
    sg_params(preset = "gp2014"),
    sg_params(female_receptive_age_dd = 20, female_background_mortality = 0.05)
  )
  for (pp in settings) {
    for (temp in c(22, 29)) {
      m <- sg_model(pp, cube3, sg_forcing(const_weather(temp), pp), area = 10)
      x <- sg_init(m, "equilibrium", temperature = temp)
      expect_s3_class(x, "sg_state")
      expect_gt(sum(x), 0)
      r <- ns$.sg_rates_r(20, as.numeric(x), ns$.sg_reference_setup(m))
      expect_lt(max(abs(r$deriv)) / max(abs(r$counters)), 1e-10)
      expect_identical(sum(x[!grepl("AA", names(x))]), 0) # wild type only
      e <- attr(x, "equilibrium")
      expect_gt(e$R0, 1)
      expect_equal(e$eggs_per_day, r$counters[["eggs_laid[AA]"]], tolerance = 1e-10)
    }
  }
})

test_that("an ODE started at the equilibrium stays there", {
  pp <- sg_params()
  m <- sg_model(pp, cube3, sg_forcing(const_weather(27, 60), pp))
  x <- sg_init(m, "equilibrium", temperature = 27)
  sim <- simulate_ode(m, init = x, rtol = 1e-10, atol = 1e-12, record = "state")
  last <- sim$state[nrow(sim$state), ]
  expect_lt(max(abs(last - as.numeric(x))) / max(x), 1e-6)
})

test_that("no equilibrium where the population cannot persist", {
  m <- sg_model(sg_params(), cube3, sg_forcing(const_weather(16), sg_params()))
  expect_message(x <- sg_init(m, "equilibrium", temperature = 16), "No positive equilibrium")
  expect_identical(sum(x), 0)
  expect_lt(attr(x, "equilibrium")$R0, 1)
  expect_message(sg_init(m, "equilibrium", temperature = 10), "R0 = 0") # no development
  expect_true(any(grepl("equilibrium at 16 degC", capture.output(print(x)))))
})

test_that("the equilibrium needs a single wild type that breeds true", {
  m <- sg_model(sg_params(), cube3, sg_forcing(const_weather(25), sg_params()))
  # constant forcing: the temperature defaults to the forcing's
  expect_identical(sg_init(m, "equilibrium"), sg_init(m, "equilibrium", temperature = 25))
  expect_error(sg_init(m, "equilibrium", temperature = c(20, 25)), "single `temperature`")
  skip_if_not_installed("MGDrivE")
  wx25 <- sg_forcing(const_weather(25), sg_params())
  ms <- sg_model(sg_params(), sg_cube(MGDrivE::cubeXShredderY()), wx25)
  expect_error(sg_init(ms, "equilibrium", temperature = 25), "single wild-type genotype")
})

test_that("a burn-in forgets its starting state and uses the right weather", {
  p <- sg_params()
  m <- sg_model(p, cube3, sg_forcing(wx_all, p, "1971-01-01", "1971-12-31"))
  b <- sg_init(m, "burn_in")
  info <- attr(b, "burn_in")
  expect_identical(info$years, 1L)
  expect_identical(info$days, 365L)
  expect_match(info$source, "repeated")
  expect_lt(info$start_effect, 1e-6)
  big <- sg_init(m, "burn_in", values = c(E = 25, L = 25, P = 25, M = 12, F = 12))
  expect_lt(max(abs(big - b)) / max(b), 1e-6)
  # the repeated first year: two years give the same state as one
  expect_lt(max(abs(sg_init(m, "burn_in", years = 2) - b)) / max(b), 1e-6)
  # with weather: the calendar year 1970 before the start
  bw <- sg_init(m, "burn_in", weather = wx_all)
  expect_match(attr(bw, "burn_in")$source, "weather before")
  m70 <- sg_model(p, cube3, sg_forcing(wx_all, p, "1970-01-01", "1970-12-31"))
  direct <- simulate_ode(m70, init = sg_init(m70, "burn_in", weather = wx_all), record = "state")
  expect_lt(max(abs(direct$state[366, ] - bw)) / max(bw), 1e-6)
  expect_true(any(grepl("1-year burn-in", capture.output(print(bw)))))
  # Hardy-Weinberg start: genotype frequencies carry through a neutral burn-in
  hw <- sg_init(m, "burn_in", alleles = c(A = 0.5, a = 0.5))
  expect_gt(sum(hw[grepl("aa", names(hw))]), 0)
})

test_that("burn-in inputs are checked", {
  p <- sg_params()
  short <- sg_model(p, cube3, sg_forcing(wx_all, p, "1971-01-01", "1971-03-31"))
  expect_error(sg_init(short, "burn_in"), "only 90 days")
  expect_error(sg_init(short, "burn_in", years = 0, weather = wx_all), "whole number >= 1")
  expect_error(sg_init(short, "burn_in", years = 1.5, weather = wx_all), "whole number >= 1")
  expect_s3_class(sg_init(short, "burn_in", weather = wx_all), "sg_state")
})

# ---- population size from adult_females (D-028) ------------------------------

adult_females_of <- function(m, x) sum(x[c(sg_state_index(m, "U"), sg_state_index(m, "F"))])

test_that("sg_constant_weather() gives one temperature on consecutive days", {
  wx <- sg_constant_weather(20, days = 10, start = "2005-03-01")
  expect_named(wx, c("date", "tmean"))
  expect_identical(wx$date, as.Date("2005-03-01") + 0:9)
  expect_true(all(wx$tmean == 20))
  f <- sg_forcing(wx, sg_params())
  expect_identical(ns$.sg_forcing_constant_temp(f, sg_params()), 20)
  expect_null(ns$.sg_forcing_constant_temp(sg_forcing(wx_all[1:40, ], sg_params()), sg_params()))
  expect_error(sg_constant_weather(c(20, 21)), "single finite number")
  expect_error(sg_constant_weather(20, days = 0), "whole number >= 1")
  expect_error(sg_constant_weather(20, start = "not a date"), "single date")
})

test_that("adult_females sets the area so the equilibrium holds that many females", {
  p <- sg_params()
  f <- sg_forcing(sg_constant_weather(20, days = 60), p)
  m <- sg_model(p, cube3, f, adult_females = 1000)
  e <- sg_init(m, "equilibrium") # temperature taken from the forcing
  expect_identical(attr(e, "equilibrium")$temperature, 20)
  expect_equal(adult_females_of(m, e), 1000, tolerance = 1e-10)
  # densities stay as given; only the area (and so the wound count) changes
  unit <- sg_model(p, cube3, f)
  density <- adult_females_of(unit, sg_init(unit, "equilibrium"))
  expect_equal(m$area, 1000 / density, tolerance = 1e-12)
  expect_identical(m$params$wound_density, p$wound_density)
  expect_equal(m$counts$H, p$wound_density * m$area, tolerance = 1e-12)
  expect_identical(m$meta$adult_females$temperature, 20)
  expect_true(any(grepl("area chosen for 1000 adult females", capture.output(print(m)))))
  # the ODE stays at the equilibrium
  fem <- sg_females(simulate_ode(m, e, record = "state"))$females
  expect_lt(max(abs(fem - 1000)), 1e-6)
})

test_that("adult_females is exact with larval crowding and an Allee effect", {
  p <- sg_params(
    larval_crowding_mode = "per_wound", larval_crowding_coef = 0.05,
    male_density_half_mating = 0.5
  )
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(22), p), adult_females = 250)
  expect_equal(adult_females_of(m, sg_init(m, "equilibrium")), 250, tolerance = 1e-8)
})

test_that("adult_females and the equilibrium default are checked", {
  p <- sg_params()
  f20 <- sg_forcing(sg_constant_weather(20, days = 30), p)
  varying <- sg_forcing(wx_all[1:40, ], p)
  expect_error(sg_model(p, cube3, f20, area = 2, adult_females = 10), "not both")
  expect_error(sg_model(p, cube3, f20, adult_females = -1), "single positive number")
  expect_error(sg_model(p, cube3, varying, adult_females = 10), "one constant temperature")
  expect_error(
    sg_model(p, cube3, sg_forcing(sg_constant_weather(15, days = 30), p), adult_females = 10),
    "no wild-type population persists at 15 degC"
  )
  expect_error(sg_init(sg_model(p, cube3, varying), "equilibrium"), "needs a `temperature`")
  # drought stress the equilibrium would miss makes the forcing not constant
  wet <- transform(sg_constant_weather(20, days = 30), rain = 0)
  pd <- sg_params(pupa_drought_mortality = 0.1)
  expect_error(sg_model(pd, cube3, sg_forcing(wet, pd), adult_females = 10), "one constant")
})

test_that("calibrating the wound density drops the adult_females note", {
  p <- sg_params()
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(20), p), adult_females = 1000)
  expect_null(ns$.sg_model_with_h(m, 10)$meta$adult_females)
})

# ---- Allee threshold (D-029) ---------------------------------------------------

test_that("the Allee threshold is the lower, unstable equilibrium", {
  p <- sg_params(male_density_half_mating = 20)
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 730), p), adult_females = 1000)
  thr <- sg_init(m, "threshold") # temperature from the forcing
  info <- attr(thr, "threshold")
  expect_identical(info$temperature, 20)
  expect_equal(info$adult_females, adult_females_of(m, thr), tolerance = 1e-12)
  expect_gt(info$adult_females, 0)
  expect_lt(info$adult_females, 1000)
  expect_true(any(grepl("Allee threshold", capture.output(print(thr)))))
  # an equilibrium: the ODE stays there over a few weeks
  short <- sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 20), p), area = m$area)
  fem <- sg_females(simulate_ode(short, init = sg_init(short, "threshold")))$females
  expect_lt(max(abs(fem / info$adult_females - 1)), 1e-6)
  # unstable: just below it the population declines, just above it grows
  end_from <- function(k) tail(sg_females(simulate_ode(m, init = thr * k))$females, 1)
  expect_lt(end_from(0.9), 1e-3 * info$adult_females)
  expect_equal(end_from(1.1), 1000, tolerance = 1e-4)
})

test_that("the Allee threshold rises in proportion to h_M", {
  thr_of <- function(h) {
    p <- sg_params(male_density_half_mating = h)
    m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 30), p), area = 100)
    attr(sg_init(m, "threshold"), "threshold")$adult_females
  }
  t1 <- thr_of(1)
  expect_equal(thr_of(5) / t1, 5, tolerance = 1e-3)
  expect_equal(thr_of(20) / t1, 20, tolerance = 1e-3)
})

test_that("there is no Allee threshold without an Allee effect", {
  f <- sg_forcing(sg_constant_weather(20, days = 30), sg_params())
  m <- sg_model(sg_params(), cube3, f)
  expect_message(x <- sg_init(m, "threshold"), "No Allee threshold")
  expect_identical(sum(x), 0)
  expect_identical(attr(x, "threshold")$adult_females, 0)
  pi <- sg_params(male_density_half_mating = 5, mating_mode = "instant")
  mi <- sg_model(pi, cube3, sg_forcing(sg_constant_weather(20, days = 30), pi))
  expect_message(sg_init(mi, "threshold"), "No Allee threshold")
  # the stable equilibrium is unchanged by the lower-root option
  p <- sg_params(male_density_half_mating = 5)
  m5 <- sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 30), p))
  expect_equal(
    attr(sg_init(m5, "equilibrium"), "equilibrium")$adult_females,
    attr(sg_init(m, "equilibrium"), "equilibrium")$adult_females,
    tolerance = 1e-3
  )
})

test_that("host-model equilibria below the oviposition threshold are empty, not an error", {
  # found in phase 1 WP7: no egg laying made the fertile batches per egg 0 / 0
  p <- sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  wx <- sg_forcing(sg_constant_weather(13, days = 30), p)
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx)
  expect_message(x <- sg_init(m, "equilibrium"), "No positive equilibrium at 13 degC")
  fly <- unlist(lapply(c("E", "L", "P", "M", "U", "F"), function(b) sg_state_index(m, b)))
  expect_identical(sum(x[fly]), 0)
  expect_equal(sum(x[-fly]), m$hosts$n)
  expect_identical(nrow(sg_equilibria(m)), 1L)
})
