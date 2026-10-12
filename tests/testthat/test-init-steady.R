# Initial states from the model's dynamics (D-025): burn-in and the
# closed-form equilibrium at constant temperature.

skip_if_not_installed("deSolve")

ns <- asNamespace("SGDrivE")
wx_all <- sg_example_weather()
cube3 <- sg_cube(sg_test_cube_mendelian())
const_weather <- function(temp, days = 40) {
  data.frame(
    date = as.Date("2001-01-01") + seq_len(days) - 1, tmax = temp, tmin = temp, tmean = temp
  )
}

test_that("the equilibrium is a fixed point of the ODE at that temperature", {
  settings <- list(
    sg_params(), sg_params(larval_dd = "per_wound", gamma_L = 0.05), sg_params(h_M = 3),
    sg_params(mating = "instant", search_denominator = "all_adults"),
    sg_params(preset = "gp2014"), sg_params(X_rec = 20, mu_b_F = 0.05)
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
  sim <- simulate_ode(m, init = x, rtol = 1e-10, atol = 1e-12)
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
  expect_error(sg_init(m, "equilibrium"), "single `temperature`")
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
  direct <- simulate_ode(m70, init = sg_init(m70, "burn_in", weather = wx_all))
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
