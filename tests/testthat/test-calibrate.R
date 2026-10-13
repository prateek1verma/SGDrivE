# calibrate_H() (D-026): closed form from the exact proportionality to H, and
# root search with an Allee effect or releases.

skip_if_not_installed("deSolve")

wx_all <- sg_example_weather()
cube1 <- sg_cube(sg_test_cube_wild()) # calibration concerns the wild population
year <- function(p) sg_forcing(wx_all, p, "1971-01-01", "1971-12-31")
mean_adults <- function(m) {
  sim <- simulate_ode(m, init = sg_init(m, "burn_in"))
  mean(sg_totals(sim, by = character(0), stages = "adults", per = "km2")$value[-1])
}

test_that("the closed form hits the target in two runs, whatever the node area", {
  p <- sg_params()
  m <- sg_model(p, cube1, year(p))
  cal <- calibrate_H(m, target = 100)
  expect_s3_class(cal, "sg_calibration")
  expect_identical(cal$method, "proportional")
  expect_identical(cal$evaluations, 2L)
  expect_equal(cal$achieved, 100, tolerance = 1e-8)
  expect_equal(cal$params$H, cal$H)
  expect_equal(cal$model$counts$H, cal$H) # area 1 km^2
  expect_false(identical(cal$model$meta$hash, m$meta$hash))
  expect_equal(mean_adults(cal$model), 100, tolerance = 1e-6) # the calibrated model
  big <- calibrate_H(sg_model(p, cube1, year(p), area = 1000), target = 100)
  expect_equal(big$H, cal$H, tolerance = 1e-10)
  expect_true(any(grepl("H = .* wounds per km\\^2 \\(proportional", capture.output(print(cal)))))
})

test_that("other summaries and reference periods", {
  p <- sg_params()
  m <- sg_model(p, cube1, year(p))
  summer <- c("1971-06-01", "1971-08-31")
  for (st in c("peak", "min", "quantile")) {
    cal <- calibrate_H(m, target = 50, stat = st, period = summer, prob = 0.9)
    expect_equal(cal$achieved, 50, tolerance = 1e-8)
  }
  juv <- calibrate_H(m, target = 200, stages = "juveniles")
  expect_equal(juv$achieved, 200, tolerance = 1e-8)
})

test_that("root search with an Allee effect or releases agrees with the closed form", {
  skip_on_cran()
  p <- sg_params()
  m <- sg_model(p, cube1, year(p))
  closed <- calibrate_H(m, target = 100)
  # a negligible Allee effect forces the root search on the same problem
  pa <- sg_params(h_M = 1e-9)
  ra <- calibrate_H(sg_model(pa, cube1, year(pa)), target = 100)
  expect_identical(ra$method, "root search")
  expect_equal(ra$achieved, 100, tolerance = 1e-8)
  expect_equal(ra$H, closed$H, tolerance = 1e-6)
  rel <- sg_releases("W", "adult", "male", times = 100, number = 50)
  expect_warning(rr <- calibrate_H(sg_model(p, cube1, year(p), rel), target = 100), "releases")
  expect_identical(rr$method, "root search")
  expect_equal(rr$achieved, 100, tolerance = 1e-8)
})

test_that("calibrate_H() checks its input and impossible targets", {
  p <- sg_params()
  m <- sg_model(p, cube1, year(p))
  expect_error(calibrate_H(list(), 10), "sg_model object")
  expect_error(calibrate_H(m, -1), "positive density")
  expect_error(calibrate_H(m, 10, stat = "median"), "must be one of")
  expect_error(calibrate_H(m, 10, stat = "quantile", prob = 2), "probability")
  expect_error(calibrate_H(m, 10, stages = "adult"), "Unknown stage")
  expect_error(calibrate_H(m, 10, period = c("1972-01-01", "1972-02-01")), "outside")
  expect_error(calibrate_H(m, 10, period = "1971-01-01"), "two dates")
  expect_error(calibrate_H(m, 10, years = 0), "whole number >= 1")
  m0 <- sg_model(sg_params(H = 0), cube1, year(p))
  expect_error(calibrate_H(m0, 10), "no wounds")
  cold <- data.frame(date = as.Date("2001-01-01") + 0:399, tmax = 12, tmin = 12, tmean = 12)
  mc <- sg_model(p, cube1, sg_forcing(cold, p))
  expect_error(calibrate_H(mc, 10), "population is zero")
})
