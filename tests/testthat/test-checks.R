# Model checks (D-026): mating in the largest-rate diagnostic (O-007), R0 by
# temperature, check_model(), flow balance and runtime guards.

skip_if_not_installed("deSolve")

ns <- asNamespace("ScrewDrivE")
wx_all <- sg_example_weather()
cube3 <- sg_cube(sg_test_cube_mendelian())
year <- function(p) sg_forcing(wx_all, p, "1971-01-01", "1971-12-31")
month <- function(p) sg_forcing(wx_all, p, "1971-05-01", "1971-05-31")
const_weather <- function(temp, days = 40) {
  data.frame(
    date = as.Date("2001-01-01") + seq_len(days) - 1, tmax = temp, tmin = temp, tmean = temp
  )
}

test_that("the largest-rate diagnostic includes mating of unmated females", {
  pp <- sg_params(mating_rate = 100)
  m <- sg_model(pp, cube3, month(pp))
  x <- as.numeric(sg_init(m, values = c(E = 1, L = 1, P = 1, M = 1, U = 1, F = 1)))
  s <- ns$.sg_reference_setup(m)
  v <- ns$.sg_model_view(m)
  r <- ns$.sg_max_rate(5, x, s)
  expect_gt(r, 100)
  expect_equal(ns$sg_max_rate_cpp(v, 5L, x), r, tolerance = 1e-14)
  # no males: nobody can mate, so the mating rate does not count
  x0 <- x
  x0[sg_state_index(m, "M")] <- 0
  expect_lt(ns$.sg_max_rate(5, x0, s), 100)
  expect_equal(ns$sg_max_rate_cpp(v, 5L, x0), ns$.sg_max_rate(5, x0, s), tolerance = 1e-14)
  sim <- simulate_ode(m, init = x)
  expect_true(all(sim$diagnostics$max_rate > 100))
})

test_that("R0 by temperature: thermal window and agreement with the equilibrium", {
  p <- sg_params()
  m <- sg_model(p, cube3, month(p))
  temps <- c(10, 14, 20, 27, 35, 45)
  r0 <- sg_r0(m, temps)
  expect_length(r0, 6)
  expect_identical(r0[1], 0) # below the developmental threshold
  expect_lt(r0[2], 1)
  expect_gt(r0[4], max(r0[c(3, 5)])) # highest near the optimum
  expect_lt(r0[6], r0[5])
  mc <- sg_model(p, cube3, sg_forcing(const_weather(27), p))
  eq <- sg_init(mc, "equilibrium", temperature = 27)
  expect_equal(sg_r0(mc, 27), attr(eq, "equilibrium")$R0, tolerance = 1e-12)
  expect_error(sg_r0(m, NA), "finite numbers")
  expect_error(sg_r0(list(), 20), "sg_model object")
  skip_if_not_installed("MGDrivE")
  xs <- sg_model(p, sg_cube(MGDrivE::cubeXShredderY()), month(p))
  expect_error(sg_r0(xs, 25), "single wild-type genotype")
})

test_that("check_model() reports size, rates, R0 season and the population scale", {
  p <- sg_params()
  ck <- check_model(sg_model(p, cube3, year(p)))
  expect_s3_class(ck, "sg_check")
  expect_identical(ck$size$genotypes, 3L)
  expect_true(ck$r0$available)
  expect_gt(ck$r0$days_above_1, 200)
  expect_lte(ck$r0$longest_run, ck$r0$days_above_1)
  expect_equal(ck$rates$mating_rate, p$mating_rate)
  expect_lt(ck$scale$min_wild, 10) # 1 km^2: a few individuals in winter
  expect_true(any(grepl("stochastic runs may go extinct", ck$warnings)))
  out <- capture.output(print(ck))
  expect_match(out[1], "<sg_check>")
  expect_true(any(grepl("R0 \\(eggs per egg at low density\\)", out)))
  big <- check_model(sg_model(p, cube3, year(p), area = 1000))
  expect_false(any(grepl("extinct", big$warnings)))
  expect_equal(big$scale$min_wild, 1000 * ck$scale$min_wild, tolerance = 1e-6) # scales with area
})

test_that("check_model() warns about stiffness, inert genotypes and cold climates", {
  raw <- sg_test_cube_sit(residual = 0)
  pp <- sg_params(mating_rate = 2000)
  m <- sg_model(pp, sg_cube(raw, lifespan_fraction = c(W = 1, S = 0.25)), year(pp),
    sg_releases("S", "adult", "male", times = 10, number = 100),
    area = 1000
  )
  ck <- check_model(m, scale = FALSE)
  expect_true(any(grepl("may be stiff", ck$warnings)))
  expect_true(any(grepl("S: females die before they can lay eggs", ck$warnings)))
  expect_identical(ck$size$sterile_releases, "S")
  expect_true(any(grepl("sterile genotype\\(s\\): S", capture.output(print(ck)))))
  cold <- sg_model(sg_params(), cube3, sg_forcing(const_weather(12), sg_params()))
  ckc <- check_model(cold)
  expect_true(any(grepl("R0 never exceeds 1", ckc$warnings)))
  expect_identical(ckc$scale$note, "needs at least 365 days of forcing")
  skip_if_not_installed("MGDrivE")
  xs <- check_model(sg_model(sg_params(), sg_cube(MGDrivE::cubeXShredderY()), month(sg_params())))
  expect_false(xs$r0$available)
  expect_true(any(grepl("R0: not available", capture.output(print(xs)))))
})

test_that("flow balance: round-off for the ODE, exactly zero for stochastic runs", {
  p <- sg_params()
  m <- sg_model(p, cube3, month(p), sg_releases("aa", "adult", "male", times = 3, number = 40),
    area = 100
  )
  sim <- simulate_ode(m)
  fb <- sg_flow_balance(sim)
  expect_identical(names(fb), c("day", "change", "releases", "eggs_laid", "deaths", "residual"))
  expect_lt(max(abs(fb$residual)) / max(rowSums(sim$state)), 1e-10)
  expect_identical(fb$releases[3], 40)
  init <- sg_init(m, values = c(E = 5, L = 5, P = 5, M = 3, F = 3), per = "km2")
  for (keep in c("state", "totals")) {
    cb <- simulate_cb(m, init, n_rep = 2, seed = 1, output_every = 7, keep = keep)
    fc <- sg_flow_balance(cb)
    expect_identical(fc$day, rep(c(7L, 14L, 21L, 28L, 31L), 2))
    expect_identical(max(abs(fc$residual)), 0)
  }
})

test_that("runtime guards stop runs whose numbers blow up", {
  p <- sg_params(eggs_per_day_max = 1e300)
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), month(p), area = 10)
  for (b in c("cpp", "R")) {
    expect_error(suppressWarnings(simulate_ode(m, backend = b)), "failed on day 1")
  }
  p2 <- sg_params(eggs_per_day_max = 1e308, clutch_size = 1e-300)
  m2 <- sg_model(p2, sg_cube(sg_test_cube_wild()), month(p2), area = 10)
  expect_error(suppressWarnings(simulate_cb(m2, seed = 1)), "non-finite or negative on day 1")
})
