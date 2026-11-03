wx <- sg_example_weather()
legacy <- readRDS(test_path("fixtures", "gp2014_mcallen_forcing.rds"))

test_that("McAllen forcing matches the legacy daily terms exactly (gp2014 preset)", {
  f <- sg_forcing(wx, sg_params(preset = "gp2014"))
  expect_s3_class(f, "sg_forcing")
  expect_identical(f$n_days, nrow(legacy))
  expect_identical(f$dates, as.Date(legacy$date))
  expect_equal(f$tmean, legacy$tmean, tolerance = sg_tol_exact)
  expect_equal(f$tmean5, legacy$aver5, tolerance = sg_tol_exact)
  expect_equal(f$dx, legacy$dd, tolerance = sg_tol_exact)
  expect_equal(f$p_A, legacy$temp_mort, tolerance = sg_tol_exact)
  expect_equal(exp(-f$mu_T), legacy$temp_surv, tolerance = sg_tol_exact)
  expect_equal(f$egg_viability, pmax(0, legacy$egg_viab), tolerance = sg_tol_exact)
  expect_equal(f$mu_E, -log(pmax(1e-5, legacy$egg_viab)), tolerance = sg_tol_exact)
  expect_equal(f$phi_ov, legacy$phi_temp, tolerance = sg_tol_exact)
})

test_that("the forcing records the parameter values it used", {
  f <- sg_forcing(wx[1:30, ], sg_params(optimal_survival_temp = 27))
  expect_identical(f$param_values$optimal_survival_temp, 27)
  expect_named(f$param_values, ScrewDrivE:::.sg_forcing_param_names())
})

test_that("start/end select the simulation period; the 5-day mean uses earlier days", {
  full <- sg_forcing(wx)
  f <- sg_forcing(wx, start = "1961-01-01", end = "1961-12-31")
  expect_identical(f$n_days, 365L)
  expect_identical(f$day, 1:365)
  expect_identical(f$dates[1], as.Date("1961-01-01"))
  i <- which(full$dates == as.Date("1961-01-01"))
  expect_identical(f$tmean5[1], full$tmean5[i])
  expect_error(sg_forcing(wx, start = "1950-01-01"), "not covered")
  expect_error(sg_forcing(wx, start = "1962-01-01", end = "1961-01-01"), "before `start`")
  expect_error(sg_forcing(wx, start = "not a date"), "single date")
})

test_that("tmean is computed from tmax and tmin when absent", {
  w <- wx[1:20, c("date", "tmax", "tmin")]
  f <- sg_forcing(w)
  expect_equal(f$tmean, (w$tmax + w$tmin) / 2)
  w$tmin[3] <- w$tmax[3] + 1
  expect_warning(sg_forcing(w), "tmin > tmax")
})

test_that("host egg temperature gives constant egg degree-days (D-014)", {
  f <- sg_forcing(wx[1:50, ], sg_params(egg_temp_source = "host"))
  expect_identical(unique(f$dx_egg), sg_degree_days(38.5))
  f2 <- sg_forcing(wx[1:50, ])
  expect_identical(f2$dx_egg, f2$dx)
})

test_that("moisture index follows annual rainfall and is 1 without rain", {
  f <- sg_forcing(wx)
  yr <- format(f$dates, "%Y")
  i1983 <- yr == "1983"
  rain_1983 <- sum(wx$rain[format(wx$date, "%Y") == "1983"])
  expect_equal(unique(f$moisture[i1983]), min(1, rain_1983 / 450))
  expect_true(all(f$moisture >= 0 & f$moisture <= 1))
  dry <- wx[1:400, c("date", "tmean")]
  expect_identical(unique(sg_forcing(dry)$moisture), 1)
  expect_false(sg_forcing(dry)$has_rain)
  # an incomplete year is scaled to a full year
  half <- data.frame(date = as.Date("2001-01-01") + 0:181, tmean = 25, rain = 1)
  m_half <- sg_forcing(half, sg_params(drought_free_rainfall = 400))$moisture
  expect_equal(unique(m_half), min(1, 365 / 400))
})

test_that("invalid weather gives informative errors", {
  w <- wx[1:30, ]
  expect_error(sg_forcing(as.list(w)), "must be a data frame")
  expect_error(sg_forcing(w[, c("tmean", "rain")]), "needs a `date` column")
  expect_error(sg_forcing(w[, c("date", "rain")]), "needs a `tmean` column")
  expect_error(sg_forcing(w[0, ]), "no rows")
  expect_error(sg_forcing(w[c(2, 1, 3:30), ]), "sorted by date")
  expect_error(sg_forcing(w[c(1, 1:30), ]), "repeated dates")
  expect_error(sg_forcing(w[-5, ]), "gap after 1960-01-04")
  bad <- w
  bad$tmean[7] <- NA
  expect_error(sg_forcing(bad), "missing temperatures, e.g. on 1960-01-07")
  bad <- w
  bad$tmean[2] <- 300
  expect_error(sg_forcing(bad), "Implausible")
  bad <- w
  bad$rain[3] <- -1
  expect_error(sg_forcing(bad), "rain` must be numeric, >= 0")
  bad <- w
  bad$date <- "x"
  expect_error(sg_forcing(bad), "must be dates")
  bad <- w
  bad$tmean <- as.character(bad$tmean)
  expect_error(sg_forcing(bad), "must be numeric")
  expect_error(sg_forcing(w, params = list()), "sg_params object")
})

test_that("print() summarises the forcing", {
  expect_snapshot(print(sg_forcing(wx, start = "1971-01-01", end = "1971-12-31")))
})
