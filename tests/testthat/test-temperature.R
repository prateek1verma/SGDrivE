grid <- utils::read.csv(test_path("fixtures", "gp2014_temperature_grid.csv"))
gp <- sg_params(preset = "gp2014")

test_that("degree-days match the legacy Gutierrez & Ponti code exactly (gp2014 preset)", {
  expect_equal(sg_degree_days(grid$temp, gp), grid$degree_days, tolerance = sg_tol_exact)
})

test_that("temperature mortality matches the legacy code (proportion and survival)", {
  expect_equal(sg_temp_mortality(grid$temp, gp, type = "proportion"), grid$temp_mortality,
    tolerance = sg_tol_exact
  )
  expect_equal(exp(-sg_temp_mortality(grid$temp, gp)), grid$temp_survival,
    tolerance = sg_tol_exact
  )
})

test_that("egg viability matches the legacy code where it is non-negative", {
  # The legacy code keeps negative 'viabilities' and clamps them when applied;
  # ScrewDrivE clamps at 0 directly. The hazards are identical.
  expect_equal(sg_egg_viability(grid$temp, gp), pmax(0, grid$egg_viability),
    tolerance = sg_tol_exact
  )
  expect_equal(sg_egg_viability(grid$temp, gp, type = "hazard"),
    -log(pmax(1e-5, grid$egg_viability)),
    tolerance = sg_tol_exact
  )
})

test_that("the oviposition scalar matches the legacy code", {
  expect_equal(sg_ovip_scalar(grid$temp, gp), grid$ovip_scalar, tolerance = sg_tol_exact)
})

test_that("known values from the paper are reproduced", {
  p <- sg_params()
  expect_equal(round(sg_degree_days(c(27.2, 38.5), p), 2), c(12.08, 17.22))
  expect_equal(sg_temp_mortality(27.2, p, type = "proportion"), 0.0035)
  expect_equal(sg_egg_viability(27.2, p), 0.7)
  expect_equal(sg_ovip_scalar(29, p), 1)
  expect_identical(sg_ovip_scalar(c(14.5, 43.5, 10, 45), p), c(0, 0, 0, 0))
})

test_that("development stops below T_L with the default floor and not with gp2014", {
  expect_identical(sg_degree_days(c(5, 14.5), sg_params()), c(0, 0))
  expect_equal(sg_degree_days(5, gp), 0.61)
})

test_that("functions use the supplied parameters", {
  p <- sg_params(dd_scale_preoviposition = 122)
  expect_equal(sg_degree_days(27.2, p), 2 * sg_degree_days(27.2, sg_params()))
  expect_equal(
    sg_temp_mortality(27.2, sg_params(adult_mortality_min = 0.01), type = "proportion"),
    0.01
  )
})

test_that("functions are vectorised and validate input", {
  expect_length(sg_degree_days(1:10), 10)
  expect_error(sg_degree_days("20"), "must be numeric")
  expect_error(sg_degree_days(20, params = list()), "sg_params object")
  expect_error(sg_temp_mortality(20, type = "rate"), "should be one of")
})
