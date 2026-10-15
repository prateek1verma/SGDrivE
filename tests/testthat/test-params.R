test_that("sg_params() returns validated defaults in table order", {
  p <- sg_params()
  expect_s3_class(p, "sg_params")
  expect_identical(attr(p, "preset"), "default")
  expect_identical(names(p), sg_param_table()$name)
  expect_identical(unclass(p)[names(p)], SGDrivE:::.sg_param_defaults())
})

test_that("integer parameters are stored as integers", {
  p <- sg_params(adult_substages = 8)
  expect_identical(p$adult_substages, 8L)
})

test_that("overrides change only the named values", {
  p <- sg_params(eggs_per_day_max = 67, mating_rate = 1)
  expect_identical(p$eggs_per_day_max, 67)
  expect_identical(p$mating_rate, 1)
  expect_identical(p$adult_lifespan_dd, 202.8)
})

test_that("unknown, unnamed, duplicated and NULL parameters are rejected", {
  expect_error(
    sg_params(dd_scale_preovipositon = 60),
    "Unknown parameter.*Did you mean 'dd_scale_preoviposition'"
  )
  expect_error(sg_params(R = 60), "'R' \\(renamed to 'eggs_per_day_max'\\)") # D-027
  expect_error(sg_params(n_A = 8, X_A = 200), "'n_A' \\(renamed to 'adult_substages'\\)")
  expect_error(sg_params(foo_bar_baz = 1), "Unknown parameter.*sg_param_table")
  expect_error(sg_params(8), "must be named")
  expect_error(sg_params(eggs_per_day_max = 1, eggs_per_day_max = 2), "more than once")
  expect_error(sg_params(eggs_per_day_max = NULL), "is NULL")
})

test_that("type, length and range errors are informative", {
  expect_error(sg_params(adult_substages = 2.5), "whole number")
  expect_error(sg_params(adult_substages = "ten"), "finite number")
  expect_error(sg_params(eggs_per_day_max = c(1, 2)), "single value")
  expect_error(sg_params(eggs_per_day_max = Inf), "finite number")
  expect_error(sg_params(eggs_per_day_max = NA), "cannot be NA")
  expect_error(sg_params(larva_survival_daily = 1.2), "outside its allowed range")
  expect_error(sg_params(adult_lifespan_dd = 0), "outside its allowed range \\(0") # open bound
  expect_error(sg_params(adult_substages = 0), "outside its allowed range")
  expect_error(sg_params(wound_competitors = "all"), "must be one of")
  expect_error(sg_params(mating_mode = "rates"), "Did you mean 'rate'")
})

test_that("cross-checks catch impossible combinations", {
  expect_error(sg_params(oviposition_age_dd = 300), "shorter than the adult lifespan")
  expect_silent(sg_params(dd_scale_preoviposition = 300)) # X_pre only scales degree-days (D-018)
  expect_error(sg_params(female_receptive_age_dd = 202.8), "shorter than the adult lifespan")
  expect_error(sg_params(dev_temp_min = 30, dev_temp_decline = 25), "greater than dev_temp_min")
  expect_error(sg_params(dev_temp_min = 35, oviposition_temp_max = 30), "greater than dev_temp_min")
})

test_that("ineffective settings give warnings", {
  expect_warning(sg_params(larval_crowding_coef = 0.1), "no effect")
  expect_silent(sg_params(larval_crowding_coef = 0.1, larval_crowding_mode = "per_wound"))
  expect_warning(sg_params(mating_mode = "instant", female_receptive_age_dd = 10), "ignored")
})

test_that("gp2014 preset reproduces the Gutierrez & Ponti code configuration", {
  gp <- sg_params(preset = "gp2014")
  expect_identical(attr(gp, "preset"), "gp2014")
  expect_identical(
    unlist(unclass(gp)[c("egg_substages", "larva_substages", "pupa_substages", "adult_substages")]),
    c(egg_substages = 25L, larva_substages = 25L, pupa_substages = 25L, adult_substages = 25L)
  )
  expect_equal(gp$oviposition_age_dd, 56.784)
  expect_identical(gp$dev_rate_floor, 0.01)
  expect_identical(gp$wound_competitors, "all_adults")
  expect_identical(gp$mating_mode, "instant")
  # overrides apply on top of the preset
  expect_identical(sg_params(preset = "gp2014", eggs_per_day_max = 67)$eggs_per_day_max, 67)
  expect_identical(sg_params(preset = "gp2014", eggs_per_day_max = 67)$adult_substages, 25L)
})

test_that("update() changes values, re-validates and keeps the preset", {
  gp <- sg_params(preset = "gp2014")
  gp2 <- update(gp, eggs_per_day_max = 50)
  expect_s3_class(gp2, "sg_params")
  expect_identical(gp2$eggs_per_day_max, 50)
  expect_identical(attr(gp2, "preset"), "gp2014")
  expect_error(update(gp, oviposition_age_dd = 300), "shorter than the adult lifespan")
})

test_that("as.data.frame() lists every parameter and flags changes", {
  df <- as.data.frame(sg_params(eggs_per_day_max = 67))
  expect_identical(nrow(df), nrow(sg_param_table()))
  expect_identical(df$name[df$changed], "eggs_per_day_max")
  expect_identical(df$value[df$name == "eggs_per_day_max"], "67")
})

test_that("print() shows groups and marks changed values", {
  expect_snapshot(print(sg_params(adult_substages = 8)))
  expect_invisible(print(sg_params()))
})
