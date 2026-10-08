test_that("sg_params() returns validated defaults in table order", {
  p <- sg_params()
  expect_s3_class(p, "sg_params")
  expect_identical(attr(p, "preset"), "default")
  expect_identical(names(p), sg_param_table()$name)
  expect_identical(unclass(p)[names(p)], SGDrivE:::.sg_param_defaults())
})

test_that("integer parameters are stored as integers", {
  p <- sg_params(n_A = 8)
  expect_identical(p$n_A, 8L)
  expect_identical(sg_params()$i0, NA_integer_)
  expect_identical(sg_params(i0 = 4)$i0, 4L)
})

test_that("overrides change only the named values", {
  p <- sg_params(R = 67, nu = 1)
  expect_identical(p$R, 67)
  expect_identical(p$nu, 1)
  expect_identical(p$X_A, 202.8)
})

test_that("unknown, unnamed, duplicated and NULL parameters are rejected", {
  expect_error(sg_params(X_pree = 60), "Unknown parameter.*Did you mean 'X_pre'")
  expect_error(sg_params(foo_bar_baz = 1), "Unknown parameter.*sg_param_table")
  expect_error(sg_params(8), "must be named")
  expect_error(sg_params(R = 1, R = 2), "more than once")
  expect_error(sg_params(R = NULL), "is NULL")
})

test_that("type, length and range errors are informative", {
  expect_error(sg_params(n_A = 2.5), "whole number")
  expect_error(sg_params(n_A = "ten"), "finite number")
  expect_error(sg_params(R = c(1, 2)), "single value")
  expect_error(sg_params(R = Inf), "finite number")
  expect_error(sg_params(R = NA), "cannot be NA")
  expect_error(sg_params(s_L = 1.2), "outside its allowed range")
  expect_error(sg_params(X_A = 0), "outside its allowed range \\(0") # open bound
  expect_error(sg_params(n_A = 0), "outside its allowed range")
  expect_error(sg_params(search_denominator = "all"), "must be one of")
  expect_error(sg_params(mating = "rates"), "Did you mean 'rate'")
})

test_that("cross-checks catch impossible combinations", {
  expect_error(sg_params(X_pre = 300), "shorter than the adult lifespan")
  expect_error(sg_params(X_rec = 202.8), "shorter than the adult lifespan")
  expect_error(sg_params(i0 = 11), "exceeds the number of adult classes")
  expect_error(sg_params(T_L = 30, T_U = 25), "greater than T_L")
  expect_error(sg_params(T_L = 35, TU_ovip = 30), "greater than T_L")
})

test_that("ineffective settings give warnings", {
  expect_warning(sg_params(gamma_L = 0.1), "no effect")
  expect_silent(sg_params(gamma_L = 0.1, larval_dd = "per_wound"))
  expect_warning(sg_params(mating = "instant", X_rec = 10), "ignored")
})

test_that("gp2014 preset reproduces the Gutierrez & Ponti code configuration", {
  gp <- sg_params(preset = "gp2014")
  expect_identical(attr(gp, "preset"), "gp2014")
  expect_identical(
    unlist(unclass(gp)[c("n_E", "n_L", "n_P", "n_A")]),
    c(n_E = 25L, n_L = 25L, n_P = 25L, n_A = 25L)
  )
  expect_identical(gp$i0, 8L)
  expect_identical(gp$dx_floor, 0.01)
  expect_identical(gp$search_denominator, "all_adults")
  expect_identical(gp$mating, "instant")
  # overrides apply on top of the preset
  expect_identical(sg_params(preset = "gp2014", R = 67)$R, 67)
  expect_identical(sg_params(preset = "gp2014", R = 67)$n_A, 25L)
})

test_that("update() changes values, re-validates and keeps the preset", {
  gp <- sg_params(preset = "gp2014")
  gp2 <- update(gp, R = 50)
  expect_s3_class(gp2, "sg_params")
  expect_identical(gp2$R, 50)
  expect_identical(attr(gp2, "preset"), "gp2014")
  expect_error(update(gp, i0 = 30), "exceeds the number of adult classes")
  expect_error(update(gp, n_A = 5), "exceeds the number of adult classes") # preset i0 is 8, above 5
})

test_that("as.data.frame() lists every parameter and flags changes", {
  df <- as.data.frame(sg_params(R = 67))
  expect_identical(nrow(df), nrow(sg_param_table()))
  expect_identical(df$name[df$changed], "R")
  expect_identical(df$value[df$name == "R"], "67")
})

test_that("print() shows groups and marks changed values", {
  expect_snapshot(print(sg_params(n_A = 8)))
  expect_invisible(print(sg_params()))
})
