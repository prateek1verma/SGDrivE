# M4 milestone check: the C++ chain-binomial step and day loop
# (src/cb_step.h, src/cb_step.cpp) reproduce the R reference step draw for
# draw with the same seed, and its deterministic map to round-off.

ns <- asNamespace("SGDrivE")
run_r <- ns$.sg_cb_run_r
run_cpp <- ns$.sg_cb_run_cpp
wx_all <- sg_example_weather()
sit_raw <- sg_test_cube_sit(residual = 0.1)
sit_raw$eta[, "S"] <- 0.5

mk <- function(p, cube, end = "1971-05-15", releases = NULL, area = 50, start = "1971-05-01") {
  sg_model(p, cube, sg_forcing(wx_all, p, start, end), releases, area = area)
}
flow_error <- function(r) {
  k <- r$counters
  births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
  deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
  change <- rowSums(r$releases) + births - deaths - k[, "emergence_deaths"]
  max(abs(diff(rowSums(r$state)) - change))
}
init_of <- function(m) sg_init(m, values = c(E = 5, L = 5, P = 5, M = 3, U = 1, F = 3))

cases <- list(
  "defaults (batches, rate mating)" = list(sg_params(), sg_cube(sg_test_cube_mendelian())),
  "single eggs, instant mating, Allee" = list(
    sg_params(oviposition = "eggs", mating = "instant", h_M = 1),
    sg_cube(sg_test_cube_mendelian())
  ),
  "crowding, clutch size, SIT cube" = list(
    sg_params(larval_dd = "per_wound", gamma_L = 0.05, clutch_size = 200), sg_cube(sit_raw)
  ),
  "short-lived genotype, receptivity, all adults search" = list(
    sg_params(X_rec = 20, search_denominator = "all_adults"),
    sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25))
  ),
  "gp2014 preset" = list(sg_params(preset = "gp2014"), sg_cube(sg_test_cube_wild()))
)

for (nm in names(cases)) {
  test_that(paste("C++ runs equal the R reference draw for draw:", nm), {
    p <- cases[[nm]][[1]]
    cb <- cases[[nm]][[2]]
    last_g <- cb$genotypes[cb$n_genotypes]
    rel <- sg_releases(last_g, "adult", "male", times = c(3, 9), number = 40.5)
    m <- mk(p, cb, releases = rel)
    set.seed(11)
    a <- run_r(m, init_of(m))
    set.seed(11)
    b <- run_cpp(m, init_of(m))
    expect_identical(b$state, a$state)
    expect_identical(b$counters, a$counters)
    expect_identical(b$releases, a$releases)
    expect_identical(flow_error(b), 0)
    da <- run_r(m, init_of(m), stochastic = FALSE)
    db <- run_cpp(m, init_of(m), stochastic = FALSE)
    expect_lt(max(abs(da$state - db$state)) / max(da$state), 1e-12)
    expect_lt(max(abs(da$counters - db$counters)) / max(abs(da$counters)), 1e-12)
  })
}

test_that("a single C++ step equals the R step", {
  skip_if_not_installed("MGDrivE")
  for (cb in list(sg_cube(MGDrivE::cubeHomingDrive()), sg_cube(MGDrivE::cubeXShredderY()))) {
    m <- mk(sg_params(), cb)
    s <- ns$.sg_reference_setup(m)
    v <- ns$.sg_model_view(m)
    # whole counts in every block, adults only of genotypes of their sex
    x <- round(as.numeric(sg_init(m, values = c(E = 20, L = 20, P = 20, M = 8, U = 2, F = 8))))
    info <- ns$.sg_cb_day_r(7, s, 0.1)
    det_r <- ns$.sg_cb_step_r(x, info, s, 0.1, ns$.sg_cb_draws(FALSE))
    det_c <- ns$sg_cb_step_cpp(v, 7L, 0.1, x, FALSE)
    expect_lt(max(abs(det_r$x - det_c$x)) / max(det_r$x), 1e-12)
    expect_lt(max(abs(det_r$flows - det_c$flows)) / max(abs(det_r$flows)), 1e-12)
    set.seed(5)
    st_r <- ns$.sg_cb_step_r(x, info, s, 0.1, ns$.sg_cb_draws(TRUE))
    set.seed(5)
    st_c <- ns$sg_cb_step_cpp(v, 7L, 0.1, x, TRUE)
    expect_identical(st_c$x, st_r$x)
    expect_identical(st_c$flows, unname(st_r$flows))
  }
})

test_that("a year of C++ runs: whole counts, exact flow accounting, cost independent of size", {
  sit <- sg_test_cube_sit(residual = 0)
  sit$eta[, "S"] <- 0.35
  p <- sg_params()
  for (area in c(1000, 1e6)) {
    m <- mk(p, sg_cube(sit),
      start = "1971-01-01", end = "1971-12-31", area = area,
      releases = sg_releases("S", "adult", "male", times = seq(1, 365, 7), number = 30, per = "km2")
    )
    set.seed(1)
    r <- run_cpp(m, sg_init(m, per = "km2"))
    expect_true(all(r$state == round(r$state)) && all(r$state >= 0))
    expect_identical(flow_error(r), 0)
    expect_identical(r$extinct_day, 0L)
    if (area == 1e6) expect_gt(max(r$state), 1e7) # large counts are handled
  }
})

test_that("runs stop early once extinct with no releases left, and continue otherwise", {
  p <- sg_params()
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), sg_forcing(wx_all, p, "1971-01-01", "1971-12-31"))
  set.seed(3)
  r <- run_cpp(m, sg_init(m))
  expect_gt(r$extinct_day, 0) # the 1 km^2 population dies out in winter
  expect_identical(r$days_run, r$extinct_day)
  expect_true(all(r$state[(r$extinct_day + 1):366, ] == 0))
  expect_true(all(r$counters[(r$extinct_day + 1):365, ] == 0))
  # a release after the extinction: the run goes on and the population returns
  m2 <- sg_model(
    p, sg_cube(sg_test_cube_wild()), sg_forcing(wx_all, p, "1971-01-01", "1971-12-31"),
    sg_releases("W", "pupa", times = 200, number = 500)
  )
  set.seed(3)
  r2 <- run_cpp(m2, sg_init(m2))
  expect_identical(r2$extinct_day, r$extinct_day)
  expect_identical(r2$days_run, 365L)
  expect_gt(sum(r2$state[202, ]), 0)
  expect_identical(flow_error(r2), 0)
})

test_that("the C++ entry points check their input", {
  m <- mk(sg_params(), sg_cube(sg_test_cube_wild()))
  v <- ns$.sg_model_view(m)
  x <- round(as.numeric(sg_init(m, values = c(E = 4))))
  rec <- rep(TRUE, 15)
  expect_error(ns$sg_cb_run_cpp(v, x[-1], numeric(0), 0.1, TRUE, 15L, rec), "wrong length")
  expect_error(ns$sg_cb_run_cpp(v, x, numeric(0), 0.3, TRUE, 15L, rec), "whole number of steps")
  expect_error(ns$sg_cb_run_cpp(v, x + 0.5, numeric(0), 0.1, TRUE, 15L, rec), "whole, non-negative")
  expect_error(ns$sg_cb_run_cpp(v, x, numeric(0), 0.1, TRUE, 16L, rec), "n_days")
  expect_error(ns$sg_cb_run_cpp(v, x, numeric(0), 0.1, TRUE, 15L, rec[-1]), "one value per day")
  expect_error(ns$sg_cb_run_cpp(v, x, 1, 0.1, TRUE, 15L, rec), "one value per release")
  expect_error(ns$sg_cb_step_cpp(v, 0L, 0.1, x, TRUE), "outside the forcing")
  expect_error(ns$sg_cb_step_cpp(v, 1L, 0.1, x[-1], TRUE), "wrong length")
  expect_error(ns$sg_cb_step_cpp(v, 1L, 0.1, -x - 1, TRUE), "whole, non-negative")
  # deterministic runs accept fractional counts
  expect_silent(ns$sg_cb_step_cpp(v, 1L, 0.1, x + 0.5, FALSE))
})
