# M3 milestone check: the C++ rate kernel equals the R reference rate function
# (.sg_rates_r, itself checked against hand calculations and the legacy model)
# to floating-point round-off, for every option, cube type and preset.

rates_r <- SGDrivE:::.sg_rates_r
setup <- SGDrivE:::.sg_reference_setup
view <- SGDrivE:::.sg_model_view
rates_cpp <- SGDrivE:::sg_rates_cpp
max_rate_r <- SGDrivE:::.sg_max_rate
max_rate_cpp <- SGDrivE:::sg_max_rate_cpp

wx_all <- sg_example_weather()
days <- c(1, 45, 120, 200, 300) # winter, spring, summer, autumn (1971)

expect_rates_equal <- function(m, x, label) {
  s <- setup(m)
  v <- view(m)
  for (d in days) {
    r <- rates_r(d, x, s)
    cc <- rates_cpp(v, d, x)
    scale_d <- max(abs(r$deriv), 1e-300)
    scale_c <- max(abs(r$counters), 1e-300)
    testthat::expect_lt(max(abs(r$deriv - cc$deriv)) / scale_d, 1e-12,
      label = paste(label, "deriv, day", d)
    )
    testthat::expect_lt(max(abs(r$counters - cc$counters)) / scale_c, 1e-12,
      label = paste(label, "counters, day", d)
    )
    testthat::expect_equal(max_rate_cpp(v, d, x), max_rate_r(d, x, s), tolerance = 1e-14)
  }
}

rand_state <- function(m, seed) {
  set.seed(seed)
  stats::runif(m$n_states, 0, 10) * stats::rbinom(m$n_states, 1, 0.9) # some zeros
}

mk <- function(p, cube, ...) {
  wx <- sg_forcing(wx_all, p, start = "1971-01-01", end = "1971-12-31")
  sg_model(p, cube, wx, ...)
}

cubes <- list(
  wild = sg_cube(sg_test_cube_wild()),
  sit = sg_cube(sg_test_cube_sit(residual = 1.2e-8), lifespan_fraction = c(W = 1, S = 0.52)),
  mendelian = sg_cube(sg_test_cube_mendelian()),
  short_lived = sg_cube(sg_test_cube_sit(), lifespan_fraction = c(W = 1, S = 0.25)),
  infertile_females = sg_cube(sg_test_cube_sit(residual = 0, s_s = 0)),
  partly_fertile = sg_cube(sg_test_cube_sit(residual = 0.3, s_s = 0.4))
)

param_sets <- list(
  default = sg_params(),
  instant = sg_params(mating_mode = "instant"),
  all_adults = sg_params(wound_competitors = "all_adults"),
  crowding = sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.05),
  allee = sg_params(male_density_half_mating = 2),
  thresholds = sg_params(female_receptive_age_dd = 24, male_mature_age_dd = 15),
  background = sg_params(female_background_mortality = 0.15, male_background_mortality = 0.08),
  moisture = sg_params(pupa_drought_mortality = 0.2, drought_free_rainfall = 2000),
  sigma0 = sg_params(wound_finding_max = 0.05),
  host_eggs = sg_params(egg_temp_source = "host"),
  gp2014 = sg_params(preset = "gp2014"),
  custom_stages = sg_params(
    egg_substages = 1, larva_substages = 2, pupa_substages = 3, adult_substages = 4
  )
)

test_that("C++ rates equal R rates for every option and hand-built cube", {
  for (pn in names(param_sets)) {
    for (cn in names(cubes)) {
      m <- mk(param_sets[[pn]], cubes[[cn]], area = 3)
      expect_rates_equal(m, rand_state(m, seed = nchar(pn) + nchar(cn)), paste(pn, cn))
    }
  }
})

test_that("C++ rates equal R rates for MGDrivE cubes, including sex-linked ones", {
  skip_if_not_installed("MGDrivE")
  for (raw in list(
    MGDrivE::cubeHomingDrive(), MGDrivE::cubeXShredderY(), MGDrivE::cubeRIDL(),
    MGDrivE::cubeReciprocalTranslocations()
  )) {
    cb <- sg_cube(raw)
    for (p in list(sg_params(), sg_params(mating_mode = "instant", male_density_half_mating = 1))) {
      m <- mk(p, cb)
      expect_rates_equal(m, rand_state(m, seed = cb$n_genotypes), paste(cb$genotypes[1]))
    }
  }
})

test_that("C++ rates equal R rates in edge cases", {
  m <- mk(sg_params(), cubes$mendelian)
  # empty population
  z <- numeric(m$n_states)
  cc <- rates_cpp(view(m), 100, z)
  expect_identical(cc$deriv, numeric(m$n_states))
  expect_true(all(cc$counters == 0))
  # no males at all: no mating, instant mating sends females to U
  for (p in list(sg_params(), sg_params(mating_mode = "instant"))) {
    mm <- mk(p, cubes$mendelian)
    x <- rand_state(mm, 3)
    x[sg_state_index(mm, "M")] <- 0
    expect_rates_equal(mm, x, paste("no males", p$mating_mode))
  }
  # males present but rejected by every female (eta = 0)
  raw <- sg_test_cube_sit()
  raw$eta["W", ] <- 0
  ms <- mk(sg_params(), sg_cube(raw))
  expect_rates_equal(ms, rand_state(ms, 4), "eta zero")
  # minimal chains: one juvenile sub-stage, two adult sub-stages (one per segment)
  m1 <- mk(
    sg_params(
      egg_substages = 1, larva_substages = 1,
      pupa_substages = 1, adult_substages = 2
    ),
    cubes$wild
  )
  expect_rates_equal(m1, rand_state(m1, 5), "minimal chains")
})

test_that("the C++ test entry points check their input", {
  m <- mk(sg_params(), cubes$wild)
  v <- view(m)
  expect_error(rates_cpp(v, 1, numeric(3)), "wrong length")
  expect_error(rates_cpp(v, 0, numeric(m$n_states)), "outside the forcing")
  expect_error(max_rate_cpp(v, 366, numeric(m$n_states)), "outside the forcing")
})

test_that("C++ day rates and state aggregates equal the R reference", {
  agg_cpp <- SGDrivE:::sg_aggregates_cpp
  for (pn in c("default", "instant", "all_adults", "crowding", "allee", "gp2014")) {
    for (cn in c("sit", "mendelian")) {
      m <- mk(param_sets[[pn]], cubes[[cn]], area = 3)
      s <- setup(m)
      v <- view(m)
      x <- rand_state(m, seed = 7)
      for (d in days) {
        a <- rates_r(d, x, s, aux = TRUE)$aux
        cc <- agg_cpp(v, d, x)
        lab <- paste(pn, cn, d)
        expect_equal(cc$q, unname(a$q), tolerance = 1e-14, label = lab)
        expect_equal(cc$mu, unname(a$mu[c("E", "L", "P", "F", "M")]), tolerance = 1e-14)
        expect_equal(cc$w_lay, as.vector(a$w_lay), tolerance = 1e-14) # mate fastest
        expect_equal(cc$n_s, a$n_s, tolerance = 1e-14)
        expect_equal(cc$psi, a$psi, tolerance = 1e-14)
        expect_equal(cc$births, unname(a$births), tolerance = 1e-13)
        expect_equal(cc$m_bar, unname(a$m_bar), tolerance = 1e-14)
        expect_equal(cc$g_find, a$g_find, tolerance = 1e-14)
        expect_equal(cc$pi, as.vector(t(a$pi)), tolerance = 1e-14) # [g, m], m fastest
        expect_identical(as.logical(cc$mate_ok), unname(rowSums(a$pi) > 0))
      }
    }
  }
  expect_error(agg_cpp(view(m), 0, x), "outside the forcing")
  expect_error(agg_cpp(view(m), 1, x[-1]), "wrong length")
})
