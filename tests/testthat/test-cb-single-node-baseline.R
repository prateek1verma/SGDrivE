# One-node chain-binomial runs must reproduce the results of the code before
# phase 2 of the multi-node work (dev/MOVEMENT_PHASE2_PLAN.md, WP1). The
# fixture is made by data-raw/make_cb_baseline.R.
#
# With SCREWDRIVE_BITWISE_BASELINE=true (set by dev/check.R; the machine and R
# version that made the fixture) stochastic runs and steps must be identical
# draw for draw. Elsewhere only the deterministic map and the deterministic
# step are compared, with tolerances: on another platform a last-bit
# difference in a probability can change a single draw, after which a
# stochastic run takes a different path.

cbb_bitwise <- function() identical(Sys.getenv("SCREWDRIVE_BITWISE_BASELINE"), "true")

cbb_expect <- function(actual, expected, info, tolerance = 1e-9) {
  if (cbb_bitwise()) {
    testthat::expect_identical(actual, expected, info = info)
  } else {
    testthat::expect_equal(actual, expected, tolerance = tolerance, info = info)
  }
}

cbb_cases <- function() {
  wx <- sg_example_weather()
  mk <- function(p, cube, releases = NULL, area = 20) {
    sg_model(p, cube, sg_forcing(wx, p, start = "1971-04-01", end = "1971-06-29"), releases,
      area = area
    )
  }
  list(
    default_mendelian = mk(
      sg_params(), sg_cube(MGDrivE::cubeMendelian()),
      sg_releases("aa",
        stage = "adult", sex = "male", times = seq(20, by = 7, length.out = 6), number = 50
      )
    ),
    gp2014_homing = mk(
      sg_params(preset = "gp2014"), sg_cube(MGDrivE::cubeHomingDrive()),
      sg_releases("HH", stage = "adult", sex = "male", times = 15, number = 20)
    ),
    hosts_crowding_allee = mk(
      sg_params(
        wound_model = "hosts", host_density = 50, wounded_fraction = 0.1,
        larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
        male_density_half_mating = 2
      ),
      sg_cube(MGDrivE::cubeMendelian()),
      area = 10
    )
  )
}

test_that("one-node chain-binomial runs reproduce the pre-phase-2 baseline", {
  skip_if_not_installed("MGDrivE")
  ref <- readRDS(test_path("fixtures", "cb_single_node_baseline.rds"))
  cases <- cbb_cases()
  expect_identical(names(cases), names(ref))
  for (nm in names(cases)) {
    m <- cases[[nm]]
    r <- ref[[nm]]
    x0 <- sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10))
    det <- .sg_cb_run_cpp(m, init = x0, dt = 0.1, stochastic = FALSE)
    cbb_expect(unname(det$state[seq(1L, nrow(det$state), by = 10L), ]), r$deterministic$state, nm)
    cbb_expect(unname(det$counters), r$deterministic$counters, nm)
    view <- .sg_model_view(m)
    cbb_expect(sg_cb_step_cpp(view, 30L, 0.1, r$step$x, FALSE), r$step$deterministic, nm, 1e-12)
    if (cbb_bitwise()) {
      expect_identical(m$meta$hash, r$hash, info = nm)
      sto <- simulate_cb(m, init = x0, n_rep = 2L, dt = 0.1, output_every = 10L, seed = 11L)
      expect_identical(unname(sto$state), r$stochastic$state, info = nm)
      expect_identical(unname(sto$counters), r$stochastic$counters, info = nm)
      expect_identical(unname(sto$extinct_day), r$stochastic$extinct_day, info = nm)
      expect_identical(unname(sto$wild_extinct_day), r$stochastic$wild_extinct_day, info = nm)
      set.seed(5L)
      expect_identical(sg_cb_step_cpp(view, 30L, 0.1, r$step$x, TRUE), r$step$stochastic, info = nm)
    }
  }
})
