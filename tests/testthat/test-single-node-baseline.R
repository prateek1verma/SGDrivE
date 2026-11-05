# One-node models must reproduce the results of the code before the
# multi-node work (dev/MOVEMENT_PHASE1_PLAN.md, section 6). The fixture is
# made by data-raw/make_single_node_baseline.R.
#
# Bitwise equality only holds on the machine, compiler and R version that
# made the fixture: other platforms differ in the last bits of floating-point
# results (and so in the model hash, computed from 17-digit numbers). With
# SCREWDRIVE_BITWISE_BASELINE=true (set by dev/check.R) the test demands
# identical results; otherwise (e.g. on CI) it compares with tight
# tolerances, far below any model change.
baseline_bitwise <- function() identical(Sys.getenv("SCREWDRIVE_BITWISE_BASELINE"), "true")

expect_baseline <- function(actual, expected, info, tolerance) {
  if (baseline_bitwise()) {
    expect_identical(actual, expected, info = info)
  } else {
    expect_equal(actual, expected, tolerance = tolerance, info = info)
  }
}

baseline_cases <- function() {
  list(
    default_mendelian = list(
      params = sg_params(),
      cube = sg_cube(MGDrivE::cubeMendelian()),
      releases = sg_releases("aa",
        stage = "adult", sex = "male",
        times = seq(60, by = 7, length.out = 8), number = 200
      )
    ),
    gp2014_homing = list(
      params = sg_params(preset = "gp2014"),
      cube = sg_cube(MGDrivE::cubeHomingDrive()),
      releases = sg_releases("HH", stage = "adult", sex = "male", times = 30, number = 50)
    ),
    hosts_crowding_mendelian = list(
      params = sg_params(
        wound_model = "hosts", host_density = 50, wounded_fraction = 0.1,
        larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
        male_density_half_mating = 5
      ),
      cube = sg_cube(MGDrivE::cubeMendelian()),
      releases = NULL
    )
  )
}

test_that("one-node models reproduce the pre-phase-1 baseline", {
  skip_if_not_installed("MGDrivE")
  ref <- readRDS(test_path("fixtures", "single_node_baseline.rds"))
  cases <- baseline_cases()
  expect_identical(names(cases), names(ref))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    r <- ref[[nm]]
    wx <- sg_forcing(sg_example_weather(), cs$params, start = "1971-01-01", end = "1971-12-31")
    mod <- sg_model(cs$params, cs$cube, wx, cs$releases, area = 2)
    if (baseline_bitwise()) expect_identical(mod$meta$hash, r$hash, info = nm)
    expect_identical(sg_state_names(mod), r$state_names, info = nm)
    expect_identical(mod$layout, r$layout, info = nm)
    expect_baseline(mod$counts, r$counts, nm, 1e-14)
    expect_baseline(mod$releases, r$releases, nm, 1e-14)
    sim <- simulate_ode(mod, output_every = 30L)
    expect_identical(sim$day, r$day, info = nm)
    expect_baseline(unname(sim$state), r$state, nm, 1e-9)
    expect_baseline(sim$counters, r$counters, nm, 1e-9)
    expect_baseline(sim$diagnostics$max_rate, r$max_rate, nm, 1e-9)
    view <- .sg_model_view(mod)
    for (pr in r$probe) {
      expect_baseline(sg_rates_cpp(view, pr$day, pr$state), pr$rates, nm, 1e-12)
      expect_baseline(sg_max_rate_cpp(view, pr$day, pr$state), pr$max_rate, nm, 1e-12)
    }
  }
})
