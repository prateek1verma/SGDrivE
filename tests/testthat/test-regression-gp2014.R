# M2 milestone check: with one wild genotype and the "gp2014" preset, the
# reference ODE reproduces the legacy Gutierrez & Ponti (2014) ODE translation
# (run without population floors at rtol 1e-10; data-raw/make_reference_tight.R)
# to a relative error of 1e-6 in every daily total. Milestone M3 check: so does
# the C++ backend over the full 24 years, routinely (D-021).

skip_if_not_installed("deSolve")

gp <- sg_params(preset = "gp2014")
wild <- sg_cube(sg_test_cube_wild())

# Daily totals of an sg_sim in the column layout of the legacy output.
legacy_totals <- function(sim) {
  tot <- SGDrivE:::.sg_sim_block_totals(sim)[-1, , drop = FALSE]
  k <- sim$counters
  data.frame(
    eggs = tot[, "E"], larvae = tot[, "L"], pupae = tot[, "P"],
    females = tot[, "U"] + tot[, "F"], adults = tot[, "M"] + tot[, "U"] + tot[, "F"],
    eggs_laid = k[, "eggs_laid[W]"], egg_batches = k[, "egg_batches"],
    pupae_in = k[, "pupae_in"], adults_in = k[, "emerged_male"] + k[, "emerged_female"]
  )
}
max_rel_error <- function(ours, ref) {
  vapply(names(ours), function(v) max(abs(ours[[v]] - ref[[v]])) / max(abs(ref[[v]])), 0)
}

test_that("R backend: gp2014 matches the legacy ODE over 1960 (relative error < 1e-6)", {
  ref <- readRDS(test_path("fixtures", "gp2014_mcallen_ode_tight_3yr.rds"))
  wx <- sg_forcing(sg_example_weather(), gp, start = "1960-01-01", end = "1960-12-31")
  m <- sg_model(gp, wild, wx)
  sim <- simulate_ode(m, backend = "R", rtol = 1e-7, atol = 1e-9)
  ref <- ref[seq_len(wx$n_days), ]
  expect_identical(sim$date[-1], as.Date(ref$date))
  err <- max_rel_error(legacy_totals(sim), ref)
  expect_true(all(err < 1e-6), label = paste(names(err), signif(err, 3), collapse = ", "))
})

test_that("C++ backend: gp2014 matches the legacy ODE over 1960-1983 (relative error < 1e-6)", {
  ref <- readRDS(test_path("fixtures", "gp2014_mcallen_ode_tight_24yr_weekly.rds"))
  wx <- sg_forcing(sg_example_weather(), gp)
  m <- sg_model(gp, wild, wx)
  sim <- simulate_ode(m, backend = "cpp", rtol = 1e-8, atol = 1e-10)
  rows <- match(as.Date(ref$date), sim$date[-1])
  expect_false(anyNA(rows))
  err <- max_rel_error(legacy_totals(sim)[rows, ], ref)
  expect_true(all(err < 1e-6), label = paste(names(err), signif(err, 3), collapse = ", "))
  # the explicit integrator never needed the stiff fallback over 24 years (O-006)
  expect_identical(sim$meta$fallback_days, 0L)
})

test_that("R backend: gp2014 matches the legacy ODE over 1960-1983 (slow)", {
  skip_if_not(
    identical(Sys.getenv("SGDRIVE_SLOW_TESTS"), "true"),
    "slow test: set SGDRIVE_SLOW_TESTS=true to run"
  )
  ref_file <- test_path(
    "..", "..", "data-raw", "reference", "gp2014_mcallen_ode_no_floors_tight.rds"
  )
  skip_if_not(file.exists(ref_file), "legacy reference not available (run from the source tree)")
  ref <- readRDS(ref_file)
  wx <- sg_forcing(sg_example_weather(), gp)
  m <- sg_model(gp, wild, wx)
  sim <- simulate_ode(m, backend = "R", rtol = 1e-8, atol = 1e-10)
  expect_identical(nrow(ref), wx$n_days)
  err <- max_rel_error(legacy_totals(sim), ref)
  expect_true(all(err < 1e-6), label = paste(names(err), signif(err, 3), collapse = ", "))
  # O-006: record whether lsoda needed its stiff method over 24 years
  message(sprintf(
    "24-year gp2014 run: %.0f steps/day, stiff method on %d days, largest rate %.1f/day",
    mean(sim$diagnostics$steps), sum(sim$diagnostics$stiff), max(sim$diagnostics$max_rate)
  ))
})
