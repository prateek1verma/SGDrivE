# Single-node baseline for the multi-node work (dev/MOVEMENT_PHASE1_PLAN.md,
# section 6): states, counters, names, hashes and kernel rates of one-node
# models computed with the code *before* phase 1. The test
# tests/testthat/test-single-node-baseline.R checks that every phase 1 change
# reproduces them exactly (identical(), i.e. bitwise).
#
# Run from the package root, on the commit before phase 1 code (branch
# movement-phase1a at f4a7e10):
#   Rscript data-raw/make_single_node_baseline.R
# Do not regenerate it after phase 1 code has changed, unless a deliberate
# model change is recorded in DECISIONS.md.

stopifnot(file.exists("DESCRIPTION"))
devtools::load_all(quiet = TRUE)

weather <- sg_example_weather()
cases <- list(
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

baseline <- lapply(cases, function(cs) {
  wx <- sg_forcing(weather, cs$params, start = "1971-01-01", end = "1971-12-31")
  mod <- sg_model(cs$params, cs$cube, wx, cs$releases, area = 2)
  sim <- simulate_ode(mod, output_every = 30L)
  view <- .sg_model_view(mod)
  # kernel rates at three states of the run (start, middle, end)
  probe_rows <- c(1L, ceiling(nrow(sim$state) / 2), nrow(sim$state))
  probe <- lapply(probe_rows, function(r) {
    day <- max(1L, sim$day[r])
    list(
      day = day, state = unname(sim$state[r, ]),
      rates = sg_rates_cpp(view, day, unname(sim$state[r, ])),
      max_rate = sg_max_rate_cpp(view, day, unname(sim$state[r, ]))
    )
  })
  list(
    hash = mod$meta$hash,
    state_names = sg_state_names(mod),
    layout = mod$layout,
    counts = mod$counts,
    releases = mod$releases,
    day = sim$day,
    state = unname(sim$state),
    counters = sim$counters,
    max_rate = sim$diagnostics$max_rate,
    probe = probe
  )
})

out <- file.path("tests", "testthat", "fixtures", "single_node_baseline.rds")
saveRDS(baseline, out, version = 3)
cat("wrote", out, format(file.size(out), big.mark = ","), "bytes\n")
