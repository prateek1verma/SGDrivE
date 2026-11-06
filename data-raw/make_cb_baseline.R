# Single-node chain-binomial baseline for the multi-node work, phase 2
# (dev/MOVEMENT_PHASE2_PLAN.md, WP1): seeded stochastic runs, deterministic-map
# runs and single steps of one-node models, computed with the code before
# phase 2. tests/testthat/test-cb-single-node-baseline.R checks that every
# phase 2 change reproduces them.
#
# Run from the package root, on the commit before phase 2 code (branch
# movement-phase2a at 3220f94 plus the plan document):
#   Rscript data-raw/make_cb_baseline.R
# Do not regenerate it after phase 2 code has changed, unless a deliberate
# model change is recorded in DECISIONS.md.

stopifnot(file.exists("DESCRIPTION"))
devtools::load_all(quiet = TRUE)

cb_baseline_cases <- function() {
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

init_of <- function(m) sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10))

baseline <- lapply(cb_baseline_cases(), function(m) {
  x0 <- init_of(m)
  sto <- simulate_cb(m, init = x0, n_rep = 2L, dt = 0.1, output_every = 10L, seed = 11L, record = "state")
  det <- .sg_cb_run_cpp(m, init = x0, dt = 0.1, stochastic = FALSE)
  view <- .sg_model_view(m)
  xs <- unname(sto$state[3L, , 1L]) # a whole-number state from the run
  set.seed(5L)
  step_sto <- sg_cb_step_cpp(view, 30L, 0.1, xs, TRUE)
  step_det <- sg_cb_step_cpp(view, 30L, 0.1, xs, FALSE)
  list(
    hash = m$meta$hash,
    stochastic = list(
      state = unname(sto$state), counters = unname(sto$counters),
      extinct_day = unname(sto$extinct_day), wild_extinct_day = unname(sto$wild_extinct_day)
    ),
    deterministic = list(
      state = unname(det$state[seq(1L, nrow(det$state), by = 10L), ]),
      counters = unname(det$counters)
    ),
    step = list(x = xs, stochastic = step_sto, deterministic = step_det)
  )
})

out <- file.path("tests", "testthat", "fixtures", "cb_single_node_baseline.rds")
saveRDS(baseline, out, version = 3)
cat("wrote", out, format(file.size(out), big.mark = ","), "bytes\n")
