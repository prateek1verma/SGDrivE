# Speed of the chain-binomial solver (milestone M4): seconds per simulated
# year and replicate of simulate_cb() for the scenarios of the M3 benchmark,
# started from a one-year burn-in, at two node areas (population sizes) and
# two time steps; and the speed-up of the C++ step over the R reference step
# (.sg_cb_run_r, 20 days). M4 target: at most about 0.2 s per simulated year
# and replicate for the two-genotype SIT model at dt = 0.1.
#
# Run from the package root:  Rscript bench/m4_cb_speed.R
# Writes bench/results/m4_cb_speed.csv (about two minutes). Run with the
# computer plugged in and awake.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-cubes.R")

n_rep <- 5L
wx_all <- sg_example_weather()
sit <- sg_test_cube_sit(residual = 0)
sit$eta[, "S"] <- 0.35

scenario <- function(name, p, cube, area, releases = NULL, dts = c(0.1, 0.05), ref = TRUE) {
  wx <- sg_forcing(wx_all, p, start = "1971-01-01", end = "1971-12-31")
  m <- sg_model(p, sg_cube(cube), wx, releases, area = area)
  init <- sg_init(m, "burn_in", weather = wx_all)
  rows <- lapply(dts, function(dt) {
    t0 <- proc.time()[["elapsed"]]
    sim <- simulate_cb(m, init, n_rep = n_rep, dt = dt, seed = 1, keep = "totals")
    el <- proc.time()[["elapsed"]] - t0
    r_per_day <- NA_real_
    if (ref && dt == 0.1) {
      m20 <- sg_model(p, sg_cube(cube), sg_forcing(wx_all, p, "1971-05-01", "1971-05-20"),
        area = area
      )
      x20 <- sg_init(m20, values = c(E = 20, L = 20, P = 20, M = 10, F = 10), per = "km2")
      set.seed(1)
      t1 <- proc.time()[["elapsed"]]
      .sg_cb_run_r(m20, x20, dt = dt)
      r_per_day <- (proc.time()[["elapsed"]] - t1) / 20
      set.seed(1)
      t2 <- proc.time()[["elapsed"]]
      .sg_cb_run_cpp(m20, x20, dt = dt)
      cpp_per_day <- (proc.time()[["elapsed"]] - t2) / 20
    }
    tot <- .sg_cb_block_totals(sim)
    data.frame(
      scenario = name, genotypes = m$dims$n_genotypes, states = m$n_states, area = area,
      dt = dt, mean_max_adults = round(mean(apply(tot[, c("M", "U", "F"), , drop = FALSE], 3L,
        function(a) max(rowSums(a))
      ))),
      seconds_per_rep_year = round(el / n_rep / wx$n_days * 365, 3),
      r_reference_seconds_per_day = round(r_per_day, 3),
      speedup_vs_r = if (is.na(r_per_day)) NA else round(r_per_day / cpp_per_day)
    )
  })
  do.call(rbind, rows)
}

weekly <- function(g, n) sg_releases(g, "adult", "male", times = seq(1, 365, 7), number = n, per = "km2")
res <- list(
  scenario("wild, defaults", sg_params(), sg_test_cube_wild(), 1000),
  scenario("wild, defaults", sg_params(), sg_test_cube_wild(), 1e5, dts = 0.1, ref = FALSE),
  scenario("wild, gp2014 preset", sg_params(preset = "gp2014"), sg_test_cube_wild(), 1000),
  scenario("SIT W/S, weekly releases", sg_params(), sit, 1000, weekly("S", 30)),
  scenario("SIT W/S, weekly releases", sg_params(), sit, 1e5, weekly("S", 30), dts = 0.1, ref = FALSE),
  scenario("SIT, single eggs", sg_params(oviposition = "eggs"), sit, 1000, weekly("S", 30),
    dts = 0.1, ref = FALSE
  )
)
if (requireNamespace("MGDrivE", quietly = TRUE)) {
  res <- c(res, list(
    scenario("Mendelian (3)", sg_params(), MGDrivE::cubeMendelian(), 1000),
    scenario("homing drive (10)", sg_params(), MGDrivE::cubeHomingDrive(), 1000,
      weekly("HH", 5),
      dts = 0.1
    )
  ))
}
res <- do.call(rbind, res)
info <- Sys.info()
res$machine <- paste(info[["sysname"]], info[["machine"]])
res$r_version <- R.version.string
res$date <- format(Sys.Date())
print(res[, c(
  "scenario", "states", "area", "dt", "mean_max_adults", "seconds_per_rep_year",
  "speedup_vs_r"
)], row.names = FALSE)
utils::write.csv(res, "bench/results/m4_cb_speed.csv", row.names = FALSE)
