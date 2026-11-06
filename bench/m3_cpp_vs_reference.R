# Speed of the C++ backend (milestone M3) against the pure-R reference backend
# (M2) on the scenarios of m2_reference_baseline.R, at the same tolerances
# (rtol 1e-8, atol 1e-10). M3 target: at least 10 times faster.
#
# For each scenario: run time of the R backend (lsoda, one run), the C++
# backend with its own Dormand-Prince day loop and with C++ rates driven by
# lsoda (median of `reps` runs each), the speed-up and the largest relative
# difference between the R and C++ states. The last two scenarios are full
# 24-year runs (1960-1983), timed for the C++ backend only (the R backend
# takes about 8 minutes for gp2014; see the slow regression test).
#
# Run from the package root:  Rscript bench/m3_cpp_vs_reference.R
# Writes bench/results/m3_cpp_vs_reference.csv. Run with the computer plugged
# in and awake: on battery, background processes can be throttled heavily.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-cubes.R")

reps <- 3L
wx_all <- sg_example_weather()
sit <- sg_test_cube_sit(residual = 0)
sit$eta[, "S"] <- 0.35

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, seconds = proc.time()[["elapsed"]] - t0)
}
median_run <- function(m, ...) {
  runs <- lapply(seq_len(reps), function(i) timed(simulate_ode(m, ..., record = "state")))
  secs <- vapply(runs, `[[`, 0, "seconds")
  list(value = runs[[1L]]$value, seconds = stats::median(secs))
}

scenario <- function(name, p, cube, start = NULL, end = NULL, releases = NULL, with_r = TRUE) {
  wx <- sg_forcing(wx_all, p, start = start, end = end)
  m <- sg_model(p, sg_cube(cube), wx, releases)
  dp <- median_run(m, backend = "cpp")
  ls <- median_run(m, backend = "cpp", method = "lsoda")
  r <- if (with_r) timed(simulate_ode(m, backend = "R", record = "state")) else list(value = NULL, seconds = NA)
  years <- wx$n_days / 365.25
  d <- dp$value$diagnostics
  data.frame(
    scenario = name, genotypes = m$dims$n_genotypes, states = m$n_states, days = wx$n_days,
    r_seconds = round(r$seconds, 3), cpp_seconds = round(dp$seconds, 3),
    cpp_lsoda_seconds = round(ls$seconds, 3),
    speedup = round(r$seconds / dp$seconds, 1),
    speedup_cpp_lsoda = round(r$seconds / ls$seconds, 1),
    cpp_seconds_per_year = round(dp$seconds / years, 3),
    cpp_steps_per_day = round(mean(d$steps), 1),
    cpp_us_per_evaluation = round(1e6 * dp$seconds / sum(d$evaluations), 2),
    fallback_days = dp$value$meta$fallback_days,
    max_rel_diff = if (with_r) {
      signif(max(abs(r$value$state - dp$value$state)) / max(abs(r$value$state)), 3)
    } else {
      NA
    },
    rtol = dp$value$meta$rtol, atol = dp$value$meta$atol
  )
}

res <- list(
  scenario("wild, gp2014 preset, 1960", sg_params(preset = "gp2014"), sg_test_cube_wild(),
           "1960-01-01", "1960-12-31"),
  scenario("wild, defaults, 1971", sg_params(), sg_test_cube_wild(), "1971-01-01", "1971-12-31"),
  scenario("SIT W/S, defaults, weekly releases, 1971", sg_params(), sit,
           "1971-01-01", "1971-12-31",
           sg_releases("S", "adult", "male", times = seq(1, 365, by = 7), number = 1000))
)
if (requireNamespace("MGDrivE", quietly = TRUE)) {
  res <- c(res, list(
    scenario("Mendelian (3), defaults, 1971", sg_params(), MGDrivE::cubeMendelian(),
             "1971-01-01", "1971-12-31"),
    scenario("homing drive (10), defaults, 90 days", sg_params(), MGDrivE::cubeHomingDrive(),
             "1971-05-01", "1971-07-29")
  ))
}
res <- c(res, list(
  scenario("wild, gp2014 preset, 1960-1983", sg_params(preset = "gp2014"), sg_test_cube_wild(),
           with_r = FALSE),
  scenario("SIT W/S, defaults, weekly releases, 1960-1983", sg_params(), sit,
           releases = sg_releases("S", "adult", "male", times = seq(1, 8766, by = 7),
                                  number = 1000),
           with_r = FALSE)
))
res <- do.call(rbind, res)
info <- Sys.info()
res$machine <- paste(info[["sysname"]], info[["machine"]])
res$r_version <- R.version.string
res$date <- format(Sys.Date())
print(res[, c("scenario", "states", "days", "r_seconds", "cpp_seconds", "speedup",
              "speedup_cpp_lsoda", "cpp_us_per_evaluation", "max_rel_diff")], row.names = FALSE)
utils::write.csv(res, "bench/results/m3_cpp_vs_reference.csv", row.names = FALSE)
