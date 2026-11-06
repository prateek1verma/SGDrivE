# Runtime baseline of the pure-R reference solver (milestone M2), against which
# the C++ backend of milestone M3 is measured (target: at least 10 times faster
# at the same tolerances).
#
# Run from the package root:  Rscript bench/m2_reference_baseline.R
# Writes bench/results/m2_reference_baseline.csv. Run with the computer plugged
# in and awake: on battery, background processes can be throttled heavily.

devtools::load_all(".", quiet = TRUE)

wild <- list(
  genotypesID = "W", wildType = "W", ih = array(1, c(1, 1, 1)), tau = array(1, c(1, 1, 1)),
  eta = matrix(1, 1, 1), phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1),
  omega = c(W = 1)
)
sit <- local({
  g <- c("W", "S")
  ih <- array(0, c(2, 2, 2), dimnames = list(g, g, g))
  ih[, , "W"] <- 1
  tau <- array(1, c(2, 2, 2), dimnames = list(g, g, g))
  tau[, "S", ] <- 0
  list(
    genotypesID = g, wildType = "W", releaseType = "S", ih = ih, tau = tau,
    eta = matrix(c(1, 1, 0.35, 0.35), 2, dimnames = list(g, g)),
    phi = c(W = .5, S = .5), xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1),
    s = c(W = 1, S = 1), omega = c(W = 1, S = 1)
  )
})
wx_all <- sg_example_weather()

scenario <- function(name, p, cube, start, end, releases = NULL) {
  wx <- sg_forcing(wx_all, p, start = start, end = end)
  m <- sg_model(p, sg_cube(cube), wx, releases)
  t0 <- proc.time()
  sim <- simulate_ode(m, record = "state")
  el <- (proc.time() - t0)[["elapsed"]]
  data.frame(
    scenario = name, genotypes = m$dims$n_genotypes, states = m$n_states,
    days = wx$n_days, seconds = round(el, 2),
    seconds_per_year = round(el / wx$n_days * 365.25, 2),
    steps_per_day = round(mean(sim$diagnostics$steps), 1),
    evaluations = sum(sim$diagnostics$evaluations),
    us_per_evaluation = round(1e6 * el / sum(sim$diagnostics$evaluations), 1),
    rtol = sim$meta$rtol, atol = sim$meta$atol
  )
}

res <- list(
  scenario("wild, gp2014 preset, 1960", sg_params(preset = "gp2014"), wild,
           "1960-01-01", "1960-12-31"),
  scenario("wild, defaults, 1971", sg_params(), wild, "1971-01-01", "1971-12-31"),
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
res <- do.call(rbind, res)
info <- Sys.info()
res$machine <- paste(info[["sysname"]], info[["machine"]])
res$r_version <- R.version.string
res$date <- format(Sys.Date())
print(res[, c("scenario", "states", "days", "seconds", "seconds_per_year", "steps_per_day",
              "us_per_evaluation")], row.names = FALSE)
utils::write.csv(res, "bench/results/m2_reference_baseline.csv", row.names = FALSE)
