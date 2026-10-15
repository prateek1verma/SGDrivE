# Time-step study for the chain-binomial solver (milestone M4, D-022): the
# deterministic map of the reference chain-binomial step (every draw replaced
# by its expectation) against the ODE solved at tight tolerance, for several
# step sizes dt over one year of McAllen weather. Chain moves are exact, so
# the remaining error comes from quantities held fixed within a step (wound
# search, male availability, larval crowding) and from arrivals in a chain
# during a step; it is used to choose the default dt.
#
# Run from the package root:  Rscript bench/m4_dt_convergence.R
# Writes bench/results/m4_dt_convergence.csv (about two minutes).

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-cubes.R")

wx_all <- sg_example_weather()
period <- c("1971-01-01", "1971-12-31")
dts <- c(1, 0.5, 0.25, 0.1, 0.05)
sit <- sg_test_cube_sit(residual = 0)
sit$eta[, "S"] <- 0.35

mk <- function(p, cube, releases = NULL) {
  sg_model(p, sg_cube(cube), sg_forcing(wx_all, p, period[1], period[2]), releases)
}
scenarios <- list(
  "wild, defaults" = mk(sg_params(), sg_test_cube_wild()),
  "wild, gp2014 preset" = mk(sg_params(preset = "gp2014"), sg_test_cube_wild()),
  "Mendelian, male releases" = mk(
    sg_params(), sg_test_cube_mendelian(),
    sg_releases("aa", "adult", "male", times = seq(30, 300, by = 7), number = 20)
  ),
  "SIT, weekly releases" = mk(
    sg_params(), sit, sg_releases("S", "adult", "male", times = seq(1, 365, by = 7), number = 200)
  ),
  "larval crowding" = mk(sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.05), sg_test_cube_wild()),
  "instant mating, Allee" = mk(sg_params(mating_mode = "instant", male_density_half_mating = 2), sg_test_cube_mendelian())
)

rel_err <- function(a, b) if (max(abs(b)) > 0) max(abs(a - b)) / max(abs(b)) else max(abs(a))
res <- list()
for (nm in names(scenarios)) {
  m <- scenarios[[nm]]
  ode <- simulate_ode(m, rtol = 1e-10, atol = 1e-12)
  tot_ode <- .sg_sim_block_totals(ode)
  adults_ode <- rowSums(tot_ode[, c("M", "U", "F")])
  eggs_ode <- sum(ode$counters[, grepl("^eggs_laid", colnames(ode$counters))])
  for (dt in dts) {
    t0 <- proc.time()[["elapsed"]]
    cb <- .sg_cb_run_r(m, dt = dt, stochastic = FALSE)
    el <- proc.time()[["elapsed"]] - t0
    tot <- .sg_sim_block_totals(cb)
    res[[length(res) + 1L]] <- data.frame(
      scenario = nm, dt = dt,
      block_error = signif(max(vapply(colnames(tot), function(b) {
        rel_err(tot[, b], tot_ode[, b])
      }, 0)), 3),
      adults_error = signif(rel_err(rowSums(tot[, c("M", "U", "F")]), adults_ode), 3),
      eggs_laid_error = signif(
        sum(cb$counters[, grepl("^eggs_laid", colnames(cb$counters))]) / eggs_ode - 1, 3
      ),
      r_seconds = round(el, 2)
    )
    message(nm, ", dt ", dt, ": done")
  }
}
res <- do.call(rbind, res)
res$date <- format(Sys.Date())
print(res[, c("scenario", "dt", "block_error", "adults_error", "eggs_laid_error")],
  row.names = FALSE
)
utils::write.csv(res, "bench/results/m4_dt_convergence.csv", row.names = FALSE)
