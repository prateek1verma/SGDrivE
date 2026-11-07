# Run time of both solvers against the number of nodes (multi-node work
# phase 3, WP7; design note dev/MOVEMENT_DESIGN.md section 8a.4). Regular
# grids of 1, 10, 100 and 254 cells of 50 km (254 = the Texas counties),
# built with sg_movement_diffusion(); the SIT model (2 genotypes, wild type W
# and sterile S); one year of McAllen weather (1971), shifted by up to
# +-1 degC across the nodes so that every node has its own forcing series
# ("distinct"), and at 254 nodes also one series for all ("shared"). ODE:
# Dormand-Prince in C++ at the default tolerances; chain-binomial: one
# replicate at dt = 0.1. Both keep totals (the default output).
#
# Targets (8a.4): linear growth in the number of nodes (about 18 ns per state
# and rate evaluation for the ODE); 254-node SIT ODE year under about 30 s;
# a chain-binomial target set from this first measurement.
#
# Run from the package root:  Rscript bench/m6_multinode_scaling.R
# Writes bench/results/m6_multinode_scaling.csv (a few minutes). Run with the
# computer plugged in and awake.

devtools::load_all(".", quiet = TRUE)
source("bench/m6_grid.R")

run_case <- function(n, weather) {
  t_build <- proc.time()[["elapsed"]]
  m <- m6_grid_model(n, weather)
  x0 <- m6_init(m)
  t_build <- proc.time()[["elapsed"]] - t_build
  t0 <- proc.time()[["elapsed"]]
  ode <- simulate_ode(m, x0)
  t_ode <- proc.time()[["elapsed"]] - t0
  evals <- sum(ode$diagnostics$evaluations)
  t1 <- proc.time()[["elapsed"]]
  simulate_cb(m, x0, n_rep = 1, seed = 1)
  t_cb <- proc.time()[["elapsed"]] - t1
  row <- data.frame(
    nodes = n, weather = weather, states = m$n_states,
    edges = if (is.null(m$movement)) 0L else length(m$movement$groups[[1L]]$col),
    build_s = t_build, ode_s = t_ode, ode_steps_per_day = mean(ode$diagnostics$steps),
    ode_evaluations = evals, ode_ns_per_state_eval = 1e9 * t_ode / (evals * m$n_states),
    cb_s = t_cb, cb_ms_per_node = 1e3 * t_cb / n
  )
  cat(sprintf(
    "%4d nodes (%s): %6d states; ODE %6.2f s (%.0f steps/day, %.1f ns/state/eval); CB %6.2f s\n",
    n, weather, m$n_states, t_ode, row$ode_steps_per_day, row$ode_ns_per_state_eval, t_cb
  ))
  row
}

cases <- rbind(
  data.frame(n = c(1L, 10L, 100L, 254L), weather = "distinct"),
  data.frame(n = 254L, weather = "shared")
)
res <- do.call(rbind, Map(run_case, cases$n, cases$weather))
res$r_version <- paste(R.version$major, R.version$minor, sep = ".")
res$package_version <- as.character(utils::packageVersion("ScrewDrivE"))
res$date <- format(Sys.Date())
dir.create("bench/results", showWarnings = FALSE)
utils::write.csv(res, "bench/results/m6_multinode_scaling.csv", row.names = FALSE)
print(res[c("nodes", "weather", "states", "ode_s", "ode_ns_per_state_eval", "cb_s", "cb_ms_per_node")])
