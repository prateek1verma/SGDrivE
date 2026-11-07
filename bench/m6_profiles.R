# Profiles for the deferred performance items (multi-node work phase 3, WP7;
# design note dev/MOVEMENT_DESIGN.md section 8a.2 items 4 and 5):
#  1. ODE: time per rate evaluation on the 100-node grid with one forcing
#     series per node or one shared series (would sharing the daily
#     preparation across nodes help?), and with and without movement (the
#     share of the movement pass);
#  2. chain-binomial: time per replicate-year on 100 nodes with and without
#     larval crowding, with distinct and shared forcing (the larval chain is
#     rebuilt per node and step under crowding; without it, nodes on one
#     series share their chains).
#
# Run from the package root:  Rscript bench/m6_profiles.R
# Writes bench/results/m6_profiles.csv (a few minutes).

devtools::load_all(".", quiet = TRUE)
source("bench/m6_grid.R")

eval_time <- function(m) {
  x0 <- as.numeric(m6_init(m))
  ptr <- sg_kernel_new_cpp(.sg_model_view(m))
  y <- c(x0, numeric(length(attr(.sg_model_view(m), "counter_names"))))
  b <- bench::mark(sg_kernel_rates_cpp(ptr, 180L, y), min_iterations = 200, check = FALSE)
  as.numeric(b$median) # seconds
}

rows <- list()
for (case in list(
  list(label = "distinct, movement", weather = "distinct", movement = TRUE),
  list(label = "shared, movement", weather = "shared", movement = TRUE),
  list(label = "distinct, no movement", weather = "distinct", movement = FALSE)
)) {
  m <- m6_grid_model(100, case$weather, movement = case$movement)
  s <- eval_time(m)
  rows[[length(rows) + 1L]] <- data.frame(
    profile = "ode_rate_evaluation", case = case$label, value = s * 1e6, unit = "us",
    per_state_ns = 1e9 * s / m$n_states
  )
  cat(sprintf("ODE evaluation, 100 nodes, %-22s %8.1f us (%.1f ns per state)\n",
    case$label, s * 1e6, 1e9 * s / m$n_states))
}

crowd <- sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.01)
for (weather in c("distinct", "shared")) {
  for (crowding in c(FALSE, TRUE)) {
    m <- m6_grid_model(100, weather, params = if (crowding) crowd else m6_params)
    x0 <- m6_init(m)
    t0 <- proc.time()[["elapsed"]]
    simulate_cb(m, x0, n_rep = 1, seed = 1)
    el <- proc.time()[["elapsed"]] - t0
    label <- sprintf("%s, crowding %s", weather, if (crowding) "on" else "off")
    rows[[length(rows) + 1L]] <- data.frame(
      profile = "cb_replicate_year", case = label, value = el, unit = "s", per_state_ns = NA_real_
    )
    cat(sprintf("CB replicate-year, 100 nodes, %-24s %6.2f s\n", label, el))
  }
}
res <- do.call(rbind, rows)
res$date <- format(Sys.Date())
utils::write.csv(res, "bench/results/m6_profiles.csv", row.names = FALSE)
print(res)
