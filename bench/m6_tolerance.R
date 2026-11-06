# Tolerances of the ODE solver on landscapes (multi-node work phase 3, WP7;
# design note dev/MOVEMENT_DESIGN.md section 8a item 7): accuracy against
# speed for rtol and atol (atol in individuals) on the 100-node grid of
# bench/m6_grid.R, with flies in every node ("all") and with flies starting
# in one corner and spreading ("corner", D = 50 km^2 per day, so most nodes
# hold tiny numbers for a while). Reference: rtol = 1e-11, atol = 1e-13.
# Errors: the largest difference in a node's daily adults, relative to that
# node's largest value, over nodes and days ("node_adults"), and the
# relative difference in eggs laid over the year in the whole landscape.
#
# Run from the package root:  Rscript bench/m6_tolerance.R
# Writes bench/results/m6_tolerance.csv (a few minutes).

devtools::load_all(".", quiet = TRUE)
source("bench/m6_grid.R")

tols <- data.frame(
  rtol = c(1e-6, 1e-6, 1e-7, 1e-8),
  atol = c(1e-6, 1e-8, 1e-9, 1e-10)
)

adults_by_node <- function(sim) {
  t <- sg_totals(sim, by = "node", stages = "adults", per = "node", nodes = seq_len(100))
  matrix(t$value, nrow = length(sim$day))
}

rows <- list()
for (start in c("all", "corner")) {
  m <- m6_grid_model(100, "distinct", d_coef = if (start == "corner") 50 else 0.4)
  x0 <- m6_init(m, start)
  ref <- simulate_ode(m, x0, rtol = 1e-11, atol = 1e-13)
  a_ref <- adults_by_node(ref)
  eggs_ref <- sum(ref$counters[, grepl("^eggs_laid", colnames(ref$counters))])
  for (k in seq_len(nrow(tols))) {
    t0 <- proc.time()[["elapsed"]]
    sim <- simulate_ode(m, x0, rtol = tols$rtol[k], atol = tols$atol[k])
    el <- proc.time()[["elapsed"]] - t0
    a <- adults_by_node(sim)
    node_err <- max(sweep(abs(a - a_ref), 2L, pmax(apply(a_ref, 2L, max), 1e-300), "/"))
    eggs <- sum(sim$counters[, grepl("^eggs_laid", colnames(sim$counters))])
    rows[[length(rows) + 1L]] <- data.frame(
      start = start, rtol = tols$rtol[k], atol = tols$atol[k], seconds = el,
      steps_per_day = mean(sim$diagnostics$steps), node_adults_error = node_err,
      eggs_laid_error = abs(eggs - eggs_ref) / eggs_ref
    )
    cat(sprintf(
      "%-6s rtol %.0e atol %.0e: %6.2f s, %5.1f steps/day, node adults error %.2e, eggs laid error %.2e\n",
      start, tols$rtol[k], tols$atol[k], el, mean(sim$diagnostics$steps), node_err,
      abs(eggs - eggs_ref) / eggs_ref
    ))
  }
}
res <- do.call(rbind, rows)
res$date <- format(Sys.Date())
utils::write.csv(res, "bench/results/m6_tolerance.csv", row.names = FALSE)
print(res)
