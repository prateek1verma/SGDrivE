# Speed of the chain-binomial solver with several nodes (multi-node work
# phase 2, WP6; dev/MOVEMENT_PHASE2_PLAN.md): seconds per simulated year and
# replicate of simulate_cb() for 1, 10 and 50 nodes in a row, with and
# without adult movement, starting with flies in every node or in the first
# node only (empty nodes skip their fly draws until flies arrive). Expected:
# time grows about linearly with the number of occupied nodes.
#
# Run from the package root:  Rscript bench/phase2_cb_nodes.R
# Writes bench/results/phase2_cb_nodes.csv (about two minutes). Run with the
# computer plugged in and awake.

devtools::load_all(".", quiet = TRUE)

n_rep <- 3L
p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, start = "1971-03-01", end = "1972-02-29")
cube <- sg_cube(MGDrivE::cubeMendelian())

# n nodes in a row, 50 km^2 each; flies leave to their neighbours (closed ends)
row_model <- function(n, movement) {
  land <- sg_landscape(data.frame(id = paste0("n", seq_len(n)), area = 50), wx)
  if (n == 1L || !movement) {
    return(sg_model(p, cube, landscape = land))
  }
  dest <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nb <- c(i - 1L, i + 1L)
    nb <- nb[nb >= 1L & nb <= n]
    dest[i, nb] <- 1 / length(nb)
  }
  sg_model(p, cube,
    landscape = land,
    movement = sg_movement_layer(list(male = 0.05, virgin = 0.1, mated = 0.08), dest)
  )
}

rows <- list()
for (n in c(1L, 10L, 50L)) {
  for (movement in c(FALSE, TRUE)) {
    if (n == 1L && movement) next
    m <- row_model(n, movement)
    for (start in c("all", "first")) {
      if (n == 1L && start == "first") next
      nodes <- if (start == "all") NULL else "n1"
      x0 <- sg_init(m,
        values = c(E = 200, L = 200, P = 200, M = 100, U = 20, F = 100), nodes = nodes
      )
      t0 <- proc.time()[["elapsed"]]
      sim <- simulate_cb(m, x0, n_rep = n_rep, seed = 1, record = "totals", output_every = 30)
      el <- (proc.time()[["elapsed"]] - t0) / n_rep
      fly <- !grepl("^H", rownames(sim$final_state))
      nm <- rownames(sim$final_state)[fly]
      node <- if (n == 1L) rep("n1", length(nm)) else sub("^[^@]*@", "", nm)
      occupied <- mean(apply(sim$final_state[fly, , drop = FALSE], 2L, function(v) {
        sum(tapply(v, node, sum) > 0)
      }))
      rows[[length(rows) + 1L]] <- data.frame(
        nodes = n, movement = movement, start = start, states = m$n_states,
        seconds_per_year = el, seconds_per_node_year = el / n,
        occupied_at_end = occupied
      )
      cat(sprintf(
        "%2d nodes, movement %-5s, start %-5s: %.3f s per year (%.4f per node), %.1f occupied\n",
        n, movement, start, el, el / n, occupied
      ))
    }
  }
}
res <- do.call(rbind, rows)
res$r_version <- paste(R.version$major, R.version$minor, sep = ".")
res$package_version <- as.character(utils::packageVersion("ScrewDrivE"))
dir.create("bench/results", showWarnings = FALSE)
utils::write.csv(res, "bench/results/phase2_cb_nodes.csv", row.names = FALSE)
print(res)
