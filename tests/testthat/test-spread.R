# sg_spread() and the node and front plots (phase 3 WP6): metrics against
# direct calculations, the order of arrival along a chain, sources of data
# (states, totals, a recorded quantity), positions, plot data.

sp_cube <- sg_cube(sg_test_cube_mendelian())

# n cells of 10 km in a row, flies in the first; the last node is a source
sp_model <- function(n = 6, days = 60, d_coef = 5, ...) {
  p <- sg_params(...)
  wx <- sg_forcing(sg_constant_weather(27, days = days), p)
  nd <- data.frame(
    id = paste0("n", seq_len(n)), area = 100, x = (seq_len(n) - 1) * 10,
    role = c(rep("study", n - 1L), "source")
  )
  sg_model(p, sp_cube,
    landscape = sg_landscape(nd, wx),
    movement = sg_movement_diffusion(stats::setNames(rep(100, n), nd$id),
      data.frame(from = nd$id[-n], to = nd$id[-1L], border = 10, dist = 10),
      D = d_coef
    )
  )
}

sp_init <- function(m) {
  sg_init(m, values = c(E = 50, L = 50, P = 50, M = 20, U = 5, F = 20), nodes = "n1")
}

# flies per node [time, node] of a full-state run (replicate r)
sp_direct <- function(sim, r = 1L) {
  st <- if (length(dim(sim$state)) == 3L) sim$state[, , r] else sim$state
  nm <- colnames(st)
  fly <- !.sg_is_host_name(nm)
  node <- .sg_name_node(nm)
  ids <- unique(node)
  vapply(ids, function(id) rowSums(st[, fly & node == id, drop = FALSE]), numeric(nrow(st)))
}

test_that("metrics equal direct calculations on states", {
  m <- sp_model()
  x0 <- sp_init(m)
  for (solver in c("ode", "cb")) {
    sim <- if (solver == "ode") {
      simulate_ode(m, x0, node_extinction = TRUE, record = "state")
    } else {
      simulate_cb(m, x0, n_rep = 2, seed = 5, record = "state")
    }
    sp <- sg_spread(sim, nodes = 1:6)
    expect_s3_class(sp, "sg_spread")
    expect_identical(sp$threshold, 1)
    for (r in if (solver == "cb") 1:2 else 1L) {
      flies <- sp_direct(sim, r)
      occ <- flies >= 1
      nd <- if (solver == "cb") sp$nodes[sp$nodes$replicate == r, ] else sp$nodes
      tt <- if (solver == "cb") sp$time[sp$time$replicate == r, ] else sp$time
      first <- apply(occ, 2L, function(v) if (any(v)) sim$day[which(v)[1L]] else NA)
      expect_identical(as.numeric(nd$first_arrival), as.numeric(first), info = solver)
      expect_identical(nd$days_occupied, as.integer(colSums(occ[-1L, ])))
      expect_identical(tt$occupied, as.integer(rowSums(occ)))
      x <- (0:5) * 10
      expect_identical(tt$front, apply(occ, 1L, function(v) if (any(v)) max(x[v]) else NA_real_))
    }
    expect_identical(sp$nodes$first_arrival[sp$nodes$node == "n1"][1], 0L)
  }
})

test_that("flies reach the nodes of a chain in order", {
  m <- sp_model(n = 8, days = 90, d_coef = 2)
  ode <- simulate_ode(m, sp_init(m), node_extinction = TRUE)
  sp <- sg_spread(ode, nodes = 1:8)
  arr <- sp$nodes$first_arrival
  expect_false(anyNA(arr))
  expect_true(all(diff(arr) >= 0))
  expect_gt(arr[8], arr[2])
  # the front only moves forward while every node stays occupied
  expect_true(all(diff(sp$time$front) >= 0))
})

test_that("states, totals and a recorded per-node quantity give the same metrics", {
  m <- sp_model()
  x0 <- sp_init(m)
  rec <- sg_record("totals", sg_quantity("flies", stage = "all", by = "node", node = 1:6))
  st <- simulate_cb(m, x0, n_rep = 2, seed = 9, record = "state")
  tt <- simulate_cb(m, x0, n_rep = 2, seed = 9, record = rec)
  expect_identical(sg_spread(tt, nodes = 1:6), sg_spread(st, nodes = 1:6))
  expect_identical(sg_spread(tt, nodes = 1:6, quantity = "flies"), sg_spread(st, nodes = 1:6))
  # the default: the study nodes (the source node n6 left out)
  expect_identical(unique(sg_spread(tt)$nodes$node), paste0("n", 1:5))
  only_q <- simulate_cb(m, x0, n_rep = 1, seed = 9, record = sg_record(
    sg_quantity("flies", stage = "all", by = "node", node = 1:6)
  ))
  expect_error(sg_spread(only_q), "states or totals")
  expect_identical(
    sg_spread(only_q, nodes = 1:6, quantity = "flies")$nodes,
    sg_spread(st, nodes = 1:6)$nodes[sg_spread(st, nodes = 1:6)$nodes$replicate == 1, ]
  )
  expect_error(sg_spread(only_q, quantity = "adults"), "No recorded quantity \"adults\"")
})

test_that("thresholds and positions", {
  m <- sp_model(extinction_threshold = 20)
  sim <- simulate_ode(m, sp_init(m), node_extinction = TRUE)
  expect_identical(sg_spread(sim)$threshold, 20)
  low <- sg_spread(sim, threshold = 0.5, nodes = 1:6)
  high <- sg_spread(sim, threshold = 50, nodes = 1:6)
  expect_true(all(high$nodes$first_arrival >= low$nodes$first_arrival, na.rm = TRUE))
  # positions given in node order or by id
  pos <- c(n1 = 0, n2 = 3, n3 = 7, n4 = 20, n5 = 21, n6 = 50)
  a <- sg_spread(sim, nodes = 1:6, position = unname(pos))
  b <- sg_spread(sim, nodes = 1:6, position = rev(pos))
  expect_identical(a$time$front, b$time$front)
  expect_true(all(a$time$front %in% pos))
  expect_error(sg_spread(sim, position = 1:3), "one value per node")
  expect_error(sg_spread(sim, position = c(n1 = 0)), "no value for node")
  expect_error(sg_spread(sim, threshold = 0), "positive number")
  # no coordinates: no front
  p <- sg_params()
  wx <- sg_forcing(sg_constant_weather(27, days = 20), p)
  plain <- sg_model(p, sp_cube,
    landscape = sg_landscape(data.frame(id = c("a", "b"), area = 50), wx),
    movement = sg_movement_layer(
      list(male = 0.1, virgin = 0.1, mated = 0.1), rbind(c(0, 1), c(1, 0))
    )
  )
  ps <- sg_spread(simulate_ode(plain, sg_init(plain, nodes = "a")))
  expect_true(all(is.na(ps$time$front)))
  expect_output(print(ps), "nodes ever occupied: 2 of 2")
})

test_that("plot data for the node heat map and the front", {
  m <- sp_model()
  x0 <- sp_init(m)
  cb <- simulate_cb(m, x0, n_rep = 3, seed = 2)
  hm <- sg_plot_data(cb, "nodes", nodes = 1:6)
  expect_identical(levels(hm$series), paste0("n", 1:6))
  expect_identical(unique(hm$replicate), 1L) # the first replicate by default
  expect_identical(nrow(hm), length(cb$day) * 6L)
  expect_identical(unique(sg_plot_data(cb, "nodes", replicates = 2:3)$replicate), 2:3)
  fr <- sg_plot_data(cb, "front")
  sp <- sg_spread(cb)
  expect_identical(fr$value, sp$time$front)
  expect_identical(as.character(unique(fr$series)), "front")
  ode <- simulate_ode(m, x0)
  expect_equal(
    sg_plot_data(ode, "nodes", per = "node")$value,
    sg_totals(ode, by = "node", per = "node")$value
  )
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(cb, "nodes"), "ggplot")
  expect_s3_class(plot(cb, "front"), "ggplot")
  p <- sg_params()
  wx <- sg_forcing(sg_constant_weather(27, days = 20), p)
  two <- sg_landscape(data.frame(id = c("a", "b"), area = 50), wx)
  plain <- sg_model(p, sp_cube, landscape = two)
  expect_error(plot(simulate_ode(plain), "front"), "needs node positions")
})
