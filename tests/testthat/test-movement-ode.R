# ODE simulations of models with several nodes and adult movement
# (MODEL_SPEC section 8a, D-032; phase 1 WP5).

mo_forcing <- function(p, days = 60, start = "1971-04-01") {
  sg_forcing(sg_example_weather(), p, start = start, end = as.Date(start) + days - 1)
}

mo_chain <- function(n, rate = list(male = 0.05, virgin = 0.1, mated = 0.17), open_end = FALSE,
                     clock = "physiological") {
  dest <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nb <- c(i - 1L, i + 1L)
    nb <- nb[nb >= 1L & nb <= n]
    dest[i, nb] <- 1 / length(nb)
  }
  if (open_end) dest[n, ] <- dest[n, ] * 0.5
  sg_movement_layer(rate, dest, clock = clock)
}

mo_land <- function(p, n = 3, area = 10, days = 60) {
  sg_landscape(data.frame(id = paste0("n", seq_len(n)), area = area), mo_forcing(p, days))
}

# Wild-type state of a one-node model in the first node, empty elsewhere.
mo_state_first <- function(p, cube, n, area = 10, days = 60) {
  one <- sg_model(p, cube, mo_forcing(p, days), area = area)
  x1 <- as.numeric(sg_init(one))
  c(x1, numeric((n - 1L) * length(x1)))
}

test_that("two identical nodes with symmetric movement each equal the one-node run", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  one <- sg_model(p, cube, mo_forcing(p), area = 10)
  x1 <- as.numeric(sg_init(one))
  ref <- simulate_ode(one, init = x1, output_every = 10L)
  m <- sg_model(p, cube, landscape = mo_land(p, 2), movement = mo_chain(2))
  sim <- simulate_ode(m, init = c(x1, x1), output_every = 10L)
  n_s <- length(x1)
  scale <- max(abs(ref$state))
  expect_lt(max(abs(sim$state[, seq_len(n_s)] - ref$state)) / scale, 1e-6)
  expect_lt(max(abs(sim$state[, n_s + seq_len(n_s)] - ref$state)) / scale, 1e-6)
  expect_identical(colnames(sim$state)[1L], "E[AA,1]@n1")
  expect_true(all(c("emigrants_F@n2", "lost_M@n1") %in% colnames(sim$counters)))
})

test_that("released sterile males spread over a closed chain without changing their total", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_sit())
  land <- mo_land(p, 4)
  rel <- sg_releases("S", stage = "adult", sex = "male", times = 1, number = 1000, node = "n1")
  x0 <- numeric(sg_model(p, cube, releases = rel, landscape = land)$n_states)
  still <- simulate_ode(sg_model(p, cube, releases = rel, landscape = land), init = x0)
  moving <- simulate_ode(
    sg_model(p, cube, releases = rel, landscape = land, movement = mo_chain(4)),
    init = x0
  )
  males_s <- function(sim, node = NULL) {
    nm <- colnames(sim$state)
    keep <- grepl("^M\\[S,", nm) & (is.null(node) | grepl(paste0("@", node, "$"), nm))
    rowSums(sim$state[, keep, drop = FALSE])
  }
  # same forcing in every node: movement does not change the network total
  expect_equal(males_s(moving), males_s(still), tolerance = 1e-7)
  # without movement they stay in node 1; with movement they reach node 4
  expect_identical(max(males_s(still, "n4")), 0)
  expect_gt(max(males_s(moving, "n4")), 1)
  # emigrants equal immigrants over the run (closed chain)
  cn <- colnames(moving$counters)
  expect_equal(
    sum(moving$counters[, grepl("^emigrants_M@", cn)]),
    sum(moving$counters[, grepl("^immigrants_M@", cn)]),
    tolerance = 1e-10
  )
  expect_identical(sum(moving$counters[, grepl("^lost_", cn)]), 0)
})

test_that("with an open end the network loses the counted flies", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_sit())
  land <- mo_land(p, 3)
  rel <- sg_releases("S", stage = "adult", sex = "male", times = 1, number = 1000, node = "n3")
  x0 <- numeric(sg_model(p, cube, releases = rel, landscape = land)$n_states)
  still <- simulate_ode(sg_model(p, cube, releases = rel, landscape = land), init = x0)
  open <- simulate_ode(
    sg_model(p, cube, releases = rel, landscape = land, movement = mo_chain(3, open_end = TRUE)),
    init = x0
  )
  nm <- colnames(open$state)
  ms <- function(sim) sum(sim$state[nrow(sim$state), grepl("^M\\[S,", nm)])
  lost <- sum(open$counters[, grepl("^lost_M@", colnames(open$counters))])
  expect_gt(lost, 1)
  expect_lt(ms(open), ms(still))
})

test_that("an empty node is colonised only through movement, from a source node", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  land <- sg_landscape(
    data.frame(id = c("src", "study"), area = 10, role = c("source", "study")),
    mo_forcing(p)
  )
  one_way <- sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17),
    data.frame(from = "src", to = "study", prob = 1)
  )
  x0 <- mo_state_first(p, cube, 2)
  m0 <- sg_model(p, cube, landscape = land)
  m1 <- sg_model(p, cube, landscape = land, movement = one_way)
  s0 <- simulate_ode(m0, init = x0, output_every = 10L)
  s1 <- simulate_ode(m1, init = x0, output_every = 10L)
  study <- grepl("@study$", colnames(s1$state))
  expect_identical(max(abs(s0$state[, study])), 0)
  expect_gt(sum(s1$state[nrow(s1$state), study]), 10)
  # nothing comes back to the source; the study node has no destinations, so
  # its leavers leave the landscape (a row of dest summing to 0 < 1)
  cn <- colnames(s1$counters)
  expect_identical(sum(s1$counters[, grepl("^immigrants_.@src$", cn)]), 0)
  expect_equal(
    sum(s1$counters[, grepl("^emigrants_.@study$", cn)]),
    sum(s1$counters[, grepl("^lost_.@study$", cn)])
  )
})

test_that("releases happen in their node", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_sit())
  rel <- sg_releases("S",
    stage = "adult", sex = "male", times = 5, number = 100, per = "km2",
    node = "n2"
  )
  m <- sg_model(p, cube, releases = rel, landscape = mo_land(p, 3, area = c(1, 2, 3)))
  sim <- simulate_ode(m, init = numeric(m$n_states))
  expect_identical(unname(sim$releases[5L, "M"]), 200)
  nm <- colnames(sim$state)
  day5 <- sim$state[sim$day == 5L, ]
  expect_gt(sum(day5[grepl("^M\\[S,.*@n2$", nm)]), 150)
  expect_identical(sum(day5[grepl("^M\\[S,.*@n[13]$", nm)]), 0)
})

test_that("the C++ and R backends agree on a multi-node run with hosts and movement", {
  skip_if_not_installed("deSolve")
  p <- sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  cube <- sg_cube(sg_test_cube_mendelian())
  land <- mo_land(p, 3, days = 20)
  m <- sg_model(p, cube, landscape = land, movement = mo_chain(3, open_end = TRUE))
  one <- sg_model(p, cube, mo_forcing(p, 20), area = 10)
  x1 <- as.numeric(sg_init(one))
  hosts_only <- replace(numeric(one$n_states), sg_state_index(one, "HS"), one$hosts$n)
  x0 <- c(x1, x1, hosts_only)
  a <- simulate_ode(m, init = x0, output_every = 5L, rtol = 1e-10, atol = 1e-10)
  b <- simulate_ode(m,
    init = x0, backend = "R", method = "ode45", output_every = 5L,
    rtol = 1e-10, atol = 1e-10
  )
  expect_lt(max(abs(a$state - b$state)) / max(abs(a$state)), 1e-6)
  expect_lt(max(abs(a$counters - b$counters)) / max(abs(a$counters)), 1e-6)
})

test_that("several nodes: stiff fallback, memory limits and the limiting-node diagnostic", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  m <- sg_model(p, cube, landscape = mo_land(p, 3, days = 20), movement = mo_chain(3))
  x0 <- mo_state_first(p, cube, 3, days = 20)
  expect_error(
    simulate_ode(m, init = x0, stiff_fallback = TRUE),
    "not available with several nodes"
  )
  expect_error(simulate_ode(m, init = x0, max_memory = 1e4), "more than max_memory")
  skip_if_not_installed("deSolve")
  expect_error(
    simulate_ode(m, init = x0, method = "lsoda", max_memory = 1e5),
    "needs a dense Jacobian"
  )
  sim <- simulate_ode(m, init = x0)
  expect_false(sim$meta$stiff_fallback)
  d <- sim$diagnostics
  expect_true("limiting_node" %in% names(d))
  expect_true(all(d$limiting_node %in% c("n1", "n2", "n3")))
  # one-node runs keep their diagnostics columns and the fallback default
  one <- sg_model(p, cube, mo_forcing(p, 20), area = 10)
  s1 <- simulate_ode(one, init = x0[seq_len(one$n_states)])
  expect_false("limiting_node" %in% names(s1$diagnostics))
  expect_true(s1$meta$stiff_fallback)
  expect_error(simulate_ode(one, stiff_fallback = NA), "TRUE or FALSE \\(or NULL")
  expect_error(simulate_ode(one, max_memory = 0), "single positive number")
  expect_error(simulate_ode(one, max_memory = 100), "more than max_memory")
})

test_that("print summarises nodes and movement", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  m <- sg_model(p, cube, landscape = mo_land(p, 3, days = 20), movement = mo_chain(3))
  sim <- simulate_ode(m, init = mo_state_first(p, cube, 3, days = 20))
  out <- capture.output(print(sim))
  expect_true(any(grepl("nodes: 3, [0-9] with flies at the end", out)))
  expect_true(any(grepl("movers [0-9.e+]+, lost from the landscape 0", out)))
})
