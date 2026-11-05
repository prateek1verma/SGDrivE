# Node extinction rule of the ODE (MODEL_SPEC N5, D-032; phase 1 WP6).

# built once at the top level, so helper functions below can use it
ne_cube <- sg_cube(sg_test_cube_mendelian())

ne_forcing <- function(p, start = "1971-10-01", days = 90) {
  sg_forcing(sg_example_weather(), p, start = start, end = as.Date(start) + days - 1)
}

ne_fly_idx <- function(m, node = NULL) {
  unlist(lapply(c("E", "L", "P", "M", "U", "F"), function(b) sg_state_index(m, b, node = node)))
}

# One node with a few flies (total `total`) in stages that lay no eggs (no
# females), so the node cannot grow above the threshold in a day.
ne_small_state <- function(m, total) {
  x <- numeric(m$n_states)
  idx <- unlist(lapply(c("E", "L", "P", "M"), function(b) sg_state_index(m, b)))
  x[idx] <- total / length(idx)
  x
}

test_that("node_extinction = FALSE leaves runs unchanged", {
  p <- sg_params()
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), ne_forcing(p), area = 5)
  x0 <- as.numeric(sg_init(m))
  a <- simulate_ode(m, init = x0)
  b <- simulate_ode(m, init = x0, node_extinction = FALSE)
  expect_identical(a$state, b$state)
  expect_identical(a$counters, b$counters)
  expect_false(any(grepl("^removed_", colnames(a$counters))))
  expect_false(b$meta$node_extinction)
  expect_error(simulate_ode(m, init = x0, node_extinction = NA), "TRUE or FALSE")
})

test_that("a node below the threshold is set to zero and the removal is counted", {
  p <- sg_params()
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), ne_forcing(p), area = 5)
  x0 <- ne_small_state(m, 0.5)
  free <- simulate_ode(m, init = x0)
  ruled <- simulate_ode(m, init = x0, node_extinction = TRUE)
  expect_true(ruled$meta$node_extinction)
  stages <- c("E", "L", "P", "M", "U", "F")
  expect_identical(colnames(ruled$counters)[ncol(ruled$counters) - 5:0], paste0("removed_", stages))
  # day 1: the removed flies by stage are the flies of the free run at the end of day 1
  end1 <- free$state[free$day == 1L, ]
  by_stage <- vapply(stages, function(b) sum(end1[sg_state_index(m, b)]), 0)
  expect_equal(ruled$counters[1L, paste0("removed_", stages)], by_stage, ignore_attr = TRUE)
  # afterwards the node stays empty (no releases), nothing more is removed
  expect_identical(max(ruled$state[-1L, ]), 0)
  expect_identical(sum(ruled$counters[-1L, paste0("removed_", stages)]), 0)
  # above the threshold nothing happens
  big <- as.numeric(sg_init(m))
  expect_equal(
    simulate_ode(m, init = big, node_extinction = TRUE)$state,
    simulate_ode(m, init = big)$state
  )
})

test_that("with hosts, a removed node's infested hosts return to HS and count as ended", {
  p <- sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), ne_forcing(p), area = 2)
  x0 <- ne_small_state(m, 0.5)
  x0[sg_state_index(m, "HI")] <- 2
  x0[sg_state_index(m, "HW")] <- 10
  x0[sg_state_index(m, "HS")] <- m$hosts$n - 10 - 2 * m$dims$n_I
  free <- simulate_ode(m, init = x0)
  ruled <- simulate_ode(m, init = x0, node_extinction = TRUE)
  hi1 <- sum(free$state[free$day == 1L, sg_state_index(m, "HI")])
  expect_gt(hi1, 0)
  r1 <- ruled$state[ruled$day == 1L, ]
  expect_identical(sum(r1[sg_state_index(m, "HI")]), 0)
  expect_equal(
    r1[sg_state_index(m, "HS")],
    free$state[free$day == 1L, sg_state_index(m, "HS")] + hi1,
    ignore_attr = TRUE
  )
  expect_equal(
    ruled$counters[1L, "infestations_ended"],
    free$counters[1L, "infestations_ended"] + hi1,
    ignore_attr = TRUE
  )
  # hosts are conserved every day
  host_idx <- c(sg_state_index(m, "HS"), sg_state_index(m, "HW"), sg_state_index(m, "HI"))
  expect_equal(rowSums(ruled$state[, host_idx]), rep(m$hosts$n, nrow(ruled$state)))
})

# Two nodes: warm McAllen weather and a constant cold node fed by movement.
ne_two_nodes <- function(p) {
  wx <- ne_forcing(p, days = 120)
  cold <- sg_forcing(sg_constant_weather(12, days = wx$n_days, start = "1971-10-01"), p)
  land <- sg_landscape(data.frame(id = c("warm", "cold"), area = 10), list(wx, cold))
  layer <- sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17), rbind(c(0, 0.2), c(1, 0)),
    clock = "calendar"
  )
  m <- sg_model(p, ne_cube, landscape = land, movement = layer)
  one <- sg_model(p, ne_cube, wx, area = 10)
  x1 <- as.numeric(sg_init(one))
  list(m = m, x0 = c(x1, x1))
}

test_that("a cold node keeps fractions of flies without the rule and is emptied with it", {
  p <- sg_params()
  nn <- ne_two_nodes(p)
  free <- simulate_ode(nn$m, init = nn$x0, output_every = 30L)
  ruled <- simulate_ode(nn$m, init = nn$x0, output_every = 30L, node_extinction = TRUE)
  cold <- ne_fly_idx(nn$m, "cold")
  last <- nrow(free$state)
  expect_gt(sum(free$state[last, cold]), 0)
  expect_identical(sum(ruled$state[last, cold]), 0)
  cn <- colnames(ruled$counters)
  expect_true(all(c("removed_E@cold", "removed_F@warm") %in% cn))
  expect_gt(sum(ruled$counters[, grepl("^removed_.@cold$", cn)]), 0)
  expect_identical(sum(ruled$counters[, grepl("^removed_.@warm$", cn)]), 0)
})

test_that("flies change only by births, deaths, removals and losses", {
  p <- sg_params()
  nn <- ne_two_nodes(p)
  sim <- simulate_ode(nn$m, init = nn$x0, node_extinction = TRUE)
  cn <- colnames(sim$counters)
  tot <- function(pattern) sum(sim$counters[, grepl(pattern, cn)])
  fly <- !grepl("^H", colnames(sim$state))
  change <- sum(sim$state[nrow(sim$state), fly]) - sum(sim$state[1L, fly])
  flows <- tot("^eggs_laid") - tot("^deaths_") - tot("^senescence_") -
    tot("^emergence_deaths") - tot("^removed_") - tot("^lost_")
  expect_equal(change, flows, tolerance = 1e-7)
})

test_that("the C++ day loop and the R deSolve path apply the rule alike", {
  skip_if_not_installed("deSolve")
  p <- sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  nn <- ne_two_nodes(p)
  hosts <- unlist(lapply(c("HS", "HW", "HI"), function(b) sg_state_index(nn$m, b)))
  x0 <- nn$x0
  x0[hosts] <- 0
  x0[sg_state_index(nn$m, "HS")] <- nn$m$hosts$n
  a <- simulate_ode(nn$m, init = x0, node_extinction = TRUE, output_every = 30L)
  b <- simulate_ode(nn$m,
    init = x0, node_extinction = TRUE, output_every = 30L,
    backend = "R", method = "ode45"
  )
  expect_lt(max(abs(a$state - b$state)) / max(a$state), 1e-6)
  removed <- grepl("^removed_", colnames(a$counters))
  expect_equal(sum(a$counters[, removed]), sum(b$counters[, removed]), tolerance = 1e-6)
  expect_gt(sum(a$counters[, removed]), 0)
})
