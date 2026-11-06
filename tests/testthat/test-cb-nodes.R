# Chain-binomial steps of models with several nodes, without movement yet
# (phase 2 WP2; MODEL_SPEC section 10): the C++ step against the R reference
# draw for draw, each node against its own one-node model, empty nodes.

cbn_ns <- asNamespace("ScrewDrivE")
cbn_cube <- sg_cube(sg_test_cube_mendelian())

cbn_params <- function(hosts) {
  if (hosts) {
    sg_params(
      wound_model = "hosts", host_density = 50, wounded_fraction = 0.1,
      larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
      male_density_half_mating = 2
    )
  } else {
    sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.01)
  }
}

# three nodes: two on McAllen weather (one forcing series), one at 22 degC
cbn_model <- function(p, cube = cbn_cube, days = 12) {
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  mild <- sg_forcing(sg_constant_weather(22, days = days, start = start), p)
  nd <- data.frame(id = c("a", "b", "c"), area = c(5, 8, 10))
  if (p$wound_model == "hosts") {
    nd$host_density <- c(NA, 30, 80)
  } else {
    nd$wound_density <- c(NA, 8, 3)
  }
  sg_model(p, cube, landscape = sg_landscape(nd, list(wx, wx, mild)))
}

# A start with flies in nodes a and c; node b has none (only its hosts).
cbn_init <- function(m) {
  sg_init(m, values = c(E = 8, L = 8, P = 8, M = 4, U = 1, F = 4), nodes = c("a", "c"))
}

test_that("C++ runs of several nodes equal the R reference draw for draw", {
  cubes <- list(mendelian = cbn_cube)
  if (requireNamespace("MGDrivE", quietly = TRUE)) {
    cubes$homing <- sg_cube(MGDrivE::cubeHomingDrive())
  }
  for (hosts in c(FALSE, TRUE)) {
    for (cn in names(cubes)) {
      m <- cbn_model(cbn_params(hosts), cubes[[cn]], days = 6)
      x0 <- cbn_init(m)
      info <- paste(hosts, cn)
      set.seed(42)
      r <- cbn_ns$.sg_cb_run_r(m, init = x0, dt = 0.25, stochastic = TRUE)
      set.seed(42)
      cc <- cbn_ns$.sg_cb_run_cpp(m, init = x0, dt = 0.25, stochastic = TRUE)
      expect_identical(cc$state, r$state, info = info)
      expect_identical(cc$counters, r$counters, info = info)
      # node b had no flies at the start: they stay absent without movement
      fly_b <- !grepl("^H", colnames(r$state)) & endsWith(colnames(r$state), "@b")
      expect_identical(max(r$state[, fly_b]), 0, info = info)
      # deterministic maps agree to round-off
      dr <- cbn_ns$.sg_cb_run_r(m, init = x0, dt = 0.25, stochastic = FALSE)
      dc <- cbn_ns$.sg_cb_run_cpp(m, init = x0, dt = 0.25, stochastic = FALSE)
      expect_equal(dc$state, dr$state, tolerance = 1e-12, info = info)
      expect_equal(dc$counters, dr$counters, tolerance = 1e-12, info = info)
    }
  }
})

test_that("in the deterministic map each node runs as its own one-node model", {
  for (hosts in c(FALSE, TRUE)) {
    p <- cbn_params(hosts)
    m <- cbn_model(p)
    x0 <- as.numeric(cbn_init(m))
    det <- cbn_ns$.sg_cb_run_cpp(m, init = x0, dt = 0.1, stochastic = FALSE)
    n_s <- m$n_states_node
    n_c <- ncol(det$counters) / 3L
    for (i in 1:3) {
      one <- cbn_ns$.sg_node_model(m, i)
      at <- (i - 1L) * n_s + seq_len(n_s)
      ref <- cbn_ns$.sg_cb_run_cpp(one, init = x0[at], dt = 0.1, stochastic = FALSE)
      expect_equal(unname(det$state[, at]), unname(ref$state), tolerance = 1e-13)
      expect_equal(unname(det$counters[, (i - 1L) * n_c + seq_len(n_c)]), unname(ref$counters),
        tolerance = 1e-13
      )
    }
  }
})

test_that("a one-node landscape runs like the plain one-node model", {
  p <- cbn_params(TRUE)
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + 9)
  plain <- sg_model(p, cbn_cube, wx, area = 5)
  land <- sg_model(p, cbn_cube, landscape = sg_landscape(data.frame(id = "only", area = 5), wx))
  x0 <- as.numeric(sg_init(plain, values = c(E = 8, L = 8, P = 8, M = 4, U = 1, F = 4)))
  set.seed(7)
  a <- cbn_ns$.sg_cb_run_cpp(plain, init = x0, dt = 0.1, stochastic = TRUE)
  set.seed(7)
  b <- cbn_ns$.sg_cb_run_cpp(land, init = x0, dt = 0.1, stochastic = TRUE)
  expect_identical(unname(b$state), unname(a$state))
  expect_identical(unname(b$counters), unname(a$counters))
})

test_that("extinction days refer to the whole landscape", {
  p <- sg_params()
  m <- cbn_model(p, days = 30)
  x0 <- numeric(m$n_states)
  # a few old adults in node a only: they die, and the landscape is extinct
  x0[sg_state_index(m, "M", substage = 10, node = "a")] <- 1
  set.seed(1)
  r <- cbn_ns$.sg_cb_run_cpp(m, init = x0, dt = 0.5, stochastic = TRUE)
  expect_gt(r$extinct_day, 0L)
  expect_identical(max(r$state[r$extinct_day + 1L, ]), 0)
})
