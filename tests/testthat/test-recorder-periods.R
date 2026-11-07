# Counter sums by period and the "cases" set (phase 3 WP4): period sums
# equal the sums of the daily counters for both solvers and every
# replicate, across C++ chunks and R paths; sg_cases() reads the "cases"
# set; memory estimates count period rows.

rp_cube <- sg_cube(sg_test_cube_mendelian())

rp_model <- function(hosts = FALSE, n_node = 2L, start = "1971-04-10", end = "1971-07-20") {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = end)
  if (n_node == 1L) {
    return(sg_model(p, rp_cube, wx, area = 50))
  }
  sg_model(p, rp_cube,
    landscape = sg_landscape(data.frame(id = c("a", "b"), area = c(50, 30)), wx),
    movement = sg_movement_layer(
      list(male = 0.1, virgin = 0.2, mated = 0.15), rbind(c(0, 1), c(1, 0))
    )
  )
}

rp_init <- function(m) {
  v <- c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10)
  if (m$dims$n_node == 1L) sg_init(m, values = v) else sg_init(m, values = v, nodes = "a")
}

# sums of the daily rows of `k` [day, column] by the periods of `every`
rp_by_period <- function(k, dates, every) {
  idx <- .sg_counter_periods(dates, every)$index
  unname(rowsum(k, idx, reorder = FALSE))
}

test_that("periods are 7-day blocks from the first day, or calendar months and years", {
  dates <- as.Date("1971-04-10") + 0:101
  w <- .sg_counter_periods(dates, "week")
  expect_identical(w$periods$period[1:3], as.Date(c("1971-04-10", "1971-04-17", "1971-04-24")))
  expect_identical(w$periods$days, c(rep(7L, 14), 4L))
  expect_identical(w$index[c(1, 7, 8, 102)], c(1L, 1L, 2L, 15L))
  mo <- .sg_counter_periods(dates, "month")
  firsts <- as.Date(c("1971-04-10", "1971-05-01", "1971-06-01", "1971-07-01"))
  expect_identical(mo$periods$period, firsts)
  expect_identical(mo$periods$days, c(21L, 31L, 30L, 20L))
  yr <- .sg_counter_periods(as.Date("1971-12-30") + 0:3, "year")
  expect_identical(yr$periods$days, c(2L, 2L))
  expect_identical(.sg_counter_periods(dates, "day")$index, seq_along(dates))
})

test_that("ODE counter sums by period equal the sums of the daily counters", {
  m <- rp_model()
  x0 <- rp_init(m)
  full <- simulate_ode(m, x0, record = "state")
  dates <- m$forcing$dates
  eggs_laid <- c("eggs_laid[AA]", "eggs_laid[Aa]", "eggs_laid[aa]")
  terms <- sg_counter("eggs", counter = eggs_laid, by = "node")
  for (every in c("day", "week", "month", "year")) {
    sim <- simulate_ode(m, x0, record = sg_record(terms, "cases", every = every))
    expect_null(sim$counters)
    expect_identical(sim$meta$record$every, every)
    expect_identical(sim$counter_periods, .sg_counter_periods(dates, every)$periods)
    for (col in c("fertile_batches@a", "egg_batches@b")) {
      ref <- rp_by_period(full$counters[, col, drop = FALSE], dates, every)[, 1]
      expect_equal(unname(sim$counter_sums[, col]), ref,
        tolerance = 1e-12, info = paste(every, col)
      )
    }
    eggs_a <- rowSums(full$counters[, grepl("^eggs_laid\\[.*@a$", colnames(full$counters))])
    ref_eggs <- rp_by_period(cbind(eggs_a), dates, every)[, 1]
    expect_equal(unname(sim$counter_sums[, "eggs@a"]), ref_eggs,
      tolerance = 1e-12
    )
  }
  # daily sums of one counter equal that counter
  d <- simulate_ode(m, x0, record = sg_record("cases"))
  expect_identical(
    unname(d$counter_sums[, "egg_batches@a"]), unname(full$counters[, "egg_batches@a"])
  )
})

test_that("periods continue across C++ chunks and in the R paths", {
  m <- rp_model(n_node = 1L)
  x0 <- as.numeric(rp_init(m))
  full <- simulate_ode(m, x0, record = "state")
  rec <- .sg_recorder(sg_record("cases", every = "week"), m, colnames(full$counters))
  ptr <- sg_kernel_new_cpp(.sg_model_view(m))
  n <- m$forcing$n_days
  keep <- rep(TRUE, n)
  args <- rec[c("state", "counters", "wx", "wc", "c_period")]
  # two chunks split in the middle of week 2 (days 1-10, 11-n)
  a <- sg_ode_run_cpp(ptr, x0, 1L, 10L, 1e-8, 1e-10, keep, FALSE, 1e6, recorder = args)
  b <- sg_ode_run_cpp(ptr, a$x, 11L, n, 1e-8, 1e-10, keep, FALSE, 1e6, recorder = args)
  expect_identical(a$sum_periods, 1:2)
  expect_identical(b$sum_periods[1:2], 2:3)
  sums <- matrix(0, nrow(rec$counter_periods), length(rec$names_c))
  sums[a$sum_periods, ] <- sums[a$sum_periods, ] + a$counter_sums
  sums[b$sum_periods, ] <- sums[b$sum_periods, , drop = FALSE] + b$counter_sums
  ref <- rp_by_period(full$counters[, c("fertile_batches", "egg_batches")], m$forcing$dates, "week")
  expect_equal(sums, ref, tolerance = 1e-12)
  # the R backend (deSolve path) adds days into period rows in R
  skip_if_not_installed("deSolve")
  r_full <- simulate_ode(m, x0, backend = "R", record = "state")
  r_week <- simulate_ode(m, x0, backend = "R", record = sg_record("cases", every = "week"))
  expect_equal(unname(r_week$counter_sums),
    rp_by_period(r_full$counters[, c("fertile_batches", "egg_batches")], m$forcing$dates, "week"),
    tolerance = 1e-12
  )
})

test_that("chain-binomial counter sums by period are exact in every replicate", {
  for (hosts in c(FALSE, TRUE)) {
    m <- rp_model(hosts)
    x0 <- rp_init(m)
    full <- simulate_cb(m, x0, n_rep = 3, seed = 2, record = "state")
    for (every in c("week", "month")) {
      rec <- sg_record("state", "cases", every = every)
      sim <- simulate_cb(m, x0, n_rep = 3, seed = 2, record = rec)
      expect_identical(sim$state, full$state) # recording draws no random numbers
      expect_identical(sim$counter_periods, .sg_counter_periods(m$forcing$dates, every)$periods)
      cols <- dimnames(sim$counter_sums)[[2L]]
      for (r in 1:3) {
        expect_identical(unname(sim$counter_sums[, , r]),
          rp_by_period(full$counters[, cols, r], m$forcing$dates, every),
          info = paste(hosts, every, r)
        )
      }
    }
  }
})

test_that("an early stop leaves zero periods; a partial last period is summed", {
  m <- rp_model(n_node = 1L, end = "1971-06-08") # 60 days: 8 weeks and 4 days
  x0 <- numeric(m$n_states)
  x0[sg_state_index(m, "M", genotype = "AA", substage = 10)] <- 2
  sim <- simulate_cb(m, x0, n_rep = 2, seed = 1, record = sg_record("cases", every = "week"))
  expect_identical(dim(sim$counter_sums), c(9L, 2L, 2L))
  expect_identical(sim$counter_periods$days, c(rep(7L, 8), 4L))
  expect_identical(max(sim$counter_sums), 0) # males only: no eggs at all
})

test_that("the cases set names its sums like the counters and drops full counters", {
  m1 <- rp_model(n_node = 1L)
  r1 <- .sg_recorder(sg_record("cases"), m1, .sg_model_counter_names(m1))
  expect_identical(r1$names_c, c("fertile_batches", "egg_batches"))
  expect_false(r1$counters)
  mh <- rp_model(hosts = TRUE)
  rh <- .sg_recorder(sg_record("cases", "counters"), mh, .sg_model_counter_names(mh))
  host_cases <- c("infestations", "egg_batches", "superinfestations")
  expect_identical(rh$names_c, paste0(rep(host_cases, each = 2), "@", c("a", "b")))
  expect_true(rh$counters)
  # the cases set covers every node, also source nodes
  wx <- sg_forcing(sg_constant_weather(25, days = 10), sg_params())
  land <- sg_landscape(data.frame(id = c("s", "t"), area = 10, role = c("source", "study")), wx)
  ms <- sg_model(sg_params(), rp_cube, landscape = land)
  rs <- .sg_recorder(sg_record("cases"), ms, .sg_model_counter_names(ms))
  expect_true("fertile_batches@s" %in% rs$names_c)
  expect_error(sg_record("cases", every = "fortnight"), "must be one of")
  expect_output(
    print(sg_record("totals", "cases", every = "week")), "sets: totals, cases; counter sums by week"
  )
})

test_that("sg_cases() reads the cases set as it reads the counters", {
  for (hosts in c(FALSE, TRUE)) {
    for (n_node in 1:2) {
      m <- rp_model(hosts, n_node)
      x0 <- rp_init(m)
      info <- paste(hosts, n_node)
      full <- simulate_ode(m, x0, record = "state")
      daily <- simulate_ode(m, x0, record = sg_record("state", "cases"))
      weekly <- simulate_ode(m, x0, record = sg_record("state", "cases", every = "week"))
      for (period in c("day", "week", "month")) {
        expect_equal(sg_cases(daily, period), sg_cases(full, period),
          tolerance = 1e-12, info = info
        )
      }
      expect_equal(sg_cases(weekly, "week"), sg_cases(full, "week"), tolerance = 1e-12, info = info)
      expect_error(sg_cases(weekly, "month"), "recorded by week")
      expect_error(sg_cases(weekly, "day"), "needs daily sums \\(")
      cb_full <- simulate_cb(m, x0, n_rep = 2, seed = 4, record = "state")
      weekly_rec <- sg_record("state", "cases", every = "week")
      cb_week <- simulate_cb(m, x0, n_rep = 2, seed = 4, record = weekly_rec)
      expect_identical(sg_cases(cb_week, "week"), sg_cases(cb_full, "week"), info = info)
    }
  }
  m <- rp_model()
  no_cases <- sg_record("totals", sg_counter("x", counter = "egg_batches"))
  only_q <- simulate_ode(m, rp_init(m), record = no_cases)
  expect_error(sg_cases(only_q), "\"cases\" set")
})

test_that("memory estimates count period rows", {
  m <- rp_model(end = "1971-05-09") # 30 days
  rec <- sg_record("cases", every = "week")
  # 5 periods x 4 sums, 30 days x 8 release columns, the final state; 2 replicates
  bytes <- 8 * 2 * (5 * 4 + 30 * 8 + m$n_states)
  expect_no_error(simulate_cb(m, rp_init(m), n_rep = 2, seed = 1, record = rec, max_memory = bytes))
  expect_error(
    simulate_cb(m, rp_init(m), n_rep = 2, seed = 1, record = rec, max_memory = bytes - 1),
    "max_memory"
  )
  ode_bytes <- 8 * (5 * 4 + m$n_states)
  expect_no_error(simulate_ode(m, rp_init(m), record = rec, max_memory = ode_bytes))
  expect_error(simulate_ode(m, rp_init(m), record = rec, max_memory = ode_bytes - 1), "max_memory")
})
