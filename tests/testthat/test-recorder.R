# The recorder of simulate_ode() (phase 1 WP9, D-032, design note 7.1).

rc_cube <- sg_cube(sg_test_cube_mendelian())

rc_one <- function(hosts = FALSE, days = 60) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-04-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  sg_model(p, rc_cube, wx, area = 5)
}

rc_nodes <- function(hosts = FALSE, days = 60) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-04-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  land <- sg_landscape(
    data.frame(id = c("src", "a", "b"), area = c(10, 20, 30), role = c("source", "study", "study")),
    wx
  )
  layer <- sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17),
    rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 0.6, 0))
  )
  sg_model(p, rc_cube, landscape = land, movement = layer)
}

# Weighted sum of state columns selected by a predicate on the state names.
rc_sum <- function(sim, keep) rowSums(sim$state[, keep(colnames(sim$state)), drop = FALSE])

test_that("the default record keeps states and counters as before", {
  m <- rc_one()
  x0 <- as.numeric(sg_init(m))
  a <- simulate_ode(m, init = x0, output_every = 7L, record = "state")
  b <- simulate_ode(m, init = x0, output_every = 7L, record = "state")
  expect_identical(a$state, b$state)
  expect_identical(a$counters, b$counters)
  expect_null(a$totals)
  expect_null(a$recorded)
  expect_null(a$counter_sums)
  expect_identical(unname(a$final_state), unname(a$state[nrow(a$state), ]))
  expect_identical(names(a$final_state), sg_state_names(m))
  expect_identical(a$meta$record$sets, c("state", "counters"))
  expect_error(simulate_ode(m, record = "everything"), "must be \"state\", \"totals\"")
})

test_that("recorded quantities equal the same sums of the full state", {
  for (hosts in c(FALSE, TRUE)) {
    m <- rc_nodes(hosts)
    x0 <- sg_init(m)
    qs <- list(
      sg_quantity("adults", stage = "adults", by = "node"),
      sg_quantity("females", sex = "female"),
      sg_quantity("wildF", stage = "F", genotype = "AA", node = c("a", "b")),
      sg_quantity("mated", stage = "F", by = c("genotype", "mate"), node = "a"),
      sg_quantity("young", stage = c("M", "U"), substage = 1:2, by = "stage"),
      sg_quantity("larvae", stage = "L", by = c("node", "genotype"), node = c("src", "a", "b"))
    )
    if (hosts) qs <- c(qs, list(sg_quantity("infested", host = "HI", by = "node")))
    rec <- do.call(sg_record, c(list("state"), qs))
    sim <- simulate_ode(m, init = x0, output_every = 10L, record = rec)
    r <- sim$recorded
    at <- function(nm, ids) sub("^.*@", "", nm) %in% ids
    blk <- function(nm, blocks) substr(nm, 1L, 1L) %in% blocks & !grepl("^H", nm)
    study <- c("a", "b")
    expect_equal(r[, "adults@b"], rc_sum(sim, function(nm) blk(nm, c("M", "U", "F")) & at(nm, "b")))
    expect_equal(r[, "females"], rc_sum(sim, function(nm) blk(nm, c("U", "F")) & at(nm, study)))
    expect_equal(r[, "wildF"], rc_sum(sim, function(nm) startsWith(nm, "F[AA|") & at(nm, study)))
    starts_at <- function(prefix, ids) function(nm) startsWith(nm, prefix) & at(nm, ids)
    expect_equal(r[, "mated[Aa,aa]"], rc_sum(sim, starts_at("F[Aa|aa,", "a")))
    expect_equal(r[, "young[M]"], rc_sum(sim, function(nm) {
      grepl("^M\\[[^,]+,[12]\\]", nm) & at(nm, study)
    }))
    expect_equal(r[, "larvae[aa]@src"], rc_sum(sim, starts_at("L[aa,", "src")))
    if (hosts) {
      expect_equal(r[, "infested@a"], rc_sum(sim, function(nm) startsWith(nm, "HI[") & at(nm, "a")))
    }
    # the R backend applies the same sums
    b <- simulate_ode(m,
      init = x0, output_every = 10L, record = rec, backend = "R", method = "ode45"
    )
    expect_lt(max(abs(b$recorded - r)) / max(abs(r)), 1e-6)
  }
})

test_that("user weights as a matrix or a table give the same sums", {
  m <- rc_one()
  x0 <- as.numeric(sg_init(m))
  nm <- sg_state_names(m)
  w <- matrix(0, 2, m$n_states)
  w[1, startsWith(nm, "M[")] <- 1
  w[2, startsWith(nm, "F[")] <- 0.5
  tab <- data.frame(
    quantity = c(rep("males", sum(startsWith(nm, "M["))), rep("half_f", sum(startsWith(nm, "F[")))),
    state = c(nm[startsWith(nm, "M[")], nm[startsWith(nm, "F[")]),
    weight = c(rep(1, sum(startsWith(nm, "M["))), rep(0.5, sum(startsWith(nm, "F["))))
  )
  sim <- simulate_ode(m, init = x0, record = sg_record(
    "state", sg_quantity_weights(w, c("m1", "f1")), sg_quantity_weights(tab)
  ))
  expect_equal(sim$recorded[, "m1"], rc_sum(sim, function(n) startsWith(n, "M[")))
  expect_equal(sim$recorded[, "f1"], 0.5 * rc_sum(sim, function(n) startsWith(n, "F[")))
  expect_equal(sim$recorded[, "males"], sim$recorded[, "m1"])
  expect_equal(sim$recorded[, "half_f"], sim$recorded[, "f1"])
  expect_error(
    simulate_ode(m, init = x0, record = sg_record(sg_quantity_weights(w[, -1], c("a", "b")))),
    "columns but the model has"
  )
  bad <- tab
  bad$state[1] <- "Z[?,1]"
  expect_error(
    simulate_ode(m, init = x0, record = sg_record(sg_quantity_weights(bad))),
    "unknown state"
  )
})

test_that("counter sums equal sums of the full counters, also removals", {
  m <- rc_nodes()
  x0 <- sg_init(m)
  rec <- sg_record(
    "counters",
    sg_counter("cases", "fertile_batches", by = "node"),
    sg_counter("adult_deaths", c("deaths_M", "deaths_U", "deaths_F"), node = c("a", "b")),
    sg_counter("removed", paste0("removed_", c("E", "L", "P", "M", "U", "F")))
  )
  sim <- simulate_ode(m, init = x0, record = rec, node_extinction = TRUE)
  k <- sim$counters
  s <- sim$counter_sums
  expect_equal(s[, "cases@b"], unname(k[, "fertile_batches@b"]), ignore_attr = TRUE)
  expect_equal(
    s[, "adult_deaths"],
    rowSums(k[, grepl("^deaths_[MUF]@(a|b)$", colnames(k))]),
    ignore_attr = TRUE
  )
  expect_equal(s[, "removed"], rowSums(k[, grepl("^removed_", colnames(k))]), ignore_attr = TRUE)
  # without "counters" only the sums are kept
  only <- simulate_ode(m,
    init = x0, record = sg_record("totals", sg_counter("cases", "fertile_batches"))
  )
  expect_null(only$counters)
  cases_ab <- rowSums(k[, grepl("^fertile_batches@(a|b)$", colnames(k))])
  expect_equal(only$counter_sums[, "cases"], cases_ab, ignore_attr = TRUE)
  expect_error(sg_cases(only), "add \"counters\" to sg_record")
  expect_error(sg_flow_balance(only), "add \"counters\" to sg_record")
  expect_error(
    simulate_ode(m, init = x0, record = sg_record(sg_counter("x", "fertile_batch"))),
    "Unknown counter \"fertile_batch\""
  )
})

test_that("totals give the same summaries as the full state", {
  for (hosts in c(FALSE, TRUE)) {
    for (m in list(rc_one(hosts), rc_nodes(hosts))) {
      x0 <- sg_init(m)
      st <- simulate_ode(m, init = x0, output_every = 10L, record = "state")
      tt <- simulate_ode(m, init = x0, output_every = 10L, record = "totals")
      expect_null(tt$state)
      expect_identical(colnames(tt$totals), colnames(.sg_totals_matrix(m)))
      expect_equal(sg_totals(tt), sg_totals(st))
      expect_equal(sg_allele_freq(tt), sg_allele_freq(st))
      f_t <- sg_females(tt)
      f_s <- sg_females(st)
      classes <- c("females", "unmated", "mated", "wild_type")
      expect_equal(f_t[classes], f_s[classes])
      expect_true(all(is.na(f_t$egg_laying)))
      expect_equal(sg_cases(tt, "week"), sg_cases(st, "week"))
      fb_t <- sg_flow_balance(tt)
      expect_lt(max(abs(fb_t$residual)), 1e-6)
      if (hosts) expect_equal(sg_hosts(tt), sg_hosts(st))
      expect_error(sg_totals(tt, by = "mate"), "need the full state")
    }
  }
  # the chain-binomial totals use the same column names
  m <- rc_one()
  expect_identical(
    colnames(simulate_ode(m, record = "totals")$totals),
    colnames(simulate_cb(m, n_rep = 1L, keep = "totals", seed = 1L)$totals)
  )
})

test_that("a run continued from its final state equals one long run", {
  p <- sg_params()
  wx <- sg_example_weather()
  period <- function(from, to) {
    sg_model(p, rc_cube, sg_forcing(wx, p, start = from, end = to), area = 5)
  }
  whole <- period("1971-04-01", "1971-05-30")
  first <- period("1971-04-01", "1971-04-30")
  second <- period("1971-05-01", "1971-05-30")
  x0 <- as.numeric(sg_init(whole))
  adults <- sg_record(sg_quantity("adults", stage = "adults"))
  long <- simulate_ode(whole, init = x0, record = adults)
  a <- simulate_ode(first, init = x0, record = adults)
  b <- simulate_ode(second, init = a$final_state, record = adults)
  expect_equal(b$final_state, long$final_state)
  expect_equal(unname(b$recorded[-1L, "adults"]), unname(long$recorded[32:61, "adults"]))
})

test_that("only quantities: summaries ask for totals, print shows the quantities", {
  m <- rc_one()
  sim <- simulate_ode(m, record = sg_record(sg_quantity("adults", stage = "adults")))
  expect_null(sim$state)
  expect_null(sim$totals)
  expect_error(sg_totals(sim), "record = \"totals\" or \"state\"")
  out <- capture.output(print(sim))
  expect_true(any(grepl("recorded \\(last output day\\): adults", out)))
  expect_false(any(grepl("final totals", out)))
})

test_that("the memory limit uses what is recorded", {
  m <- rc_nodes(days = 60)
  x0 <- sg_init(m)
  expect_error(
    simulate_ode(m, init = x0, max_memory = 2e5, record = "state"),
    "more than max_memory"
  )
  small <- sg_record(
    sg_quantity("adults", stage = "adults", by = "node"), sg_counter("cases", "fertile_batches")
  )
  expect_silent(simulate_ode(m, init = x0, max_memory = 2e5, record = small))
})

test_that("bad record requests give clear errors", {
  m <- rc_nodes()
  expect_error(sg_record(), "at least one")
  expect_error(sg_record("states"), "Unknown record set \"states\"")
  expect_error(sg_record(42), "takes \"state\"")
  expect_error(sg_quantity("a[1]"), "without")
  expect_error(sg_quantity("x", stage = "Q"), "Unknown stage")
  expect_error(sg_quantity("x", stage = "M", host = "HS"), "not both")
  expect_error(sg_quantity("x", stage = "M", by = "mate"), "needs stage = \"F\"")
  expect_error(sg_quantity("x", stage = "M", sex = "female"), "No stage of the chosen")
  expect_error(sg_quantity("x", host = "HQ"), "any of \"HS\"")
  expect_error(sg_counter("x", 3), "counter names")
  expect_error(sg_quantity_weights(matrix(1, 2, 2)), "one distinct name per row")
  expect_error(
    simulate_ode(m, record = sg_record(sg_quantity("x", genotype = "QQ"))),
    "Unknown genotype \"QQ\""
  )
  expect_error(
    simulate_ode(m, record = sg_record(sg_quantity("x", node = "z"))),
    "Unknown node \"z\""
  )
  expect_error(simulate_ode(m, record = sg_record(sg_quantity("x", host = "HI"))), "has no hosts")
  expect_error(
    simulate_ode(m, record = sg_record(sg_quantity("x"), sg_quantity("x"))),
    "distinct; repeated: x"
  )
  rec <- sg_record("totals", sg_quantity("adults", stage = "adults", by = "node"))
  out <- capture.output(print(rec))
  expect_match(out[1L], "sets: totals")
  expect_true(any(grepl("quantity \"adults\" by node", out)))
})
