# The recorder in simulate_cb() (phase 2 WP5; D-032): with the same seed the
# recorded quantities, totals and counter sums equal the sums of a full-state
# run per replicate (recording uses no random numbers), final states,
# keep as a shortcut for record, memory limit, print().

cbr_cube <- sg_cube(sg_test_cube_mendelian())

cbr_model <- function(n_node = 3L, hosts = FALSE, days = 25) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  if (n_node == 1L) {
    return(sg_model(p, cbr_cube, wx, area = 40))
  }
  nd <- data.frame(
    id = c("a", "b", "c"), area = c(40, 60, 80), role = c("study", "study", "source")
  )
  lay <- sg_movement_layer(
    list(male = 0.1, virgin = 0.2, mated = 0.15),
    rbind(c(0, 1, 0), c(0.5, 0, 0.3), c(0, 1, 0))
  )
  sg_model(p, cbr_cube, landscape = sg_landscape(nd, wx), movement = lay)
}

cbr_init <- function(m) {
  v <- c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10)
  if (m$dims$n_node == 1L) sg_init(m, values = v) else sg_init(m, values = v, nodes = "a")
}

cbr_record <- function(m) {
  multi <- m$dims$n_node > 1L
  sg_record(
    "totals",
    sg_quantity("adults", stage = "adults", by = if (multi) "node" else character(0)),
    sg_quantity("wild_mated", stage = "F", genotype = "AA"),
    sg_quantity_weights(matrix(seq_len(m$n_states) / 7, 1L), names = "weighted"),
    sg_counter("cases", counter = "fertile_batches", by = if (multi) "node" else character(0)),
    sg_counter("deaths", counter = c("deaths_M", "deaths_F"))
  )
}

test_that("recorded sums equal the sums of a full-state run with the same seed", {
  for (n_node in c(1L, 3L)) {
    for (hosts in c(FALSE, TRUE)) {
      m <- cbr_model(n_node, hosts)
      x0 <- cbr_init(m)
      info <- paste(n_node, hosts)
      full <- simulate_cb(m, x0, n_rep = 3, seed = 6, output_every = 4)
      rr <- cbr_record(m)
      rec <- simulate_cb(m, x0, n_rep = 3, seed = 6, output_every = 4, record = rr)
      cn <- dimnames(full$counters)[[2L]]
      r_full <- .sg_recorder(rr, m, cn)
      expect_null(rec$state)
      expect_null(rec$counters) # only the counter sums (no "counters" in the record)
      expect_identical(rec$extinct_day, full$extinct_day, info = info)
      for (r in 1:3) {
        x <- full$state[, , r]
        # totals: the sums by block and genotype, exactly (whole numbers)
        expect_identical(unname(rec$totals[, , r]), unname(x %*% .sg_totals_matrix(m)), info = info)
        wx <- vapply(
          seq_len(nrow(x)), function(t) .sg_csr_multiply(r_full$wx, x[t, ]),
          numeric(length(r_full$names_x))
        )
        q <- !seq_along(r_full$names_x) %in% seq_len(r_full$n_totals)
        expect_equal(unname(rec$recorded[, , r]), unname(t(wx[q, , drop = FALSE])),
          tolerance = 1e-14, info = info
        )
        cs <- vapply(seq_len(dim(full$counters)[1L]), function(d) {
          .sg_csr_multiply(r_full$wc, full$counters[d, , r])
        }, numeric(length(r_full$names_c)))
        cs <- t(matrix(cs, ncol = dim(full$counters)[1L]))
        expect_identical(unname(rec$counter_sums[, , r]), unname(cs), info = info)
        expect_identical(unname(rec$final_state[, r]), unname(full$state[dim(full$state)[1L], , r]))
      }
      expect_identical(dimnames(rec$recorded)[[2L]], r_full$names_x[q])
      expect_identical(dimnames(rec$counter_sums)[[2L]], r_full$names_c)
      expect_identical(rownames(rec$final_state), sg_state_names(m))
      expect_identical(rec$meta$record$counter_sums, r_full$names_c)
    }
  }
})

test_that("keep is a shortcut for record; both together are an error", {
  m <- cbr_model()
  x0 <- cbr_init(m)
  a <- simulate_cb(m, x0, n_rep = 2, seed = 2, keep = "totals")
  b <- simulate_cb(m, x0, n_rep = 2, seed = 2, record = "totals")
  st <- simulate_cb(m, x0, n_rep = 2, seed = 2)
  expect_identical(a$totals, b$totals)
  expect_identical(a$counters, b$counters)
  expect_null(a$state)
  # the former keep = "totals": the full states summed by block and genotype
  agg <- .sg_totals_matrix(m)
  for (r in 1:2) expect_identical(unname(a$totals[, , r]), unname(st$state[, , r] %*% agg))
  expect_identical(simulate_cb(m, x0, n_rep = 2, seed = 2, keep = "state")$state, st$state)
  # record = "state" keeps the states and all counters, as by default
  expect_identical(simulate_cb(m, x0, n_rep = 2, seed = 2, record = "state")$state, st$state)
  expect_error(simulate_cb(m, x0, keep = "totals", record = "totals"), "not both")
  expect_error(simulate_cb(m, x0, keep = "all"), "must be one of")
  expect_error(simulate_cb(m, x0, record = "everything"), "sg_record")
})

test_that("states and counters together with recorded quantities", {
  m <- cbr_model()
  x0 <- cbr_init(m)
  st <- simulate_cb(m, x0, n_rep = 2, seed = 5)
  both <- simulate_cb(m, x0, n_rep = 2, seed = 5, record = sg_record(
    "state", "counters", sg_counter("cases", counter = "fertile_batches")
  ))
  expect_identical(both$state, st$state)
  expect_identical(both$counters, st$counters)
  study <- paste0("fertile_batches@", c("a", "b"))
  cases <- apply(st$counters[, study, , drop = FALSE], c(1L, 3L), sum)
  expect_identical(unname(both$counter_sums[, "cases", ]), unname(cases))
})

test_that("an early stop leaves zero rows and a zero final state", {
  m <- cbr_model(1L, days = 40)
  x0 <- numeric(m$n_states)
  x0[sg_state_index(m, "M", genotype = "AA", substage = 10)] <- 2
  rec <- simulate_cb(m, x0, n_rep = 2, seed = 1, record = cbr_record(m))
  expect_true(all(rec$extinct_day < 40L))
  expect_identical(max(rec$final_state), 0)
  d <- rec$extinct_day[1L]
  expect_identical(max(rec$totals[rec$day > d, , 1L]), 0)
  expect_identical(max(rec$counter_sums[seq.int(d + 1L, 40L), , 1L]), 0)
})

test_that("summaries say what a recorded run is missing", {
  m <- cbr_model()
  # quantities only: no states or totals (all counters are kept without
  # sg_counter() terms, as in simulate_ode())
  q_only <- simulate_cb(m, cbr_init(m), seed = 1, record = sg_record(
    sg_quantity("adults", stage = "adults")
  ))
  expect_false(is.null(q_only$counters))
  expect_error(sg_totals(q_only), "states or totals")
  # totals and a counter sum: no full counters
  tot <- simulate_cb(m, cbr_init(m), seed = 1, record = sg_record(
    "totals", sg_counter("cases", counter = "fertile_batches")
  ))
  expect_s3_class(sg_totals(tot), "data.frame")
  expect_error(sg_cases(tot), "needs all daily counters")
  expect_error(sg_flow_balance(tot), "needs all daily counters")
})

test_that("parallel replicates record the same", {
  skip_if_not_installed("future.apply")
  old <- future::plan("sequential")
  on.exit(future::plan(old))
  m <- cbr_model()
  rr <- cbr_record(m)
  a <- simulate_cb(m, cbr_init(m), n_rep = 2, seed = 8, record = rr)
  b <- simulate_cb(m, cbr_init(m), n_rep = 2, seed = 8, record = rr, parallel = TRUE)
  expect_identical(b$recorded, a$recorded)
  expect_identical(b$counter_sums, a$counter_sums)
  expect_identical(b$final_state, a$final_state)
})

test_that("the memory limit counts what is recorded", {
  m <- cbr_model(days = 10)
  rr <- cbr_record(m)
  r <- .sg_recorder(rr, m, attr(.sg_model_view(m), "counter_names"))
  # 11 output days of W_x rows, 10 days of counter sums and 12 release
  # columns, the final state; 2 replicates
  bytes <- 8 * 2 * (11 * length(r$names_x) + 10 * (length(r$names_c) + 12) + m$n_states)
  expect_no_error(simulate_cb(m, cbr_init(m), n_rep = 2, seed = 1, record = rr, max_memory = bytes))
  expect_error(
    simulate_cb(m, cbr_init(m), n_rep = 2, seed = 1, record = rr, max_memory = bytes - 1),
    "max_memory"
  )
})

test_that("print() shows recorded quantities and counter sums", {
  m <- cbr_model()
  out <- capture.output(print(simulate_cb(m, cbr_init(m), n_rep = 2, seed = 1,
    record = sg_record(
      sg_quantity("adults", stage = "adults", by = "node"),
      sg_counter("cases", counter = "fertile_batches")
    )
  )))
  expect_true(any(grepl("recorded (last output day, mean): adults@a", out, fixed = TRUE)))
  expect_true(any(grepl("counter sums (whole run, mean): cases", out, fixed = TRUE)))
  expect_false(any(grepl("final totals", out)))
  expect_true(any(grepl("nodes: 3, with flies at the end", out)))
})
