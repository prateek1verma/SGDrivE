# Totals as the default output (decision 13) and summaries from recorded
# quantities (phase 3 WP5): every summary gives the same result from states,
# totals and the recorded quantities it accepts; sg_alleles(); clear errors.

sr_cube <- sg_cube(sg_test_cube_mendelian())

sr_model <- function(hosts = FALSE, n_node = 3L) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-06-19")
  rel <- function(node) {
    sg_releases("aa", stage = "adult", sex = "male", times = 5, number = 300, node = node)
  }
  if (n_node == 1L) {
    return(sg_model(p, sr_cube, wx, rel(NULL), area = 40))
  }
  nd <- data.frame(
    id = c("a", "b", "c"), area = c(40, 60, 30), role = c("study", "study", "source")
  )
  sg_model(p, sr_cube,
    landscape = sg_landscape(nd, wx), releases = rel("a"),
    movement = sg_movement_layer(
      list(male = 0.1, virgin = 0.2, mated = 0.15),
      rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0))
    )
  )
}

sr_init <- function(m) sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10))

test_that("both solvers keep totals, counters and the final state by default", {
  m <- sr_model()
  x0 <- sr_init(m)
  ode <- simulate_ode(m, x0)
  cb <- simulate_cb(m, x0, n_rep = 2, seed = 1)
  for (sim in list(ode, cb)) {
    expect_null(sim$state)
    expect_false(is.null(sim$totals))
    expect_false(is.null(sim$counters))
    expect_false(is.null(sim$final_state))
    expect_identical(sim$meta$record$sets, c("totals", "counters"))
  }
  expect_identical(colnames(ode$totals), colnames(.sg_totals_matrix(m)))
  # the default totals are the summed states of a full run
  full <- simulate_ode(m, x0, record = "state")
  expect_equal(unname(ode$totals), unname(full$state %*% .sg_totals_matrix(m)), tolerance = 1e-12)
  cb_full <- simulate_cb(m, x0, n_rep = 2, seed = 1, record = "state")
  expect_identical(unname(cb$totals[, , 2]), unname(cb_full$state[, , 2] %*% .sg_totals_matrix(m)))
  expect_identical(cb$final_state, cb_full$final_state)
})

test_that("summaries agree from states, totals and recorded quantities", {
  for (hosts in c(FALSE, TRUE)) {
    for (n_node in c(1L, 3L)) {
      m <- sr_model(hosts, n_node)
      x0 <- sr_init(m)
      info <- paste(hosts, n_node)
      st <- simulate_ode(m, x0, record = "state")
      tt <- simulate_ode(m, x0)
      expect_equal(sg_totals(tt), sg_totals(st), tolerance = 1e-12, info = info)
      expect_equal(sg_totals(tt, by = c("sex", "genotype"), per = "km2"),
        sg_totals(st, by = c("sex", "genotype"), per = "km2"),
        tolerance = 1e-12, info = info
      )
      expect_equal(sg_allele_freq(tt), sg_allele_freq(st), tolerance = 1e-12, info = info)
      ft <- sg_females(tt)
      fs <- sg_females(st)
      cols <- c("females", "unmated", "mated", "wild_type")
      expect_equal(ft[cols], fs[cols], tolerance = 1e-12, info = info)
      expect_true(all(is.na(ft$fertile_equivalent)))
      if (hosts) expect_equal(sg_hosts(tt), sg_hosts(st), tolerance = 1e-12, info = info)
      expect_equal(sg_cases(tt, "week"), sg_cases(st, "week"), tolerance = 1e-12, info = info)
      if (n_node == 3L) {
        by_node_t <- sg_totals(tt, by = "node", nodes = 1:3)
        expect_equal(by_node_t, sg_totals(st, by = "node", nodes = 1:3),
          tolerance = 1e-12, info = info
        )
      }
      # allele copies only (no states, no totals)
      by <- if (n_node == 3L) "node" else character(0)
      al <- simulate_ode(m, x0, record = sg_record(sg_alleles("ad", by = by)))
      expect_equal(sg_allele_freq(al), sg_allele_freq(st), tolerance = 1e-10, info = info)
      # stochastic: totals and allele copies give exactly the full-state results
      cs <- simulate_cb(m, x0, n_rep = 2, seed = 3, record = "state")
      ct <- simulate_cb(m, x0, n_rep = 2, seed = 3)
      ca <- simulate_cb(m, x0, n_rep = 2, seed = 3, record = sg_record(sg_alleles("ad", by = by)))
      expect_identical(sg_totals(ct), sg_totals(cs), info = info)
      expect_identical(sg_allele_freq(ct), sg_allele_freq(cs), info = info)
      expect_identical(sg_allele_freq(ca), sg_allele_freq(cs), info = info)
    }
  }
})

test_that("sg_alleles() weights each state by its genotype's copies", {
  m <- sr_model(n_node = 3L)
  x0 <- as.numeric(sr_init(m))
  st <- simulate_ode(m, x0, record = "state")
  rec <- sg_record(
    sg_alleles("ad", by = "node"),
    sg_alleles("fem", sex = "female", node = "b"),
    sg_alleles("custom", alleles = cbind(drive = c(AA = 0, Aa = 1, aa = 2)), stage = "M")
  )
  sim <- simulate_ode(m, x0, record = rec)
  expect_identical(
    colnames(sim$recorded),
    c("ad[A]@a", "ad[a]@a", "ad[A]@b", "ad[a]@b", "fem[A]", "fem[a]", "custom[drive]")
  )
  nm <- sg_state_names(m)
  geno <- sub("^.\\[([^],|]+).*$", "\\1", nm)
  copies_a <- c(AA = 0, Aa = 1, aa = 2)[geno]
  last <- st$state[nrow(st$state), ]
  pick <- function(blocks, node) substr(nm, 1L, 1L) %in% blocks & endsWith(nm, paste0("@", node))
  expect_equal(unname(sim$recorded[nrow(sim$recorded), "ad[a]@b"]),
    sum(last[pick(c("M", "U", "F"), "b")] * copies_a[pick(c("M", "U", "F"), "b")]),
    tolerance = 1e-12
  )
  expect_equal(unname(sim$recorded[nrow(sim$recorded), "fem[A]"]),
    sum(last[pick(c("U", "F"), "b")] * (2 - copies_a)[pick(c("U", "F"), "b")]),
    tolerance = 1e-12
  )
  # "custom": males of the study nodes a and b, drive copies as given
  male_ab <- substr(nm, 1L, 1L) == "M" & (endsWith(nm, "@a") | endsWith(nm, "@b"))
  drive_males <- sum(last[male_ab] * copies_a[male_ab])
  expect_equal(unname(sim$recorded[nrow(sim$recorded), "custom[drive]"]), drive_males,
    tolerance = 1e-12
  )
  expect_output(print(rec), "allele copies \"ad\" by node")
  expect_error(sg_alleles("a[1]"), "record name")
  expect_error(sg_alleles("x", by = "genotype"), "`by` must contain only")
  expect_error(sg_alleles("x", stage = "eggs_and_more"), "stage")
  expect_error(sg_alleles("x", alleles = c(A = 1)), "numeric matrix")
})

test_that("summaries name what a run is missing", {
  m <- sr_model()
  x0 <- sr_init(m)
  tt <- simulate_ode(m, x0)
  expect_error(sg_totals(tt, by = "mate"), "record = \"state\"")
  al <- simulate_ode(m, x0, record = sg_record(sg_alleles("ad")))
  # pooled over the study nodes: other node selections need a split by node
  expect_error(sg_allele_freq(al, nodes = "c"), "sg_alleles\\(\\)")
  expect_error(sg_allele_freq(al, stages = "eggs"), "sg_alleles\\(\\)")
  expect_error(sg_totals(al), "states or totals")
  nothing <- simulate_ode(m, x0, record = sg_record(sg_quantity("adults", stage = "adults")))
  expect_error(sg_allele_freq(nothing), "add \"totals\" or sg_alleles")
})

test_that("check_model() finds the same smallest wild population from totals", {
  p <- sg_params()
  wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
  m <- sg_model(p, sr_cube, wx, area = 40)
  ck <- check_model(m)
  expect_true(ck$scale$available)
  # the same burn-in and run with full states
  full <- simulate_ode(m, suppressWarnings(sg_init(m, "burn_in")), record = "state")
  geno <- sub("^.\\[([^],|]+).*$", "\\1", colnames(full$state))
  wild <- rowSums(full$state[, geno %in% m$cube$wild_type, drop = FALSE])
  expect_equal(ck$scale$min_wild, min(wild), tolerance = 1e-10)
})
