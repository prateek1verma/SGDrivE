# simulate_cb() for landscapes (phase 2 WP4; D-032): results per node,
# releases per node, summaries against direct sums of states and counters,
# flow balance per node, landscape-wide extinction, parallel replicates and
# print().

cbl_cube <- sg_cube(sg_test_cube_mendelian())

cbl_default_releases <- function() {
  sg_releases("aa", stage = "adult", sex = "male", times = c(5, 12), number = 40, node = "b")
}

cbl_model <- function(hosts = FALSE, days = 30, releases = cbl_default_releases()) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  mild <- sg_forcing(sg_constant_weather(22, days = days, start = start), p)
  nd <- data.frame(
    id = c("a", "b", "c"), area = c(40, 60, 80), role = c("study", "study", "source")
  )
  lay <- sg_movement_layer(
    list(male = 0.1, virgin = 0.2, mated = 0.15),
    rbind(c(0, 1, 0), c(0.5, 0, 0.3), c(0, 1, 0)) # node b open
  )
  sg_model(p, cbl_cube,
    landscape = sg_landscape(nd, list(wx, wx, mild)), movement = lay, releases = releases
  )
}

cbl_init <- function(m) {
  sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10), nodes = "a")
}

cbl_sim <- local({
  cache <- list()
  function(hosts = FALSE, keep = "state") {
    key <- paste(hosts, keep)
    if (is.null(cache[[key]])) {
      m <- cbl_model(hosts)
      cache[[key]] <<- simulate_cb(m, cbl_init(m),
        n_rep = 3, seed = 4, output_every = 5, keep = keep
      )
    }
    cache[[key]]
  }
})

test_that("simulate_cb() runs landscapes with results per node", {
  m <- cbl_model()
  sim <- cbl_sim()
  expect_s3_class(sim, "sg_sim_cb")
  expect_identical(dimnames(sim$state)[[2L]], sg_state_names(m))
  expect_identical(dim(sim$state), c(7L, m$n_states, 3L))
  expect_identical(dimnames(sim$counters)[[2L]], attr(.sg_model_view(m), "counter_names"))
  move_cols <- c("emigrants_M@a", "immigrants_F@c", "lost_U@b")
  expect_true(all(move_cols %in% dimnames(sim$counters)[[2L]]))
  # releases per block and node, as in simulate_ode()
  expect_identical(dimnames(sim$releases)[[2L]], colnames(.sg_release_matrix(m)))
  expect_identical(unname(apply(sim$releases, c(2L, 3L), sum)["M@b", ]), rep(80, 3))
  expect_identical(sum(sim$releases), 240)
  # flies started in node a only and reached c (through b) by movement
  fly <- !.sg_is_host_name(sg_state_names(m))
  at_c <- fly & endsWith(sg_state_names(m), "@c")
  expect_identical(sum(sim$state[1L, at_c, ]), 0)
  expect_gt(min(apply(sim$state[7L, at_c, , drop = FALSE], 3L, sum)), 0)
  # flies are lost only from the open node b
  k <- sim$counters
  expect_gt(sum(k[, c("lost_M@b", "lost_U@b", "lost_F@b"), ]), 0)
  expect_identical(sum(k[, grepl("^lost_.@[ac]$", dimnames(k)[[2L]]), ]), 0)
})

test_that("summaries of stochastic landscapes equal direct sums", {
  for (hosts in c(FALSE, TRUE)) {
    m <- cbl_model(hosts)
    sim <- cbl_sim(hosts)
    nm <- sg_state_names(m)
    node <- .sg_name_node(nm)
    block <- substr(nm, 1L, 1L)
    direct <- function(b, ids, r) {
      rowSums(matrix(sim$state[, block == b & !.sg_is_host_name(nm) & node %in% ids, r], 7L))
    }
    tot <- sg_totals(sim, by = "stage")
    for (r in 1:3) {
      # default: the study nodes a and b
      v <- tot$value[tot$replicate == r & tot$stage == "pupae"]
      expect_identical(v, direct("P", c("a", "b"), r))
    }
    tc <- sg_totals(sim, by = "stage", nodes = "c")
    expect_identical(tc$value[tc$replicate == 2L & tc$stage == "eggs"], direct("E", "c", 2L))
    bn <- sg_totals(sim, by = "node", stages = "adults", nodes = 1:3)
    adults <- function(id, r) direct("M", id, r) + direct("U", id, r) + direct("F", id, r)
    expect_identical(bn$value[bn$replicate == 3L & bn$node == "c"], adults("c", 3L))
    # sg_cases(): daily counters summed over the chosen nodes
    cs <- sg_cases(sim, per = "node", nodes = c("b", "c"))
    col <- if (hosts) "infestations" else "fertile_batches"
    for (r in 1:3) {
      expect_identical(
        cs$cases[cs$replicate == r],
        unname(rowSums(sim$counters[, paste0(col, c("@b", "@c")), r]))
      )
    }
    if (hosts) {
      hs <- sg_hosts(sim, nodes = "c")
      hi <- block == "H" & substr(nm, 2L, 2L) == "I" & node == "c"
      expect_identical(hs$infested[hs$replicate == 1L], rowSums(sim$state[, hi, 1L]))
    }
    # as.data.frame(): all nodes, by node
    d <- as.data.frame(sim)
    expect_setequal(unique(d$node), c("a", "b", "c"))
    last <- d$replicate == 1L & d$day == 30
    expect_identical(sum(d$value[last]), sum(sim$state[7L, !.sg_is_host_name(nm), 1L]))
  }
})

test_that("the flow balance of every node closes exactly", {
  for (hosts in c(FALSE, TRUE)) {
    for (keep in c("state", "totals")) {
      fb <- sg_flow_balance(cbl_sim(hosts, keep))
      expect_setequal(unique(fb$node), c("a", "b", "c"))
      expect_identical(max(abs(fb$residual)), 0)
      expect_true(any(fb$movement != 0))
      if (hosts) expect_identical(max(fb$host_residual), 0)
    }
  }
})

test_that("totals of a landscape equal the summed states of the same replicates", {
  m <- cbl_model()
  a <- cbl_sim()
  b <- cbl_sim(keep = "totals")
  agg <- .sg_totals_matrix(m)
  expect_true(all(grepl("@[abc]$", colnames(agg))))
  for (r in 1:3) expect_identical(unname(b$totals[, , r]), unname(a$state[, , r] %*% agg))
  expect_identical(b$counters, a$counters)
})

test_that("extinction refers to the whole landscape; releases elsewhere keep a run going", {
  old_males <- function(m) {
    x0 <- numeric(m$n_states)
    x0[sg_state_index(m, "M", genotype = "AA", substage = 10, node = "a")] <- 2
    x0
  }
  # two old males in node a: they die and the landscape is extinct
  m <- cbl_model(days = 40, releases = NULL)
  sim <- simulate_cb(m, old_males(m), n_rep = 2, seed = 1, record = "state")
  expect_true(all(sim$extinct_day > 0))
  expect_identical(sum(sim$state[41L, , ][!.sg_is_host_name(sg_state_names(m)), ]), 0)
  # a release in node c on day 20 keeps the run going after the extinction
  m2 <- cbl_model(days = 40, releases = sg_releases("AA",
    stage = "adult", sex = "male", times = 20, number = 5, node = "c"
  ))
  s2 <- simulate_cb(m2, old_males(m2), n_rep = 2, seed = 1, record = "state")
  expect_true(all(s2$extinct_day < 20L))
  at_c <- endsWith(sg_state_names(m2), "@c") & !.sg_is_host_name(sg_state_names(m2))
  expect_gt(min(apply(s2$state[21L, at_c, , drop = FALSE], 3L, sum)), 0)
})

test_that("parallel replicates of a landscape equal sequential ones", {
  skip_if_not_installed("future.apply")
  old <- future::plan("sequential")
  on.exit(future::plan(old))
  m <- cbl_model(days = 15)
  a <- simulate_cb(m, cbl_init(m), n_rep = 2, seed = 8, record = "state")
  b <- simulate_cb(m, cbl_init(m), n_rep = 2, seed = 8, parallel = TRUE, record = "state")
  expect_identical(b$state, a$state)
  expect_identical(b$counters, a$counters)
  expect_identical(b$releases, a$releases)
})

test_that("the memory estimate counts every node", {
  m <- cbl_model(days = 15)
  n_rel <- 12L # 4 release blocks x 3 nodes
  n_c <- length(attr(.sg_model_view(m), "counter_names"))
  bytes <- 8 * 2 * (16 * m$n_states + 15 * (n_c + n_rel) + m$n_states) # and the final state
  expect_no_error(simulate_cb(m, cbl_init(m),
    n_rep = 2, seed = 1, max_memory = bytes,
    record = "state"
  ))
  expect_error(
    simulate_cb(m, cbl_init(m), n_rep = 2, seed = 1, max_memory = bytes - 1, record = "state"),
    "max_memory"
  )
})

test_that("print() shows the study nodes and the movers", {
  out <- capture.output(print(cbl_sim()))
  expect_match(out[2L], "final totals (study nodes)", fixed = TRUE)
  expect_true(any(grepl("nodes: 3, with flies at the end", out)))
  expect_true(any(grepl("movers [0-9,]+, lost from the landscape", out)))
  expect_true(any(grepl("all nodes", out, fixed = TRUE)))
  outh <- capture.output(print(cbl_sim(TRUE)))
  expect_true(any(grepl("hosts infested at the end", outh)))
})
