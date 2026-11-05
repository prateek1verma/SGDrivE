# Summaries, flow balance and checks of models with several nodes (D-032;
# phase 1 WP8).

on_cube <- sg_cube(sg_test_cube_mendelian())

# Three nodes (a source and two study nodes) on McAllen weather with a chain
# of movement, an open end, a release in node "a" and the node extinction rule.
on_sim <- function(hosts = FALSE, days = 120) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-03-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  nd <- data.frame(
    id = c("src", "a", "b"), area = c(10, 20, 30), role = c("source", "study", "study")
  )
  if (!hosts) nd$wound_density <- c(NA, 8, 3)
  land <- sg_landscape(nd, wx)
  layer <- sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17),
    rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 0.6, 0))
  )
  rel <- sg_releases("aa",
    stage = "adult", sex = "male", times = c(20, 50), number = 100, node = "a"
  )
  m <- sg_model(p, on_cube, releases = rel, landscape = land, movement = layer)
  sim <- simulate_ode(m,
    init = sg_init(m, nodes = "src"), output_every = 10L, node_extinction = TRUE
  )
  list(m = m, sim = sim)
}

# Direct sums of state columns of the given nodes and blocks at every output day.
on_sum <- function(sim, nodes, blocks) {
  nm <- colnames(sim$state)
  keep <- substr(nm, 1L, 1L) %in% blocks & !grepl("^H", nm) & sub("^.*@", "", nm) %in% nodes
  rowSums(sim$state[, keep, drop = FALSE])
}

test_that("sg_totals() sums the study nodes by default, or the nodes asked for", {
  r <- on_sim()
  sim <- r$sim
  tot <- sg_totals(sim, by = character(0))
  expect_equal(tot$value, on_sum(sim, c("a", "b"), c("E", "L", "P", "M", "U", "F")))
  all3 <- sg_totals(sim, by = character(0), nodes = c("src", "a", "b"))
  expect_equal(all3$value, on_sum(sim, c("src", "a", "b"), c("E", "L", "P", "M", "U", "F")))
  all_stages <- c("E", "L", "P", "M", "U", "F")
  expect_equal(sg_totals(sim, by = character(0), nodes = 1)$value, on_sum(sim, "src", all_stages))
  # per km^2: divided by the summed area of the chosen nodes
  expect_equal(sg_totals(sim, by = character(0), per = "km2")$value, tot$value / 50)
  expect_error(sg_totals(sim, nodes = "c"), "Unknown node \"c\"")
  expect_error(sg_totals(sim, nodes = 4), "between 1 and 3")
})

test_that("sg_totals() can group by node, per km^2 by each node's area", {
  r <- on_sim()
  sim <- r$sim
  tn <- sg_totals(sim, by = c("stage", "node"), stages = "adults", nodes = c("src", "a", "b"))
  expect_setequal(unique(tn$node), c("src", "a", "b"))
  males_b <- tn$value[tn$stage == "males" & tn$node == "b"]
  expect_equal(males_b, on_sum(sim, "b", "M"))
  tk <- sg_totals(sim, by = "node", per = "km2", nodes = c("src", "a", "b"))
  expect_equal(tk$value[tk$node == "a"], on_sum(sim, "a", c("E", "L", "P", "M", "U", "F")) / 20)
  # as.data.frame lists every node by node
  df <- as.data.frame(sim)
  expect_true(all(c("stage", "sex", "genotype", "mate", "node") %in% names(df)))
  expect_setequal(unique(df$node), c("src", "a", "b"))
  last <- df[df$day == max(df$day), ]
  expect_equal(sum(last$value), sum(sim$state[nrow(sim$state), !grepl("^H", colnames(sim$state))]))
})

test_that("females, alleles and hosts are summed over the chosen nodes", {
  r <- on_sim(hosts = TRUE)
  sim <- r$sim
  f <- sg_females(sim)
  expect_equal(f$females, on_sum(sim, c("a", "b"), c("U", "F")))
  expect_equal(sg_females(sim, per = "km2", nodes = "b")$mated, on_sum(sim, "b", "F") / 30)
  af <- sg_allele_freq(sim, nodes = "a")
  last <- af[af$day == max(af$day), ]
  expect_equal(sum(last$freq), 1)
  h <- sg_hosts(sim)
  expect_equal(h$hosts, rep(sum(r$m$hosts$n[2:3]), nrow(h)))
  expect_equal(sg_hosts(sim, nodes = "src")$hosts[1L], r$m$hosts$n[1L])
  expect_equal(sg_hosts(sim, per = "km2", nodes = "src")$hosts[1L], 50)
})

test_that("cases are summed over the chosen nodes", {
  r <- on_sim()
  sim <- r$sim
  k <- sim$counters
  cs <- sg_cases(sim, period = "day", per = "node")
  expect_equal(cs$cases, unname(k[, "fertile_batches@a"] + k[, "fertile_batches@b"]))
  expect_equal(cs$events, unname(k[, "egg_batches@a"] + k[, "egg_batches@b"]))
  # wound fraction: per wound of the chosen nodes (8 x 20 + 3 x 30)
  expect_equal(cs$wound_fraction, cs$cases / 250)
  src <- sg_cases(sim, period = "month", per = "km2", nodes = "src")
  expect_equal(sum(src$cases), sum(k[, "fertile_batches@src"]) / 10)
  rh <- on_sim(hosts = TRUE)
  kh <- rh$sim$counters
  ch <- sg_cases(rh$sim, period = "day", per = "node", nodes = c("a", "b"))
  expect_equal(ch$cases, unname(kh[, "infestations@a"] + kh[, "infestations@b"]))
  expect_equal(ch$cases_per_host, ch$cases / sum(rh$m$hosts$n[2:3]))
})

test_that("the flow balance closes node by node with movement, losses and removals", {
  for (hosts in c(FALSE, TRUE)) {
    r <- on_sim(hosts)
    fb <- sg_flow_balance(r$sim)
    expect_identical(names(fb)[1:2], c("day", "node"))
    expect_true(all(c("movement", "removed") %in% names(fb)))
    expect_setequal(unique(fb$node), c("src", "a", "b"))
    scale <- max(abs(fb$change), abs(fb$eggs_laid))
    expect_lt(max(abs(fb$residual)) / scale, 1e-9)
    if (hosts) expect_lt(max(abs(fb$host_residual)), 1e-8)
    # releases land in node a only
    expect_identical(sum(fb$releases[fb$node != "a"]), 0)
    expect_equal(sum(fb$releases[fb$node == "a"]), 200)
  }
  # the release matrix has one column per stage and node
  expect_true(all(c("M@a", "U@src") %in% colnames(on_sim()$sim$releases)))
})

test_that("check_model() checks every node and the movement", {
  r <- on_sim(days = 60)
  ck <- check_model(r$m, scale = FALSE)
  expect_s3_class(ck, "sg_check")
  expect_identical(ck$nodes$node, c("src", "a", "b"))
  expect_identical(ck$nodes$wounds, r$m$counts$H)
  expect_identical(ck$movement$groups[[1L]]$open_nodes, "b")
  expect_equal(unname(ck$movement$groups[[1L]]$max_rate), c(0.05, 0.1, 0.17))
  expect_true(any(grepl("fractions of a fly", ck$warnings)))
  out <- capture.output(print(ck))
  expect_true(any(grepl("states in 3 nodes", out)))
  expect_true(any(grepl("open edges in b", out)))
  # the population scale is checked with movement over a year
  r1 <- on_sim(days = 365)
  ck1 <- suppressWarnings(check_model(r1$m))
  expect_true("min_wild" %in% names(ck1$nodes))
  expect_true(isTRUE(ck1$scale$available))
})

test_that("plots and print use the study nodes unless told otherwise", {
  r <- on_sim()
  d <- sg_plot_data(r$sim, "stages", per = "node")
  d_src <- sg_plot_data(r$sim, "stages", per = "node", nodes = "src")
  expect_false(isTRUE(all.equal(d$value, d_src$value)))
  out <- capture.output(print(r$sim))
  expect_true(any(grepl("final totals \\(study nodes\\)", out)))
})
