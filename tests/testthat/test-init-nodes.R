# Initial states of models with several nodes (D-032; phase 1 WP7).

in_cube <- sg_cube(sg_test_cube_mendelian())

in_host_p <- function(...) {
  sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1, ...)
}

# three nodes at constant temperatures (27, 22, 13 degC) with different areas
in_landscape <- function(p, temps = c(27, 22, 13), days = 400, extra = list()) {
  nd <- data.frame(id = c("hot", "mild", "cold"), area = c(10, 20, 30))
  for (nm in names(extra)) nd[[nm]] <- extra[[nm]]
  sg_landscape(nd, lapply(temps, function(t) sg_forcing(sg_constant_weather(t, days = days), p)))
}

in_chain <- function() {
  sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17),
    rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0))
  )
}

# The one-node model of node i built from scratch.
in_one_node <- function(p, land, i, density_name) {
  nd <- land$nodes[i, ]
  q <- p
  if (!is.na(nd[[density_name]])) {
    q <- do.call(update, c(list(p), stats::setNames(list(nd[[density_name]]), density_name)))
  }
  sg_model(q, in_cube, land$forcing[[nd$forcing_series]], area = nd$area)
}

test_that("the node helper gives each node as its own one-node model", {
  for (hosts in c(FALSE, TRUE)) {
    p <- if (hosts) in_host_p() else sg_params()
    dn <- if (hosts) "host_density" else "wound_density"
    land <- in_landscape(p, extra = stats::setNames(list(c(NA, 20, 8)), dn))
    m <- sg_model(p, in_cube, landscape = land, movement = in_chain())
    for (i in 1:3) {
      nm <- .sg_node_model(m, i)
      ref <- in_one_node(p, land, i, dn)
      expect_identical(nm$counts, ref$counts)
      expect_identical(nm$hosts, ref$hosts)
      expect_identical(nm$area, ref$area)
      expect_identical(nm$forcing, ref$forcing)
      expect_identical(nm$layout, ref$layout)
      expect_identical(nm$n_states, ref$n_states)
      expect_null(nm$landscape)
      expect_null(nm$movement)
      expect_identical(sg_state_names(nm), sg_state_names(ref))
    }
  }
  one <- sg_model(sg_params(), in_cube, sg_forcing(sg_constant_weather(25, days = 5), sg_params()))
  expect_identical(.sg_node_model(one, 1L), one)
})

test_that("constant starts go into the chosen nodes; every node gets its hosts", {
  p <- in_host_p()
  land <- in_landscape(p, extra = list(host_density = c(NA, 20, 80)))
  m <- sg_model(p, in_cube, landscape = land)
  x <- sg_init(m)
  expect_s3_class(x, "sg_state")
  expect_identical(attr(x, "model_hash"), m$meta$hash)
  for (i in 1:3) {
    ref <- as.numeric(sg_init(in_one_node(p, land, i, "host_density")))
    expect_identical(unname(as.numeric(x)[.sg_node_positions(m, i)]), ref)
  }
  # only some nodes; values per km^2 scale with each node's area
  y <- sg_init(m, values = c(M = 1, F = 2), per = "km2", nodes = c("cold", "hot"))
  fly <- function(id) {
    unlist(lapply(c("E", "L", "P", "M", "U", "F"), function(b) sg_state_index(m, b, node = id)))
  }
  expect_identical(sum(y[fly("mild")]), 0)
  expect_equal(sum(y[sg_state_index(m, "M", node = "cold")]), 30 * 10)
  expect_equal(sum(y[sg_state_index(m, "M", node = "hot")]), 10 * 10)
  hosts_of <- function(id) {
    sum(y[unlist(lapply(c("HS", "HW", "HI"), function(b) sg_state_index(m, b, node = id)))])
  }
  expect_equal(vapply(c("hot", "mild", "cold"), hosts_of, 0), m$hosts$n, ignore_attr = TRUE)
  expect_equal(sum(y[sg_state_index(m, "HW", node = "mild")]), 0.1 * m$hosts$n[2])
  expect_error(sg_init(m, nodes = "north"), "Unknown node \"north\"")
  # Hardy-Weinberg per node
  z <- sg_init(m, "hardy_weinberg", alleles = c(A = 0.5, a = 0.5), nodes = 2)
  expect_gt(sum(z[sg_state_index(m, "M", genotype = "Aa", node = "mild")]), 0)
  expect_identical(sum(z[sg_state_index(m, "M", node = "hot")]), 0)
})

test_that("a state can be given as a list of per-node vectors", {
  p <- sg_params()
  m <- sg_model(p, in_cube, landscape = in_landscape(p))
  x <- as.numeric(sg_init(m))
  parts <- lapply(1:3, function(i) x[.sg_node_positions(m, i)])
  names(parts) <- c("hot", "mild", "cold")
  expect_identical(as.numeric(sg_init(m, "state", state = rev(parts))), x)
  expect_error(sg_init(m, "state", state = parts[1:2]), "one element per node")
  bad <- parts
  bad$mild <- bad$mild[-1]
  expect_error(sg_init(m, "state", state = bad), "state\\$mild` must be a numeric vector of length")
})

test_that("equilibria are computed node by node", {
  p <- in_host_p()
  land <- in_landscape(p, extra = list(host_density = c(NA, 20, 80)))
  m <- sg_model(p, in_cube, landscape = land, movement = in_chain())
  expect_message(x <- sg_init(m, "equilibrium"), "node cold: No positive equilibrium at 13 degC")
  info <- attr(x, "equilibrium")
  expect_identical(info$node, c("hot", "mild", "cold"))
  expect_identical(info$temperature, c(27, 22, 13))
  for (i in 1:2) {
    ref <- sg_init(in_one_node(p, land, i, "host_density"), "equilibrium")
    expect_identical(unname(as.numeric(x)[.sg_node_positions(m, i)]), as.numeric(ref))
    expect_identical(info$adult_females[i], attr(ref, "equilibrium")$adult_females)
  }
  # the cold node starts without flies, with its hosts
  cold <- as.numeric(x)[.sg_node_positions(m, 3L)]
  expect_equal(sum(cold), m$hosts$n[3])
  # one temperature for all nodes, or one per node
  x25 <- sg_init(m, "equilibrium", temperature = 25)
  expect_identical(attr(x25, "equilibrium")$temperature, rep(25, 3))
  expect_identical(
    attr(sg_init(m, "equilibrium", temperature = c(25, 26, 27)), "equilibrium")$temperature,
    c(25, 26, 27)
  )
  expect_error(sg_init(m, "equilibrium", temperature = c(25, 26)), "one per node \\(3\\)")
  # a varying forcing needs a temperature
  wx <- sg_forcing(sg_example_weather()[1:400, ], p)
  flat <- sg_forcing(sg_constant_weather(27, days = 400, start = wx$dates[1L]), p)
  land_v <- sg_landscape(data.frame(id = c("a", "b"), area = 1), list(wx, flat))
  expect_error(
    sg_init(sg_model(p, in_cube, landscape = land_v), "equilibrium"),
    "forcing of node\\(s\\) a does not have one constant temperature"
  )
})

test_that("Allee thresholds are computed node by node", {
  p <- sg_params(male_density_half_mating = 5)
  land <- in_landscape(p, temps = c(27, 22, 20))
  m <- sg_model(p, in_cube, landscape = land)
  x <- sg_init(m, "threshold")
  for (i in 1:3) {
    ref <- sg_init(in_one_node(p, land, i, "wound_density"), "threshold")
    expect_identical(unname(as.numeric(x)[.sg_node_positions(m, i)]), as.numeric(ref))
  }
  expect_true(all(attr(x, "threshold")$adult_females > 0))
})

test_that("sg_equilibria() lists the equilibria of every node", {
  p <- sg_params(male_density_half_mating = 5)
  land <- in_landscape(p, temps = c(27, 22, 20))
  m <- sg_model(p, in_cube, landscape = land)
  tab <- sg_equilibria(m)
  expect_identical(names(tab)[1L], "node")
  expect_identical(unique(tab$node), c("hot", "mild", "cold"))
  for (i in 1:3) {
    ref <- sg_equilibria(in_one_node(p, land, i, "wound_density"))
    got <- tab[tab$node == land$nodes$id[i], -1L]
    rownames(got) <- NULL
    expect_equal(got, ref, ignore_attr = TRUE)
  }
  expect_named(attr(tab, "R0"), c("hot", "mild", "cold"))
})

test_that("a burn-in runs all nodes with movement", {
  p <- sg_params()
  land <- in_landscape(p, temps = c(27, 27, 27))
  m <- sg_model(p, in_cube, landscape = land, movement = in_chain())
  fly <- function(x, id) {
    stages <- c("E", "L", "P", "M", "U", "F")
    sum(x[unlist(lapply(stages, function(b) sg_state_index(m, b, node = id)))])
  }
  # flies started in one node reach the others through movement
  x <- sg_init(m, "burn_in", nodes = "hot")
  expect_gt(fly(x, "cold"), 1)
  expect_match(attr(x, "burn_in")$source, "each node's first year")
  # without movement they stay where they started
  x0 <- sg_init(sg_model(p, in_cube, landscape = land), "burn_in", nodes = "hot")
  expect_identical(fly(x0, "cold"), 0)
})

test_that("two identical nodes with symmetric movement burn in like one node", {
  p <- sg_params()
  wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
  land <- sg_landscape(data.frame(id = c("a", "b"), area = 5), wx)
  layer <- sg_movement_layer(list(male = 0.05, virgin = 0.1, mated = 0.17), rbind(c(0, 1), c(1, 0)))
  m <- sg_model(p, in_cube, landscape = land, movement = layer)
  one <- sg_model(p, in_cube, wx, area = 5)
  x <- as.numeric(sg_init(m, "burn_in"))
  ref <- as.numeric(sg_init(one, "burn_in"))
  n_s <- length(ref)
  expect_lt(max(abs(x[seq_len(n_s)] - ref)) / max(ref), 1e-6)
  expect_lt(max(abs(x[n_s + seq_len(n_s)] - ref)) / max(ref), 1e-6)
})

test_that("burn-in weather is one data frame for all nodes or a list per node", {
  p <- sg_params()
  wx_all <- sg_example_weather()
  wx <- sg_forcing(wx_all, p, start = "1972-01-01", end = "1972-03-31")
  land <- sg_landscape(data.frame(id = c("a", "b"), area = 5), wx)
  m <- sg_model(p, in_cube, landscape = land, movement = sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = 0.17), rbind(c(0, 1), c(1, 0))
  ))
  burn_land <- function(weather) .sg_burn_in_landscape(m, 1L, weather)$landscape
  one <- burn_land(wx_all)
  expect_length(one$forcing, 1L)
  expect_identical(one$forcing[[1L]]$dates[1L], as.Date("1971-01-01"))
  cool <- wx_all
  cool$tmean <- cool$tmean - 5
  cool$tmax <- cool$tmax - 5
  cool$tmin <- cool$tmin - 5
  two <- burn_land(list(b = cool, a = wx_all))
  expect_length(two$forcing, 2L)
  tm <- vapply(two$nodes$forcing_series, function(k) mean(two$forcing[[k]]$tmean), 0)
  expect_lt(tm[2L], tm[1L])
  expect_error(burn_land(list(wx_all)), "one per node \\(2\\)")
  expect_error(burn_land(list(a = wx_all, c = wx_all)), "node ids")
  expect_error(.sg_burn_in_landscape(m, 1L, NULL), "repeats each node's first year")
  x <- sg_init(m, "burn_in", weather = wx_all)
  expect_match(attr(x, "burn_in")$source, "weather before the start date")
})

test_that("the burn-in can apply the node extinction rule", {
  p <- sg_params()
  # a trickle of mated females into a node too cold for development: it
  # holds a fraction of a fly without the rule and none with it
  land <- sg_landscape(
    data.frame(id = c("hot", "cold"), area = 10),
    lapply(c(27, 12), function(t) sg_forcing(sg_constant_weather(t, days = 400), p))
  )
  trickle <- sg_movement_layer(
    list(male = 0, virgin = 0, mated = c(1e-4, 0)),
    data.frame(from = "hot", to = "cold", prob = 1),
    clock = "calendar"
  )
  m <- sg_model(p, in_cube, landscape = land, movement = trickle)
  cold <- function(x) {
    sum(x[unlist(lapply(c("E", "L", "P", "M", "U", "F"), function(b) {
      sg_state_index(m, b, node = "cold")
    }))])
  }
  free <- sg_init(m, "burn_in", nodes = "hot")
  ruled <- sg_init(m, "burn_in", nodes = "hot", node_extinction = TRUE)
  expect_gt(cold(free), 0)
  expect_lt(cold(free), 1)
  expect_identical(cold(ruled), 0)
  expect_error(sg_init(m, "burn_in", node_extinction = NA), "TRUE or FALSE")
})

test_that("the default start lets simulate_ode() run a landscape", {
  p <- sg_params()
  m <- sg_model(p, in_cube, landscape = in_landscape(p, days = 30), movement = in_chain())
  sim <- simulate_ode(m)
  expect_identical(unname(sim$state[1L, ]), as.numeric(sg_init(m)))
})

test_that("print shows totals by node", {
  p <- in_host_p()
  m <- sg_model(p, in_cube, landscape = in_landscape(p), movement = in_chain())
  out <- capture.output(print(suppressMessages(sg_init(m, "equilibrium"))))
  expect_match(out[1L], "645 states in 3 nodes")
  expect_true(any(grepl("equilibrium per node", out)))
  expect_true(any(grepl("^ *hot +mild +cold", out)))
  expect_true(any(grepl("without wound", out)))
})
