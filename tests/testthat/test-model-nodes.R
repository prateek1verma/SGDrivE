# Models with several nodes: node-major layout, per-node quantities,
# releases and movement setup (MODEL_SPEC section 8a, D-032).

nodes_forcing <- function(p, temp = 25, days = 20) {
  sg_forcing(sg_constant_weather(temp, days = days), p)
}

nodes_landscape <- function(p, extra = list()) {
  nd <- data.frame(
    id = c("s", "m", "n"), area = c(100, 50, 150), role = c("source", "study", "study")
  )
  for (nm in names(extra)) nd[[nm]] <- extra[[nm]]
  sg_landscape(nd, list(nodes_forcing(p, 27), nodes_forcing(p, 27), nodes_forcing(p, 20)))
}

chain_layer <- function(mated = 0.17) {
  sg_movement_layer(
    list(male = 0.05, virgin = 0.1, mated = mated),
    rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 0.6, 0))
  )
}

host_p <- function(...) {
  sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1, ...)
}

test_that("each node holds a copy of the one-node layout (node slowest)", {
  p <- host_p()
  cube <- sg_cube(sg_test_cube_mendelian())
  one <- sg_model(p, cube, nodes_forcing(p))
  m <- sg_model(p, cube, landscape = nodes_landscape(p))
  s_node <- one$n_states
  expect_identical(m$dims$n_node, 3L)
  expect_identical(m$n_states_node, s_node)
  expect_identical(m$n_states, 3L * s_node)
  expect_identical(m$layout$block, one$layout$block)
  expect_identical(m$layout$start, one$layout$start)
  expect_identical(m$layout$length, one$layout$length)
  expect_match(m$layout$dims[1L], "^per node: host 1")
  for (b in c("E", "M", "F", "HS", "HI")) {
    for (k in 1:3) {
      expect_identical(
        sg_state_index(m, b, node = k),
        sg_state_index(one, b) + (k - 1L) * s_node
      )
    }
  }
  expect_identical(
    sg_state_index(m, "F", genotype = "Aa", mate = "aa", substage = 4, node = "n"),
    sg_state_index(one, "F", genotype = "Aa", mate = "aa", substage = 4) + 2L * s_node
  )
  # all nodes, node slowest
  expect_identical(sg_state_index(m, "HW"), sg_state_index(one, "HW") + c(0L, 1L, 2L) * s_node)
  expect_error(sg_state_index(m, "M", node = "x"), "Unknown node \"x\"")
  expect_error(sg_state_index(m, "M", node = 4), "between 1 and 3")
})

test_that("state names carry the node id with more than one node", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  one <- sg_model(p, cube, nodes_forcing(p))
  m <- sg_model(p, cube, landscape = nodes_landscape(p))
  nm <- sg_state_names(m)
  expect_length(nm, m$n_states)
  expect_false(anyDuplicated(nm) > 0)
  expect_identical(nm[seq_len(one$n_states)], paste0(sg_state_names(one), "@s"))
  i_n <- sg_state_index(m, "M", genotype = "aa", substage = 2, node = "n")
  expect_identical(nm[i_n], "M[aa,2]@n")
  # a one-node landscape keeps the plain names
  land1 <- sg_landscape(data.frame(id = "only", area = 1), nodes_forcing(p))
  expect_identical(sg_state_names(sg_model(p, cube, landscape = land1)), sg_state_names(one))
})

test_that("wounds, hosts and mate finding use each node's area and densities", {
  p <- sg_params(wound_density = 5, male_density_half_mating = 2)
  cube <- sg_cube(sg_test_cube_mendelian())
  m <- sg_model(p, cube, landscape = nodes_landscape(p, list(wound_density = c(NA, 8, NA))))
  expect_identical(m$area, c(100, 50, 150))
  expect_identical(m$counts$H, c(500, 400, 750))
  expect_identical(m$counts$h_M, c(200, 100, 300))
  hp <- host_p()
  mh <- sg_model(hp, cube, landscape = nodes_landscape(hp, list(host_density = c(NA, 20, NA))))
  expect_identical(mh$hosts$n, c(5000, 1000, 7500))
  expect_identical(mh$counts$H, rep(NA_real_, 3))
  expect_warning(
    sg_model(hp, cube, landscape = nodes_landscape(hp, list(wound_density = c(1, NA, NA)))),
    "wound_density has no effect"
  )
  expect_warning(
    sg_model(p, cube, landscape = nodes_landscape(p, list(host_density = c(1, NA, NA)))),
    "host_density has no effect"
  )
})

test_that("the movement layer is merged and stored with the clock settings", {
  p <- sg_params(movement_kappa_max = Inf, movement_t_ref = 25)
  cube <- sg_cube(sg_test_cube_mendelian())
  m <- sg_model(p, cube, landscape = nodes_landscape(p), movement = chain_layer())
  mv <- m$movement
  expect_length(mv$groups, 1L)
  expect_identical(mv$t_ref, 25)
  expect_identical(mv$kappa_max, Inf)
  expect_identical(mv$dx_ref, sg_degree_days(25, p))
  ref <- .sg_movement_groups(list(chain_layer()), c("s", "m", "n"))[[1L]]
  expect_identical(mv$groups[[1L]], ref)
  # a list with one layer is accepted
  m2 <- sg_model(p, cube, landscape = nodes_landscape(p), movement = list(chain_layer()))
  expect_identical(m2$movement, mv)
  expect_null(sg_model(p, cube, landscape = nodes_landscape(p))$movement)
})

test_that("invalid combinations give clear errors", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  land <- nodes_landscape(p)
  wx <- nodes_forcing(p)
  expect_error(sg_model(p, cube, wx, landscape = land), "Give `forcing` or `landscape`")
  expect_error(sg_model(p, cube, landscape = land, area = 2), "Give `area` or `landscape`")
  expect_error(sg_model(p, cube, landscape = land, adult_females = 10), "only for one-node")
  expect_error(sg_model(p, cube, wx, movement = chain_layer()), "needs a `landscape`")
  expect_error(sg_model(p, cube, landscape = "land"), "sg_landscape object")
  expect_error(sg_model(p, cube, landscape = land, movement = "chain"), "sg_movement_layer object")
  expect_error(
    sg_model(p, cube, landscape = land, movement = list(chain_layer(), chain_layer())),
    "one movement layer per model"
  )
  expect_error(sg_model(p, cube), "`forcing` must be an sg_forcing object")
  # forcings built with other parameters
  land_other <- nodes_landscape(sg_params(dev_temp_min = 14))
  expect_error(sg_model(p, cube, landscape = land_other), "different parameter values")
})

test_that("releases go to the named nodes and use each node's area", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  rel <- c(
    sg_releases("aa",
      stage = "adult", sex = "male", times = c(5, 12), number = 10, per = "km2",
      node = c("m", "n")
    ),
    sg_releases("aa", stage = "pupa", times = 3, number = 7, node = "s")
  )
  m <- sg_model(p, cube, releases = rel, landscape = nodes_landscape(p))
  r <- m$releases
  expect_named(r, c("day", "index", "count", "genotype", "stage", "sex", "node"))
  expect_identical(r$day, c(3L, 5L, 5L, 12L, 12L))
  expect_identical(r$node, c("s", "m", "n", "m", "n"))
  expect_identical(r$count, c(7, 500, 1500, 500, 1500))
  expect_identical(r$index[1L], sg_state_index(m, "P", genotype = "aa", substage = 1, node = "s"))
  expect_identical(r$index[2L], sg_state_index(m, "M", genotype = "aa", substage = 1, node = "m"))
  expect_identical(r$index[3L], sg_state_index(m, "M", genotype = "aa", substage = 1, node = "n"))
  # every release needs a known node
  expect_error(
    sg_model(p, cube,
      releases = sg_releases("aa", stage = "pupa", times = 3, number = 7),
      landscape = nodes_landscape(p)
    ),
    "every release needs a node"
  )
  expect_error(
    sg_model(p, cube,
      releases = sg_releases("aa", stage = "pupa", times = 3, number = 7, node = "north"),
      landscape = nodes_landscape(p)
    ),
    "\"north\" is not a node"
  )
  # one-node models ignore the node
  m1 <- sg_model(p, cube, nodes_forcing(p), sg_releases("aa",
    stage = "pupa", times = 3, number = 7, node = "anywhere"
  ))
  expect_named(m1$releases, c("day", "index", "count", "genotype", "stage", "sex"))
  expect_identical(m1$releases$count, 7)
})

test_that("the hash sees the landscape and movement", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  h0 <- sg_model(p, cube, landscape = nodes_landscape(p))$meta$hash
  expect_identical(sg_model(p, cube, landscape = nodes_landscape(p))$meta$hash, h0)
  land2 <- nodes_landscape(p)
  land2$nodes$area[2L] <- 51
  expect_false(identical(sg_model(p, cube, landscape = land2)$meta$hash, h0))
  hm <- sg_model(p, cube, landscape = nodes_landscape(p), movement = chain_layer())$meta$hash
  expect_false(identical(hm, h0))
  expect_false(identical(
    sg_model(p, cube, landscape = nodes_landscape(p), movement = chain_layer(0.2))$meta$hash, hm
  ))
  # movement parameters count only with movement
  p2 <- sg_params(movement_kappa_max = 2)
  expect_identical(sg_model(p2, cube, landscape = nodes_landscape(p2))$meta$hash, h0)
  expect_false(identical(
    sg_model(p2, cube, landscape = nodes_landscape(p2), movement = chain_layer())$meta$hash, hm
  ))
})

test_that("print summarises nodes and movement", {
  p <- host_p()
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()),
    landscape = nodes_landscape(p), movement = chain_layer()
  )
  out <- capture.output(print(m))
  expect_true(any(grepl(
    "3 nodes \\(2 study, 1 source\\), total area 300 km\\^2; 2 forcing series", out
  )))
  expect_true(any(grepl("hosts N_H per node: 2500-7500", out)))
  expect_true(any(grepl(
    "movement \\(physiological clock\\): 4 edges, 1 node with open edges", out
  )))
  expect_true(any(grepl("kappa = min\\(1, dx / ", out)))
})

test_that("functions without multi-node support yet stop clearly", {
  p <- sg_params()
  cube <- sg_cube(sg_test_cube_mendelian())
  m <- sg_model(p, cube, landscape = nodes_landscape(p), movement = chain_layer())
  x0 <- numeric(m$n_states)
  expect_error(simulate_cb(m, init = x0), "simulate_cb\\(\\) is not available")
  expect_error(calibrate_population(m, 10), "calibrate_population\\(\\) is not available")
  # (sg_init() and sg_equilibria() handle several nodes since WP7 and
  # check_model() since WP8); a full state vector is accepted
  x <- sg_init(m, "state", state = x0)
  expect_length(x, m$n_states)
  expect_identical(names(x), sg_state_names(m))
})

test_that("host states given for several nodes are checked node by node", {
  p <- host_p()
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), landscape = nodes_landscape(p))
  x <- numeric(m$n_states)
  x[sg_state_index(m, "HS")] <- m$hosts$n
  expect_silent(sg_init(m, "state", state = x))
  x[sg_state_index(m, "HS", node = "n")] <- 10
  expect_warning(sg_init(m, "state", state = x), "host states of node n sum to 10")
})
