# Rates of models with several nodes and adult movement (MODEL_SPEC section
# 8a, N1-N3, D-032): the C++ kernel against the R reference, each node
# against the one-node model, conservation, losses and the movement clock.

# built once at the top level, so helper functions below can use it
rn_cube <- sg_cube(sg_test_cube_mendelian())

rn_params <- function(hosts = FALSE, ...) {
  if (hosts) {
    sg_params(
      wound_model = "hosts", host_density = 50, wounded_fraction = 0.1,
      larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
      male_density_half_mating = 3, ...
    )
  } else {
    sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.01, ...)
  }
}

# three nodes: two share McAllen weather, the third is a constant 18 degC
rn_landscape <- function(p, days = 15) {
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  cool <- sg_forcing(sg_constant_weather(18, days = days, start = "1971-05-01"), p)
  nd <- data.frame(id = c("s", "m", "n"), area = c(100, 50, 150))
  if (p$wound_model == "hosts") {
    nd$host_density <- c(NA, 20, 80)
  } else {
    nd$wound_density <- c(NA, 8, 3)
  }
  sg_landscape(nd, list(wx, wx, cool))
}

rn_layer <- function(open = FALSE, clock = "physiological") {
  dest <- rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, if (open) 0.6 else 1, 0))
  sg_movement_layer(list(male = 0.05, virgin = 0.1, mated = c(0.17, 0.3, 0.2)), dest, clock = clock)
}

# A positive random state with host totals that match each node's hosts.
rn_state <- function(m, seed = 1) {
  set.seed(seed)
  x <- stats::runif(m$n_states, 0, 50)
  if (!is.null(m$hosts)) {
    for (i in seq_len(m$dims$n_node)) {
      busy <- c(sg_state_index(m, "HW", node = i), sg_state_index(m, "HI", node = i))
      x[busy] <- x[busy] * m$hosts$n[i] / (4 * sum(x[busy]))
      x[sg_state_index(m, "HS", node = i)] <- m$hosts$n[i] - sum(x[busy])
    }
  }
  x
}

test_that("C++ rates equal the R reference for several nodes, with and without movement", {
  cubes <- list(mendelian = sg_cube(sg_test_cube_mendelian()))
  if (requireNamespace("MGDrivE", quietly = TRUE)) {
    cubes$homing <- sg_cube(MGDrivE::cubeHomingDrive())
  }
  layers <- list(
    none = NULL, physiological = rn_layer(), calendar = rn_layer(clock = "calendar"),
    open = rn_layer(open = TRUE)
  )
  for (hosts in c(FALSE, TRUE)) {
    p <- rn_params(hosts)
    land <- rn_landscape(p)
    for (cn in names(cubes)) {
      for (ln in names(layers)) {
        m <- sg_model(p, cubes[[cn]], landscape = land, movement = layers[[ln]])
        x <- rn_state(m)
        s <- .sg_reference_setup(m)
        v <- .sg_model_view(m)
        expect_identical(attr(v, "counter_names"), s$counter_names)
        for (d in c(1L, 9L)) {
          r <- .sg_rates_r(d, x, s)
          cc <- sg_rates_cpp(v, d, x)
          info <- paste(hosts, cn, ln, d)
          expect_equal(cc$deriv, r$deriv, tolerance = 1e-12, info = info)
          expect_equal(cc$counters, unname(r$counters), tolerance = 1e-12, info = info)
          expect_equal(sg_max_rate_cpp(v, d, x), .sg_max_rate(d, x, s),
            tolerance = 1e-12, info = info
          )
        }
      }
    }
  }
})

test_that("without movement each node has the rates of its own one-node model", {
  for (hosts in c(FALSE, TRUE)) {
    p <- rn_params(hosts)
    cube <- sg_cube(sg_test_cube_mendelian())
    land <- rn_landscape(p)
    m <- sg_model(p, cube, landscape = land)
    x <- rn_state(m)
    v <- .sg_model_view(m)
    cc <- sg_rates_cpp(v, 5L, x)
    r <- .sg_rates_r(5L, x, .sg_reference_setup(m))
    n_s <- m$n_states_node
    n_c <- length(cc$counters) / 3L
    for (i in 1:3) {
      nd <- land$nodes[i, ]
      pi <- if (hosts) {
        update(p, host_density = c(50, 20, 80)[i])
      } else {
        update(p, wound_density = c(5, 8, 3)[i])
      }
      one <- sg_model(pi, cube, land$forcing[[nd$forcing_series]], area = nd$area)
      at <- (i - 1L) * n_s + seq_len(n_s)
      ref <- sg_rates_cpp(.sg_model_view(one), 5L, x[at])
      expect_identical(cc$deriv[at], ref$deriv)
      expect_identical(cc$counters[(i - 1L) * n_c + seq_len(n_c)], ref$counters)
      expect_identical(unname(r$deriv[at]), .sg_rates_r(5L, x[at], .sg_reference_setup(one))$deriv)
    }
  }
})

# Movement part of the derivatives and its counters (C++ kernel).
rn_movement_part <- function(p, layer, x = NULL, day = 3L,
                             cube = rn_cube) {
  land <- rn_landscape(p)
  m0 <- sg_model(p, cube, landscape = land)
  m1 <- sg_model(p, cube, landscape = land, movement = layer)
  if (is.null(x)) x <- rn_state(m1)
  r0 <- sg_rates_cpp(.sg_model_view(m0), day, x)
  r1 <- sg_rates_cpp(.sg_model_view(m1), day, x)
  cn <- attr(.sg_model_view(m1), "counter_names")
  list(
    m = m1, x = x, d = r1$deriv - r0$deriv,
    counters = stats::setNames(r1$counters, cn)[grepl("^(emigrants|immigrants|lost)_", cn)]
  )
}

test_that("with closed edges movement only redistributes flies", {
  mp <- rn_movement_part(rn_params(), rn_layer())
  n_s <- mp$m$n_states_node
  by_position <- rowSums(matrix(mp$d, n_s, 3L))
  expect_lt(max(abs(by_position)), 1e-12 * max(abs(mp$d)))
  # only adult blocks change
  adult <- unlist(lapply(c("M", "U", "F"), function(b) sg_state_index(mp$m, b)))
  expect_identical(max(abs(mp$d[-adult])), 0)
  ct <- mp$counters
  expect_identical(sum(ct[grepl("^lost_", names(ct))]), 0)
  expect_equal(sum(ct[grepl("^emigrants_", names(ct))]), sum(ct[grepl("^immigrants_", names(ct))]))
})

test_that("with open edges the network loses exactly the counted flies", {
  mp <- rn_movement_part(rn_params(), rn_layer(open = TRUE))
  ct <- mp$counters
  for (b in c("M", "U", "F")) {
    change <- sum(mp$d[sg_state_index(mp$m, b)])
    lost <- sum(ct[grepl(paste0("^lost_", b, "@"), names(ct))])
    expect_equal(change, -lost, tolerance = 1e-12)
  }
  expect_gt(sum(ct[grepl("^lost_", names(ct))]), 0)
  expect_equal(
    sum(ct[grepl("^emigrants_", names(ct))]),
    sum(ct[grepl("^immigrants_", names(ct))]) + sum(ct[grepl("^lost_", names(ct))])
  )
})

test_that("movers keep their genotype, mate genotype and sub-stage", {
  p <- rn_params()
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()),
    landscape = rn_landscape(p), movement = rn_layer()
  )
  x <- numeric(m$n_states)
  src <- sg_state_index(m, "F", genotype = "Aa", mate = "aa", substage = 6, node = "m")
  x[src] <- 100
  mp <- rn_movement_part(p, rn_layer(), x = x)
  changed <- which(mp$d != 0)
  same_slot <- vapply(c("s", "m", "n"), function(nd) {
    sg_state_index(m, "F", genotype = "Aa", mate = "aa", substage = 6, node = nd)
  }, 1L)
  expect_setequal(changed, same_slot)
  # node m sends half to s and half to n
  expect_equal(mp$d[same_slot[["s"]]], mp$d[same_slot[["n"]]])
  expect_equal(mp$d[same_slot[["s"]]], -mp$d[src] / 2)
})

# Leave rate of mated females in each node on `day`, from the C++ emigrant
# counter: emigrants_F = lambda * sum(F).
rn_lambda <- function(temps, layer_clock = "physiological", ...) {
  p <- sg_params(...)
  wx <- lapply(temps, function(t) sg_forcing(sg_constant_weather(t, days = 3), p))
  land <- sg_landscape(data.frame(id = paste0("n", seq_along(temps)), area = 1), wx)
  n <- length(temps)
  dest <- matrix(1 / (n - 1), n, n)
  diag(dest) <- 0
  layer <- sg_movement_layer(list(male = 0, virgin = 0, mated = 0.2), dest, clock = layer_clock)
  m <- sg_model(p, rn_cube, landscape = land, movement = layer)
  x <- numeric(m$n_states)
  x[sg_state_index(m, "F")] <- 1
  v <- .sg_model_view(m)
  cc <- sg_rates_cpp(v, 2L, x)
  em <- cc$counters[grepl("^emigrants_F@", attr(v, "counter_names"))]
  em / sum(x[sg_state_index(m, "F", node = 1)])
}

test_that("the movement clock follows N3", {
  p <- sg_params()
  dx <- function(t, q = p) sg_degree_days(t, q)
  temps <- c(10, 14.5, 20, 27, 35)
  # calendar: the given rate every day
  expect_equal(rn_lambda(temps, "calendar"), rep(0.2, 5), tolerance = 1e-14)
  # physiological, capped at 1: zero at and below T_L, dx ratio below T_ref, 1 above
  expect_equal(rn_lambda(temps), 0.2 * pmin(1, dx(temps) / dx(27)), tolerance = 1e-14)
  expect_identical(rn_lambda(temps)[1:2], c(0, 0))
  expect_equal(rn_lambda(temps)[5], 0.2, tolerance = 1e-14)
  # uncapped: the dx ratio also above T_ref
  unc <- rn_lambda(temps, movement_kappa_max = Inf)
  expect_equal(unc, 0.2 * dx(temps) / dx(27), tolerance = 1e-14)
  expect_gt(unc[5], 0.3)
  # an intermediate cap
  expect_equal(rn_lambda(temps, movement_kappa_max = 1.3)[5], 0.2 * 1.3, tolerance = 1e-14)
  # with a development-rate floor flies keep moving a little in the cold
  q <- sg_params(dev_rate_floor = 0.01)
  fl <- rn_lambda(temps, dev_rate_floor = 0.01)
  expect_equal(fl, 0.2 * pmin(1, dx(temps, q) / dx(27, q)), tolerance = 1e-14)
  expect_gt(fl[1], 0)
  # another reference temperature
  expect_equal(rn_lambda(temps, movement_t_ref = 25)[3], 0.2 * dx(20) / dx(25), tolerance = 1e-14)
})

test_that("the multi-node model view is checked by C++", {
  p <- rn_params(TRUE)
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()),
    landscape = rn_landscape(p), movement = rn_layer()
  )
  v <- unclass(.sg_model_view(m))
  info <- sg_model_view_info_cpp(v)
  expect_identical(info$dims[6], m$n_states)
  bad <- function(field, value) {
    w <- v
    w[[field]] <- value
    expect_error(sg_model_view_info_cpp(w), "model view")
  }
  bad("forcing_index", c(0L, 0L, 5L))
  bad("n_node", 2L)
  bad("group_col", c(v$group_col[-1L], 7L))
  bad("group_prob", replace(v$group_prob, 1L, 1.5))
  bad("group_rate", replace(v$group_rate, 2L, -1))
  bad("group_row_ptr", rev(v$group_row_ptr))
  bad("node_N_H", c(v$node_N_H[1L], 0, v$node_N_H[3L]))
  bad("kappa_max", NA_real_)
  bad("dx_ref", 0)
  bad("n_counters", v$n_counters - 9L)
  w <- v
  w$kappa_max <- Inf
  expect_silent(sg_model_view_info_cpp(w))
})
