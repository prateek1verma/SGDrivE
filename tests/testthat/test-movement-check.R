# check_movement() and the movement part of check_model() (phase 3 WP3;
# design note 6.2 and section 8): implied D on a chain and a grid, daily and
# lifetime leave probabilities, the table of design note 5.3, Matlock &
# Skoda's daily leave fractions, warnings.

mc_cube <- sg_cube(sg_test_cube_mendelian())

mc_model <- function(layer, nodes, temp = 27, days = 30, ...) {
  p <- sg_params(...)
  wx <- sg_forcing(sg_constant_weather(temp, days = days), p)
  sg_model(p, mc_cube, landscape = sg_landscape(nodes, wx), movement = layer)
}

# n cells of spacing h in a row, with coordinates
mc_chain <- function(n = 5, h = 10, d_coef = 0.4, ...) {
  edges <- data.frame(from = seq_len(n - 1L), to = 2:n, border = h, dist = h)
  layer <- sg_movement_diffusion(rep(h^2, n), edges, D = d_coef, ...)
  mc_model(layer, data.frame(id = seq_len(n), area = h^2, x = (seq_len(n) - 1) * h))
}

test_that("on a chain the implied D is the D of sg_movement_diffusion()", {
  d_coef <- list(male = 0.25, virgin = 0.4, mated = 0.6)
  m <- mc_chain(d_coef = d_coef)
  ck <- check_movement(m)
  expect_s3_class(ck, "sg_movement_check")
  expect_identical(ck$dimension, 1L)
  for (cl in names(d_coef)) {
    d <- ck$nodes[ck$nodes$class == cl, ]
    expect_equal(d$D[2:4], rep(d_coef[[cl]], 3)) # interior nodes
    expect_equal(d$D[c(1, 5)], rep(d_coef[[cl]] / 2, 2)) # closed ends: half (one neighbour)
    life <- .sg_adult_life(m, cl, 27)
    expect_equal(d$rate, c(1, 2, 2, 2, 1) * d_coef[[cl]] / 100)
    expect_equal(d$daily_leave, 1 - exp(-d$rate))
    expect_equal(d$lifetime_leave, .sg_lifetime_leave(d$rate, life))
    expect_equal(d$moves_per_life, d$rate * .sg_expected_life(life))
    expect_equal(d$path_per_life, d$moves_per_life * 10) # every move is 10 km
    expect_identical(d$lost_share, rep(0, 5))
  }
  life_f <- .sg_expected_life(.sg_adult_life(m, "mated", 27))
  expect_equal(ck$life$life[ck$life$class == "mated"], life_f)
})

test_that("on a grid the implied D is the D of sg_movement_diffusion()", {
  h <- 5
  k <- 4
  xy <- expand.grid(x = (seq_len(k) - 1) * h, y = (seq_len(k) - 1) * h)
  dd <- as.matrix(stats::dist(xy))
  border <- (abs(dd - h) < 1e-9) * h
  layer <- sg_movement_diffusion(rep(h^2, k^2), border, dd, D = 0.3)
  m <- mc_model(layer, data.frame(id = seq_len(k^2), area = h^2, x = xy$x, y = xy$y))
  ck <- check_movement(m, temperatures = NULL)
  expect_identical(ck$dimension, 2L)
  expect_null(ck$temperatures)
  inner <- which(xy$x %in% c(h, 2 * h) & xy$y %in% c(h, 2 * h))
  d <- ck$nodes[ck$nodes$class == "mated", ]
  expect_equal(d$D[inner], rep(0.3, length(inner)))
})

test_that("Matlock & Skoda's daily leave fractions come back from their rates", {
  rate <- as.list(sg_rate_from_fraction(c(male = 0.10, virgin = 0.16, mated = 0.16)))
  layer <- sg_movement_layer(rate, rbind(c(0, 1), c(1, 0)))
  ck <- check_movement(mc_model(layer, data.frame(id = c("a", "b"), area = 10)))
  expect_equal(ck$nodes$daily_leave[ck$nodes$class == "male"], c(0.10, 0.10))
  expect_equal(ck$nodes$daily_leave[ck$nodes$class == "mated"], c(0.16, 0.16))
  # no coordinates: no path or D
  expect_true(all(is.na(ck$nodes$D)))
  expect_true(all(is.na(ck$nodes$path_per_life)))
})

test_that("the temperature table reproduces design note 5.3 for all three clocks", {
  nd <- data.frame(id = c("a", "b"), area = 10)
  dest <- rbind(c(0, 1), c(1, 0))
  rate <- list(male = 0.1, virgin = 0.17, mated = 0.17)
  expected <- list(
    calendar = c(4.19, 3.96, 2.81, 2.24, 1.59),
    capped = c(1.17, 2.38, 2.81, 2.24, 1.59),
    uncapped = c(1.17, 2.38, 2.81, 2.77, 2.50)
  )
  for (clock in names(expected)) {
    layer_clock <- if (clock == "calendar") "calendar" else "physiological"
    layer <- sg_movement_layer(rate, dest, clock = layer_clock)
    m <- mc_model(layer, nd, movement_kappa_max = if (clock == "uncapped") Inf else 1)
    tt <- check_movement(m)$temperatures
    mated <- tt[tt$class == "mated", ]
    expect_identical(mated$temperature, c(18, 22, 27, 30, 35))
    expect_equal(round(mated$moves_per_life, 2), expected[[clock]], label = clock)
    expect_equal(round(mated$life, 1), c(24.6, 23.3, 16.5, 13.1, 9.4))
    if (clock == "calendar") expect_equal(mated$kappa, rep(1, 5))
  }
  # kappa below T_ref follows dx(T) / dx(T_ref)
  m <- mc_model(sg_movement_layer(rate, dest), nd)
  tt <- check_movement(m, temperatures = c(14, 20))$temperatures
  p <- m$params
  expect_equal(tt$kappa[tt$class == "male"], c(0, sg_degree_days(20, p) / sg_degree_days(27, p)))
})

test_that("open edges give the lost share; implausible rates are warned about", {
  m <- mc_chain(d_coef = 0.4, outer_border = c(10, 0, 0, 0, 0))
  d <- check_movement(m)$nodes
  expect_equal(d$lost_share[d$class == "male" & d$node == 1], 0.5)
  expect_identical(d$lost_share[d$class == "male" & d$node == 3], 0)
  expect_length(check_movement(m)$warnings, 0L)
  # a 1 km-scale rate on a 10 km grid: D of 0.4 km^2/d per 1 km^2 cell
  fast <- mc_chain(d_coef = 40)
  ck <- check_movement(fast)
  expect_match(ck$warnings, "above 0.99")
  out <- capture.output(print(ck))
  expect_true(any(grepl("by class, range over nodes", out)))
  expect_true(any(grepl("warnings:", out)))
  one_node <- sg_model(sg_params(), mc_cube, sg_forcing(sg_constant_weather(27), sg_params()))
  expect_error(check_movement(one_node), "no movement")
})

test_that("check_movement() uses the same shares as the chain-binomial step", {
  # two groups (one per clock), put together by hand as in test-cb-movement.R
  nd <- data.frame(id = c("a", "b", "c"), area = 10)
  m <- mc_model(sg_movement_layer(
    list(male = 0.1, virgin = 0.3, mated = 0.2),
    rbind(c(0, 0.6, 0.2), c(0.5, 0, 0.5), c(0, 1, 0))
  ), nd, temp = 22)
  other <- mc_model(sg_movement_layer(
    list(male = 0.05, virgin = 0.05, mated = 0.1),
    rbind(c(0, 0, 0.5), c(0, 0, 0), c(0.25, 0.25, 0)),
    clock = "calendar"
  ), nd)
  m$movement$groups <- c(m$movement$groups, other$movement$groups)
  s <- .sg_reference_setup(m)
  cb <- .sg_cb_move_shares(m$movement, 3L, s, 0.1)
  kappa_22 <- function(g) .sg_clock_kappa(g$clock, 22, m$params, m$movement)
  sh <- .sg_movement_shares(m$movement, 3L, kappa_22)
  expect_equal(-expm1(-sh$lambda * 0.1), cb$emigrate)
  expect_equal(sh$lost, cb$lost)
  for (i in 1:3) {
    expect_identical(sh$prob[[i]]$to, cb$dest[[i]])
    expect_equal(sh$prob[[i]]$p, cb$prob[[i]])
  }
})

test_that("check_model() reports movement and warns about large steps", {
  wx <- sg_forcing(sg_example_weather(), sg_params(), start = "1971-05-01", end = "1971-06-30")
  h <- 10
  edges <- data.frame(from = 1:4, to = 2:5, border = h, dist = h)
  land <- sg_landscape(data.frame(id = 1:5, area = h^2, x = (0:4) * h), wx)
  diffusion <- function(d_coef) sg_movement_diffusion(rep(h^2, 5), edges, D = d_coef)
  m <- sg_model(sg_params(), mc_cube, landscape = land, movement = diffusion(0.4))
  ck <- check_model(m, scale = FALSE)
  expect_s3_class(ck$movement$check, "sg_movement_check")
  # the warmest day's kappa (capped at 1) times the leave rate
  dx_hi <- max(wx$dx) / m$movement$dx_ref
  expect_equal(ck$nodes$max_move_rate, min(1, dx_hi) * c(1, 2, 2, 2, 1) * 0.4 / 100)
  expect_equal(ck$movement$D[, "mated"], c(0.2, 0.4))
  expect_false(any(grepl("chain-binomial step", ck$warnings)))
  out <- capture.output(print(ck))
  expect_true(any(grepl("lifetime leave probability at 27 degC", out)))
  expect_true(any(grepl("implied D", out)))
  # fast movement: lifetime and step warnings; the largest rate includes movement
  fast <- sg_model(sg_params(), mc_cube, landscape = land, movement = diffusion(400))
  ckf <- check_model(fast, scale = FALSE, dt = 0.1)
  expect_true(any(grepl("above 0.99", ckf$warnings)))
  expect_true(any(grepl("within one chain-binomial step", ckf$warnings)))
  expect_equal(ckf$rates$max_rate, max(ckf$nodes$max_move_rate))
  small_dt <- check_model(fast, scale = FALSE, dt = 0.01)
  expect_false(any(grepl("chain-binomial step", small_dt$warnings)))
})
