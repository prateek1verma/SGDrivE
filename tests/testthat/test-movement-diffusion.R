# sg_movement_diffusion() (phase 3 WP1; design note 6.1 and section 8):
# finite-volume rates on a chain, a grid and an irregular network, open
# edges, edge tables, the diffusion variance of a marked cohort, errors.

# A chain of n cells of spacing h and width w (matrices).
md_chain <- function(n, h = 10, w = 10) {
  border <- matrix(0, n, n)
  border[cbind(seq_len(n - 1L), 2:n)] <- w
  border <- border + t(border)
  list(area = rep(h * w, n), border = border, dist = matrix(h, n, n))
}

# The generator Q [n, n] of a layer's class (rates i -> j, rows summing to 0
# minus losses).
md_generator <- function(layer, n, class = "mated") {
  q <- matrix(0, n, n)
  e <- layer$edges
  rate <- rep_len(layer$rate[[class]], n)
  from <- as.integer(e$from)
  q[cbind(from, as.integer(e$to))] <- rate[from] * e$prob
  diag(q) <- -rate
  q
}

test_that("a chain gets D/h^2 to each neighbour and the end rates of the design note", {
  h <- 10
  d_coef <- 0.5
  g <- md_chain(5, h)
  closed <- sg_movement_diffusion(g$area, g$border, g$dist, D = d_coef)
  expect_s3_class(closed, "sg_movement_layer")
  expect_identical(closed$n_nodes, 5L)
  q <- md_generator(closed, 5)
  expect_equal(q[2, 1], d_coef / h^2)
  expect_equal(q[2, 3], d_coef / h^2)
  expect_equal(closed$rate$mated, c(1, 2, 2, 2, 1) * d_coef / h^2)
  expect_equal(rowSums(q), rep(0, 5)) # closed: nothing lost
  # open ends: 2D/h^2, half of it lost
  open <- sg_movement_diffusion(g$area, g$border, g$dist,
    D = d_coef, outer_border = c(10, 0, 0, 0, 10)
  )
  expect_equal(open$rate$male, rep(2 * d_coef / h^2, 5))
  qo <- md_generator(open, 5)
  expect_equal(-rowSums(qo)[c(1, 5)], c(1, 1) * d_coef / h^2) # lost rate at each end
  expect_equal(sum(open$edges$prob[open$edges$from == 1]), 0.5)
  # the width does not matter on a transect
  wide <- md_chain(5, h, w = 37)
  wide_layer <- sg_movement_diffusion(wide$area, wide$border, wide$dist, D = d_coef)
  expect_equal(wide_layer$rate, closed$rate)
})

test_that("per-class D scales the rates and leaves the destinations alone", {
  g <- md_chain(4)
  per_class <- list(male = 0.1, virgin = 0.4, mated = 0.2)
  a <- sg_movement_diffusion(g$area, g$border, g$dist, D = per_class)
  b <- sg_movement_diffusion(g$area, g$border, g$dist, D = 0.2)
  expect_equal(a$rate$male, b$rate$mated / 2)
  expect_equal(a$rate$virgin, b$rate$mated * 2)
  expect_equal(a$rate$mated, b$rate$mated)
  expect_identical(a$edges, b$edges)
  expect_identical(a$clock, "physiological")
  expect_identical(
    sg_movement_diffusion(g$area, g$border, g$dist, D = 0.2, clock = "calendar")$clock,
    "calendar"
  )
  # a named numeric vector works too
  v <- sg_movement_diffusion(g$area, g$border, g$dist, D = c(male = 0.1, virgin = 0.4, mated = 0.2))
  expect_equal(v$rate, a$rate)
})

test_that("an interior cell of a regular grid leaves at 4D/h^2", {
  h <- 5
  k <- 4 # 4 x 4 grid
  n <- k * k
  xy <- expand.grid(x = seq_len(k), y = seq_len(k))
  dd <- as.matrix(stats::dist(xy)) * h
  border <- (abs(dd - h) < 1e-9) * h
  layer <- sg_movement_diffusion(rep(h * h, n), border, dd, D = 0.3)
  inner <- which(xy$x %in% 2:3 & xy$y %in% 2:3)
  corner <- which(xy$x %in% c(1, k) & xy$y %in% c(1, k))
  expect_equal(layer$rate$mated[inner], rep(4 * 0.3 / h^2, length(inner)))
  expect_equal(layer$rate$mated[corner], rep(2 * 0.3 / h^2, length(corner)))
  expect_equal(nrow(layer$edges), 2L * 2L * k * (k - 1L)) # each shared side, both ways
})

test_that("an irregular network matches a hand calculation", {
  # three polygons: a touches b (12 km border) and c (5 km); b touches c (8 km)
  area <- c(a = 150, b = 90, c = 60)
  border <- matrix(c(0, 12, 5, 12, 0, 8, 5, 8, 0), 3, 3, dimnames = list(names(area), names(area)))
  dist <- matrix(c(0, 11, 13, 11, 0, 9, 13, 9, 0), 3, 3, dimnames = dimnames(border))
  d_coef <- 0.7
  layer <- sg_movement_diffusion(area, border, dist,
    D = d_coef, outer_border = c(a = 20, b = 0, c = 0)
  )
  r_ab <- d_coef * 12 / (11 * 150)
  r_ac <- d_coef * 5 / (13 * 150)
  r_out <- d_coef * 20 / ((150 / 20) * 150) # default outer distance area / outer border
  r_ba <- d_coef * 12 / (11 * 90)
  r_bc <- d_coef * 8 / (9 * 90)
  r_ca <- d_coef * 5 / (13 * 60)
  r_cb <- d_coef * 8 / (9 * 60)
  expect_equal(layer$rate$mated, c(r_ab + r_ac + r_out, r_ba + r_bc, r_ca + r_cb))
  e <- layer$edges
  p <- function(f, t) e$prob[e$from == f & e$to == t]
  expect_equal(p("a", "b"), r_ab / (r_ab + r_ac + r_out))
  expect_equal(p("c", "b"), r_cb / (r_ca + r_cb))
  expect_equal(sum(e$prob[e$from == "a"]), (r_ab + r_ac) / (r_ab + r_ac + r_out))
  # an explicit outer distance
  l2 <- sg_movement_diffusion(area, border, dist,
    D = d_coef, outer_border = c(20, 0, 0),
    outer_dist = c(4, 1, 1)
  )
  expect_equal(l2$rate$mated[1], r_ab + r_ac + d_coef * 20 / (4 * 150))
})

test_that("an edge table gives the same layer as matrices", {
  area <- c(a = 150, b = 90, c = 60)
  border <- matrix(c(0, 12, 5, 12, 0, 8, 5, 8, 0), 3, 3, dimnames = list(names(area), names(area)))
  dist <- matrix(c(0, 11, 13, 11, 0, 9, 13, 9, 0), 3, 3, dimnames = dimnames(border))
  m <- sg_movement_diffusion(area, border, dist, D = 0.7)
  once <- data.frame(
    from = c("a", "a", "b"), to = c("b", "c", "c"), border = c(12, 5, 8),
    dist = c(11, 13, 9)
  )
  twice <- rbind(once, data.frame(
    from = once$to, to = once$from, border = once$border,
    dist = once$dist
  ))
  expect_equal(sg_movement_diffusion(area, once, D = 0.7), m)
  expect_equal(sg_movement_diffusion(area, twice, D = 0.7), m)
  # node numbers with unnamed areas
  nums <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3), border = c(12, 5, 8), dist = c(11, 13, 9))
  ln <- sg_movement_diffusion(unname(area), nums, D = 0.7)
  expect_equal(ln$rate, m$rate)
  expect_identical(as.integer(ln$edges$from), match(m$edges$from, names(area)))
})

test_that("a marked cohort on a long chain spreads with variance 2Dt", {
  n <- 121
  h <- 1
  d_coef <- 0.5
  g <- md_chain(n, h)
  q <- md_generator(sg_movement_diffusion(g$area, g$border, g$dist, D = d_coef), n)
  # p(t) = p(0) exp(Q t) by uniformisation (no deaths, no births)
  t_end <- 20
  lam <- max(-diag(q))
  s <- diag(n) + q / lam
  p0 <- numeric(n)
  p0[61] <- 1
  p <- numeric(n)
  v <- p0
  k_max <- stats::qpois(1 - 1e-15, lam * t_end) + 10
  for (k in 0:k_max) {
    p <- p + stats::dpois(k, lam * t_end) * v
    v <- as.vector(v %*% s)
  }
  x <- (seq_len(n) - 61) * h
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_equal(sum(p * x), 0, tolerance = 1e-10)
  expect_equal(sum(p * x^2), 2 * d_coef * t_end, tolerance = 1e-8)
})

test_that("the layer runs in a model", {
  g <- md_chain(3, h = 20, w = 20)
  p <- sg_params()
  wx <- sg_forcing(sg_constant_weather(25, days = 10), p)
  land <- sg_landscape(data.frame(id = c("s", "m", "n"), area = g$area), wx)
  layer <- sg_movement_diffusion(c(s = 400, m = 400, n = 400),
    data.frame(from = c("s", "m"), to = c("m", "n"), border = 20, dist = 20),
    D = 50, outer_border = c(0, 0, 20)
  )
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), landscape = land, movement = layer)
  sim <- simulate_ode(m, sg_init(m, values = c(M = 100, U = 10, F = 100), nodes = "s"))
  expect_gt(sum(sim$counters[, "immigrants_F@n"]), 0)
  expect_gt(sum(sim$counters[, "lost_F@n"]), 0)
  expect_identical(sum(sim$counters[, "lost_F@s"]), 0)
})

test_that("invalid inputs give informative errors", {
  g <- md_chain(3)
  expect_error(sg_movement_diffusion(c(1, -1, 1), g$border, g$dist, D = 1), "positive node areas")
  expect_error(sg_movement_diffusion(g$area, g$border[1:2, ], g$dist, D = 1), "3 x 3")
  bad <- g$border
  bad[1, 2] <- 3
  expect_error(sg_movement_diffusion(g$area, bad, g$dist, D = 1), "symmetric")
  d0 <- g$dist
  d0[1, 2] <- d0[2, 1] <- 0
  expect_error(sg_movement_diffusion(g$area, g$border, d0, D = 1), "positive, finite")
  expect_error(sg_movement_diffusion(g$area, g$border, g$dist, D = -1), "diffusion coefficient")
  expect_error(
    sg_movement_diffusion(g$area, g$border, g$dist, D = list(male = 1, mated = 1)),
    "named list"
  )
  expect_error(
    sg_movement_diffusion(g$area, g$border, g$dist, D = 1, outer_border = c(1, 1)),
    "3 lengths"
  )
  expect_error(
    sg_movement_diffusion(g$area, g$border, g$dist, D = 1, outer_dist = c(1, 1, 1)),
    "needs `outer_border`"
  )
  e <- data.frame(from = c("a", "b"), to = c("b", "z"), border = 1, dist = 1)
  expect_error(sg_movement_diffusion(c(a = 1, b = 1, c = 1), e, D = 1), "Unknown node \"z\"")
  expect_error(sg_movement_diffusion(c(a = 1, b = 1, c = 1), e[1, ], g$dist, D = 1), "leave `dist`")
  conflict <- data.frame(from = c(1, 2), to = c(2, 1), border = c(1, 2), dist = 1)
  expect_error(sg_movement_diffusion(c(1, 1), conflict, D = 1), "different")
})
