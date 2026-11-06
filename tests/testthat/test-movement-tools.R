# sg_movement_from_prob() and the rate conversions (phase 3 WP2; design note
# 6.1): MGDrivE2's conversion, losses, daily leave probabilities, Matlock &
# Skoda fractions, the exact lifetime leave probability and its inverse.

mt_cube <- sg_cube(sg_test_cube_mendelian())

mt_model <- function(temp = 27, ...) {
  p <- sg_params(...)
  sg_model(p, mt_cube, sg_forcing(sg_constant_weather(temp), p))
}

test_that("sg_movement_from_prob() converts daily probabilities as MGDrivE2 does", {
  P <- rbind( # nolint: object_name_linter.
    c(0.90, 0.10, 0.00),
    c(0.05, 0.90, 0.05),
    c(0.00, 0.20, 0.80)
  )
  layer <- sg_movement_from_prob(P, clock = "calendar")
  expect_s3_class(layer, "sg_movement_layer")
  expect_identical(layer$clock, "calendar")
  expect_identical(layer$n_nodes, 3L)
  expect_equal(layer$rate$male, -log(diag(P)))
  expect_identical(layer$rate$virgin, layer$rate$male)
  expect_identical(layer$rate$mated, layer$rate$male)
  # the daily leave probability is reproduced exactly
  expect_equal(1 - exp(-layer$rate$mated), 1 - diag(P))
  e <- layer$edges
  expect_equal(e$prob[e$from == 2 & e$to == 3], 0.5)
  expect_equal(e$prob[e$from == 3 & e$to == 2], 1)
  # MGDrivE2 2.1.1 movement_prob2rate(P) for this matrix (computed with that
  # package once; not a dependency): gamma and mat
  expect_equal(layer$rate$male, c(0.10536051565782628, 0.10536051565782628, 0.22314355131420971))
  m <- matrix(0, 3, 3)
  m[cbind(as.integer(e$from), as.integer(e$to))] <- e$prob
  expect_equal(m, rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0)))
})

test_that("rows below 1 become losses; node ids come from dimnames", {
  ids <- c("s", "m", "n")
  P <- matrix(c( # nolint: object_name_linter.
    0.95, 0.03, 0,
    0.02, 0.96, 0.02,
    0, 0.01, 0.97
  ), 3, 3, byrow = TRUE, dimnames = list(ids, ids))
  layer <- sg_movement_from_prob(P)
  e <- layer$edges
  expect_identical(sort(unique(e$from)), sort(ids))
  # node s: 5% leave a day, 3% to m and 2% out of the landscape
  expect_equal(e$prob[e$from == "s"], 0.03 / 0.05)
  expect_equal(e$prob[e$from == "n"], 0.01 / 0.03)
  # a node that keeps all its flies has rate 0 and no edges
  p_keep <- P
  p_keep["n", ] <- c(0, 0, 1)
  l2 <- sg_movement_from_prob(p_keep)
  expect_identical(l2$rate$male[3], 0)
  expect_false("n" %in% l2$edges$from)
})

test_that("sg_movement_from_prob() rejects invalid matrices", {
  expect_error(sg_movement_from_prob(matrix(0.5, 2, 3)), "square matrix")
  expect_error(sg_movement_from_prob(diag(2) * 1.2), "probabilities in")
  expect_error(sg_movement_from_prob(matrix(c(0, 1, 0.5, 0.5), 2)), "P\\[i, i\\] > 0")
  expect_error(sg_movement_from_prob(matrix(c(0.9, 0.2, 0.2, 0.9), 2)), "at most 1")
  named <- matrix(c(1, 0, 0, 1), 2, dimnames = list(c("a", "b"), c("b", "a")))
  expect_error(sg_movement_from_prob(named), "same node ids")
})

test_that("daily leave fractions convert to rates (Matlock & Skoda 2009)", {
  r <- sg_rate_from_fraction(c(female = 0.16, male = 0.10))
  expect_equal(unname(r), c(-log(0.84), -log(0.9)))
  expect_identical(names(r), c("female", "male"))
  expect_equal(1 - exp(-r), c(female = 0.16, male = 0.10))
  expect_identical(sg_rate_from_fraction(0), 0)
  expect_error(sg_rate_from_fraction(1), "\\[0, 1\\)")
  expect_error(sg_rate_from_fraction(-0.1), "\\[0, 1\\)")
})

test_that("the lifetime leave probability is exact for the adult chain", {
  m <- mt_model(female_background_mortality = 0.02)
  life <- .sg_adult_life(m, "mated", 27)
  # independent calculation with the sub-generator S of the transient
  # classes: E[exp(-delta tau)] = alpha (delta I - S)^-1 s0, E[tau] = alpha (-S)^-1 1
  n <- length(life$q)
  s <- diag(-(life$q + life$mu), n)
  s[cbind(seq_len(n - 1L), 2:n)] <- life$q[-n]
  s0 <- -rowSums(s)
  alpha <- c(1, numeric(n - 1L))
  for (delta in c(0, 0.01, 0.17, 2)) {
    lt <- 1 - sum(alpha * solve(delta * diag(n) - s, s0))
    expect_equal(.sg_lifetime_leave(delta, life), lt, tolerance = 1e-12)
  }
  expect_equal(.sg_expected_life(life), sum(alpha * solve(-s, rep(1, n))), tolerance = 1e-12)
  expect_identical(.sg_lifetime_leave(0, life), 0)
  # background mortality shortens life
  expect_lt(.sg_expected_life(life), .sg_expected_life(.sg_adult_life(mt_model(), "mated", 27)))
})

test_that("adult lifespans and moves per lifetime reproduce the design note table", {
  # design note 5.3: female rate 0.17 at 27 degC, default chain
  m <- mt_model()
  p <- m$params
  tab <- data.frame(
    temp = c(18, 22, 27, 30, 35), life = c(24.6, 23.3, 16.5, 13.1, 9.4),
    calendar = c(4.19, 3.96, 2.81, 2.24, 1.59), capped = c(1.17, 2.38, 2.81, 2.24, 1.59),
    uncapped = c(1.17, 2.38, 2.81, 2.77, 2.50)
  )
  for (k in seq_len(nrow(tab))) {
    el <- .sg_expected_life(.sg_adult_life(m, "mated", tab$temp[k]))
    ratio <- sg_degree_days(tab$temp[k], p) / sg_degree_days(27, p)
    expect_equal(round(el, 1), tab$life[k])
    expect_equal(round(0.17 * el, 2), tab$calendar[k])
    expect_equal(round(0.17 * min(1, ratio) * el, 2), tab$capped[k])
    expect_equal(round(0.17 * ratio * el, 2), tab$uncapped[k])
  }
})

test_that("sg_rate_from_lifetime() inverts the lifetime leave probability", {
  m <- mt_model()
  life_f <- .sg_adult_life(m, "mated", 27)
  for (target in c(0.01, 0.3, 0.9, 0.999)) {
    r <- sg_rate_from_lifetime(target, m)
    expect_equal(.sg_lifetime_leave(r, life_f), target, tolerance = 1e-10)
  }
  # vectorised, names kept, zero
  r <- sg_rate_from_lifetime(c(a = 0, b = 0.5), m)
  expect_identical(names(r), c("a", "b"))
  expect_identical(r[["a"]], 0)
  # round trip from a rate
  expect_equal(sg_rate_from_lifetime(.sg_lifetime_leave(0.17, life_f), m), 0.17, tolerance = 1e-9)
  # males use the male chain and background mortality
  mm <- mt_model(male_background_mortality = 0.1)
  r_m <- sg_rate_from_lifetime(0.5, mm, class = "male")
  expect_equal(.sg_lifetime_leave(r_m, .sg_adult_life(mm, "male", 27)), 0.5, tolerance = 1e-10)
  expect_gt(r_m, sg_rate_from_lifetime(0.5, mm, class = "mated"))
  # the physiological clock at T_ref: kappa = min(kappa_max, 1)
  half <- mt_model(movement_kappa_max = 0.5)
  expect_equal(
    sg_rate_from_lifetime(0.5, half),
    2 * sg_rate_from_lifetime(0.5, half, clock = "calendar")
  )
  expect_equal(sg_rate_from_lifetime(0.5, mt_model(movement_kappa_max = Inf)), r[["b"]])
})

test_that("sg_rate_from_lifetime() checks its inputs", {
  m <- mt_model()
  expect_error(sg_rate_from_lifetime(1, m), "\\[0, 1\\)")
  expect_error(sg_rate_from_lifetime(0.5, list()), "sg_model")
  expect_error(sg_rate_from_lifetime(0.5, m, class = "female"), "must be one of")
  expect_error(sg_rate_from_lifetime(0.5, m, genotype = "XX"), "Unknown genotype")
})
