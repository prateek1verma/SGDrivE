# Adult movement in the chain-binomial step (phase 2 WP3; MODEL_SPEC N4,
# D-032): the C++ step against the R reference draw for draw, conservation
# of flies, movers keeping their compartment, the order of draws (no fly
# moves twice in a step) and the binomial and multinomial statistics.

cbm_ns <- asNamespace("ScrewDrivE")
cbm_cube <- sg_cube(sg_test_cube_mendelian())

cbm_params <- function(hosts) {
  if (hosts) {
    sg_params(
      wound_model = "hosts", host_density = 50, wounded_fraction = 0.1,
      larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
      male_density_half_mating = 2
    )
  } else {
    sg_params(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.01)
  }
}

# Destinations of three nodes: closed rows, or node a open (40% of its
# emigrants leave the landscape).
cbm_dest <- function(open = FALSE) {
  d <- rbind(c(0, 0.7, 0.3), c(0.5, 0, 0.5), c(0, 1, 0))
  if (open) d[1L, ] <- c(0, 0.4, 0.2)
  d
}

cbm_layer <- function(clock = "physiological", open = FALSE, dest = cbm_dest(open)) {
  sg_movement_layer(list(male = 0.3, virgin = 0.5, mated = c(0.4, 0.2, 0.6)), dest, clock = clock)
}

# three nodes: two on McAllen weather (one forcing series), one at 22 degC
cbm_model <- function(p, cube = cbm_cube, movement = NULL, days = 6) {
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  mild <- sg_forcing(sg_constant_weather(22, days = days, start = start), p)
  nd <- data.frame(id = c("a", "b", "c"), area = c(5, 8, 10))
  if (p$wound_model == "hosts") {
    nd$host_density <- c(NA, 30, 80)
  } else {
    nd$wound_density <- c(NA, 8, 3)
  }
  sg_model(p, cube, landscape = sg_landscape(nd, list(wx, wx, mild)), movement = movement)
}

cbm_init <- function(m) {
  sg_init(m, values = c(E = 8, L = 8, P = 8, M = 30, U = 6, F = 30), nodes = c("a", "c"))
}

# C++ and R runs with the same seed, and the two deterministic maps.
cbm_expect_same_runs <- function(m, x0, info, dt = 0.25) {
  set.seed(42)
  r <- cbm_ns$.sg_cb_run_r(m, init = x0, dt = dt, stochastic = TRUE)
  set.seed(42)
  cc <- cbm_ns$.sg_cb_run_cpp(m, init = x0, dt = dt, stochastic = TRUE)
  testthat::expect_identical(cc$state, r$state, info = info)
  testthat::expect_identical(cc$counters, r$counters, info = info)
  dr <- cbm_ns$.sg_cb_run_r(m, init = x0, dt = dt, stochastic = FALSE)
  dc <- cbm_ns$.sg_cb_run_cpp(m, init = x0, dt = dt, stochastic = FALSE)
  testthat::expect_equal(dc$state, dr$state, tolerance = 1e-12, info = info)
  testthat::expect_equal(dc$counters, dr$counters, tolerance = 1e-12, info = info)
  invisible(r)
}

test_that("C++ steps with movement equal the R reference draw for draw", {
  cubes <- list(mendelian = cbm_cube)
  if (requireNamespace("MGDrivE", quietly = TRUE)) {
    cubes$homing <- sg_cube(MGDrivE::cubeHomingDrive())
  }
  layers <- list(
    none = NULL, physiological = cbm_layer(), calendar = cbm_layer("calendar"),
    open = cbm_layer(open = TRUE)
  )
  for (hosts in c(FALSE, TRUE)) {
    for (cn in names(cubes)) {
      for (ln in names(layers)) {
        m <- cbm_model(cbm_params(hosts), cubes[[cn]], layers[[ln]], days = 4)
        r <- cbm_expect_same_runs(m, cbm_init(m), paste(hosts, cn, ln))
        moved <- sum(r$counters[, grepl("^emigrants_", colnames(r$counters))])
        if (ln == "none") {
          expect_false(any(grepl("^emigrants_", colnames(r$counters))))
        } else {
          expect_gt(moved, 0)
        }
      }
    }
  }
})

test_that("destinations are the union of the groups' edges", {
  # one layer per model for now; two groups (one per clock) are put together
  # by hand to cover the union of edges and the weights w_G
  p <- cbm_params(FALSE)
  m <- cbm_model(p, movement = cbm_layer())
  other <- cbm_model(p, movement = cbm_layer("calendar", dest = rbind(
    c(0, 0, 0.5), c(0, 0, 0), c(0.25, 0.25, 0)
  )))
  m$movement$groups <- c(m$movement$groups, other$movement$groups)
  cbm_expect_same_runs(m, cbm_init(m), "two groups")
})

test_that("closed edges keep every fly; open edges lose exactly the lost ones", {
  for (open in c(FALSE, TRUE)) {
    m <- cbm_model(cbm_params(FALSE), movement = cbm_layer(open = open))
    set.seed(3)
    r <- cbm_ns$.sg_cb_run_cpp(m, init = cbm_init(m), dt = 0.1, stochastic = TRUE)
    cn <- colnames(r$counters)
    for (cl in c("M", "U", "F")) {
      sum_of <- function(what) rowSums(r$counters[, startsWith(cn, paste0(what, "_", cl, "@"))])
      expect_identical(sum_of("emigrants"), sum_of("immigrants") + sum_of("lost"))
      if (open) {
        expect_gt(sum(sum_of("lost")), 0)
        expect_identical(sum(r$counters[, paste0("lost_", cl, c("@b", "@c"))]), 0)
      } else {
        expect_identical(sum(sum_of("lost")), 0)
      }
    }
  }
})

# A work area of three nodes (only the adult blocks the movement step uses)
# with random whole counts.
cbm_work <- function(n_m = 4L, n_f = 8L) {
  lapply(1:3, function(i) {
    list(
      new_m = as.numeric(stats::rpois(n_m, 20)), new_u = as.numeric(stats::rpois(n_m, 5)),
      new_f = as.numeric(stats::rpois(n_f, 20))
    )
  })
}

test_that("movers keep their compartment; the network total is unchanged", {
  m <- cbm_model(cbm_params(FALSE), movement = cbm_layer("calendar"))
  s <- cbm_ns$.sg_reference_setup(m)
  draw <- cbm_ns$.sg_cb_draws(TRUE)
  set.seed(9)
  for (rep in 1:20) {
    work <- cbm_work()
    moved <- cbm_ns$.sg_cb_move_r(work, 2L, s, 0.5, draw)
    for (b in c("new_m", "new_u", "new_f")) {
      before <- Reduce(`+`, lapply(work, `[[`, b))
      after <- Reduce(`+`, lapply(moved$work, `[[`, b))
      expect_identical(after, before) # each compartment (genotype, mate, sub-stage) across nodes
    }
    cnt <- matrix(moved$counters, 9L)
    expect_identical(sum(cnt[1:3, ]), sum(cnt[4:6, ]))
    expect_identical(sum(cnt[7:9, ]), 0)
  }
})

test_that("no fly moves twice in a step (ring a -> b -> c -> a)", {
  ring <- rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
  lay <- sg_movement_layer(list(male = 5, virgin = 5, mated = 5), ring, clock = "calendar")
  m <- cbm_model(cbm_params(FALSE), movement = lay)
  x0 <- as.numeric(sg_init(m, values = c(M = 50, U = 50, F = 50), nodes = "a"))
  view <- cbm_ns$.sg_model_view(m)
  fly_c <- !grepl("^H", sg_state_names(m)) & endsWith(sg_state_names(m), "@c")
  fly_b <- !grepl("^H", sg_state_names(m)) & endsWith(sg_state_names(m), "@b")
  set.seed(1)
  for (rep in 1:10) {
    y <- sg_cb_step_cpp(view, 2L, 0.1, x0, TRUE)
    expect_identical(max(y$x[fly_c]), 0)
    expect_gt(sum(y$x[fly_b]), 0)
  }
})

test_that("emigrant counts are binomial, destinations multinomial", {
  m <- cbm_model(cbm_params(FALSE), movement = cbm_layer("calendar", open = TRUE))
  s <- cbm_ns$.sg_reference_setup(m)
  draw <- cbm_ns$.sg_cb_draws(TRUE)
  dt <- 0.5
  pe <- -expm1(-0.4 * dt) # mated females of node a
  prob <- c(0.4, 0.2, 0.4) # to b, to c, lost
  x <- 60
  zero <- list(new_m = 0, new_u = 0, new_f = 0)
  work <- list(list(new_m = 0, new_u = 0, new_f = x), zero, zero)
  n_rep <- 4000L
  set.seed(11)
  out <- t(vapply(seq_len(n_rep), function(k) {
    cnt <- matrix(cbm_ns$.sg_cb_move_r(work, 1L, s, dt, draw)$counters, 9L)
    c(cnt[3L, 1L], cnt[6L, 2L], cnt[6L, 3L], cnt[9L, 1L]) # emigrants F@a; F to b, c; lost
  }, numeric(4)))
  z <- function(v, mean, var) (mean(v) - mean) / sqrt(var / n_rep)
  expect_lt(abs(z(out[, 1L], x * pe, x * pe * (1 - pe))), 4)
  expect_equal(var(out[, 1L]), x * pe * (1 - pe), tolerance = 0.1)
  q <- pe * prob # each count unconditionally binomial with x and pe times p_j
  for (j in 1:3) {
    expect_lt(abs(z(out[, j + 1L], x * q[j], x * q[j] * (1 - q[j]))), 4)
    expect_equal(var(out[, j + 1L]), x * q[j] * (1 - q[j]), tolerance = 0.1)
  }
  expect_identical(out[, 1L], rowSums(out[, 2:4]))
  # covariance of two destinations: -x q_b q_c
  expect_equal(cov(out[, 2L], out[, 3L]), -x * q[1L] * q[2L], tolerance = 0.25)
})

test_that("C++ splits emigrants over destinations with the right shares", {
  m <- cbm_model(cbm_params(FALSE), movement = cbm_layer("calendar", open = TRUE))
  x0 <- as.numeric(sg_init(m, values = c(M = 40, U = 20, F = 40), nodes = "a"))
  view <- cbm_ns$.sg_model_view(m)
  set.seed(12)
  out <- sg_cb_step_rep_cpp(view, 2L, 0.25, x0, 3000L)
  cn <- attr(view, "counter_names")
  cnt <- out[, m$n_states + seq_along(cn)]
  colnames(cnt) <- cn
  # b and c start empty: their immigrants are all from a; given the
  # emigrants E of a class, the count at b is binomial(E, p_ab)
  for (cl in c("M", "U", "F")) {
    e <- cnt[, paste0("emigrants_", cl, "@a")]
    for (to in list(c("b", 0.4), c("c", 0.2))) {
      p <- as.numeric(to[2L])
      imm <- cnt[, paste0("immigrants_", cl, "@", to[1L])]
      expect_lt(abs(sum(imm - e * p) / sqrt(sum(e * p * (1 - p)))), 4)
    }
    split <- paste0(c("immigrants_", "immigrants_", "lost_"), cl, c("@b", "@c", "@a"))
    expect_identical(e, rowSums(cnt[, split]))
  }
})
