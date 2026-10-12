# Exact advancement along a chain over one time step (D-022, S16): the C++
# uniformization (src/chains.h) and the R matrix exponential
# (.sg_chain_advance_r) against closed forms and each other.

adv_cpp <- SGDrivE:::sg_chain_advance_cpp
adv_r <- SGDrivE:::.sg_chain_advance_r
avg_r <- SGDrivE:::.sg_chain_average_r

test_that("equal rates give Poisson advancement (Erlang chain)", {
  for (dt in c(0.05, 0.4, 3)) {
    q <- 2.5
    n <- 6
    for (adv in list(adv_cpp, adv_r)) {
      a <- adv(rep(q, n), n, dt)
      for (i in seq_len(n)) {
        k <- 0:(n - i)
        expect_equal(a[i, i + k], dpois(k, q * dt), tolerance = 1e-13)
        expect_equal(a[i, n + 1], ppois(n - i, q * dt, lower.tail = FALSE), tolerance = 1e-12)
        if (i > 1) expect_true(all(a[i, seq_len(i - 1)] == 0)) # no moving backwards
      }
    }
  }
})

test_that("two different rates match the hypoexponential formula", {
  q1 <- 1.3
  q2 <- 4.1
  t <- 0.7
  a <- adv_cpp(c(q1, q2), 2, t)
  expect_equal(a[1, 1], exp(-q1 * t), tolerance = 1e-14)
  expect_equal(a[1, 2], q1 / (q2 - q1) * (exp(-q1 * t) - exp(-q2 * t)), tolerance = 1e-14)
  expect_equal(a[2, 2], exp(-q2 * t), tolerance = 1e-14)
  expect_equal(a[1, 3], 1 - a[1, 1] - a[1, 2], tolerance = 1e-14)
  expect_equal(a[2, 3], 1 - exp(-q2 * t), tolerance = 1e-14)
})

test_that("C++ and R agree for segmented chains that end early, at any step size", {
  set.seed(11)
  for (rep in 1:20) {
    n <- sample(2:10, 1)
    last <- sample(seq_len(n), 1)
    q <- rep(stats::runif(3, 0.1, 15), length.out = n)[order(sample(n))]
    q[-seq_len(last)] <- 0 # classes never reached have no rate
    for (dt in c(1e-3, 0.1, 1, 30, 120)) { # 120 days at up to 15/day: squaring
      a <- adv_cpp(q, last, dt)
      b <- adv_r(q, last, dt)
      expect_lt(max(abs(a - b)), 1e-12)
      expect_equal(rowSums(a)[seq_len(last)], rep(1, last), tolerance = 1e-14) # no death
      expect_true(all(a >= 0) && all(a <= 1))
      if (last < n) expect_true(all(a[(last + 1):n, ] == 0))
      expect_true(all(a[lower.tri(a)] == 0))
    }
  }
})

test_that("death in the chain: closed forms, and C++ against R", {
  q <- 2
  mu <- 0.7
  t <- 0.4
  a <- adv_cpp(c(q, q), 2, t, mu = mu)
  expect_equal(a[1, 1], exp(-(q + mu) * t), tolerance = 1e-14)
  expect_equal(a[2, 2], exp(-(q + mu) * t), tolerance = 1e-14)
  # leaving the chain alive from the last class: the hazard stops at exit
  expect_equal(a[2, 3], q / (q + mu) * (1 - exp(-(q + mu) * t)), tolerance = 1e-14)
  expect_lt(sum(a[1, ]), 1)
  # time average for one class: mean over a uniform arrival time
  r <- (q + mu) * t
  avg <- adv_cpp(q, 1, t, average = TRUE, mu = mu)
  expect_equal(avg[1, 1], (1 - exp(-r)) / r, tolerance = 1e-14)
  expect_equal(avg[1, 2], q / (q + mu) * (1 - (1 - exp(-r)) / r), tolerance = 1e-14)
  # hazard only: survival exp(-mu t) when nobody can move
  expect_equal(adv_cpp(c(0, 0), 2, t, mu = mu), diag(exp(-mu * t), 2, 3), tolerance = 1e-14)
  set.seed(4)
  for (rep in 1:20) {
    n <- sample(2:10, 1)
    last <- sample(seq_len(n), 1)
    qq <- stats::runif(n, 0.1, 15)
    qq[-seq_len(last)] <- 0
    m <- stats::runif(1, 0, 3)
    for (dt in c(0.01, 0.1, 1, 60)) {
      expect_lt(max(abs(adv_cpp(qq, last, dt, mu = m) - adv_r(qq, last, dt, m))), 1e-13)
      expect_lt(max(abs(adv_cpp(qq, last, dt, TRUE, m) - avg_r(qq, last, dt, m))), 1e-13)
    }
  }
})

test_that("time-averaged advancement against its R counterpart and limits", {
  set.seed(2)
  for (rep in 1:20) {
    n <- sample(2:10, 1)
    last <- sample(seq_len(n), 1)
    q <- stats::runif(n, 0.1, 15)
    q[-seq_len(last)] <- 0
    for (dt in c(1e-3, 0.1, 1, 100)) {
      a <- adv_cpp(q, last, dt, average = TRUE)
      expect_lt(max(abs(a - avg_r(q, last, dt))), 1e-13)
      expect_equal(rowSums(a)[seq_len(last)], rep(1, last), tolerance = 1e-14)
    }
  }
  # a short step: arrivals have moved on by about half of A(dt)
  q <- c(4, 4, 4)
  h <- 1e-6
  expect_equal(adv_cpp(q, 3, h, average = TRUE)[1, 2], q[1] * h / 2, tolerance = 1e-5)
})

test_that("advancement has the semigroup property and the right short-step limit", {
  q <- c(3, 3, 3, 0.8, 0.8, 0.8, 0.8)
  p1 <- adv_cpp(q, 7, 0.3)
  p2 <- adv_cpp(q, 7, 0.45)
  p12 <- adv_cpp(q, 7, 0.75)
  pad <- function(p) rbind(p, c(rep(0, ncol(p) - 1), 1)) # add the absorbing exit row
  expect_lt(max(abs(pad(p1) %*% pad(p2) - pad(p12))), 1e-14)
  h <- 1e-7
  ph <- adv_cpp(q, 7, h)
  expect_equal(diag(ph), 1 - q * h, tolerance = 1e-12)
  expect_equal(ph[cbind(1:6, 2:7)], q[1:6] * h, tolerance = 1e-6)
  expect_equal(ph[7, 8], q[7] * h, tolerance = 1e-6) # last class exits
  long <- adv_cpp(q, 7, 1e3)
  expect_equal(long[, 8], rep(1, 7), tolerance = 1e-14) # everyone has left
})

test_that("degenerate chains are handled", {
  expect_identical(adv_cpp(c(2, 2), 2, 0), diag(1, 2, 3)) # no time: nobody moves
  expect_identical(adv_cpp(c(0, 0), 2, 1), diag(1, 2, 3)) # no development
  expect_identical(adv_cpp(c(1, 1), 0, 1), matrix(0, 2, 3)) # nothing reached
  expect_identical(adv_r(c(1, 1), 0, 1), matrix(0, 2, 3))
  expect_error(adv_cpp(numeric(0), 0, 1), "at least one class")
  expect_error(adv_cpp(1, 2, 1), "between 0")
  expect_error(adv_cpp(1, 1, -1), ">= 0")
  expect_error(adv_cpp(-1, 1, 1), "finite and >= 0")
  expect_error(adv_cpp(1, 1, 1, mu = -1), "`mu` must be")
  expect_identical(avg_r(c(1, 1), 0, 1), matrix(0, 2, 3))
})

test_that("chains of a real model: adult chains of a short-lived genotype", {
  p <- sg_params()
  # S lives 0.25 x 202.8 = 50.7 dd, ending in the pre-oviposition segment
  cube <- sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25))
  m <- sg_model(p, cube, sg_forcing(sg_example_weather(), p, "1971-06-01", "1971-06-30"))
  v <- SGDrivE:::.sg_model_view(m)
  dx <- v$dx[1]
  for (g in 1:2) {
    iw <- v$inv_w_f[(g - 1) * v$n_A + seq_len(v$n_A)]
    last <- v$last_f[g]
    a <- adv_cpp(dx * iw, last, 0.1)
    expect_lt(max(abs(a - adv_r(dx * iw, last, 0.1))), 1e-14)
    # the last reached class exits at its own rate
    expect_equal(a[last, v$n_A + 1], 1 - exp(-dx * iw[last] * 0.1), tolerance = 1e-14)
  }
  expect_lt(v$last_f[2], v$last_f[1]) # the short-lived genotype stops earlier
})
