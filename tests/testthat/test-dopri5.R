# Dormand-Prince 5(4) integrator (src/dopri5.h), tested on problems with known
# solutions independently of the screwworm model.

dp <- function(problem, y0, t1, rtol = 1e-10, atol = 1e-12, fixed_h = 0, par = 1,
               max_steps = 1e6, t0 = 0) {
  ScrewDrivE:::sg_dopri5_test_cpp(problem, y0, t0, t1, rtol, atol, fixed_h, max_steps, par)
}

test_that("exponential decay is solved to the requested accuracy", {
  for (tol in c(1e-6, 1e-9, 1e-12)) {
    r <- dp("decay", c(1, 2, 0.5), 2, rtol = tol, atol = tol * 1e-2)
    expect_identical(r$status, 0L)
    expect_lt(max(abs(r$y - c(1, 2, 0.5) * exp(-2)) / (c(1, 2, 0.5) * exp(-2))), 50 * tol)
  }
})

test_that("a harmonic oscillator stays accurate over ten periods", {
  r <- dp("oscillator", c(1, 0), 20 * pi)
  expect_lt(max(abs(r$y - c(1, 0))), 1e-8)
})

test_that("an Erlang chain matches its closed-form solution", {
  q <- 5
  n <- 25
  t1 <- 3
  r <- dp("erlang", c(1, rep(0, n - 1)), t1, par = q)
  exact <- exp(-q * t1) * (q * t1)^(0:(n - 1)) / factorial(0:(n - 1))
  expect_lt(max(abs(r$y - exact)), 1e-10)
  expect_equal(sum(r$y), sum(exact), tolerance = 1e-10)
})

test_that("fixed steps show 5th-order convergence", {
  errs <- vapply(c(0.1, 0.05, 0.025), function(h) {
    abs(dp("oscillator", c(1, 0), 2 * pi, fixed_h = h)$y[1] - 1)
  }, 0)
  order <- log2(errs[1:2] / errs[2:3])
  expect_true(all(order > 4.7 & order < 5.3), label = paste(signif(order, 3), collapse = ", "))
})

test_that("the evaluation count reflects first-same-as-last (6 per step)", {
  r <- dp("decay", 1, 1)
  # one initial evaluation, one for the starting step estimate, 6 per step
  expect_identical(r$evaluations, 2 + 6 * r$steps)
  expect_identical(r$steps, r$accepted + r$rejected)
})

test_that("stiffness is detected on a stiff problem and not on a non-stiff one", {
  stiff <- dp("prothero_robinson", 1, 10, rtol = 1e-6, atol = 1e-8, par = 1e4)
  expect_true(stiff$stiff)
  expect_identical(stiff$status, 0L)
  expect_lt(abs(stiff$y - cos(10)), 1e-6) # still accurate, just expensive
  mild <- dp("prothero_robinson", 1, 10, rtol = 1e-6, atol = 1e-8, par = 1)
  expect_false(mild$stiff)
  expect_lt(mild$steps, stiff$steps / 100)
})

test_that("limits and edge cases are handled", {
  r <- dp("decay", 1, 10, max_steps = 3)
  expect_identical(r$status, 1L) # too many steps
  r0 <- dp("decay", 1, 0)
  expect_identical(r0$y, 1) # empty interval: unchanged
  expect_identical(r0$steps, 0)
  expect_identical(dp("decay", 0, 5)$y, 0) # zero solution stays zero
  expect_error(dp("unknown", 1, 1), "unknown test problem")
  expect_error(dp("oscillator", 1, 1), "needs 2 states")
})

test_that("results do not depend on the starting time", {
  a <- dp("prothero_robinson", 1, 2, rtol = 1e-10, atol = 1e-12, par = 5)
  b <- dp("prothero_robinson", cos(1), 3, rtol = 1e-10, atol = 1e-12, par = 5, t0 = 1)
  # both start on the smooth solution y = cos(t), so both end on it
  expect_lt(abs(a$y - cos(2)), 1e-9)
  expect_lt(abs(b$y - cos(3)), 1e-9)
})
