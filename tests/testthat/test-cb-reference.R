# Reference chain-binomial step (R/cb-reference.R; D-022, S15-S19): its
# deterministic map against the ODE, flow accounting, and the stochastic step
# against its deterministic map.

skip_if_not_installed("deSolve")

ns <- asNamespace("SGDrivE")
setup <- ns$.sg_reference_setup
step <- ns$.sg_cb_step_r
day_info <- ns$.sg_cb_day_r
draws <- ns$.sg_cb_draws
run_cb <- ns$.sg_cb_run_r
wx_all <- sg_example_weather()

sit_raw <- sg_test_cube_sit(residual = 0.1)
sit_raw$eta[, "S"] <- 0.5
mk <- function(p, cube, releases = NULL, start = "1971-05-01", end = "1971-06-29") {
  sg_model(p, cube, sg_forcing(wx_all, p, start, end), releases)
}
flow_error <- function(r) {
  k <- r$counters
  births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
  deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
  change <- rowSums(r$releases) + births - deaths - k[, "emergence_deaths"]
  max(abs(diff(rowSums(r$state)) - change))
}

test_that("a very short deterministic step reproduces the ODE derivatives and flows", {
  h <- 1e-6
  cases <- list(
    list(sg_params(), sg_cube(sg_test_cube_mendelian())),
    list(sg_params(mating = "instant", h_M = 2), sg_cube(sg_test_cube_mendelian())),
    list(sg_params(larval_dd = "per_wound", gamma_L = 0.05), sg_cube(sit_raw)),
    list(sg_params(search_denominator = "all_adults", X_rec = 20), sg_cube(sit_raw)),
    list(sg_params(), sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25)))
  )
  for (cs in cases) {
    m <- mk(cs[[1]], cs[[2]])
    s <- setup(m)
    set.seed(1)
    x <- stats::runif(m$n_states, 0, 50)
    for (d in c(5, 40)) {
      st <- step(x, day_info(d, s, h), s, h, draws(FALSE))
      ode <- ns$.sg_rates_r(d, x, s)
      expect_lt(max(abs((st$x - x) / h - ode$deriv)) / max(abs(ode$deriv)), 1e-4)
      expect_lt(max(abs(st$flows / h - ode$counters)) / max(abs(ode$counters)), 1e-4)
    }
  }
})

test_that("the deterministic map converges to the ODE as dt shrinks", {
  m <- mk(
    sg_params(), sg_cube(sg_test_cube_mendelian()),
    sg_releases("aa", "adult", "male", times = c(5, 20), number = 40)
  )
  ode <- ns$.sg_sim_block_totals(simulate_ode(m, rtol = 1e-10, atol = 1e-12))
  err <- vapply(c(0.25, 0.1, 0.05), function(dt) {
    r <- run_cb(m, dt = dt, stochastic = FALSE)
    expect_lt(flow_error(r) / max(rowSums(r$state)), 1e-9)
    tot <- ns$.sg_sim_block_totals(r)
    max(vapply(colnames(tot), function(b) max(abs(tot[, b] - ode[, b])) / max(ode[, b]), 0))
  }, 0)
  expect_true(all(diff(err) < 0))
  expect_gt(err[2] / err[3], 1.6) # at least first order
  expect_lt(err[2], 0.01) # default dt = 0.1: within 1% in every block
})

test_that("stochastic runs give whole, non-negative counts with exact flow accounting", {
  m <- mk(
    sg_params(), sg_cube(sit_raw),
    sg_releases("S", "adult", "male", times = c(2, 9), number = 60.5),
    end = "1971-05-15"
  )
  init <- sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, U = 2, F = 10))
  set.seed(7)
  r <- run_cb(m, init, dt = 0.1)
  expect_true(all(r$state == round(r$state)) && all(r$state >= 0))
  expect_true(all(r$counters[, !grepl("egg_batches", colnames(r$counters))] %% 1 == 0))
  expect_identical(flow_error(r), 0)
  expect_true(all(r$releases[c(2, 9), "M"] %in% c(60, 61))) # rounded without bias
  set.seed(7)
  expect_identical(run_cb(m, init, dt = 0.1), r) # reproducible
  set.seed(8)
  expect_false(identical(run_cb(m, init, dt = 0.1)$state, r$state))
})

test_that("the mean of a stochastic step is the deterministic step", {
  for (mating in c("rate", "instant")) {
    p <- sg_params(mating = mating, larval_dd = "per_wound", gamma_L = 0.01)
    m <- mk(p, sg_cube(sit_raw))
    s <- setup(m)
    set.seed(5)
    x <- round(stats::runif(m$n_states, 0, 200))
    info <- day_info(20, s, 0.1)
    det <- unlist(step(x, info, s, 0.1, draws(FALSE)))
    set.seed(3)
    n_rep <- 200
    runs <- replicate(n_rep, unlist(step(x, info, s, 0.1, draws(TRUE))))
    se <- apply(runs, 1, stats::sd) / sqrt(n_rep)
    ok <- se > 0
    z <- (rowMeans(runs)[ok] - det[ok]) / se[ok]
    expect_lt(max(abs(z)), 5, label = paste(mating, "max |z|"))
    expect_lt(mean(abs(z) > 2), 0.15)
    # entries that never vary are exact: egg_batches (an expectation) or rare
    expect_true(all(abs(rowMeans(runs)[!ok] - det[!ok]) < 0.05 * pmax(1, det[!ok])))
  }
})

test_that("unbiased rounding, empty populations and sex-linked cubes", {
  set.seed(1)
  x <- c(0.25, 3, 7.5, 0)
  rr <- replicate(4000, ns$.sg_cb_round(x))
  expect_equal(rowMeans(rr), x, tolerance = 0.03)
  expect_true(all(rr[2, ] == 3) && all(rr[4, ] == 0))
  set.seed(1)
  ns$.sg_cb_round(c(1, 2, 3)) # whole numbers use no random numbers
  a <- stats::runif(1)
  set.seed(1)
  expect_identical(stats::runif(1), a)

  m <- mk(sg_params(), sg_cube(sg_test_cube_mendelian()), end = "1971-05-10")
  zero <- sg_init(m, values = c(E = 0))
  expect_identical(sum(run_cb(m, zero, stochastic = FALSE)$state), 0)
  expect_identical(sum(run_cb(m, zero)$state), 0)

  skip_if_not_installed("MGDrivE")
  cube <- sg_cube(MGDrivE::cubeXShredderY())
  ms <- mk(sg_params(), cube, end = "1971-05-20")
  set.seed(2)
  r <- run_cb(ms, sg_init(ms, values = c(E = 30, L = 30, P = 30, M = 20, U = 0, F = 20)))
  female_only <- cube$genotypes[cube$phi == 1]
  expect_identical(sum(r$state[, sg_state_index(ms, "M", genotype = female_only)]), 0)
  expect_identical(flow_error(r), 0)
})

test_that("individuals in sub-stages a genotype never reaches only die", {
  m <- mk(sg_params(), sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25)))
  s <- setup(m)
  # the female chain of S ends before oviposition (D-018); place virgins
  # beyond it, with no males present so that nobody mates
  last <- s$last_f[2]
  expect_lt(last, s$n_A)
  x <- numeric(m$n_states)
  i <- sg_state_index(m, "U", genotype = "S", substage = last + 1)
  x[i] <- 100
  info <- day_info(10, s, 0.1)
  st <- step(x, info, s, 0.1, draws(FALSE))
  expect_equal(st$x[i], 100 * exp(-info$mu[["F"]] * 0.1), tolerance = 1e-12)
  expect_identical(sum(st$x[-i]), 0)
})

test_that("invalid time steps are rejected", {
  m <- mk(sg_params(), sg_cube(sg_test_cube_wild()), end = "1971-05-03")
  expect_error(run_cb(m, dt = 0.3), "whole number of steps")
  expect_error(run_cb(m, dt = 2), "whole number of steps")
})
