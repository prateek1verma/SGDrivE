skip_if_not_installed("deSolve")

p <- sg_params()
wx60 <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-06-29")
cube3 <- sg_cube(sg_test_cube_mendelian())
m3 <- sg_model(p, cube3, wx60)

flow_error <- function(sim) {
  tot <- rowSums(sim$state)
  k <- sim$counters
  births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
  deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
  pred <- rowSums(sim$releases) + births - deaths - k[, "emergence_deaths"]
  max(abs(diff(tot) - pred)) / max(tot)
}

test_that("simulate_ode() returns daily states, counters, releases and diagnostics", {
  sim <- simulate_ode(m3)
  expect_s3_class(sim, "sg_sim")
  expect_identical(sim$day, 0:60)
  expect_identical(sim$date[1], as.Date("1971-04-30"))
  expect_identical(sim$date[61], as.Date("1971-06-29"))
  expect_identical(dim(sim$state), c(61L, m3$n_states))
  expect_identical(colnames(sim$state), sg_state_names(m3))
  expect_identical(unname(sim$state[1, ]), as.numeric(sg_init(m3)))
  expect_identical(dim(sim$counters), c(60L, length(SGDrivE:::.sg_counter_names(cube3$genotypes))))
  expect_named(sim$diagnostics, c(
    "day", "steps", "rejected", "evaluations", "stiff", "fallback", "max_rate", "min_state"
  ))
  expect_true(all(is.na(sim$diagnostics$rejected)) && !any(sim$diagnostics$fallback))
  expect_identical(sim$meta$model_hash, m3$meta$hash)
  expect_identical(sim$meta$backend, "R")
})

test_that("flow accounting closes over a simulation, with and without releases", {
  expect_lt(flow_error(simulate_ode(m3)), 1e-12)
  rel <- c(
    sg_releases("aa", "adult", "male", times = c(5, 20, 35), number = 40),
    sg_releases("aa", "adult", "female", times = 10, number = 7),
    sg_releases("Aa", "pupa", times = 15, number = 3),
    sg_releases("aa", "egg", times = 25, number = 11)
  )
  sim <- simulate_ode(sg_model(p, cube3, wx60, rel))
  expect_lt(flow_error(sim), 1e-12)
  expect_identical(colSums(sim$releases), c(E = 11, P = 3, M = 120, U = 7))
  expect_identical(unname(sim$releases[5, "M"]), 40)
})

test_that("releases are added at the start of their day", {
  rel <- sg_releases("aa", "pupa", times = 3, number = 100)
  m <- sg_model(p, cube3, wx60, rel)
  sim <- simulate_ode(m, init = sg_init(m, values = c(E = 0)), output_every = 1)
  pupae_aa <- rowSums(sim$state[, sg_state_index(m, "P", genotype = "aa")])
  expect_identical(pupae_aa[1:3], c(0, 0, 0)) # days 0 to 2
  # end of day 3: 100 released at its start, minus that day's mortality and ageing out
  expect_gt(pupae_aa[4], 90)
  expect_lt(pupae_aa[4], 100)
})

test_that("an empty population stays empty", {
  sim <- simulate_ode(m3, init = sg_init(m3, values = c(E = 0)))
  expect_identical(sum(abs(sim$state)), 0)
  expect_identical(sum(abs(sim$counters)), 0)
})

test_that("results do not depend on the integration method", {
  a <- simulate_ode(m3, method = "lsoda", rtol = 1e-10, atol = 1e-12)
  b <- simulate_ode(m3, method = "radau", rtol = 1e-10, atol = 1e-12)
  rel <- max(abs(a$state - b$state)) / max(abs(a$state))
  expect_lt(rel, 1e-7)
  expect_true(is.logical(a$diagnostics$stiff) && !anyNA(a$diagnostics$stiff))
  expect_true(all(is.na(b$diagnostics$stiff)))
})

test_that("states stay non-negative up to round-off", {
  sim <- simulate_ode(m3)
  expect_gt(min(sim$diagnostics$min_state), -1e-8)
  expect_true(all(sim$diagnostics$max_rate > 0))
  expect_true(all(sim$diagnostics$steps > 0))
})

test_that("output_every thins the stored states but keeps daily counters", {
  sim <- simulate_ode(m3, output_every = 7)
  expect_identical(sim$day, c(0L, seq(7L, 56L, by = 7L), 60L))
  expect_identical(nrow(sim$counters), 60L)
  full <- simulate_ode(m3)
  expect_equal(sim$state[nrow(sim$state), ], full$state[61, ], tolerance = sg_tol_exact)
})

test_that("init is checked against the model", {
  m_other <- sg_model(sg_params(R = 50), cube3, wx60)
  expect_error(simulate_ode(m3, init = sg_init(m_other)), "different model")
  x <- as.numeric(sg_init(m3))
  sim <- simulate_ode(m3, init = x)
  expect_identical(unname(sim$state[1, ]), x)
  expect_error(simulate_ode(m3, init = x[-1]), "length")
})

test_that("invalid arguments give informative errors", {
  expect_error(simulate_ode(list()), "sg_model object")
  expect_error(simulate_ode(m3, method = "dopri5"), "C\\+\\+ backend's own integrator")
  expect_error(simulate_ode(m3, method = 1), "single string")
  expect_error(simulate_ode(m3, stiff_fallback = NA), "TRUE or FALSE")
  expect_error(simulate_ode(m3, backend = "julia"), "must be one of")
  expect_error(simulate_ode(m3, rtol = -1), "positive")
  expect_error(simulate_ode(m3, output_every = 0), "whole number >= 1")
})

test_that("print() summarises the simulation", {
  sim <- simulate_ode(m3)
  out <- capture.output(print(sim))
  expect_match(out[1], "^<sg_sim> 60 days \\(1971-05-01 to 1971-06-29\\), backend \"R\"")
  expect_true(any(grepl("final totals: E", out)))
  expect_true(any(grepl("steps/day", out)))
})
