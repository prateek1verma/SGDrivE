# C++ backend of simulate_ode() (D-021): the C++ day loop with Dormand-Prince,
# its lsoda fallback on stiff days and C++ rates driven by deSolve methods,
# all checked against the R reference backend.

skip_if_not_installed("deSolve")

p <- sg_params()
wx60 <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-06-29")
cube3 <- sg_cube(sg_test_cube_mendelian())
rel <- c(
  sg_releases("aa", "adult", "male", times = c(5, 20, 35), number = 40),
  sg_releases("aa", "adult", "female", times = 10, number = 7),
  sg_releases("Aa", "pupa", times = 15, number = 3),
  sg_releases("aa", "egg", times = 25, number = 11)
)
m_rel <- sg_model(p, cube3, wx60, rel)

rel_diff <- function(a, b) max(abs(a$state - b$state)) / max(abs(a$state))

flow_error <- function(sim) {
  tot <- rowSums(sim$state)
  k <- sim$counters
  births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
  deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
  pred <- rowSums(sim$releases) + births - deaths - k[, "emergence_deaths"]
  max(abs(diff(tot) - pred)) / max(tot)
}

ref <- simulate_ode(m_rel, backend = "R", rtol = 1e-10, atol = 1e-12)
cpp <- simulate_ode(m_rel, backend = "cpp", rtol = 1e-10, atol = 1e-12)

test_that("the C++ backend returns the same structure as the R backend", {
  expect_s3_class(cpp, "sg_sim")
  expect_identical(cpp$day, ref$day)
  expect_identical(cpp$date, ref$date)
  expect_identical(dimnames(cpp$state), dimnames(ref$state))
  expect_identical(dimnames(cpp$counters), dimnames(ref$counters))
  expect_identical(cpp$releases, ref$releases)
  expect_identical(names(cpp$diagnostics), names(ref$diagnostics))
  expect_identical(cpp$state[1, ], ref$state[1, ])
  expect_identical(cpp$meta$backend, "cpp")
  expect_identical(cpp$meta$method, "dopri5")
  expect_true(cpp$meta$stiff_fallback)
  expect_identical(cpp$meta$fallback_days, 0L)
  expect_false(ref$meta$stiff_fallback) # not used by deSolve methods
})

test_that("the C++ backend agrees with the R backend, with releases", {
  expect_lt(rel_diff(ref, cpp), 1e-8)
  expect_lt(max(abs(ref$counters - cpp$counters)) / max(abs(ref$counters)), 1e-8)
  expect_lt(flow_error(cpp), 1e-12)
  d <- cpp$diagnostics
  expect_true(all(d$steps > 0) && all(d$rejected >= 0) && all(d$evaluations > d$steps))
  expect_false(any(d$stiff) || any(d$fallback))
  expect_equal(d$max_rate, ref$diagnostics$max_rate, tolerance = 1e-8)
  expect_gt(min(d$min_state), -1e-8)
})

test_that("C++ rates driven by a deSolve method reproduce the R backend", {
  a <- simulate_ode(m_rel, backend = "R", method = "lsoda")
  b <- simulate_ode(m_rel, backend = "cpp", method = "lsoda")
  # same integrator and rates equal to round-off, but lsoda's step choices
  # react to round-off, so agreement is at the level of the integration error
  expect_lt(rel_diff(a, b), 1e-7)
  expect_identical(b$diagnostics$stiff, a$diagnostics$stiff)
  r <- simulate_ode(m_rel, backend = "cpp", method = "radau", rtol = 1e-10, atol = 1e-12)
  expect_lt(rel_diff(ref, r), 1e-7)
  expect_true(all(is.na(r$diagnostics$stiff)))
  expect_true(any(grepl("stiffness not reported", capture.output(print(r)))))
})

test_that("default tolerances keep the backends within the integration error", {
  a <- simulate_ode(m_rel, backend = "R")
  b <- simulate_ode(m_rel, backend = "cpp")
  expect_lt(rel_diff(a, b), 1e-6)
  expect_lt(rel_diff(ref, b), 1e-6)
})

test_that("output_every and an empty population work in the C++ backend", {
  sim <- simulate_ode(m_rel, backend = "cpp", output_every = 7, rtol = 1e-10, atol = 1e-12)
  expect_identical(sim$day, c(0L, seq(7L, 56L, by = 7L), 60L))
  expect_identical(sim$state, cpp$state[sim$day + 1L, ])
  expect_identical(sim$counters, cpp$counters)
  m <- sg_model(p, cube3, wx60)
  empty <- simulate_ode(m, init = sg_init(m, values = c(E = 0)), backend = "cpp")
  expect_identical(sum(abs(empty$state)), 0)
  expect_identical(sum(abs(empty$counters)), 0)
})

test_that("stiff days fall back to lsoda and give the same result", {
  # near-instant mating (nu = 1e4 per day) makes every day stiff
  ps <- sg_params(mating_rate = 1e4)
  ms <- sg_model(ps, cube3, sg_forcing(sg_example_weather(), ps, "1971-05-01", "1971-05-10"))
  a <- simulate_ode(ms, backend = "R", rtol = 1e-9, atol = 1e-12)
  b <- simulate_ode(ms, backend = "cpp", rtol = 1e-9, atol = 1e-12)
  expect_identical(b$meta$fallback_days, 10L)
  expect_true(all(b$diagnostics$fallback & b$diagnostics$stiff))
  expect_true(all(is.na(b$diagnostics$rejected)))
  expect_lt(rel_diff(a, b), 1e-8)
  out <- capture.output(print(b))
  expect_true(any(grepl("stiffness detected on 10 of 10 days, 10 integrated by the lsoda", out)))
  # without the fallback, Dormand-Prince gets there too, with many more steps
  nf <- simulate_ode(ms, backend = "cpp", rtol = 1e-9, atol = 1e-12, stiff_fallback = FALSE)
  expect_identical(nf$meta$fallback_days, 0L)
  expect_false(nf$meta$stiff_fallback)
  expect_true(all(nf$diagnostics$stiff))
  expect_lt(rel_diff(a, nf), 1e-7)
  expect_gt(mean(nf$diagnostics$steps), 5 * mean(b$diagnostics$steps))
  expect_lt(flow_error(b), 1e-12)
})

test_that("a stiff day in the middle of a run is handed over and the loop resumes", {
  # with an Allee effect, the effective mating rate nu * M / (h_M + M) is only
  # large (stiff) while a big release of males is alive and mature
  pc <- sg_params(mating_rate = 1e5, male_density_half_mating = 1e5)
  wx <- sg_forcing(sg_example_weather(), pc, "1971-05-01", "1971-05-30")
  mc <- sg_model(pc, cube3, wx, sg_releases("aa", "adult", "male", times = 10, number = 1e6))
  a <- simulate_ode(mc, backend = "R", rtol = 1e-9, atol = 1e-12)
  b <- simulate_ode(mc, backend = "cpp", rtol = 1e-9, atol = 1e-12, output_every = 3)
  fb <- which(b$diagnostics$fallback)
  expect_gt(length(fb), 0L)
  expect_lt(length(fb), 20L)
  expect_true(all(fb > 10))
  expect_identical(b$diagnostics$stiff, b$diagnostics$fallback)
  expect_lt(max(abs(a$state[b$day + 1L, ] - b$state)) / max(a$state), 1e-7)
  expect_lt(max(abs(a$counters - b$counters)) / max(abs(a$counters)), 1e-7)
})

test_that("a Dormand-Prince failure without fallback is an informative error", {
  local_mocked_bindings(.sg_dopri5_max_steps = 2)
  expect_error(
    simulate_ode(m_rel, backend = "cpp", stiff_fallback = FALSE),
    "failed on day 1 \\(1971-05-01\\): too many steps"
  )
  # with the fallback, lsoda takes over every day
  sim <- simulate_ode(m_rel, backend = "cpp")
  expect_identical(sim$meta$fallback_days, 60L)
  expect_false(any(sim$diagnostics$stiff))
  expect_lt(rel_diff(simulate_ode(m_rel, backend = "R"), sim), 1e-7)
})

test_that("the C++ backend needs deSolve only for a fallback day or a deSolve method", {
  local_mocked_bindings(.sg_has_desolve = function() FALSE)
  sim <- simulate_ode(m_rel)
  expect_identical(sim$meta$fallback_days, 0L)
  expect_identical(sim$meta$desolve_version, NA_character_)
  expect_error(simulate_ode(m_rel, method = "lsoda"), "Method \"lsoda\" needs the 'deSolve'")
  expect_error(simulate_ode(m_rel, backend = "R"), "Method \"lsoda\" needs the 'deSolve'")
  ps <- sg_params(mating_rate = 1e4)
  ms <- sg_model(ps, cube3, sg_forcing(sg_example_weather(), ps, "1971-05-01", "1971-05-03"))
  expect_error(
    simulate_ode(ms),
    "Day 1 \\(1971-05-01\\) is stiff; its fallback needs the 'deSolve' .*stiff_fallback = FALSE"
  )
  expect_identical(simulate_ode(ms, stiff_fallback = FALSE)$meta$fallback_days, 0L)
})

test_that("the C++ day loop can stop and resume without changing the result", {
  v <- SGDrivE:::.sg_model_view(m_rel)
  ptr <- SGDrivE:::sg_kernel_new_cpp(v)
  run <- SGDrivE:::sg_ode_run_cpp
  x0 <- as.numeric(sg_init(m_rel))
  rec <- rep(TRUE, 60)
  full <- run(ptr, x0, 1L, 60L, 1e-8, 1e-10, rec, TRUE, 1e5)
  expect_identical(full$completed, 60L)
  expect_identical(full$stopped_day, 0L)
  expect_identical(full$state_days, 1:60)
  first <- run(ptr, x0, 1L, 30L, 1e-8, 1e-10, rec, TRUE, 1e5)
  second <- run(ptr, first$x, 31L, 60L, 1e-8, 1e-10, rec, TRUE, 1e5)
  expect_identical(rbind(first$states, second$states), full$states)
  expect_identical(rbind(first$counters, second$counters), full$counters)
  # a stopped run returns the state at the start of the failed day, after its releases
  stop5 <- run(ptr, full$states[4, ], 5L, 60L, 1e-8, 1e-10, rec, TRUE, 2)
  expect_identical(stop5$stopped_day, 5L)
  expect_identical(stop5$stop_status, 1L)
  expect_identical(stop5$completed, 0L)
  expect_identical(nrow(stop5$counters), 0L)
  m_idx <- sg_state_index(m_rel, "M", genotype = "aa", substage = 1)
  expect_equal(stop5$x[m_idx] - full$states[4, m_idx], 40)
  # single rate evaluations for deSolve
  y <- c(x0, numeric(length(attr(v, "counter_names"))))
  r_cpp <- SGDrivE:::sg_kernel_rates_cpp(ptr, 3L, y)
  r_ref <- SGDrivE:::sg_rates_cpp(v, 3L, x0)
  expect_identical(r_cpp, c(r_ref$deriv, r_ref$counters))
  expect_identical(
    SGDrivE:::sg_kernel_max_rate_cpp(ptr, 3L, x0), SGDrivE:::sg_max_rate_cpp(v, 3L, x0)
  )
})

test_that("the C++ entry points check their input", {
  v <- SGDrivE:::.sg_model_view(m_rel)
  ptr <- SGDrivE:::sg_kernel_new_cpp(v)
  run <- SGDrivE:::sg_ode_run_cpp
  x0 <- as.numeric(sg_init(m_rel))
  rec <- rep(TRUE, 60)
  expect_error(run(ptr, x0[-1], 1L, 60L, 1e-8, 1e-10, rec, TRUE, 1e5), "wrong length")
  expect_error(run(ptr, x0, 1L, 60L, 1e-8, 1e-10, rec[-1], TRUE, 1e5), "one value per day")
  expect_error(run(ptr, x0, 0L, 60L, 1e-8, 1e-10, rec, TRUE, 1e5), "outside the forcing")
  expect_error(run(ptr, x0, 1L, 61L, 1e-8, 1e-10, rec, TRUE, 1e5), "outside the forcing")
  expect_error(run(ptr, x0, 30L, 20L, 1e-8, 1e-10, rec, TRUE, 1e5), "before day_from")
  expect_error(SGDrivE:::sg_kernel_rates_cpp(ptr, 1L, numeric(2)), "wrong length")
  expect_error(SGDrivE:::sg_kernel_rates_cpp(ptr, 61L, x0), "outside the forcing")
  expect_error(SGDrivE:::sg_kernel_max_rate_cpp(ptr, 1L, numeric(2)), "wrong length")
})
