# Stress scenarios for both backends (M2 reference, M3 C++). They pin down
# accuracy and stiffness behaviour in demanding regimes:
#   * very large sterile-male releases (10x to 1000x the wild males);
#   * larval crowding with a very large release of fertile eggs (a stiffness
#     source that grows with population size);
#   * a very large mating rate (artificial stiffness);
#   * many egg sub-stages (fast transfer rates).
# In each scenario two different integrators must agree (the backend's own
# explicit one, lsoda for R and Dormand-Prince for C++, and the implicit
# radau), and states must stay non-negative up to round-off.

skip_if_not_installed("deSolve")
skip_on_cran()

period <- c("1971-05-01", "1971-07-29")
wx_all <- sg_example_weather()

sit_raw <- sg_test_cube_sit(residual = 0)
sit_raw$eta[, "S"] <- 0.35
sit_cube <- sg_cube(sit_raw)

# Largest relative difference between two runs in a set of states, ignoring
# values below `floor_rel` of their maximum (deep extinction).
rel_diff <- function(a, b, cols, floor_rel = 1e-9) {
  xa <- rowSums(a$state[, cols, drop = FALSE])
  xb <- rowSums(b$state[, cols, drop = FALSE])
  max(abs(xa - xb) / pmax(xb, floor_rel * max(xb)))
}
both_methods <- function(m, backend, init = NULL) {
  main <- if (backend == "R") "lsoda" else "dopri5"
  list(
    main = simulate_ode(m, init, backend, method = main, rtol = 1e-9, atol = 1e-12),
    radau = simulate_ode(m, init, backend, method = "radau", rtol = 1e-9, atol = 1e-12)
  )
}

sg_test_each_backend("very large sterile-male releases add no stiffness and stay accurate", {
  p <- sg_params()
  wx <- sg_forcing(wx_all, p, period[1], period[2])
  cube <- sit_cube
  m0 <- sg_model(p, cube, wx)
  control <- sim_ode(m0, rtol = 1e-9, atol = 1e-12)
  wm_idx <- sg_state_index(m0, "M", genotype = "W")
  wild_males <- mean(rowSums(control$state[, wm_idx]))
  for (k in c(10, 100, 1000)) {
    rel <- sg_releases("S", "adult", "male", times = seq(1, 90, by = 7), number = k * wild_males)
    m <- sg_model(p, cube, wx, rel)
    wild <- c(
      sg_state_index(m, "E"), sg_state_index(m, "L"), sg_state_index(m, "P"),
      sg_state_index(m, "M", genotype = "W"), sg_state_index(m, "U", genotype = "W"),
      sg_state_index(m, "F", genotype = "W", mate = "W")
    )
    s <- both_methods(m, backend)
    expect_lt(rel_diff(s$main, s$radau, wild), 1e-6)
    expect_identical(sum(s$main$diagnostics$stiff), 0L) # release size adds no stiffness
    expect_identical(sum(s$main$diagnostics$fallback), 0L)
    expect_equal(max(s$main$diagnostics$max_rate), max(control$diagnostics$max_rate))
    expect_gt(min(s$main$diagnostics$min_state), -1e-8)
    # releases suppress the wild population
    expect_lt(sum(s$main$state[91, wild]), sum(control$state[91, wild]))
  }
})

sg_test_each_backend("larval crowding with a huge fertile egg release stays accurate", {
  p <- sg_params(larval_dd = "per_wound", gamma_L = 0.01)
  wx <- sg_forcing(wx_all, p, period[1], period[2])
  m <- sg_model(
    p, sg_cube(sg_test_cube_wild()), wx,
    sg_releases("W", "egg", times = 10, number = 1e6)
  )
  s <- both_methods(m, backend)
  lar <- sg_state_index(m, "L")
  expect_lt(rel_diff(s$main, s$radau, lar), 1e-6)
  expect_lt(rel_diff(s$main, s$radau, seq_len(m$n_states)), 1e-6)
  expect_gt(min(s$main$diagnostics$min_state), -1e-8)
  # crowding makes the largest rate grow with larval density
  m0 <- sg_model(p, sg_cube(sg_test_cube_wild()), wx)
  no_release <- sim_ode(m0, rtol = 1e-9, atol = 1e-12)
  expect_gt(max(s$main$diagnostics$max_rate), 3 * max(no_release$diagnostics$max_rate))
})

sg_test_each_backend("a very large mating rate stays accurate", {
  p <- sg_params(nu = 1000)
  wx <- sg_forcing(wx_all, p, "1971-05-01", "1971-05-30")
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx)
  s <- both_methods(m, backend)
  expect_lt(rel_diff(s$main, s$radau, seq_len(m$n_states)), 1e-6)
  expect_gt(min(s$main$diagnostics$min_state), -1e-8)
})

sg_test_each_backend("many egg sub-stages (fast transfer rates) stay accurate", {
  p <- sg_params(n_E = 25, n_L = 25)
  wx <- sg_forcing(wx_all, p, period[1], period[2])
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), wx)
  s <- both_methods(m, backend)
  expect_lt(rel_diff(s$main, s$radau, seq_len(m$n_states)), 1e-6)
  expect_gt(max(s$main$diagnostics$max_rate), 25 * sg_degree_days(20) / 10)
})
