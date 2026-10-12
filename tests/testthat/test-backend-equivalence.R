# M3 milestone check: the C++ backend (Dormand-Prince day loop) and the R
# reference backend (lsoda) give the same trajectories and daily flows, to
# within the integration error, across options, cube types and releases.

skip_if_not_installed("deSolve")

wx_all <- sg_example_weather()
period <- c("1971-05-01", "1971-06-14")
rtol <- 1e-10
atol <- 1e-12

sit_raw <- sg_test_cube_sit(residual = 0)
sit_raw$eta[, "S"] <- 0.35
weekly <- seq(1, 45, by = 7)

scenarios <- list(
  "Hardy-Weinberg Mendelian" = list(
    p = sg_params(), cube = sg_cube(sg_test_cube_mendelian()),
    init = function(m) sg_init(m, "hardy_weinberg", alleles = c(A = 0.3, a = 0.7))
  ),
  "SIT, weekly sterile-male releases" = list(
    p = sg_params(), cube = sg_cube(sit_raw),
    releases = sg_releases("S", "adult", "male", times = weekly, number = 200)
  ),
  "instant mating, releases of every kind" = list(
    p = sg_params(mating = "instant"), cube = sg_cube(sg_test_cube_mendelian()),
    releases = c(
      sg_releases("aa", "adult", "both", times = c(3, 17), number = 30),
      sg_releases("Aa", "pupa", times = 9, number = 12),
      sg_releases("AA", "egg", times = 21, number = 50)
    )
  ),
  "all adults search, Allee effect" = list(
    p = sg_params(search_denominator = "all_adults", h_M = 2), cube = sg_cube(sit_raw),
    releases = sg_releases("S", "adult", "male", times = weekly, number = 100)
  ),
  "larval crowding, large egg release" = list(
    p = sg_params(larval_dd = "per_wound", gamma_L = 0.01),
    cube = sg_cube(sg_test_cube_wild()),
    releases = sg_releases("W", "egg", times = 10, number = 1e5)
  ),
  "genotype lifespan costs" = list(
    p = sg_params(X_rec = 24, X_mat = 15),
    cube = sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25)),
    init = function(m) sg_init(m, genotypes = c(W = 0.5, S = 0.5))
  ),
  "gp2014 preset" = list(p = sg_params(preset = "gp2014"), cube = sg_cube(sg_test_cube_wild())),
  "custom sub-stages, background and moisture mortality" = list(
    p = sg_params(n_E = 1, n_L = 2, n_P = 3, n_A = 4, mu_b_F = 0.1, mu_W = 0.2),
    cube = sg_cube(sg_test_cube_mendelian())
  )
)
if (requireNamespace("MGDrivE", quietly = TRUE)) {
  scenarios[["homing drive with releases"]] <- list(
    p = sg_params(), cube = sg_cube(MGDrivE::cubeHomingDrive()),
    releases = sg_releases("HH", "adult", "male", times = c(1, 15), number = 20)
  )
  scenarios[["sex-linked X-shredder"]] <- list(
    p = sg_params(egg_temperature = "host"), cube = sg_cube(MGDrivE::cubeXShredderY())
  )
}

for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  test_that(paste("C++ and R backends agree:", nm), {
    wx <- sg_forcing(wx_all, sc$p, period[1], period[2])
    m <- sg_model(sc$p, sc$cube, wx, sc$releases)
    init <- if (is.null(sc$init)) sg_init(m) else sc$init(m)
    r <- simulate_ode(m, init, backend = "R", rtol = rtol, atol = atol)
    cc <- simulate_ode(m, init, backend = "cpp", rtol = rtol, atol = atol)
    expect_lt(max(abs(r$state - cc$state)) / max(abs(r$state)), 1e-8)
    expect_lt(max(abs(r$counters - cc$counters)) / max(abs(r$counters)), 1e-8)
    expect_identical(cc$releases, r$releases)
    expect_equal(cc$diagnostics$max_rate, r$diagnostics$max_rate, tolerance = 1e-8)
    expect_identical(cc$meta$fallback_days, 0L)
  })
}
