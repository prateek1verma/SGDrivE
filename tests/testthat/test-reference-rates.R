rates <- ScrewDrivE:::.sg_rates_r
setup <- ScrewDrivE:::.sg_reference_setup

wx_all <- sg_example_weather()
day_warm <- which(wx_all$date == as.Date("1971-07-15"))

cube_mendel <- sg_cube(sg_test_cube_mendelian())

mk <- function(p = sg_params(), cube = cube_mendel, ...) {
  wx <- sg_forcing(wx_all, p)
  sg_model(p, cube, wx, ...)
}
rand_state <- function(m, seed = 1) {
  set.seed(seed)
  x <- stats::runif(m$n_states, 0, 10)
  # keep adult classes that genotypes never reach empty, as in a simulation
  for (g in m$cube$genotypes) {
    lf <- m$adult$female$last_class[[g]]
    lm <- m$adult$male$last_class[[g]]
    if (lf < m$dims$n_A) {
      unreached <- (lf + 1):m$dims$n_A
      x[sg_state_index(m, "U", genotype = g, substage = unreached)] <- 0
      x[sg_state_index(m, "F", genotype = g, substage = unreached)] <- 0
    }
    if (lm < m$dims$n_A) {
      x[sg_state_index(m, "M", genotype = g, substage = (lm + 1):m$dims$n_A)] <- 0
    }
  }
  x
}
block_sum <- function(m, v, b) sum(v[sg_state_index(m, b)])

test_that("an empty population has zero rates", {
  m <- mk()
  r <- rates(day_warm, numeric(m$n_states), setup(m))
  expect_identical(r$deriv, numeric(m$n_states))
  expect_true(all(r$counters == 0))
  expect_identical(names(r$counters), ScrewDrivE:::.sg_counter_names(m$cube$genotypes))
})

test_that("flow accounting closes for the total and for each block", {
  variants <- list(
    sg_params(), sg_params(mating_mode = "instant"),
    sg_params(female_background_mortality = 0.1, male_density_half_mating = 3)
  )
  for (p in variants) {
    m <- mk(p)
    x <- rand_state(m)
    r <- rates(day_warm, x, setup(m))
    k <- r$counters
    births <- sum(k[grepl("^eggs_laid", names(k))])
    deaths <- sum(k[grepl("^deaths_|^senescence_", names(k))])
    expect_equal(sum(r$deriv), births - deaths - k[["emergence_deaths"]],
      tolerance = sg_tol_exact
    )
    expect_equal(block_sum(m, r$deriv, "E"), births - k[["larvae_in"]] - k[["deaths_E"]],
      tolerance = sg_tol_exact
    )
    expect_equal(block_sum(m, r$deriv, "L"), k[["larvae_in"]] - k[["pupae_in"]] - k[["deaths_L"]],
      tolerance = sg_tol_exact
    )
    emerged <- k[["emerged_male"]] + k[["emerged_female"]] + k[["emergence_deaths"]]
    expect_equal(block_sum(m, r$deriv, "P"), k[["pupae_in"]] - emerged - k[["deaths_P"]],
      tolerance = sg_tol_exact
    )
    m_flow <- k[["emerged_male"]] - k[["deaths_M"]] - k[["senescence_M"]]
    expect_equal(block_sum(m, r$deriv, "M"), m_flow, tolerance = sg_tol_exact)
    matings <- sum(k[grepl("^matings", names(k))])
    expect_equal(block_sum(m, r$deriv, "U") + block_sum(m, r$deriv, "F"),
      k[["emerged_female"]] - k[["deaths_U"]] - k[["deaths_F"]] - k[["senescence_U"]] -
        k[["senescence_F"]],
      tolerance = sg_tol_exact
    )
    if (p$mating_mode == "rate") {
      expect_equal(block_sum(m, r$deriv, "F"), matings - k[["deaths_F"]] - k[["senescence_F"]],
        tolerance = sg_tol_exact
      )
    }
  }
})

test_that("egg production follows S8-S9 (wound search, egg tensor)", {
  p <- sg_params()
  m <- mk(p)
  x <- rand_state(m)
  s <- setup(m)
  r <- rates(day_warm, x, s, aux = TRUE)
  f <- m$forcing
  laying <- m$adult$female$flags$X_ovi
  w <- matrix(0, 3, 3, dimnames = list(m$cube$genotypes, m$cube$genotypes)) # [g, m]
  for (g in m$cube$genotypes) {
    for (mm in m$cube$genotypes) {
      w[g, mm] <- sum(x[sg_state_index(m, "F", genotype = g, mate = mm, substage = which(laying))])
    }
  }
  n_s <- sum(w)
  psi <- 1 - exp(-p$wound_search_rate * f$dx[day_warm] * p$wound_density / n_s)
  expect_equal(r$aux$psi, psi, tolerance = sg_tol_exact)
  births <- vapply(m$cube$genotypes, function(k) {
    p$eggs_per_day_max * f$phi_ov[day_warm] * psi *
      sum(w * m$cube$ih[, , k] * m$cube$tau[, , k] * m$cube$s)
  }, 0)
  expect_equal(unname(r$aux$births), unname(births), tolerance = sg_tol_exact)
  expect_equal(r$counters[["egg_batches"]], psi * n_s, tolerance = sg_tol_exact)
})

test_that("sterile-mated females lay no viable eggs but compete for wounds (D-006)", {
  m <- mk(cube = sg_cube(sg_test_cube_sit(residual = 0)))
  x <- numeric(m$n_states)
  lay <- which(m$adult$female$flags$X_ovi)
  x[sg_state_index(m, "F", genotype = "W", mate = "W", substage = lay)] <- 1
  r1 <- rates(day_warm, x, setup(m), aux = TRUE)
  x[sg_state_index(m, "F", genotype = "W", mate = "S", substage = lay)] <- 9
  r2 <- rates(day_warm, x, setup(m), aux = TRUE)
  expect_equal(r2$aux$n_s, 10 * r1$aux$n_s)
  expect_lt(sum(r2$aux$births), sum(r1$aux$births)) # wound competition
  expect_identical(r2$aux$births[["S"]], 0)
  # with all adults as competitors, males count too
  pa <- sg_params(wound_competitors = "all_adults")
  ma <- mk(pa, cube = sg_cube(sg_test_cube_sit()))
  x2 <- x
  x2[sg_state_index(ma, "M", genotype = "S")] <- 5
  ra <- rates(day_warm, x2, setup(ma), aux = TRUE)
  adults <- c(sg_state_index(ma, "M"), sg_state_index(ma, "U"), sg_state_index(ma, "F"))
  expect_equal(ra$aux$n_s, sum(x2[adults]))
})

test_that("females search and lay batches in proportion to their fecundity (D-031)", {
  # S females: no eggs (s = 0), so they do not compete for wounds
  m <- mk(cube = sg_cube(sg_test_cube_sit(residual = 0, s_s = 0)))
  x <- numeric(m$n_states)
  lay <- which(m$adult$female$flags$X_ovi)
  x[sg_state_index(m, "F", genotype = "W", mate = "W", substage = lay)] <- 1
  r1 <- rates(day_warm, x, setup(m), aux = TRUE)
  x[sg_state_index(m, "F", genotype = "S", mate = "W", substage = lay)] <- 9
  x[sg_state_index(m, "F", genotype = "S", mate = "S", substage = lay)] <- 4
  r2 <- rates(day_warm, x, setup(m), aux = TRUE)
  expect_equal(r2$aux$n_s, r1$aux$n_s)
  expect_equal(r2$aux$psi, r1$aux$psi)
  expect_equal(r2$aux$births, r1$aux$births)
  expect_equal(r2$counters[["egg_batches"]], r1$counters[["egg_batches"]])
  expect_equal(r2$counters[["fertile_batches"]], r1$counters[["fertile_batches"]])

  # partly fertile S females (s = 0.4): weight 0.4 in the search and in the
  # batches, each batch with viability sum_k ih tau; mean eggs unchanged (S9)
  p <- sg_params(clutch_size = 50)
  m <- mk(p, cube = sg_cube(sg_test_cube_sit(residual = 0.3, s_s = 0.4)))
  x <- numeric(m$n_states)
  x[sg_state_index(m, "F", genotype = "W", mate = "W", substage = lay)] <- 2
  x[sg_state_index(m, "F", genotype = "S", mate = "S", substage = lay)] <- 3
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  w_ww <- 2 * length(lay)
  w_ss <- 3 * length(lay)
  n_s <- w_ww + 0.4 * w_ss
  f <- m$forcing
  psi <- 1 - exp(-p$wound_search_rate * f$dx[day_warm] * m$counts$H / n_s)
  expect_equal(r$aux$n_s, n_s, tolerance = sg_tol_exact)
  expect_equal(r$aux$psi, psi, tolerance = sg_tol_exact)
  factor <- p$eggs_per_day_max * f$phi_ov[day_warm] / 50
  expect_equal(r$counters[["egg_batches"]], psi * n_s * factor, tolerance = sg_tol_exact)
  fertile <- psi * factor * (w_ww * -expm1(-50) + 0.4 * w_ss * -expm1(-50 * 0.3))
  expect_equal(r$counters[["fertile_batches"]], fertile, tolerance = sg_tol_exact)
  eggs <- p$eggs_per_day_max * f$phi_ov[day_warm] * psi * (w_ww + 0.4 * 0.3 * w_ss)
  expect_equal(sum(r$aux$births), eggs, tolerance = sg_tol_exact)
})

test_that("mating shares follow eta-weighted male availability (S12)", {
  raw <- sg_test_cube_sit()
  raw$eta[, "S"] <- 0.35
  m <- mk(cube = sg_cube(raw))
  x <- numeric(m$n_states)
  x[sg_state_index(m, "M", genotype = "W")] <- 1 # 10 wild males
  x[sg_state_index(m, "M", genotype = "S")] <- 5 # 50 sterile males
  x[sg_state_index(m, "U", genotype = "W")] <- 2
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  k <- r$counters
  share_s <- k[["matings[W|S]"]] / (k[["matings[W|S]"]] + k[["matings[W|W]"]])
  expect_equal(share_s, 0.35 * 50 / (0.35 * 50 + 10), tolerance = sg_tol_exact)
  # total mating rate = nu * virgins (h_M = 0, all receptive)
  expect_equal(k[["matings[W|S]"]] + k[["matings[W|W]"]], log(2) * 20, tolerance = sg_tol_exact)
})

test_that("no males means no mating; instant mating sends females to U then", {
  m <- mk()
  x <- numeric(m$n_states)
  x[sg_state_index(m, "U", genotype = "AA")] <- 1
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  expect_true(all(r$counters[grepl("^matings", names(r$counters))] == 0))
  expect_identical(r$aux$g_find, 0)
  pi <- sg_params(mating_mode = "instant")
  mi <- mk(pi)
  x <- numeric(mi$n_states)
  x[sg_state_index(mi, "P", genotype = "AA", substage = 10)] <- 10
  r <- rates(day_warm, x, setup(mi))
  expect_gt(block_sum(mi, r$deriv, "U"), 0)
  expect_identical(block_sum(mi, r$deriv, "F"), 0)
  x[sg_state_index(mi, "M", genotype = "AA")] <- 1
  r <- rates(day_warm, x, setup(mi))
  expect_identical(block_sum(mi, r$deriv, "U"), -sum(0)) # nobody in U
  expect_gt(r$deriv[sg_state_index(mi, "F", genotype = "AA", mate = "AA", substage = 1)], 0)
})

test_that("receptivity, male maturity and the Allee term act as specified", {
  pr <- sg_params(female_receptive_age_dd = 24, male_mature_age_dd = 30)
  m <- mk(pr)
  x <- numeric(m$n_states)
  x[sg_state_index(m, "U", genotype = "AA")] <- 1
  x[sg_state_index(m, "M", genotype = "AA")] <- 1
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  first_rec <- which(m$adult$female$flags$X_rec)[1]
  expect_true(all(r$aux$mating_flux[seq_len(first_rec - 1), ] == 0))
  expect_true(all(r$aux$mating_flux[first_rec:10, 1] > 0))
  mature <- m$adult$male$flags$X_mat
  expect_equal(unname(r$aux$m_bar[["AA"]]), sum(mature))
  ph <- sg_params(male_density_half_mating = 4)
  mh <- mk(ph)
  rh <- rates(day_warm, x, setup(mh), aux = TRUE)
  expect_equal(rh$aux$g_find, 10 / (4 + 10), tolerance = sg_tol_exact)
})

test_that("genotypes leave the adult chain at their last reached class (D-018)", {
  cb <- sg_cube(sg_test_cube_sit(), lifespan_fraction = c(W = 1, S = 0.25))
  m <- mk(cube = cb)
  x <- numeric(m$n_states)
  x[sg_state_index(m, "U", genotype = "S", substage = 3)] <- 1
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  expect_identical(r$deriv[sg_state_index(m, "U", genotype = "S", substage = 4)], 0)
  expect_equal(r$counters[["senescence_U"]], unname(r$aux$q_f[3, "S"]), tolerance = sg_tol_exact)
  expect_equal(unname(r$aux$q_f[1, "S"]), m$forcing$dx[day_warm] / (50.7 / 3),
    tolerance = sg_tol_exact
  )
})

test_that("emergence splits by phi and xi", {
  raw <- sg_test_cube_sit()
  raw$phi <- c(W = 0.5, S = 0.2)
  raw$xiF <- c(W = 1, S = 0.5)
  m <- mk(cube = sg_cube(raw))
  x <- numeric(m$n_states)
  x[sg_state_index(m, "P", genotype = "S", substage = 10)] <- 10
  r <- rates(day_warm, x, setup(m), aux = TRUE)
  em <- r$aux$q[["P"]] * 10
  expect_equal(r$deriv[sg_state_index(m, "M", genotype = "S", substage = 1)], 0.8 * em,
    tolerance = sg_tol_exact
  )
  expect_equal(r$deriv[sg_state_index(m, "U", genotype = "S", substage = 1)], 0.2 * 0.5 * em,
    tolerance = sg_tol_exact
  )
  expect_equal(r$counters[["emergence_deaths"]], 0.2 * 0.5 * em, tolerance = sg_tol_exact)
})

test_that("optional terms change rates in the expected direction", {
  base <- mk()
  x <- rand_state(base)
  r0 <- rates(day_warm, x, setup(base), aux = TRUE)
  with_p <- function(...) {
    p <- sg_params(...)
    rates(day_warm, x, setup(mk(p)), aux = TRUE)
  }
  expect_gt(
    with_p(larval_crowding_mode = "per_wound", larval_crowding_coef = 0.01)$counters[["deaths_L"]],
    r0$counters[["deaths_L"]]
  )
  expect_gt(
    with_p(female_background_mortality = 0.1)$counters[["deaths_F"]],
    r0$counters[["deaths_F"]]
  )
  expect_gt(
    with_p(male_background_mortality = 0.1)$counters[["deaths_M"]],
    r0$counters[["deaths_M"]]
  )
  expect_equal(sum(with_p(wound_finding_max = 0.5)$aux$births), 0.5 * sum(r0$aux$births),
    tolerance = sg_tol_exact
  )
  host <- with_p(egg_temp_source = "host")
  expect_equal(host$aux$q[["E"]], 3 * sg_degree_days(38.5) / 10, tolerance = sg_tol_exact)
  # moisture: pupal mortality rises when annual rain is below the threshold
  dry <- with_p(pupa_drought_mortality = 0.2, drought_free_rainfall = 5000)
  expect_gt(dry$counters[["deaths_P"]], r0$counters[["deaths_P"]])
  expect_identical(
    with_p(pupa_drought_mortality = 0.2, drought_free_rainfall = 1)$counters[["deaths_P"]],
    r0$counters[["deaths_P"]]
  )
})

test_that("days outside the forcing are rejected", {
  m <- mk()
  expect_error(rates(0, numeric(m$n_states), setup(m)), "outside the forcing")
  expect_error(rates(m$forcing$n_days + 1, numeric(m$n_states), setup(m)), "outside the forcing")
})
