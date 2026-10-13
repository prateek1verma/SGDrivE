# M4 milestone check: statistical behaviour of the chain-binomial solver.
# One-step distributions against binomial, multinomial and compound-Poisson
# theory; the law of large numbers (replicate mean -> deterministic map ->
# ODE); neutral genetics and drift; eradication by sterile males; the effect
# of egg batches in small populations; stress cases. Fixed seeds, bounds of
# about 5 standard errors.

skip_if_not_installed("deSolve")

ns <- asNamespace("SGDrivE")
wx_all <- sg_example_weather()
p <- sg_params()
step_rep <- function(m, x, n, day = 10L) {
  r <- ns$sg_cb_step_rep_cpp(ns$.sg_model_view(m), day, 0.1, x, as.integer(n))
  colnames(r) <- c(sg_state_names(m), ns$.sg_counter_names(m$cube$genotypes))
  r
}
adults <- function(tot) apply(tot[, c("M", "U", "F"), , drop = FALSE], c(1, 3), sum)

test_that("one chain compartment moves as a multinomial with the exact chain probabilities", {
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), sg_forcing(wx_all, p, "1971-06-01", "1971-06-30"))
  s <- ns$.sg_reference_setup(m)
  info <- ns$.sg_cb_day_r(10, s, 0.1)
  n <- 1000
  for (block in c("P", "U")) {
    x <- numeric(m$n_states)
    i <- sg_state_index(m, block, substage = 2)
    x[i] <- n
    set.seed(1)
    r <- step_rep(m, x, 40000)
    ch <- if (block == "P") info$ch_p else info$ch_f[[1]]
    k <- ncol(ch$dest) - 1L
    prob <- (1 - ch$die[2]) * ch$dest[2, 2:k] # alive in sub-stages 2..k (no arrivals)
    obs <- r[, sg_state_index(m, block)[2:k]]
    expect_equal(colMeans(obs), n * prob, tolerance = 0.01, ignore_attr = TRUE)
    big <- prob > 0.01
    expect_equal(apply(obs[, big, drop = FALSE], 2, stats::var), n * prob[big] * (1 - prob[big]),
      tolerance = 0.06, ignore_attr = TRUE
    )
    # negative covariance of two destinations: -n p_j p_k
    jk <- which(big)[1:2]
    expect_equal(stats::cov(obs[, jk[1]], obs[, jk[2]]), -n * prob[jk[1]] * prob[jk[2]],
      tolerance = 0.12
    )
    deaths <- r[, paste0("deaths_", block)]
    expect_equal(mean(deaths), n * ch$die[2], tolerance = 0.03)
  }
})

test_that("matings with a male genotype are binomial thinning of the virgins", {
  raw <- sg_test_cube_sit(residual = 0)
  raw$eta[, "S"] <- 0.5
  m <- sg_model(p, sg_cube(raw), sg_forcing(wx_all, p, "1971-06-01", "1971-06-30"))
  x <- numeric(m$n_states)
  x[sg_state_index(m, "U", genotype = "W", substage = 1)] <- 2000
  x[sg_state_index(m, "M", genotype = "W", substage = 1)] <- 100
  x[sg_state_index(m, "M", genotype = "S", substage = 1)] <- 300
  set.seed(2)
  r <- step_rep(m, x, 40000)
  pi_s <- 0.5 * 300 / (0.5 * 300 + 100)
  pr <- (1 - exp(-p$nu * 0.1)) * pi_s
  y <- r[, "matings[W|S]"]
  expect_equal(mean(y), 2000 * pr, tolerance = 0.01)
  expect_equal(stats::var(y), 2000 * pr * (1 - pr), tolerance = 0.06)
})

test_that("fertile batches and their eggs: the refinement keeps compound-Poisson eggs", {
  # a cross with 10% viable eggs and clutches of 20: lambda = 2 viable eggs
  # per batch, so batches without a viable egg are common (exp(-2))
  raw <- sg_test_cube_sit(residual = 0.1)
  pp <- sg_params(clutch_size = 20)
  m <- sg_model(pp, sg_cube(raw), sg_forcing(wx_all, pp, "1971-06-01", "1971-06-30"), area = 1000)
  x <- numeric(m$n_states)
  x[sg_state_index(m, "F", genotype = "W", mate = "S", substage = 5:10)] <- 2000
  set.seed(10)
  r <- step_rep(m, x, 40000)
  batches <- r[, "egg_batches"]
  fertile <- r[, "fertile_batches"]
  eggs <- r[, "eggs_laid[W]"]
  lambda <- 20 * 0.1
  expect_equal(mean(fertile) / mean(batches), 1 - exp(-lambda), tolerance = 0.01)
  expect_equal(mean(eggs), mean(batches) * lambda, tolerance = 0.01)
  expect_equal(stats::var(eggs) / mean(eggs), 1 + lambda, tolerance = 0.05)
  expect_true(all(eggs >= fertile)) # every fertile batch has a viable egg
  expect_true(all(eggs[fertile == 0] == 0))
})

test_that("the replicate mean converges to the deterministic map and the ODE", {
  skip_on_cran()
  err <- vapply(c(10, 10000), function(area) {
    wx <- sg_forcing(wx_all, p, "1971-05-01", "1971-06-29")
    m <- sg_model(p, sg_cube(sg_test_cube_wild()), wx, area = area)
    init <- sg_init(m, "state", state = as.numeric(sg_init(m, "equilibrium", temperature = 27)))
    sim <- simulate_cb(m, init, n_rep = 40, seed = 1, keep = "totals")
    a <- adults(ns$.sg_cb_block_totals(sim))
    det_run <- ns$.sg_cb_run_cpp(m, init, stochastic = FALSE)
    det <- rowSums(ns$.sg_sim_block_totals(det_run)[, c("M", "U", "F")])
    ode <- rowSums(ns$.sg_sim_block_totals(simulate_ode(m, init))[, c("M", "U", "F")])
    se <- apply(a, 1, stats::sd) / sqrt(ncol(a))
    expect_lt(max((abs(rowMeans(a) - det) / se)[-1]), 5) # mean = map, within noise
    expect_lt(max(abs(det - ode)) / max(ode), 0.005) # map = ODE at dt = 0.1
    max(abs(rowMeans(a) - det)) / max(det)
  }, 0)
  expect_lt(err[2], 0.005)
  expect_gt(err[1] / err[2], 10) # sqrt(1000) ~ 32 times smaller noise
})

test_that("neutral alleles keep their mean frequency; drift grows in small populations", {
  skip_on_cran()
  freq_sd <- vapply(c(5, 500), function(area) {
    wx <- sg_forcing(wx_all, p, "1971-05-01", "1971-06-29")
    m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx, area = area)
    init <- sg_init(m, "hardy_weinberg",
      alleles = c(A = 0.3, a = 0.7),
      values = c(E = 20, L = 20, P = 20, M = 10, F = 10), per = "km2"
    )
    sim <- simulate_cb(m, init, n_rep = 40, seed = 3, keep = "totals")
    tt <- sim$totals[dim(sim$totals)[1], , ]
    g <- function(geno) colSums(tt[paste0(c("M", "U", "F"), "[", geno, "]"), , drop = FALSE])
    fa <- (g("Aa") + 2 * g("aa")) / (2 * (g("AA") + g("Aa") + g("aa")))
    expect_lt(abs(mean(fa) - 0.7) / (stats::sd(fa) / sqrt(length(fa))), 4)
    stats::sd(fa)
  }, 0)
  expect_gt(freq_sd[1], 3 * freq_sd[2]) # drift ~ 1 / sqrt(population size)
})

test_that("sterile-male releases eliminate the wild type; released males remain", {
  skip_on_cran()
  raw <- sg_test_cube_sit(residual = 0)
  raw$eta[, "S"] <- 0.5
  wx <- sg_forcing(wx_all, p, "1971-05-01", "1971-10-31")
  rel <- sg_releases("S", "adult", "male", times = seq(1, 180, 7), number = 100, per = "km2")
  m <- sg_model(p, sg_cube(raw), wx, rel, area = 20)
  m0 <- sg_model(p, sg_cube(raw), wx, area = 20)
  init <- sg_init(m0, "burn_in", weather = wx_all)
  sit <- simulate_cb(m, as.numeric(init), n_rep = 20, seed = 4, keep = "totals")
  ctl <- simulate_cb(m0, init, n_rep = 20, seed = 4, keep = "totals")
  expect_true(all(!is.na(sit$wild_extinct_day)))
  expect_true(all(is.na(sit$extinct_day))) # the sterile males are still there
  expect_true(all(is.na(ctl$wild_extinct_day)))
  last <- dim(sit$totals)[1]
  expect_identical(sum(sit$totals[last, grepl("\\[W\\]", dimnames(sit$totals)[[2]]), ]), 0)
  expect_true(any(grepl("wild type eliminated: 20 of 20", capture.output(print(sit)))))
})

test_that("egg batches make small populations more likely to die out", {
  skip_on_cran()
  ext <- vapply(c("eggs", "batches"), function(ov) {
    pp <- sg_params(oviposition = ov)
    # 3 km^2 through the winter low (January to April)
    wx <- sg_forcing(wx_all, pp, "1971-01-01", "1971-04-30")
    m <- sg_model(pp, sg_cube(sg_test_cube_wild()), wx, area = 3)
    init <- sg_init(m, "burn_in", weather = wx_all)
    sim <- simulate_cb(m, init, n_rep = 100, seed = 5, keep = "totals", output_every = 120)
    mean(!is.na(sim$extinct_day))
  }, 0)
  expect_lt(ext[["eggs"]], 0.06)
  expect_gt(ext[["batches"]], ext[["eggs"]] + 0.1)
})

test_that("stress: huge counts, step sizes, fast mating, long chains, many genotypes", {
  skip_on_cran()
  flow_ok <- function(r) {
    k <- r$counters
    births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
    deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
    change <- rowSums(r$releases) + births - deaths - k[, "emergence_deaths"]
    identical(max(abs(diff(rowSums(r$state)) - change)), 0)
  }
  wx <- sg_forcing(wx_all, p, "1971-05-01", "1971-05-30")
  # counts above 1e9 per compartment stay exact
  big <- sg_model(p, sg_cube(sg_test_cube_wild()), wx, area = 1e9)
  x_big <- sg_init(big, "state", state = as.numeric(sg_init(big, "equilibrium", temperature = 27)))
  set.seed(6)
  rb <- ns$.sg_cb_run_cpp(big, x_big)
  expect_gt(max(rb$state), 1e9)
  expect_true(flow_ok(rb) && all(rb$state == round(rb$state)))
  # step sizes from 0.05 to 1 day give similar means
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), wx, area = 1e5)
  x <- sg_init(m, "state", state = as.numeric(sg_init(m, "equilibrium", temperature = 27)))
  ode <- rowSums(simulate_ode(m, x)$state)
  for (dt in c(0.05, 0.5, 1)) {
    set.seed(7)
    r <- ns$.sg_cb_run_cpp(m, x, dt = dt)
    expect_true(flow_ok(r))
    expect_lt(max(abs(rowSums(r$state) - ode)) / max(ode), if (dt == 1) 0.2 else 0.05)
  }
  # very fast mating, long chains, a ten-genotype homing drive, females only
  settings <- list(
    list(sg_params(nu = 1e4), sg_cube(sg_test_cube_mendelian())),
    list(sg_params(n_E = 25, n_L = 25, n_P = 25, n_A = 25), sg_cube(sg_test_cube_wild()))
  )
  if (requireNamespace("MGDrivE", quietly = TRUE)) {
    settings <- c(settings, list(list(p, sg_cube(MGDrivE::cubeHomingDrive()))))
  }
  for (st in settings) {
    ms <- sg_model(st[[1]], st[[2]], sg_forcing(wx_all, st[[1]], "1971-05-01", "1971-05-30"),
      sg_releases(st[[2]]$genotypes[st[[2]]$n_genotypes], "adult", "male", times = 5, number = 500),
      area = 100
    )
    set.seed(8)
    r <- ns$.sg_cb_run_cpp(ms, sg_init(ms, values = c(E = 50, L = 50, P = 50, M = 20, F = 20)))
    expect_true(flow_ok(r))
    expect_false(anyNA(r$state))
  }
  females <- sg_releases("W", "adult", "female", times = 1, number = 50)
  fem <- sg_model(p, sg_cube(sg_test_cube_wild()), wx, females)
  set.seed(9)
  rf <- ns$.sg_cb_run_cpp(fem, sg_init(fem, values = c(E = 0)))
  expect_identical(sum(rf$counters[, "matings[W|W]"]), 0) # no males: no mating
  expect_identical(sum(rf$state[, sg_state_index(fem, "F")]), 0)
})
