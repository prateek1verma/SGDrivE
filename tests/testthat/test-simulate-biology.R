# M2 biological checks of the reference solver over time: Hardy-Weinberg
# neutrality, the instant-mating limit, SIT mechanics, genotype lifespan
# costs (D-002, D-018) and sex-linked cubes.

skip_if_not_installed("deSolve")

p <- sg_params()
wx90 <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-07-29")

geno_totals <- function(m, sim, block, row) {
  vapply(m$cube$genotypes, function(g) {
    sum(sim$state[row, sg_state_index(m, block, genotype = g)])
  }, 0)
}

test_that("a neutral Mendelian cube stays at Hardy-Weinberg with constant allele frequency", {
  m <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx90)
  sim <- simulate_ode(m, sg_init(m, "hardy_weinberg", alleles = c(A = 0.3, a = 0.7)),
    rtol = 1e-10, atol = 1e-12
  )
  for (b in c("E", "L", "P", "M", "U", "F")) {
    for (row in c(31, 61, 91)) {
      x <- geno_totals(m, sim, b, row)
      if (sum(x) < 1e-6) next
      f <- x / sum(x)
      allele_a <- f[["Aa"]] / 2 + f[["aa"]]
      expect_equal(allele_a, 0.7, tolerance = 1e-8, label = paste("allele a,", b, row))
      expect_equal(unname(f), c(0.09, 0.42, 0.49), tolerance = 1e-8, label = paste("HW,", b, row))
    }
  }
  # mated pairs stay in random-mating (product) proportions
  hw <- c(AA = 0.09, Aa = 0.42, aa = 0.49)
  f_tot <- sum(sim$state[91, sg_state_index(m, "F")])
  for (g in names(hw)) {
    for (mm in names(hw)) {
      fm <- sum(sim$state[91, sg_state_index(m, "F", genotype = g, mate = mm)])
      expect_equal(fm / f_tot, hw[[g]] * hw[[mm]], tolerance = 1e-8)
    }
  }
})

test_that("rate mating approaches instant mating as nu grows (S14)", {
  pr <- sg_params(nu = 1000)
  pi <- sg_params(mating = "instant")
  cube <- sg_cube(sg_test_cube_mendelian())
  # nu = 1000 makes the system deliberately stiff, so keep the window short
  period <- c("1971-05-01", "1971-06-14")
  mr <- sg_model(pr, cube, sg_forcing(sg_example_weather(), pr, period[1], period[2]))
  mi <- sg_model(pi, cube, sg_forcing(sg_example_weather(), pi, period[1], period[2]))
  sr <- simulate_ode(mr, sg_init(mr))
  si <- simulate_ode(mi, sg_init(mi))
  tot_r <- SGDrivE:::.sg_sim_block_totals(sr)
  tot_i <- SGDrivE:::.sg_sim_block_totals(si)
  # unmated females are a vanishing fraction; all other blocks agree closely
  expect_lt(max(tot_r[, "U"]) / max(tot_r[, "F"]), 0.01)
  for (b in c("E", "L", "P", "M")) {
    expect_lt(max(abs(tot_r[, b] - tot_i[, b])) / max(tot_i[, b]), 5e-3)
  }
  females_r <- tot_r[, "U"] + tot_r[, "F"]
  females_i <- tot_i[, "U"] + tot_i[, "F"]
  expect_lt(max(abs(females_r - females_i)) / max(tot_i[, "F"]), 5e-3)
})

test_that("sterile male releases divert matings and suppress the wild population", {
  raw <- sg_test_cube_sit(residual = 0)
  raw$eta[, "S"] <- 0.5
  cube <- sg_cube(raw)
  control <- sg_model(p, cube, wx90)
  rel <- sg_releases("S", "adult", "male", times = seq(1, 90, by = 3), number = 50)
  sit <- sg_model(p, cube, wx90, rel)
  s0 <- simulate_ode(control)
  s1 <- simulate_ode(sit)
  k <- s1$counters
  share_sterile <- sum(k[, "matings[W|S]"]) / sum(k[, "matings[W|S]"] + k[, "matings[W|W]"])
  expect_gt(share_sterile, 0.8)
  expect_identical(sum(s0$counters[, "matings[W|S]"]), 0)
  expect_identical(sum(k[, "eggs_laid[S]"]), 0) # sterile matings give no offspring
  wild_adults <- function(m, s) {
    sum(s$state[nrow(s$state), c(
      sg_state_index(m, "M", genotype = "W"), sg_state_index(m, "U", genotype = "W"),
      sg_state_index(m, "F", genotype = "W")
    )])
  }
  expect_lt(wild_adults(sit, s1), 0.5 * wild_adults(control, s0))
})

test_that("a genotype whose lifespan ends before first oviposition never lays (D-002, D-018)", {
  cube <- sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25))
  m <- sg_model(p, cube, wx90)
  init <- sg_init(m, values = c(M = 1, F = 1), genotypes = c(W = 0.5, S = 0.5))
  sim <- simulate_ode(m, init)
  f_s <- sim$state[, sg_state_index(m, "F", genotype = "S")]
  expect_true(all(f_s[, 4:10] == 0)) # never enters the laying classes
  # eggs come only from wild-type mothers: compare with a run without S mothers
  init_w <- init
  init_w[sg_state_index(m, "F", genotype = "S")] <- 0
  sim_w <- simulate_ode(m, init_w)
  expect_equal(sim$counters[1:10, "eggs_laid[W]"], sim_w$counters[1:10, "eggs_laid[W]"],
    tolerance = 1e-6
  )
})

test_that("a shorter lifespan shortens adult survival", {
  cube <- sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.5))
  m <- sg_model(p, cube, wx90)
  init <- sg_init(m, values = c(M = 1), genotypes = c(W = 0.5, S = 0.5))
  sim <- simulate_ode(m, init)
  males <- function(g) rowSums(sim$state[, sg_state_index(m, "M", genotype = g)])
  life_w <- sum(males("W")) # area under the survival curve, cohorts of equal size
  life_s <- sum(males("S"))
  expect_lt(life_s, 0.75 * life_w)
  expect_gt(life_s, 0.3 * life_w)
})

test_that("sex-linked cubes: males and females only from genotypes of their sex", {
  skip_if_not_installed("MGDrivE")
  cube <- sg_cube(MGDrivE::cubeXShredderY())
  wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-06-14")
  m <- sg_model(p, cube, wx)
  sim <- simulate_ode(m)
  expect_false(anyNA(sim$state))
  female_only <- cube$genotypes[cube$phi == 1]
  male_only <- cube$genotypes[cube$phi == 0]
  expect_identical(sum(abs(sim$state[, sg_state_index(m, "M", genotype = female_only)])), 0)
  expect_identical(sum(abs(sim$state[, c(
    sg_state_index(m, "U", genotype = male_only), sg_state_index(m, "F", genotype = male_only)
  )])), 0)
  expect_gt(min(sim$diagnostics$min_state), -1e-8)
})
