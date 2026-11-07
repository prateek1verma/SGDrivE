# Output tables (D-026): totals, long tables, allele frequencies, female
# classes and the case proxy, for ODE and chain-binomial simulations.

ns <- asNamespace("ScrewDrivE")
p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, "1971-05-01", "1971-05-31")
cube3 <- sg_cube(sg_test_cube_mendelian())
mod <- sg_model(p, cube3, wx, area = 10)
hw <- sg_init(mod, "hardy_weinberg", alleles = c(A = 0.8, a = 0.2))
sim <- simulate_ode(mod, init = hw, record = "state")
hw20 <- sg_init(mod, "hardy_weinberg",
  alleles = c(A = 0.8, a = 0.2), values = c(E = 20, L = 20, P = 20, M = 10, F = 10)
)
cb <- simulate_cb(mod, hw20, n_rep = 2, seed = 1, record = "state")
block_sum <- function(state, b) rowSums(state[, substr(colnames(state), 1, 1) == b, drop = FALSE])

test_that("totals add up the states, by any grouping and per km^2", {
  all <- sg_totals(sim, by = character(0))
  expect_equal(all$value, unname(rowSums(sim$state)))
  by_stage <- sg_totals(sim, by = "stage")
  expect_setequal(unique(by_stage$stage), ns$.sg_stage_names)
  expect_equal(by_stage$value[by_stage$stage == "mated_females"], unname(block_sum(sim$state, "F")))
  sex <- sg_totals(sim, by = "sex", stages = "adults")
  expect_setequal(unique(sex$sex), c("male", "female"))
  expect_equal(
    sex$value[sex$sex == "female"],
    unname(block_sum(sim$state, "U") + block_sum(sim$state, "F"))
  )
  km2 <- sg_totals(sim, by = "genotype", per = "km2")
  node <- sg_totals(sim, by = "genotype")
  expect_equal(km2$value, node$value / 10)
  juv <- sg_totals(sim, by = "sex", stages = "juveniles")
  expect_identical(unique(juv$sex), "unsexed")
  expect_identical(names(node), c("day", "date", "genotype", "value"))
})

test_that("stochastic totals carry replicates and agree with keep = 'totals'", {
  tot <- sg_totals(cb, by = "genotype")
  expect_identical(names(tot), c("day", "date", "replicate", "genotype", "value"))
  for (r in 1:2) {
    expect_equal(
      sum(tot$value[tot$replicate == r & tot$day == 31]),
      sum(cb$state[32, , r])
    )
  }
  ct <- simulate_cb(mod, as.numeric(cb$meta$init), n_rep = 2, seed = 1, keep = "totals")
  by <- c("stage", "genotype")
  expect_identical(sg_totals(ct, by = by), sg_totals(cb, by = by))
  expect_error(sg_totals(ct, by = "mate"), "full state")
})

test_that("as.data.frame() gives the long table by stage, sex, genotype and mate", {
  d <- as.data.frame(sim)
  expect_identical(names(d), c("day", "date", "stage", "sex", "genotype", "mate", "value"))
  expect_equal(sum(d$value[d$day == 31]), sum(sim$state[32, ]))
  expect_setequal(unique(d$mate[d$stage == "mated_females"]), cube3$genotypes)
  expect_true(all(is.na(d$mate[d$stage != "mated_females"])))
  dc <- as.data.frame(cb)
  expect_true("replicate" %in% names(dc))
  expect_equal(sum(dc$value), sum(cb$state))
})

test_that("allele frequencies: Hardy-Weinberg, custom allele tables, empty populations", {
  af <- sg_allele_freq(sim)
  expect_equal(af$freq[af$allele == "a"], rep(0.2, 32), tolerance = 1e-8)
  expect_equal(
    af$count[af$allele == "A"] + af$count[af$allele == "a"],
    2 * unname(rowSums(sim$state[, substr(colnames(sim$state), 1, 1) %in% c("M", "U", "F")]))
  )
  custom <- matrix(c(2, 1, 0), 3, 1, dimnames = list(cube3$genotypes, "A"))
  expect_equal(sg_allele_freq(sim, alleles = custom)$freq, rep(1, 32))
  expect_error(sg_allele_freq(sim, alleles = custom[1:2, , drop = FALSE]), "row for every genotype")
  empty <- simulate_ode(mod, init = sg_init(mod, values = c(E = 0)), record = "state")
  expect_true(all(is.na(sg_allele_freq(empty)$freq)))
  skip_if_not_installed("MGDrivE")
  xs <- sg_model(p, sg_cube(MGDrivE::cubeXShredderY()), wx)
  axs <- sg_allele_freq(simulate_ode(xs, record = "state"))
  expect_true(all(c("X", "Y") %in% axs$allele))
  # wild type only: Y is one of four sex chromosomes among equal males and females
  expect_equal(axs$freq[axs$allele == "Y" & axs$day == 0], 0.25)
})

test_that("female classes: identities, wild type and fertility weights", {
  f <- sg_females(sim)
  expect_equal(f$females, f$unmated + f$mated)
  expect_equal(f$wild_type / f$females, rep(0.64, 32), tolerance = 1e-8) # 0.8 squared
  expect_equal(f$fertile_equivalent, f$mated) # every cross fertile
  laying <- ns$.sg_reference_setup(mod)$laying
  fcols <- grepl("^F", colnames(sim$state))
  sub <- as.integer(sub("^.*,([0-9]+)\\]$", "\\1", colnames(sim$state)))
  expect_equal(f$egg_laying, unname(rowSums(sim$state[, fcols & laying[pmax(sub, 1)]])))
  # sterile and part-fertile crosses
  raw <- sg_test_cube_sit(residual = 0.1)
  ms <- sg_model(p, sg_cube(raw), wx, area = 1)
  x <- numeric(ms$n_states)
  x[sg_state_index(ms, "F", genotype = "W", mate = "W", substage = 5)] <- 10
  x[sg_state_index(ms, "F", genotype = "W", mate = "S", substage = 5)] <- 20
  s1 <- simulate_ode(ms, init = x, record = "state")
  expect_equal(sg_females(s1)$fertile_equivalent[1], 10 + 0.1 * 20)
  ct <- simulate_cb(mod, n_rep = 1, seed = 2, keep = "totals")
  expect_true(all(is.na(sg_females(ct)$egg_laying)))
  expect_equal(sg_females(sim, per = "km2")$females, f$females / 10)
})

test_that("cases: fertile oviposition events by period, wound fraction, cases per host", {
  daily <- sg_cases(sim, per = "node")
  expect_equal(daily$cases, unname(sim$counters[, "fertile_batches"]))
  expect_equal(daily$events, unname(sim$counters[, "egg_batches"]))
  expect_true(all(daily$cases <= daily$events + 1e-12))
  expect_equal(daily$wound_fraction, daily$cases / (5 * 10))
  km2 <- sg_cases(sim)
  expect_equal(km2$cases, daily$cases / 10)
  weeks <- sg_cases(sim, period = "week", host_density = 20)
  expect_identical(weeks$days, c(7L, 7L, 7L, 7L, 3L))
  expect_equal(sum(weeks$cases), sum(km2$cases))
  expect_equal(weeks$cases_per_host, weeks$cases / 20)
  month <- sg_cases(sim, period = "month")
  expect_identical(nrow(month), 1L)
  expect_identical(month$period, as.Date("1971-05-01"))
  cc <- sg_cases(cb, period = "year", per = "node")
  expect_identical(cc$replicate, 1:2)
  expect_equal(cc$cases, unname(colSums(cb$counters[, "fertile_batches", ])))
  expect_true(all(cc$cases == round(cc$cases))) # whole cases
})

test_that("females mated to sterile males lay eggs but cause no cases", {
  raw <- sg_test_cube_sit(residual = 0) # crosses with S males give no viable eggs
  ms <- sg_model(p, sg_cube(raw), wx, area = 1)
  x <- numeric(ms$n_states)
  laying <- 5:10
  x[sg_state_index(ms, "F", genotype = "W", mate = "S", substage = laying)] <- 20
  only_sterile <- simulate_ode(ms, init = x, record = "state")
  day1 <- sg_cases(only_sterile, per = "node")[1, ]
  expect_gt(day1$events, 0)
  expect_identical(day1$cases, 0)
  x2 <- x
  x2[sg_state_index(ms, "F", genotype = "W", mate = "W", substage = laying)] <- 20
  half <- sg_cases(simulate_ode(ms, init = x2, record = "state"), per = "node")[1, ]
  expect_equal(half$cases / half$events, 0.5, tolerance = 1e-6) # wild-mated half
  # residual fertility 1e-8: essentially no cases; a cross with half the
  # viable eggs: still nearly every batch is a case
  for (v in c(1e-8, 0.5)) {
    rv <- sg_test_cube_sit(residual = v)
    mv <- sg_model(p, sg_cube(rv), wx, area = 1)
    xv <- numeric(mv$n_states)
    xv[sg_state_index(mv, "F", genotype = "W", mate = "S", substage = laying)] <- 20
    d1 <- sg_cases(simulate_ode(mv, init = xv, record = "state"), per = "node")[1, ]
    if (v < 1e-6) expect_lt(d1$cases / d1$events, 1e-5) else expect_gt(d1$cases / d1$events, 0.99)
  }
  # stochastic: sterile-mated females' batches are never cases
  set.seed(4)
  r <- asNamespace("ScrewDrivE")$.sg_cb_run_cpp(sg_model(p, sg_cube(raw), wx, area = 100), x * 100)
  expect_gt(sum(r$counters[1:3, "egg_batches"]), 0)
  expect_identical(sum(r$counters[1:3, "fertile_batches"]), 0)
})

test_that("output functions check their input", {
  expect_error(sg_totals(list()), "simulation from")
  expect_error(sg_totals(sim, by = "age"), "`by` must contain only")
  expect_error(sg_totals(sim, stages = "adult"), "Unknown stage")
  expect_error(sg_cases(sim, host_density = -1), "positive number")
  expect_error(sg_cases(sim, period = "fortnight"), "must be one of")
  old <- sim
  old$meta$info <- NULL
  expect_error(sg_totals(old), "older version")
})
