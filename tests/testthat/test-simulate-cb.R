# simulate_cb(): replicates, random-number streams (D-022), outputs and
# checks. The step itself is tested in test-cb-reference.R and test-cb-cpp.R.

ns <- asNamespace("ScrewDrivE")
p <- sg_params()
wx_all <- sg_example_weather()
cube3 <- sg_cube(sg_test_cube_mendelian())
mod <- sg_model(p, cube3, sg_forcing(wx_all, p, "1971-04-01", "1971-05-30"),
  sg_releases("aa", "adult", "male", times = c(5, 20), number = 200.5),
  area = 200
)
init <- sg_init(mod, values = c(E = 3, L = 3, P = 3, M = 2, F = 2), per = "km2")

flow_error <- function(sim, r) {
  k <- sim$counters[, , r]
  tot <- rowSums(sim$state[, , r])
  births <- rowSums(k[, grepl("^eggs_laid", colnames(k)), drop = FALSE])
  deaths <- rowSums(k[, grepl("^deaths_|^senescence_", colnames(k)), drop = FALSE])
  max(abs(diff(tot) - (rowSums(sim$releases[, , r]) + births - deaths - k[, "emergence_deaths"])))
}

test_that("simulate_cb() returns replicate arrays and metadata", {
  sim <- simulate_cb(mod, init, n_rep = 3, seed = 1, record = "state")
  expect_s3_class(sim, "sg_sim_cb")
  expect_identical(sim$day, 0:60)
  expect_identical(sim$date[c(1, 61)], as.Date(c("1971-03-31", "1971-05-30")))
  expect_identical(dim(sim$state), c(61L, mod$n_states, 3L))
  expect_identical(dimnames(sim$state)[[2]], sg_state_names(mod))
  expect_identical(dim(sim$counters), c(60L, length(ns$.sg_counter_names(cube3$genotypes)), 3L))
  expect_identical(dim(sim$releases), c(60L, 4L, 3L))
  expect_true(all(sim$releases[5, "M", ] %in% c(200, 201))) # rounded without bias
  expect_identical(names(sim$extinct_day), c("1", "2", "3"))
  expect_true(all(sim$state == round(sim$state)) && all(sim$state >= 0))
  for (r in 1:3) expect_identical(flow_error(sim, r), 0)
  expect_identical(sim$meta$seed, 1L)
  expect_identical(dim(sim$meta$streams), c(7L, 3L))
  expect_identical(sim$meta$oviposition_mode, "batches")
  expect_identical(sim$meta$model_hash, mod$meta$hash)
  expect_false(identical(sim$state[, , 1], sim$state[, , 2])) # replicates differ
})

test_that("each replicate is the reference chain-binomial run on its own stream", {
  sim <- simulate_cb(mod, init, n_rep = 2, seed = 42, record = "state")
  ref <- ns$.sg_with_stream(sim$meta$streams[, 2], ns$.sg_cb_run_r(mod, init, n_days = 20))
  expect_identical(unname(sim$state[1:21, , 2]), unname(ref$state))
  expect_identical(unname(sim$counters[1:20, , 2]), unname(ref$counters))
})

test_that("runs are reproducible and replicates independent of how they are run", {
  a <- simulate_cb(mod, init, n_rep = 4, seed = 7, record = "state")
  b <- simulate_cb(mod, init, n_rep = 2, seed = 7, record = "state")
  expect_identical(a$state[, , 1:2], b$state) # replicate r does not depend on n_rep
  set.seed(3)
  c1 <- simulate_cb(mod, init, n_rep = 2, record = "state")
  set.seed(3)
  c2 <- simulate_cb(mod, init, n_rep = 2, record = "state")
  expect_identical(c1$state, c2$state)
  c3 <- simulate_cb(mod, init, n_rep = 2, record = "state") # the session's generator moved on
  expect_false(identical(c2$state, c3$state))
  # an explicit seed leaves the session's generator untouched
  set.seed(1)
  u1 <- stats::runif(1)
  set.seed(1)
  simulate_cb(mod, init, seed = 5, record = "state")
  expect_identical(stats::runif(1), u1)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
})

test_that("parallel runs give the same replicates", {
  skip_if_not_installed("future.apply")
  old <- future::plan("sequential")
  on.exit(future::plan(old))
  a <- simulate_cb(mod, init, n_rep = 3, seed = 9, record = "state")
  b <- simulate_cb(mod, init, n_rep = 3, seed = 9, parallel = TRUE, record = "state")
  expect_identical(b$state, a$state)
  expect_identical(b$counters, a$counters)
  expect_true(b$meta$parallel)
})

test_that("totals and output_every thin the output consistently", {
  full <- simulate_cb(mod, init, n_rep = 2, seed = 11, record = "state")
  tot <- simulate_cb(mod, init, n_rep = 2, seed = 11, keep = "totals", output_every = 7)
  expect_identical(tot$day, c(0L, seq(7L, 56L, 7L), 60L))
  expect_null(tot$state)
  expect_identical(dimnames(tot$totals)[[2]][1:3], c("E[AA]", "E[Aa]", "E[aa]"))
  agg <- ns$.sg_totals_matrix(mod)
  for (r in 1:2) {
    expect_equal(tot$totals[, , r], full$state[tot$day + 1, , r] %*% agg, ignore_attr = TRUE)
  }
  expect_identical(tot$counters, full$counters)
  out <- capture.output(print(tot))
  expect_match(out[1], "^<sg_sim_cb> 2 replicates x 60 days")
  expect_true(any(grepl("extinct: 0 of 2", out)))
})

test_that("small populations go extinct and the extinction day is reported", {
  m1 <- sg_model(p, sg_cube(sg_test_cube_wild()), sg_forcing(wx_all, p, "1971-01-01", "1971-06-30"))
  sim <- simulate_cb(m1, n_rep = 5, seed = 2, keep = "totals")
  # 1 km^2: the population usually dies out in winter (a few survive by chance)
  ext <- which(!is.na(sim$extinct_day))
  expect_gte(length(ext), 4L)
  for (r in ext) {
    expect_identical(sum(sim$totals[(sim$extinct_day[r] + 1):nrow(sim$totals), , r]), 0)
  }
  pattern <- sprintf("extinct: %d of 5 replicates \\(median day", length(ext))
  expect_true(any(grepl(pattern, capture.output(print(sim)))))
})

test_that("invalid arguments give informative errors", {
  expect_error(simulate_cb(list(), record = "state"), "sg_model object")
  expect_error(simulate_cb(mod, n_rep = 0, record = "state"), "whole number >= 1")
  expect_error(simulate_cb(mod, dt = 0.3, record = "state"), "whole number of steps")
  expect_error(simulate_cb(mod, output_every = 1.5, record = "state"), "whole number >= 1")
  expect_error(simulate_cb(mod, keep = "all"), "must be one of")
  expect_error(simulate_cb(mod, seed = 1.5, record = "state"), "single whole number")
  expect_error(simulate_cb(mod, parallel = NA, record = "state"), "TRUE or FALSE")
  expect_error(simulate_cb(mod, max_memory = 0, record = "state"), "positive number")
  expect_error(
    simulate_cb(mod, n_rep = 1000, max_memory = 1e6, record = "state"),
    "record = \"totals\""
  )
  expect_error(simulate_cb(mod, init = numeric(3), record = "state"), "length")
})
