p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
cube3 <- sg_cube(sg_test_cube_mendelian())

test_that("sg_model() assembles a model with the documented components", {
  m <- sg_model(p, cube3, wx)
  expect_s3_class(m, "sg_model")
  expect_named(m, c(
    "params", "cube", "forcing", "area", "dims", "adult", "layout", "n_states",
    "counts", "releases", "meta"
  ))
  expect_identical(m$dims$n_genotypes, 3L)
  expect_identical(m$meta$package_version, as.character(utils::packageVersion("SGDrivE")))
  expect_match(m$meta$hash, "^[0-9a-f]{32}$")
  expect_identical(nrow(m$releases), 0L)
})

test_that("the state vector size follows n_G (n_E + n_L + n_P + 2 n_A) + n_G^2 n_A", {
  for (cb in list(sg_cube(sg_test_cube_wild()), sg_cube(sg_test_cube_sit()), cube3)) {
    m <- sg_model(p, cb, wx)
    n_g <- cb$n_genotypes
    expect_identical(m$n_states, as.integer(n_g * (3 + 6 + 10 + 2 * 10) + n_g^2 * 10))
    expect_identical(m$layout$block, c("E", "L", "P", "M", "U", "F"))
    # blocks are contiguous, without gaps or overlaps
    expect_identical(m$layout$start, cumsum(c(1L, m$layout$length[-6])))
    expect_identical(sum(m$layout$length), m$n_states)
  }
})

test_that("state indices and names agree and cover every state once", {
  m <- sg_model(p, cube3, wx)
  nm <- sg_state_names(m)
  expect_length(nm, m$n_states)
  expect_false(anyDuplicated(nm) > 0)
  all_idx <- unlist(lapply(m$layout$block, function(b) sg_state_index(m, b)))
  expect_identical(sort(all_idx), seq_len(m$n_states))
  expect_identical(nm[sg_state_index(m, "E", genotype = "Aa", substage = 2)], "E[Aa,2]")
  i <- sg_state_index(m, "F", genotype = "AA", mate = "aa", substage = 1)
  expect_identical(nm[i], "F[AA|aa,1]")
  # F index: start + ((g - 1) n_G + (m - 1)) n_A + (s - 1)
  expect_identical(i, m$layout$start[6] + (0L * 3L + 2L) * 10L + 0L)
  # sub-stage varies fastest
  expect_identical(diff(sg_state_index(m, "M", genotype = "AA")), rep(1L, 9))
})

test_that("state index arguments are checked", {
  m <- sg_model(p, cube3, wx)
  expect_error(sg_state_index(m, "X"), "must be one of")
  expect_error(sg_state_index(m, "E", genotype = "Zz"), "Unknown genotype")
  expect_error(sg_state_index(m, "M", mate = "AA"), "only to block \"F\"")
  expect_error(sg_state_index(m, "E", substage = 4), "between 1 and 3")
  expect_error(sg_state_index(list(), "E"), "sg_model object")
  expect_error(sg_state_names(list()), "sg_model object")
})

test_that("adult chains follow D-018 and genotype lifespans", {
  raw <- sg_test_cube_sit(omega_s = 1)
  cb <- sg_cube(raw, lifespan_fraction = c(W = 1, S = 0.25))
  m <- sg_model(p, cb, wx)
  expect_identical(m$adult$female$segment_counts, c(3L, 7L))
  expect_identical(m$adult$male$segment_counts, 10L)
  expect_identical(m$adult$female$last_class[["S"]], 3L)
  expect_identical(m$adult$male$last_class[["S"]], 10L) # no male threshold
  expect_equal(sum(m$adult$male$width["S", ]), 0.25 * 202.8, tolerance = sg_tol_exact)
  expect_match(SGDrivE:::.sg_never_reached(m), "S: first oviposition")
})

test_that("instant mating drops the receptivity boundary", {
  m <- sg_model(sg_params(mating_mode = "instant"), cube3, sg_forcing(
    sg_example_weather()[1:30, ],
    sg_params(mating_mode = "instant")
  ))
  expect_named(m$adult$female$thresholds, "X_ovi")
  pr <- sg_params(female_receptive_age_dd = 24)
  m2 <- sg_model(pr, cube3, sg_forcing(sg_example_weather()[1:30, ], pr))
  expect_identical(m2$adult$female$segment_counts, c(1L, 2L, 7L))
})

test_that("adult_substages overrides the allocation", {
  m <- sg_model(p, cube3, wx, adult_segment_substages = list(female = c(5, 5), male = 10))
  expect_identical(m$adult$female$segment_counts, c(5L, 5L))
  expect_error(
    sg_model(p, cube3, wx, adult_segment_substages = list(fem = 1)),
    "female` and/or `male"
  )
  expect_error(
    sg_model(p, cube3, wx, adult_segment_substages = list(female = c(5, 6))),
    "sum to adult_substages"
  )
})

test_that("counts per node use the node area", {
  m <- sg_model(sg_params(wound_density = 5, male_density_half_mating = 2), cube3, wx, area = 25)
  expect_identical(m$counts$H, 125)
  expect_identical(m$counts$h_M, 50)
  expect_error(sg_model(p, cube3, wx, area = 0), "positive number")
  expect_error(sg_model(p, cube3, wx, area = c(1, 2)), "positive number")
})

test_that("release events map to the first sub-stage of their block", {
  rel <- c(
    sg_releases("aa", "adult", "both", times = c(10, 10), number = 100, per = "km2"),
    sg_releases("aa", "pupa", times = 20, number = 7),
    sg_releases("Aa", "egg", times = 30, number = 3)
  )
  m <- sg_model(p, cube3, wx, rel, area = 2)
  ev <- m$releases
  nm <- sg_state_names(m)
  expect_identical(nm[ev$index], c("M[aa,1]", "U[aa,1]", "P[aa,1]", "E[Aa,1]"))
  expect_identical(ev$day, c(10L, 10L, 20L, 30L))
  # two releases on day 10 aggregated: 2 x 50 per km2 x 2 km2 for each sex
  expect_identical(ev$count, c(200, 200, 7, 3))
  dated <- sg_releases("aa", "adult", "male", times = "1971-02-01", number = 1)
  expect_identical(sg_model(p, cube3, wx, dated)$releases$day, 32L)
  bad_rel <- sg_releases("Z", "pupa", times = 1, number = 1)
  expect_error(sg_model(p, cube3, wx, bad_rel), "not in the cube")
})

test_that("the forcing must be built with the same parameters", {
  p2 <- sg_params(optimal_survival_temp = 27)
  expect_error(
    sg_model(p2, cube3, wx),
    "optimal_survival_temp \\(forcing 27.2, params 27\\).*sg_forcing"
  )
  # eggs_per_day_max does not affect the forcing
  expect_no_error(sg_model(sg_params(eggs_per_day_max = 50), cube3, wx))
})

test_that("input objects are checked", {
  expect_error(sg_model(list(), cube3, wx), "sg_params object")
  expect_error(sg_model(p, list(), wx), "sg_cube object")
  expect_error(sg_model(p, cube3, list()), "sg_forcing object")
})

test_that("a one-genotype model has one of each block and names without mates", {
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), wx)
  expect_identical(m$n_states, 49L)
  expect_identical(sg_state_names(m)[m$layout$start[6]], "F[W|W,1]")
})

test_that("the hash is stable for identical inputs and changes with any input", {
  rel <- sg_releases("aa", "adult", "male", times = 5, number = 10)
  h1 <- sg_model(p, cube3, wx, rel)$meta$hash
  expect_identical(sg_model(p, cube3, wx, rel)$meta$hash, h1)
  expect_false(identical(sg_model(sg_params(eggs_per_day_max = 50), cube3, wx, rel)$meta$hash, h1))
  expect_false(identical(sg_model(p, cube3, wx, rel, area = 2)$meta$hash, h1))
  expect_false(identical(sg_model(p, cube3, wx)$meta$hash, h1))
})

test_that("print() summarises the model", {
  rel <- sg_releases("aa", "adult", "male", times = seq(60, by = 7, length.out = 10), number = 500)
  cb <- sg_cube(sg_test_cube_mendelian(), lifespan_fraction = c(AA = 1, Aa = 1, aa = 0.25))
  m <- sg_model(p, cb, wx, rel)
  out <- capture.output(print(m))
  expect_match(out[1], "^<sg_model> 207 states; 3 genotypes; 365 days")
  expect_true(any(grepl("female chain: \\[0, 61\\) dd: 3; \\[61, 202.8\\) dd: 7", out)))
  expect_true(any(grepl("aa: first oviposition \\(never lay eggs\\)", out)))
  expect_true(any(grepl("releases: 10 events on 10 days", out)))
  expect_invisible(print(m))
})
