view <- SGDrivE:::.sg_model_view
info <- SGDrivE:::sg_model_view_info_cpp

p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-06-29")

check_round_trip <- function(m) {
  v <- view(m)
  out <- info(v)
  d <- m$dims
  testthat::expect_identical(
    out$dims,
    c(d$n_genotypes, d$n_E, d$n_L, d$n_P, d$n_A, m$n_states, m$forcing$n_days, v$n_counters)
  )
  testthat::expect_identical(out$offsets, as.integer(m$layout$start - 1L))
  testthat::expect_identical(out$last_f, unname(m$adult$female$last_class))
  testthat::expect_identical(out$last_m, unname(m$adult$male$last_class))
  testthat::expect_identical(out$idx_F_last, m$n_states - 1L)
  testthat::expect_identical(out$idx_M_first, m$layout$start[4] - 1L)
  testthat::expect_equal(out$sum_egg_tensor, sum(m$cube$egg_tensor))
  testthat::expect_equal(out$sum_dx, sum(m$forcing$dx))
  testthat::expect_identical(out$n_releases, nrow(m$releases))
  invisible(v)
}

test_that("views of 1-, 2- and 3-genotype models round-trip through C++", {
  for (cb in list(
    sg_cube(sg_test_cube_wild()), sg_cube(sg_test_cube_sit()),
    sg_cube(sg_test_cube_mendelian())
  )) {
    check_round_trip(sg_model(p, cb, wx))
  }
})

test_that("a 10-genotype model and a sex-linked model round-trip", {
  skip_if_not_installed("MGDrivE")
  check_round_trip(sg_model(p, sg_cube(MGDrivE::cubeHomingDrive()), wx))
  check_round_trip(sg_model(p, sg_cube(MGDrivE::cubeXShredderY()), wx))
})

test_that("view fields follow the C++ conventions", {
  cb <- sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25))
  rel <- sg_releases("S", "adult", "male", times = c(3, 9), number = 5)
  m <- sg_model(sg_params(mating = "instant", larval_dd = "per_wound", gamma_L = 0.1), cb,
    sg_forcing(sg_example_weather()[1:60, ], sg_params(mating = "instant")), rel,
    area = 2
  )
  v <- view(m)
  expect_s3_class(v, "sg_model_view")
  expect_identical(attr(v, "counter_names"), SGDrivE:::.sg_counter_names(cb$genotypes))
  expect_identical(v$offsets[1], 0L) # 0-based
  expect_identical(v$rel_index, as.integer(m$releases$index - 1L))
  expect_identical(v$rel_day, m$releases$day)
  expect_identical(c(v$mating_instant, v$search_all_adults, v$larval_per_wound), c(1L, 0L, 1L))
  expect_identical(v$mu_L0, -log(0.95))
  expect_identical(v$H, 5 * 2) # counts per node
  # inverse widths: 1 / width, 0 for adult classes a genotype never reaches
  w <- m$adult$female$width
  iw <- matrix(v$inv_w_f, m$dims$n_A, 2)
  expect_equal(iw[1:3, 2], 1 / w["S", 1:3])
  expect_identical(iw[4:10, 2], rep(0, 7))
  expect_true(all(is.finite(v$inv_w_f)))
  expect_identical(v$last_f, c(10L, 3L))
  expect_true(all(v$laying %in% 0:1) && all(v$receptive %in% 0:1) && all(v$mature %in% 0:1))
  expect_identical(info(v)$options, c(TRUE, FALSE, TRUE))
})

test_that("C++ rejects inconsistent views with informative errors", {
  v <- view(sg_model(p, sg_cube(sg_test_cube_mendelian()), wx))
  bad <- v
  bad$dx <- bad$dx[-1]
  expect_error(info(bad), "field 'dx' has length 59, expected 60")
  bad <- v
  bad$n_states <- 10L
  expect_error(info(bad), "n_states = 10, expected 207")
  bad <- v
  bad$offsets[3] <- 1L
  expect_error(info(bad), "block 2 starts at 1, expected 27")
  bad <- v
  bad$laying[1] <- 2L
  expect_error(info(bad), "'laying' must contain only 0/1")
  bad <- v
  bad$mu_T[5] <- NaN
  expect_error(info(bad), "'mu_T' contains a non-finite value")
  bad <- v
  bad$R <- NULL
  expect_error(info(bad), "missing field 'R'")
  bad <- v
  bad$last_f[1] <- 11L
  expect_error(info(bad), "last reached adult class out of range")
  bad <- v
  bad$rel_day <- 1L
  bad$rel_index <- 999L
  bad$rel_count <- 1
  expect_error(info(bad), "release 1 is out of range")
  bad <- v
  bad$n_counters <- 3L
  expect_error(info(bad), "n_counters = 3")
})

test_that(".sg_model_view() needs a model", {
  expect_error(view(list()), "sg_model object")
})
