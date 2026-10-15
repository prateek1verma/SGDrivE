alloc <- SGDrivE:::.sg_allocate_substages
chain <- SGDrivE:::.sg_adult_chain

test_that("sub-stages are allocated in proportion to segment lengths, at least one each", {
  expect_identical(alloc(10, c(61, 141.8)), c(3L, 7L))
  expect_identical(alloc(25, c(56.784, 202.8 - 56.784)), c(7L, 18L)) # exact quota 7
  expect_identical(alloc(10, c(24, 37, 141.8)), c(1L, 2L, 7L))
  expect_identical(alloc(3, c(1, 1000, 1)), c(1L, 1L, 1L))
  expect_identical(alloc(10, 202.8), 10L)
  expect_identical(sum(alloc(7, c(0.5, 3, 9, 1))), 7L)
  expect_error(alloc(2, c(1, 1, 1)), "too small")
})

ell <- c(W = 1, H = 0.5, Q = 0.25, L = 1.5)

test_that("default female chain: exact mean age at first oviposition for every genotype", {
  ch <- chain(c(X_rec = 0, X_ovi = 61), 202.8, 10L, ell, sex = "female")
  expect_identical(ch$boundaries, c(0, 61, 202.8))
  expect_identical(ch$segment_counts, c(3L, 7L))
  first_lay <- which(ch$flags$X_ovi)[1]
  expect_identical(first_lay, 4L)
  expect_true(all(ch$flags$X_rec)) # receptive from emergence
  for (g in c("W", "H", "L")) {
    # mean age on entering the first laying class = sum of earlier widths
    expect_equal(sum(ch$width[g, seq_len(first_lay - 1)]), 61, tolerance = sg_tol_exact)
  }
})

test_that("genotype lifespan shortens or lengthens only the segment in which it ends (D-002)", {
  ch <- chain(c(X_ovi = 61), 202.8, 10L, ell, sex = "female")
  life <- rowSums(ch$width, na.rm = TRUE)
  expect_equal(unname(life), unname(ell * 202.8), tolerance = sg_tol_exact)
  expect_equal(unname(ch$width["H", 4:10]), rep((101.4 - 61) / 7, 7), tolerance = sg_tol_exact)
  expect_equal(unname(ch$width["L", 4:10]), rep((304.2 - 61) / 7, 7), tolerance = sg_tol_exact)
  expect_identical(unname(ch$last_class[c("W", "H", "L")]), c(10L, 10L, 10L))
})

test_that("genotypes whose lifespan ends before a threshold never reach it", {
  ch <- chain(c(X_ovi = 61), 202.8, 10L, ell, sex = "female")
  expect_identical(ch$last_class[["Q"]], 3L) # 50.7 dd < 61 dd
  expect_equal(unname(ch$width["Q", 1:3]), rep(50.7 / 3, 3), tolerance = sg_tol_exact)
  expect_true(all(is.na(ch$width["Q", 4:10])))
  # lifespan exactly at the threshold: the next segment is not reached
  ch2 <- chain(c(X_ovi = 101.4), 202.8, 10L, c(H = 0.5), sex = "female")
  expect_identical(ch2$last_class[["H"]], ch2$segment_counts[1])
})

test_that("three female segments with receptivity and oviposition thresholds", {
  ch <- chain(c(X_rec = 24, X_ovi = 61), 202.8, 10L, c(W = 1), sex = "female")
  expect_identical(ch$segment_counts, c(1L, 2L, 7L))
  expect_identical(which(ch$flags$X_rec)[1], 2L)
  expect_identical(which(ch$flags$X_ovi)[1], 4L)
  expect_equal(unname(ch$width["W", 1]), 24, tolerance = sg_tol_exact)
  expect_equal(sum(ch$width["W", 1:3]), 61, tolerance = sg_tol_exact)
})

test_that("user-given sub-stages per segment are used and checked", {
  ch <- chain(c(X_ovi = 61), 202.8, 10L, c(W = 1), substages = c(5, 5), sex = "female")
  expect_identical(ch$segment_counts, c(5L, 5L))
  expect_equal(unname(ch$width["W", 1]), 61 / 5, tolerance = sg_tol_exact)
  expect_error(
    chain(c(X_ovi = 61), 202.8, 10L, c(W = 1), substages = 10, sex = "female"),
    "must give 2 whole"
  )
  expect_error(
    chain(c(X_ovi = 61), 202.8, 10L, c(W = 1), substages = c(4, 4), sex = "female"),
    "must sum to adult_substages = 10"
  )
  expect_error(
    chain(c(X_ovi = 61), 202.8, 10L, c(W = 1), substages = c(0, 10), sex = "female"),
    "whole number"
  )
})

test_that("gp2014 preset reproduces the equal-width Gutierrez & Ponti chain", {
  gp <- sg_params(preset = "gp2014")
  ch <- chain(c(X_ovi = gp$oviposition_age_dd), gp$adult_lifespan_dd, gp$adult_substages,
    c(W = 1),
    sex = "female"
  )
  expect_identical(ch$segment_counts, c(7L, 18L))
  expect_equal(unname(ch$width["W", ]), rep(202.8 / 25, 25), tolerance = 1e-12)
  expect_identical(which(ch$flags$X_ovi)[1], 8L) # laying from class 8, as in the Pascal code
})
