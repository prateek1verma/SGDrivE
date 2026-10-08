p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-03-31")
m3 <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx)
m1 <- sg_model(p, sg_cube(sg_test_cube_wild()), wx)

test_that("default start reproduces the Gutierrez & Ponti initial state (one genotype)", {
  x <- sg_init(m1)
  expect_s3_class(x, "sg_state")
  expect_identical(names(x), sg_state_names(m1))
  expect_identical(attr(x, "model_hash"), m1$meta$hash)
  expect_true(all(x[sg_state_index(m1, "E")] == 0.25))
  expect_true(all(x[sg_state_index(m1, "P")] == 0.25))
  expect_true(all(x[sg_state_index(m1, "M")] == 0.125))
  expect_true(all(x[sg_state_index(m1, "F")] == 0.125)) # all females mated
  expect_true(all(x[sg_state_index(m1, "U")] == 0))
})

test_that("default weights put everything in the wild type", {
  x <- sg_init(m3)
  tot <- SGDrivE:::.sg_state_block_totals(unclass(x))
  expect_identical(unname(tot[, "Aa"]), rep(0, 6))
  expect_identical(unname(tot[, "aa"]), rep(0, 6))
  expect_equal(unname(tot["E", "AA"]), 0.25 * 3)
  # mated wild-type females all carry wild-type sperm
  expect_equal(sum(x[sg_state_index(m3, "F", genotype = "AA", mate = "AA")]), 1.25)
  expect_identical(sum(x[sg_state_index(m3, "F", genotype = "AA", mate = c("Aa", "aa"))]), 0)
})

test_that("genotype weights split every block; mated females mate at random", {
  x <- sg_init(m3, genotypes = c(AA = 0.5, aa = 0.5))
  f_hom_hom <- sum(x[sg_state_index(m3, "F", genotype = "aa", mate = "aa")])
  f_hom_wild <- sum(x[sg_state_index(m3, "F", genotype = "aa", mate = "AA")])
  expect_equal(f_hom_hom, f_hom_wild) # equal male weights, eta = 1
  expect_equal(sum(x[sg_state_index(m3, "F")]), 0.125 * 10)
  expect_equal(sum(x[sg_state_index(m3, "L", genotype = "AA")]), 0.25 * 6 / 2)
})

test_that("Hardy-Weinberg weights follow p^2, 2pq, q^2 and F pairs are products", {
  x <- sg_init(m3, "hardy_weinberg", alleles = c(A = 0.7, a = 0.3))
  pupae <- vapply(c("AA", "Aa", "aa"), function(g) sum(x[sg_state_index(m3, "P", genotype = g)]), 0)
  expect_equal(unname(pupae / sum(pupae)), c(0.49, 0.42, 0.09), tolerance = sg_tol_exact)
  hw <- c(AA = 0.49, Aa = 0.42, aa = 0.09)
  for (g in names(hw)) {
    for (mm in names(hw)) {
      fm <- sum(x[sg_state_index(m3, "F", genotype = g, mate = mm)])
      expect_equal(fm, 1.25 * hw[[g]] * hw[[mm]], tolerance = sg_tol_exact)
    }
  }
  expect_error(sg_init(m3, "hardy_weinberg", alleles = c(A = 0.7, a = 0.2)), "summing to 1")
  expect_error(sg_init(m3, "hardy_weinberg", alleles = c(A = 1)), "missing from `alleles`")
})

test_that("Hardy-Weinberg needs two-letter genotype names", {
  skip_if_not_installed("MGDrivE")
  cb <- sg_cube(MGDrivE::cubeReciprocalTranslocations())
  m <- sg_model(p, cb, wx)
  expect_error(sg_init(m, "hardy_weinberg", alleles = c(A = 0.5, a = 0.5)), "two allele letters")
})

test_that("sex-linked cubes: adults only go to genotypes of the right sex", {
  skip_if_not_installed("MGDrivE")
  cb <- sg_cube(MGDrivE::cubeXShredderY())
  m <- sg_model(p, cb, wx)
  x <- sg_init(m) # wild types XX (female) and XY (male)
  expect_identical(sum(x[sg_state_index(m, "M", genotype = "XX")]), 0)
  expect_identical(sum(x[sg_state_index(m, "F", genotype = "XY")]), 0)
  expect_equal(sum(x[sg_state_index(m, "M", genotype = "XY")]), 0.125 * 10)
  expect_equal(sum(x[sg_state_index(m, "F", genotype = "XX", mate = "XY")]), 0.125 * 10)
})

test_that("adults only fill the classes their genotype reaches (D-018)", {
  cb <- sg_cube(sg_test_cube_sit(), ell = c(W = 1, S = 0.25))
  m <- sg_model(p, cb, wx)
  x <- sg_init(m, genotypes = c(W = 1, S = 1))
  expect_true(all(x[sg_state_index(m, "F", genotype = "S", substage = 4:10)] == 0))
  expect_true(all(x[sg_state_index(m, "F", genotype = "S", substage = 1:3)] > 0))
})

test_that("values per km2 are multiplied by the node area", {
  m <- sg_model(p, sg_cube(sg_test_cube_wild()), wx, area = 4)
  x <- sg_init(m, values = c(P = 1), per = "km2")
  expect_true(all(x[sg_state_index(m, "P")] == 4))
  expect_identical(sum(x[sg_state_index(m, "E")]), 0) # missing blocks are 0
})

test_that("a full state vector can be supplied", {
  s <- seq_len(m1$n_states) / 10
  x <- sg_init(m1, "state", state = s)
  expect_identical(unname(unclass(x)[seq_along(s)]), s)
  named <- stats::setNames(s, sg_state_names(m1))
  expect_identical(as.numeric(sg_init(m1, "state", state = named)), s)
  expect_error(sg_init(m1, "state"), "needs `state`")
  expect_error(sg_init(m1, "state", state = s[-1]), "length 49")
  expect_error(sg_init(m1, "state", state = -s), ">= 0")
  expect_error(
    sg_init(m1, "state", state = stats::setNames(s, rev(sg_state_names(m1)))),
    "do not match"
  )
})

test_that("invalid arguments give informative errors", {
  expect_error(sg_init(list()), "sg_model object")
  expect_error(sg_init(m3, "uniform"), "must be one of")
  expect_error(sg_init(m3, values = c(X = 1)), "Unknown block")
  expect_error(sg_init(m3, values = c(1, 2)), "named numeric")
  expect_error(sg_init(m3, values = c(E = -1)), ">= 0")
  expect_error(sg_init(m3, genotypes = c(ZZ = 1)), "not in the cube")
  expect_error(sg_init(m3, genotypes = c(AA = 0)), "positive sum")
  raw <- sg_test_cube_sit()
  raw$eta["S", ] <- 0 # S females accept no male
  m <- sg_model(p, sg_cube(raw), wx)
  expect_error(sg_init(m, genotypes = c(S = 1)), "cannot be initialised")
})

test_that("print() shows totals by block and genotype", {
  expect_snapshot(print(sg_init(m3, "hardy_weinberg", alleles = c(A = 0.5, a = 0.5))))
})
