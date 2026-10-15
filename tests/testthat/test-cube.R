# ---- valid cubes ---------------------------------------------------------------

test_that("a hand-built cube is read and its components are standardised", {
  cb <- sg_cube(sg_test_cube_sit())
  expect_s3_class(cb, "sg_cube")
  expect_identical(cb$genotypes, c("W", "S"))
  expect_identical(cb$n_genotypes, 2L)
  expect_identical(cb$wild_type, "W")
  expect_identical(cb$release_type, "S")
  expect_identical(dim(cb$ih), c(2L, 2L, 2L))
  expect_identical(dimnames(cb$eta), list(c("W", "S"), c("W", "S")))
  expect_identical(cb$ell, c(W = 1, S = 1))
  expect_identical(cb$ell_source, "none")
})

test_that("the egg tensor is s * ih * tau with female genotype varying slowest", {
  raw <- sg_test_cube_mendelian()
  raw$s["Aa"] <- 0.8
  raw$tau["aa", "Aa", "Aa"] <- 0.5
  cb <- sg_cube(raw)
  n <- cb$n_genotypes
  expect_identical(dim(cb$egg_tensor), c(n * n, n))
  expect_identical(rownames(cb$egg_tensor)[1:3], c("AA|AA", "AA|Aa", "AA|aa"))
  for (g in cb$genotypes) {
    for (m in cb$genotypes) {
      expected <- unname(cb$s[g] * cb$ih[g, m, ] * cb$tau[g, m, ])
      row <- (match(g, cb$genotypes) - 1L) * n + match(m, cb$genotypes)
      expect_equal(unname(cb$egg_tensor[row, ]), expected, tolerance = sg_tol_exact)
    }
  }
  expect_equal(unname(cb$egg_total), unname(rowSums(cb$egg_tensor)), tolerance = sg_tol_exact)
  expect_equal(unname(cb$egg_total["aa|Aa"]), 0.5 * 1 + 0.5 * 0.5, tolerance = sg_tol_exact)
})

test_that("sterile-male crosses have zero viable eggs in the SIT test cube", {
  cb <- sg_cube(sg_test_cube_sit(residual = 0))
  expect_identical(unname(cb$egg_total), c(1, 0, 1, 0)) # W|W, W|S, S|W, S|S
})

test_that("named per-genotype vectors are reordered; unnamed ones keep cube order", {
  raw <- sg_test_cube_sit()
  raw$s <- c(S = 0.2, W = 1)
  expect_identical(sg_cube(raw)$s, c(W = 1, S = 0.2))
  raw$s <- c(1, 0.2)
  expect_identical(sg_cube(raw)$s, c(W = 1, S = 0.2))
})

test_that("an eta vector over male genotypes is expanded to a matrix", {
  raw <- sg_test_cube_sit()
  raw$eta <- c(W = 1, S = 0.35)
  cb <- sg_cube(raw)
  expect_identical(unname(cb$eta), matrix(c(1, 1, 0.35, 0.35), 2, 2))
})

test_that("unnamed arrays get genotype dimnames", {
  raw <- sg_test_cube_sit()
  dimnames(raw$ih) <- NULL
  expect_identical(dimnames(sg_cube(raw)$ih)[[3]], c("W", "S"))
})

test_that("MGDrivE cubes, including sex-linked ones, load without warnings", {
  skip_if_not_installed("MGDrivE")
  cubes <- list(
    MGDrivE::cubeMendelian(), MGDrivE::cubeHomingDrive(), MGDrivE::cubeRIDL(),
    MGDrivE::cubeReciprocalTranslocations(), MGDrivE::cubeWolbachia(),
    MGDrivE::cubeXShredderY(), MGDrivE::cubeSplitDriveX()
  )
  for (raw in cubes) {
    expect_no_warning(cb <- sg_cube(raw))
    expect_identical(cb$genotypes, raw$genotypesID)
    expect_identical(cb$n_genotypes, length(raw$genotypesID))
  }
})

test_that("the hand-built Mendelian cube equals MGDrivE::cubeMendelian()", {
  skip_if_not_installed("MGDrivE")
  a <- sg_cube(sg_test_cube_mendelian())
  b <- sg_cube(MGDrivE::cubeMendelian())
  expect_equal(a$ih, b$ih, tolerance = sg_tol_exact)
  expect_equal(a$egg_tensor, b$egg_tensor, tolerance = sg_tol_exact)
})

test_that("sex-linked cubes allow zero ih slices for impossible crosses", {
  skip_if_not_installed("MGDrivE")
  raw <- MGDrivE::cubeXShredderY()
  sums <- apply(raw$ih, c(1, 2), sum)
  expect_true(any(sums == 0))
  expect_no_error(sg_cube(raw))
})

# ---- lifespan fractions (D-002, D-012) -------------------------------------

test_that("sg_omega_to_ell() inverts MGDrivE::calcOmega() exactly", {
  skip_if_not_installed("MGDrivE")
  for (mu in c(0.02, 0.055, 0.2)) {
    for (L in c(mu + 0.01, 0.3, 0.5, 0.9, 1)) {
      om <- MGDrivE::calcOmega(mu = mu, lifespanReduction = L)
      expect_equal(sg_omega_to_ell(om, mu), L, tolerance = sg_tol_exact)
    }
  }
})

test_that("sg_omega_to_ell() handles omega = 0, 1 and invalid input", {
  expect_identical(sg_omega_to_ell(c(a = 1, b = 0), 0.055), c(a = 1, b = 0.001))
  expect_error(sg_omega_to_ell(1, 0), "reference_mortality")
  expect_error(sg_omega_to_ell(1, 1), "reference_mortality")
  expect_error(sg_omega_to_ell(1, c(0.1, 0.2)), "reference_mortality")
  expect_error(sg_omega_to_ell(-1, 0.1), ">= 0")
  expect_error(sg_omega_to_ell(1.2, 0.1), "no adult ever dies")
})

test_that("cube omega is converted with reference_mortality, which is required", {
  raw <- sg_test_cube_sit(omega_s = 0.9417989)
  expect_error(sg_cube(raw), "give `reference_mortality`")
  cb <- sg_cube(raw, reference_mortality = 0.055)
  expect_identical(cb$ell_source, "omega")
  expect_identical(cb$mu_ref, 0.055)
  expect_equal(unname(cb$ell), c(1, 0.5), tolerance = 1e-6)
})

test_that("user-given ell is used as is and overrides omega with a warning", {
  expect_identical(
    sg_cube(sg_test_cube_sit(), lifespan_fraction = c(S = 0.5, W = 1))$ell,
    c(W = 1, S = 0.5)
  )
  raw <- sg_test_cube_sit(omega_s = 0.5)
  expect_warning(cb <- sg_cube(raw, lifespan_fraction = c(1, 0.7)), "omega values are not used")
  expect_identical(cb$ell_source, "user")
  expect_identical(unname(cb$ell), c(1, 0.7))
  expect_error(sg_cube(sg_test_cube_sit(), lifespan_fraction = c(1, 0)), "must be > 0")
  expect_error(sg_cube(sg_test_cube_sit(), lifespan_fraction = 1), "one value per genotype")
})

# ---- invalid cubes -----------------------------------------------------------------

test_that("structural problems give specific errors", {
  expect_error(sg_cube(1:3), "must be a list")
  raw <- sg_test_cube_sit()
  expect_error(sg_cube(raw[setdiff(names(raw), "tau")]), "missing element.*tau")
  bad <- raw
  bad$genotypesID <- c("W", "W")
  expect_error(sg_cube(bad), "duplicated")
  bad <- raw
  bad$genotypesN <- 3L
  expect_error(sg_cube(bad), "genotypesN")
  bad <- raw
  bad$ih <- bad$ih[, , 1, drop = FALSE]
  expect_error(sg_cube(bad), "dimension 2 x 2 x 2")
  bad <- raw
  dimnames(bad$ih)[[2]] <- c("S", "W")
  expect_error(sg_cube(bad), "Dimension 2 names of cube\\$ih")
  bad <- raw
  bad$eta <- matrix(1, 3, 3)
  expect_error(sg_cube(bad), "cube\\$eta must be")
  bad <- raw
  bad$phi <- c(0.5, 0.5, 0.5)
  expect_error(sg_cube(bad), "one value per genotype")
  bad <- raw
  bad$phi <- c(X = 0.5, S = 0.5)
  expect_error(sg_cube(bad), "Names of cube\\$phi")
  bad <- raw
  bad$wildType <- "Z"
  expect_error(sg_cube(bad), "wildType contains genotype")
  bad <- raw
  bad$releaseType <- "Z"
  expect_error(sg_cube(bad), "releaseType contains genotype")
  bad <- raw
  bad$genotypesID <- character(0)
  expect_error(sg_cube(bad), "non-empty character")
})

test_that("value problems give specific errors and warnings", {
  raw <- sg_test_cube_sit()
  bad <- raw
  bad$ih["W", "W", "W"] <- 0.9
  expect_error(sg_cube(bad), "slices must sum to 1.*ih\\[\"W\", \"W\", \\]")
  bad <- raw
  bad$ih["W", "W", "S"] <- -0.1
  expect_error(sg_cube(bad), "finite values >= 0")
  bad <- raw
  bad$xiF["S"] <- 1.5
  expect_error(sg_cube(bad), "xiF must contain finite values in \\[0, 1\\]")
  bad <- raw
  bad$omega["S"] <- -1
  expect_error(sg_cube(bad), "omega must contain")
  bad <- raw
  bad$eta[1, 2] <- NA
  expect_error(sg_cube(bad), "eta must contain finite")
  bad <- raw
  bad$tau["W", "W", "W"] <- 1.5
  expect_warning(sg_cube(bad), "tau .* exceed 1")
  bad <- raw
  bad$ih["W", "S", ] <- 0 # a possible cross with no offspring
  expect_warning(sg_cube(bad), "possible cross\\(es\\) produce no offspring")
})

test_that("release type is optional", {
  raw <- sg_test_cube_sit()
  raw$releaseType <- NULL
  expect_identical(sg_cube(raw)$release_type, character(0))
})

# ---- printing --------------------------------------------------------------------

test_that("print() summarises the cube", {
  raw <- sg_test_cube_sit(residual = 1.2e-8)
  raw$eta[, "S"] <- 0.35
  raw$s["S"] <- 2.3e-7
  expect_snapshot(print(sg_cube(raw, lifespan_fraction = c(W = 1, S = 0.52))))
  expect_invisible(print(sg_cube(sg_test_cube_sit())))
})
