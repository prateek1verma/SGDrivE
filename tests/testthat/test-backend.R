test_that("package loads and reports its version", {
  info <- sg_backend_info()
  expect_type(info, "list")
  expect_identical(info$package_version, as.character(utils::packageVersion("ScrewDrivE")))
})

test_that("C++ backend is compiled with at least C++17", {
  info <- sg_backend_info()
  expect_true(is.numeric(info$cplusplus))
  expect_gte(info$cplusplus, 201703)
  expect_type(info$compiler, "character")
  expect_true(nzchar(info$compiler))
})

test_that("shared test tolerances are defined and ordered", {
  expect_lt(sg_tol_exact, sg_tol_ode)
  expect_lt(sg_tol_ode, sg_tol_reference)
})
