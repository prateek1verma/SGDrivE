# SGDrivE 0.0.0.9000

## Milestone M1: inputs (in progress)

* `sg_params()` creates validated model parameters, with presets
  `"default"` and `"gp2014"` (the Gutierrez & Ponti 2014 code
  configuration); `update()`, `print()` and `as.data.frame()` methods.
* `sg_param_table()` returns the definition of every parameter: default,
  allowed range, units, source and model equation.
* `sg_cube()` reads and validates MGDrivE inheritance cubes (including
  sex-linked cubes, whose impossible crosses have zero inheritance slices),
  converts omega to lifespan fractions and precomputes the egg tensor
  s * ih * tau.
* `sg_omega_to_ell()` inverts `MGDrivE::calcOmega()`.

## Milestone M0: foundations

* Package skeleton: DESCRIPTION, GPL-3 licence, Rcpp/C++17 build setup.
* `sg_backend_info()`: toolchain check reporting how the C++ backend was
  compiled. No model code yet.
* Testing (testthat 3rd edition), linting (lintr, styler, clang-format) and
  GitHub Actions workflows (R CMD check on macOS, Windows and Ubuntu; lint).
* Design documents in `dev/`: decisions log, v0.1 model specification,
  roadmap, conventions.
* Reference outputs of the legacy Gutierrez & Ponti (2014) R translations in
  `data-raw/reference/` for later regression tests.
