# SGDrivE 0.0.0.9000

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
