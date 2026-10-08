# Local test coverage report (coverage is not uploaded anywhere; D-008).
# Run from the package root:  source("dev/coverage.R")
cov <- covr::package_coverage(".", type = "tests")
print(cov)
covr::report(cov, file = "coverage.html", browse = interactive())
