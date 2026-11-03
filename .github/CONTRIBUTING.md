# Contributing to ScrewDrivE

ScrewDrivE is in early development. These notes keep the code base consistent.

## Workflow

1. Create a branch from `main` for each change (`git switch -c my-change`).
2. Keep the model specification in `dev/MODEL_SPEC.md` and the decisions log
   in `dev/DECISIONS.md` in step with the code. A change to the model needs a
   matching change to the specification.
3. Before opening a pull request, run locally:
   ```r
   devtools::document()
   devtools::test()
   lintr::lint_package()
   devtools::check()
   ```
   `devtools::check()` must give 0 errors, 0 warnings and 0 notes.
4. Add a line to `NEWS.md` for user-visible changes.

## Conventions

See `dev/CONVENTIONS.md` for naming, units, array index order, testing
tolerances, error messages and comment style.

## Rendering README locally

If pandoc is not on the PATH, `dev/check.R` uses the copy bundled with
RStudio. Run the full local check with `Rscript dev/check.R`.

## Slow tests

Some tests (e.g. the 24-year legacy regression) are skipped unless
`SCREWDRIVE_SLOW_TESTS=true`. Run them locally from the package root with
```r
Sys.setenv(SCREWDRIVE_SLOW_TESTS = "true"); devtools::test()
```
or trigger the `slow-tests` workflow on GitHub.

## Coverage

Coverage is measured locally only:
```r
source("dev/coverage.R")
```
