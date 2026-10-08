# Local full check: document, test, lint and R CMD check.
# Run from the package root:  Rscript dev/check.R
#
# * Uses dev/Makevars.local to silence a warning that comes from R 4.5.0's own
#   headers with Apple clang 21 (see that file).
# * Uses RStudio's bundled pandoc if pandoc is not on the PATH (needed to
#   render README.Rmd and, later, vignettes).
Sys.setenv(
  R_MAKEVARS_USER = normalizePath("dev/Makevars.local"),
  # declare the local-only flag as known so R CMD check does not report it
  "_R_CHECK_COMPILATION_FLAGS_KNOWN_" = "-Wno-unknown-warning-option"
)
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
if (!nzchar(Sys.which("pandoc")) && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

devtools::document(".")
devtools::build_readme(".", quiet = TRUE)
lints <- lintr::lint_package(".")
print(lints)
res <- devtools::check(".", document = FALSE, error_on = "never")
print(res)
ok <- length(res$errors) == 0 && length(res$warnings) == 0 &&
  length(res$notes) == 0 && length(lints) == 0
if (!ok) quit(status = 1)
