# Local full check: document, render README, lint and R CMD check.
# Run from the package root:  Rscript dev/check.R
#
# Uses RStudio's bundled pandoc if pandoc is not on the PATH (needed to render
# README.Rmd and, later, vignettes).
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
if (!nzchar(Sys.which("pandoc")) && dir.exists(rstudio_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

devtools::document(".")
devtools::build_readme(".", quiet = TRUE)
# lintr resolves internal functions in the installed package, so install the
# current source first (CI does the same).
devtools::install(".", quick = TRUE, upgrade = FALSE, quiet = TRUE)
lints <- lintr::lint_package(".")
print(lints)
res <- devtools::check(".", document = FALSE, error_on = "never")
print(res)
ok <- length(res$errors) == 0 && length(res$warnings) == 0 &&
  length(res$notes) == 0 && length(lints) == 0
if (!ok) quit(status = 1)
