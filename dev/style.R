# Apply the package code style (tidyverse style via styler).
# Run from the package root:  Rscript dev/style.R
#
# data-raw/reference_impl/ is excluded: it holds unchanged copies of the
# legacy Gutierrez & Ponti R translations used as the regression baseline
# (their md5 sums are recorded in data-raw/reference/PROVENANCE.txt).
styler::style_pkg(
  ".",
  exclude_dirs = c("data-raw/reference_impl", "packrat", "renv"),
  exclude_files = "R/RcppExports.R"
)
