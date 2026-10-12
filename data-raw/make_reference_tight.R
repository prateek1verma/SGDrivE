# Tight-tolerance legacy reference runs for the M2 regression test (D-019).
#
# Runs the unchanged legacy ODE translation of Gutierrez & Ponti (2014)
# (data-raw/reference_impl/scripts/gp2014_ode.R) WITHOUT population floors and
# with tight solver tolerances, so that the comparison with SGDrivE's
# reference solver at 1e-6 is not limited by solver error. The existing
# reference files (make_reference.R) are not modified.
#
# Outputs:
#   data-raw/reference/gp2014_mcallen_ode_no_floors_tight.rds  (1960-1983)
#   data-raw/reference/PROVENANCE_tight.txt
#   tests/testthat/fixtures/gp2014_mcallen_ode_tight_3yr.rds    (1960-1962,
#     daily totals only, used by the routine regression test)
#   tests/testthat/fixtures/gp2014_mcallen_ode_tight_24yr_weekly.rds
#     (1960-1983, totals every 7th day, used by the routine 24-year regression
#     test of the C++ backend)
#
# Run from the package root:  Rscript data-raw/make_reference_tight.R
# Takes a few minutes.

stopifnot(file.exists("DESCRIPTION"))
pkg_root <- normalizePath(".")
out_dir <- file.path(pkg_root, "data-raw", "reference")
fix_dir <- file.path(pkg_root, "tests", "testthat", "fixtures")
impl_dir <- file.path(pkg_root, "data-raw", "reference_impl")

old_wd <- setwd(impl_dir)
on.exit(setwd(old_wd), add = TRUE)
source("scripts/gp2014_ode.R") # also sources scripts/gutierrez_ponti_2014_pbdm.R

weather_file <- "mcallen_1960_1983_filled.csv"
wx <- utils::read.csv(weather_file)
wx$date <- as.Date(wx$date)
p <- gp_params()

rtol_tight <- 1e-10
atol_tight <- 1e-12

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, time = proc.time()[["elapsed"]] - t0)
}

# 1. Full period, tight tolerances, no floors.
full <- timed(gp_ode_simulate(wx, p, floors = FALSE, rtol = rtol_tight, atol = atol_tight))
saveRDS(full$value, file.path(out_dir, "gp2014_mcallen_ode_no_floors_tight.rds"),
  compress = "xz"
)
message(sprintf("full 1960-1983, rtol %.0e: %.1f s", rtol_tight, full$time))

# 2. Convergence check of the legacy solution on the 3-year window: rerun with
#    tolerances 10 times tighter and record the largest relative difference.
#    (rtol 1e-12 exceeds lsoda's default step limit because the legacy code
#    integrates across the daily changes of forcing without restarting; its
#    function does not expose maxsteps.)
wx3 <- wx[wx$date <= as.Date("1962-12-31"), ]
tighter <- timed(gp_ode_simulate(wx3, p, floors = FALSE, rtol = 1e-11, atol = 1e-13))
cols <- c(
  "eggs", "larvae", "pupae", "females", "adults", "eggs_laid", "egg_batches",
  "pupae_in", "adults_in"
)
full3 <- full$value[seq_len(nrow(wx3)), ]
rel_diff <- vapply(cols, function(v) {
  max(abs(full3[[v]] - tighter$value[[v]])) / max(abs(tighter$value[[v]]))
}, 0)
message("legacy convergence check (rtol 1e-10 vs 1e-11), max relative difference:")
print(signif(rel_diff, 3))

# 3. Compact fixture for routine tests: daily totals 1960-1962.
fixture <- full3[, c("date", cols)]
saveRDS(fixture, file.path(fix_dir, "gp2014_mcallen_ode_tight_3yr.rds"), compress = "xz")

# 4. Compact 24-year fixture: totals every 7th day 1960-1983 (68 kB instead of
#    450 kB for daily values; any drift of the solution still shows).
weekly <- full$value[seq(7L, nrow(full$value), by = 7L), c("date", cols)]
saveRDS(weekly, file.path(fix_dir, "gp2014_mcallen_ode_tight_24yr_weekly.rds"), compress = "xz")

md5 <- tools::md5sum(c(weather_file, file.path("scripts", c(
  "gutierrez_ponti_2014_pbdm.R", "gp2014_ode.R"
))))
prov <- c(
  "# Provenance of the tight-tolerance legacy reference (make_reference_tight.R)",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("deSolve:", as.character(utils::packageVersion("deSolve"))),
  "run: gp_ode_simulate(floors = FALSE), lsoda over the whole period (hmax 0.25),",
  sprintf("  rtol %.0e, atol %.0e", rtol_tight, atol_tight),
  "weather: McAllen TX 1960-01-01 to 1983-12-31 (see PROVENANCE.txt)",
  "source scripts (unchanged copies of the legacy project scripts):",
  paste0("  ", names(md5), "  md5 ", md5),
  sprintf("run time: full %.1f s; 3-year check at rtol 1e-11 %.1f s", full$time, tighter$time),
  "legacy convergence, 1960-1962, max relative difference rtol 1e-10 vs 1e-11:",
  paste0("  ", names(rel_diff), "  ", signif(rel_diff, 3)),
  "fixture: tests/testthat/fixtures/gp2014_mcallen_ode_tight_3yr.rds (daily totals 1960-1962)",
  "fixture: tests/testthat/fixtures/gp2014_mcallen_ode_tight_24yr_weekly.rds",
  "  (totals every 7th day 1960-1983)"
)
writeLines(prov, file.path(out_dir, "PROVENANCE_tight.txt"))
message("Wrote ", out_dir, " and ", fix_dir)
