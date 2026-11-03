# Regenerate the reference outputs used for ScrewDrivE regression tests.
#
# The legacy Gutierrez & Ponti (2014) R translations are copied UNCHANGED into
# data-raw/reference_impl/scripts/ (they source each other via the relative
# path "scripts/...", so they are run with reference_impl/ as the working
# directory). Outputs are written to data-raw/reference/ as compressed .rds
# files, together with a provenance record.
#
# Run from the package root:  Rscript data-raw/make_reference.R
# Requires: deSolve. Takes about 1-2 minutes (ODE runs ~25 s each).

stopifnot(file.exists("DESCRIPTION"))
out_dir <- normalizePath(file.path("data-raw", "reference"), mustWork = TRUE)
impl_dir <- normalizePath(file.path("data-raw", "reference_impl"), mustWork = TRUE)

old_wd <- setwd(impl_dir)
on.exit(setwd(old_wd), add = TRUE)

source("scripts/gp2014_ode.R") # also sources scripts/gutierrez_ponti_2014_pbdm.R

weather_file <- "mcallen_1960_1983_filled.csv"
wx <- utils::read.csv(weather_file)
wx$date <- as.Date(wx$date)
p <- gp_params()

runs <- list(
  # Discrete daily-step distributed delay (line-by-line Pascal translation).
  discrete = function() gp_simulate(wx, p),
  # ODE counterpart, same constants, Pascal-style input floors on.
  ode = function() gp_ode_simulate(wx, p, floors = TRUE),
  # ODE without floors.
  ode_no_floors = function() gp_ode_simulate(wx, p, floors = FALSE),
  # ODE with the extra day of larval mortality of the discrete scheme.
  ode_match_larvae = function() {
    gp_ode_simulate(wx, p,
      floors = TRUE,
      match_discrete_larvae = TRUE
    )
  }
)

timing <- numeric(0)
for (nm in names(runs)) {
  t0 <- proc.time()[["elapsed"]]
  res <- runs[[nm]]()
  timing[nm] <- proc.time()[["elapsed"]] - t0
  saveRDS(res, file.path(out_dir, paste0("gp2014_mcallen_", nm, ".rds")),
    compress = "xz"
  )
  message(sprintf(
    "%-18s %6.1f s  (%d days, %d columns)", nm, timing[nm],
    nrow(res), ncol(res)
  ))
}
saveRDS(p, file.path(out_dir, "gp2014_params.rds"), compress = "xz")

md5 <- tools::md5sum(c(weather_file, file.path("scripts", c(
  "gutierrez_ponti_2014_pbdm.R", "gp2014_ode.R"
))))
prov <- c(
  "# Provenance of data-raw/reference/*.rds",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("deSolve:", as.character(utils::packageVersion("deSolve"))),
  "weather: McAllen TX, NOAA GHCN-Daily USC00415701, gaps filled from",
  "  USW00012959 then interpolation (Gutierrez et al. 2019 Appendix S1 method),",
  "  1960-01-01 to 1983-12-31; first year is burn-in in analyses.",
  "source scripts (unchanged copies of the legacy project scripts):",
  paste0("  ", names(md5), "  md5 ", md5),
  "runs:",
  paste0("  ", names(timing), sprintf("  %.1f s", timing))
)
writeLines(prov, file.path(out_dir, "PROVENANCE.txt"))
message("Wrote ", out_dir)
